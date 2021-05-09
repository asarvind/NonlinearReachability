# include "linearization.cpp"
# include "newsetreps.cpp"

//****************************************************************************************************
// Intersection of Unions reach set class
//****************************************************************************************************
class ioureach: public nonlinear{
public:
  //======================================================================
  // attributes
  //======================================================================

  // union storage capacity of iou array
  static const int MaxDivs = 128;
  
  // initial set
  IvVectorNd InitState;
  
  // current bounds on the reachable set
  IvVectorNd bounds;
  
  // maximum bounds on the reachable set in the simulated time horizon
  IvVectorNd MaxBounds;
  
  // time interval;
  Interval SimTime;
    
  // log of number of elements in union
  int LogDivs;
  
  // number of elements in union
  int intrs;

  // number of elements in union
  int ivintrs;  

  // order of zonotope
  int zonOrder;

  // parallelotopic template and inverse
  IvMatrixNNd ptemp, invptemp;

  // parallelotopic bounds
  IvVectorNd pbounds;

  // discrete state action matrix at origin and its inverse
  IvMatrixNNd stmatorigin, invmatorigin;

  // iou of parallelotope bounds
  IvVectorNd piou[StateDim][MaxDivs];
    
  // array of optimal divisions vectors
  IntVectorNd DivVecs[StateDim];
      
  /* iou array is 2-dimensional array of zonotopes.
     First dimension contains intersections. Second dimension
     contains unions. */
  ZonNd iou[StateDim][MaxDivs];

  // iou of interval vectors
  vector<IvVectorNd> iviou[StateDim];

  //----------------------------------------------------------------------
  // Flowpipe
  //----------------------------------------------------------------------
  struct FlowElem{
    IvVectorNd bounds;
    Interval time;
  };
  vector<FlowElem> FlowPipe;

  // current iterator pointer of flowpipe
  vector<FlowElem>::iterator FlowItr;

  // number of simulations completed
  int SimNum;

  // boolean variable for doing bloating
  bool doBloat;

  //----------------------------------------------------------------------
  // constructor
  //----------------------------------------------------------------------
  ioureach(const IvVectorNd &State, const IvVectorMd &Input, const IvVectorKd &parvals, double tstep, int k)
    :nonlinear(Input, parvals){
    InitState = State;
    TimeStep = tstep;
    zonOrder = 1;
    LogDivs = k;
    SimNum = 0;
    doBloat = true;

    // resize iviou union size and set parallelotope template
    {
      ptemp *= 0; invptemp *= 0;
      int divs = pow(2,LogDivs);
#pragma omp parallel for
      for(int i=0; i<N; ++i){
	iviou[i].resize(divs);
	ptemp(i,i) = 1; invptemp(i,i) = 1;
      }
    }
    // set discrete time state action matrix at origin
    LinVals L;
    L.state = State;
    DisLin( L, false );
    stmatorigin = L.StMatDis;
    stmatorigin /= stmatorigin.determinant();
    invmatorigin = stmatorigin.inverse();
  }

  //======================================================================
  // Methods
  //======================================================================


  //----------------------------------------------------------------------
  /* Compute optimal division vectors along different directions 
     and the number of intersections */
  //----------------------------------------------------------------------
  void SetDivVecs(){
    intrs = 0;
    bool unique;

    IntVectorNd allDivVecs[StateDim];
    
    {
      #pragma omp parallel for
      for( int i=0; i<N; ++i){
	OptErr E = OptDivision(bounds,EvRe[i],EvIm[i],LogDivs);
	allDivVecs[i] = E.divs;
      }
    }
    
    for(int i=0; i<N; ++i){
      //OptErr E = OptDivision(bounds,EvRe[i],EvIm[i],LogDivs);
      DivVecs[intrs] = allDivVecs[i];
      // retain unique vectors in the list of optimum division vectors
      // vectors should not be ones
      unique = true;
      unique = unique && (DivVecs[intrs].lpNorm<Eigen::Infinity>() > 1);
      for(int l=0; l<intrs; ++l){
	unique = unique && ((DivVecs[intrs]-DivVecs[l]).lpNorm<Eigen::Infinity>() > 0);
      }
      if(unique){
	intrs += 1;
      }
    }
    if(intrs==0){
      for(int i=0; i<N; i++){
	DivVecs[intrs](i) = 1;
	intrs = 1;
      }
    }
  }

  void SetNewDivVecs(){
    ivintrs = 0;
    bool unique;

    IntVectorNd allDivVecs[StateDim];
    
    {
      #pragma omp parallel for
      for( int i=0; i<N; ++i){
	VectorNd vect;
	vect *= 0;
	vect(i) += 1;
	VectorNd imvect = vect*0;
	OptErr E = OptDivision(bounds,vect,imvect,LogDivs);
	// OptErr E = OptDivision(bounds,EvRe[i],EvIm[i],LogDivs);
	allDivVecs[i] = E.divs;
      }
    }
    
    for(int i=0; i<N; ++i){
      //OptErr E = OptDivision(bounds,EvRe[i],EvIm[i],LogDivs);
      DivVecs[ivintrs] = allDivVecs[i];
      // retain unique vectors in the list of optimum division vectors
      // vectors should not be ones
      unique = true;
      unique = unique && (DivVecs[ivintrs].lpNorm<Eigen::Infinity>() > 1);
      for(int l=0; l<ivintrs; ++l){
	unique = unique && ((DivVecs[ivintrs]-DivVecs[l]).lpNorm<Eigen::Infinity>() > 0);
      }
      if(unique){
	ivintrs += 1;
      }
    }
    if(ivintrs==0){
      for(int i=0; i<N; i++){
	DivVecs[ivintrs](i) = 1;
	ivintrs = 1;
      }
    }
  }

  
  //----------------------------------------------------------------------
  // Method: compute iou representation
  //----------------------------------------------------------------------
  void SetIou(){
    SetDivVecs();
    // divide state into iou
    int divs = pow(2,LogDivs);
    
# pragma omp parallel for collapse(2)
    for(int i=0; i<intrs; i++){
      for(int j=0; j<divs; j++){
	int q, d;
	q = 0;
	d = 0;
	Interval iv, delta;
	IvVectorNd addvect;
	// reset zonotope
	iou[i][j] = ZonNd();
	q = j;	
	for(int k=0; k<N; k++){
	  d = DivVecs[i](k);
	  iv = bounds(k);
	  delta = (Interval(iv.upper(),iv.upper())-Interval(iv.lower(),iv.lower()))/d;
	  addvect(k) = hull(iv.lower() + (q%d)*delta,iv.lower() + ((q%d) + 1)*delta);
	  q /= d;
	}
	// set zonotope value
	iou[i][j].MinSum(addvect);
	iou[i][j].reduceOrder(1);
	iou[i][j].setBounds();
      }
    }
  }

  //----------------------------------------------------------------------
  // Method: compute iou representation
  //----------------------------------------------------------------------
  void SetIvIou(){
    SetNewDivVecs();
    // divide state into iou
    int divs = pow(2,LogDivs);
#pragma parallel omp for collapse(2);
    for(int i=0; i<ivintrs; i++){
      for(int j=0; j<divs; j++){
	int q, d;
	q = 0;
	d = 0;
	Interval iv, delta;
	q = j;
	IvVectorNd addvect;
	for(int k=0; k<N; k++){
	  d = DivVecs[i](k);
	  iv = bounds(k);
	  delta = (Interval(iv.upper(),iv.upper())-Interval(iv.lower(),iv.lower()))/d;
	  addvect(k) = hull(iv.lower() + (q%d)*delta,iv.lower() + ((q%d) + 1)*delta);
	  q /= d;
	}
	// assign interval vector value after finding linearization region
	LinVals L;
	L.state = addvect;
	LinRegion(L);
	iviou[i][j] = L.region;
      }      
    }
  }
  

  //----------------------------------------------------------------------
  // Compute bounds from iou
  //----------------------------------------------------------------------
  void SetBounds(){
    IvVectorNd U, V;
    int divs = pow(2,LogDivs);
    // initialize bounds
    for( int i=0; i<N; ++i ){
      bounds(i) = numeric_limits<double>::infinity()*Interval(-1,1);
      pbounds(i) = numeric_limits<double>::infinity()*Interval(-1,1);
    }
    // recursively get U and intersect with bounds
    for(int i=0; i<intrs; i++){
      U = iou[i][0].bounds;
      V = piou[i][0];
      for(int j=1; j<divs; j++){
	U = join(U, iou[i][j].bounds);
	V = join( V, piou[i][j] );
      }
      bounds = meet(bounds,U);
      pbounds = meet( pbounds, V );
    }
    MaxBounds = join( MaxBounds, bounds );
  }
  

  //----------------------------------------------------------------------
  // bloating initial state
  //----------------------------------------------------------------------
  void bloat(){
    LinVals L;
    L.state = InitState;
    bounds = InitState;
    LinRegion( L );
    Interval delta = Interval(0,TimeStep);
    bounds = L.region;
    pbounds = L.region;
    MaxBounds = bounds;
    SimTime = delta;
  }

  //----------------------------------------------------------------------
  // next zonotope obtained from multiple piecewise linearization
  //----------------------------------------------------------------------
  void multilin( LinVals &L ){
    int divs = pow(2,LogDivs); // no. of divisions
    
    // declare array of error vectors after linearization and abstraction
    vector<IvVectorNd> errvect[StateDim];

    {
    // calculate error array
      //#pragma omp parallel for
    for(int i=0; i<ivintrs; ++i){
      errvect[i].resize(0);
      //#pragma omp parallel for
      for(int j=0; j<divs; ++j){
	if( is_overlap( iviou[i][j], L.region ) ){
	  LinVals newL;
	  newL.state = meet( iviou[i][j], L.region );
	  for(int r=0; r<N; ++r){
	  }
	  DisLin( newL, false );
	  errvect[i].push_back( (newL.StMatDis - L.StMatDis)*newL.state + (newL.InpMatDis - L.InpMatDis)*Inp + newL.ErrDis );
	}
      }
    }
    }
    
    // compute overall linearization error
    IvVectorNd linerr = L.ErrDis;
    for(int i=0; i<ivintrs; ++i){
      IvVectorNd U = errvect[i][0];
      for(int j=1; j<errvect[i].size(); ++j){
    	U = join( U, errvect[i][j] );
      }
      linerr = meet( U, linerr );
    }

    // reset linearization error
    L.ErrDis = linerr;
  }
  

  void nextIou(){
    // SetIvIou();
    int divs = pow(2,LogDivs);
#pragma omp parallel for collapse(2)
    for(int j=0; (j<intrs); ++j){
      for(int k=0; k<divs; ++k){
	LinVals L;
	L.state = iou[j][k].bounds;
	L.state = meet(L.state,bounds);
	DisLin(L,true);
	iou[j][k].prod(L.StMatDis);
	IvVectorNd addvect = L.InpMatDis*Inp + L.ErrDis;
	iou[j][k].MinSum( addvect );
	iou[j][k].reduceOrder(zonOrder);
	iou[j][k].setBounds();
	piou[j][k] = iou[j][k].getBounds( ptemp );
      }
    }
    SetBounds();
  }

  //----------------------------------------------------------------------
  // transform paralellotope templates
  //----------------------------------------------------------------------
  void setptemp(){
    // set invptemp
    invptemp = stmatorigin*invptemp;
    // set ptemp
    ptemp = ptemp*invmatorigin;
  }

  //----------------------------------------------------------------------
  // refine iou of zonotopes
  //----------------------------------------------------------------------
  void refine(){
#pragma omp parallel for collapse(2)
    for(int i=0; i<intrs; ++i){
      for(int j=0; j<intrs; ++j){
	iou[i][j].refine( piou[i][j], pbounds, ptemp, invptemp );
      }
    }
  }
  
  //----------------------------------------------------------------------
  // flowpipe computation for a time period
  //----------------------------------------------------------------------
  void flow(double T, double MaxFlowTime){
    // set iou
    SetIou();
    
    double FlowTime = 0;
    bool do_iter = true;
    
    while(do_iter){
      // set parallelotope template
      //setptemp();
      
      // compute next iou
      nextIou();

      // refine iou
      //refine();
      
      // update clocks
      SimTime += TimeStep;
      cout << SimTime.upper() << "\n";
      FlowTime += TimeStep;

      // update flowpipe and bounds
      ++FlowItr;
      (*(FlowItr)).bounds = bounds;
      MaxBounds = join(bounds,MaxBounds);
      
      // update do_iter
      do_iter = do_iter && SimTime.upper()<T;
      do_iter = do_iter && FlowTime<MaxFlowTime;
    }        
  }


  //----------------------------------------------------------------------
  // Simulation
  //----------------------------------------------------------------------
  void simulate(const double T, double IouResetTime){

    // initialize flowpipe on first simulation
    if(SimNum==0){
      IvVectorNd InfSt;
      for(int i=0; i<N; ++i){
	InfSt(i) = numeric_limits<double>::infinity()*Interval(-1,1);
      }
      Interval flowInitTime = Interval(0,TimeStep);
      FlowElem S = {InfSt, SimTime};
      FlowPipe.push_back(S);
      while(flowInitTime.upper()<T){
	flowInitTime += TimeStep;
	S.time = flowInitTime;
	FlowPipe.push_back(S);
      }
    }
    
    // begin flowpipe 
    FlowItr = FlowPipe.begin();
    if ( doBloat ){
      bloat();
    }
    (*(FlowItr)).bounds = bounds;

    /* perform recursive flowpipe computation with intermediate
       iou resets */
    while(SimTime.upper()<T){
      // compute flowpipe until IOU reset
      flow(T,IouResetTime);
    }
  }

  void simulate(double T){
    omp_set_num_threads(64);
    omp_set_nested(3);
    simulate(T,T);
  }

  //----------------------------------------------------------------------
  // Save traces
  //----------------------------------------------------------------------
  void SaveTraces(){
    ofstream myfile;
    // save lower bounds
    myfile.open("results/lb.txt");
    for(vector<FlowElem>::iterator itr=FlowPipe.begin(); itr!=FlowPipe.end(); ++itr){
      for(int i=0; i<N; ++i){
	if(i==N-1){
	  myfile << ((*itr).bounds)(i).lower() << "\n";
	}
	else{
	  myfile << ((*itr).bounds)(i).lower() << ",";
	}
      }
    }
    myfile.close();
    // save upper bounds
    myfile.open("results/ub.txt");
    for(vector<FlowElem>::iterator itr=FlowPipe.begin(); itr!=FlowPipe.end(); ++itr){
      for(int i=0; i<N; ++i){
	if(i==N-1){
	  myfile << ((*itr).bounds)(i).upper() << "\n";
	}
	else{
	  myfile << ((*itr).bounds)(i).upper() << ",";
	}
      }
    }
    myfile.close();
    // save time interval trace
    myfile.open("results/times.txt");
    for(vector<FlowElem>::iterator itr=FlowPipe.begin(); itr!=FlowPipe.end(); ++itr){
      myfile << ((*itr).time).lower() << "," << ((*itr).time).upper() << "\n";
    }
    myfile.close();
  }


  
  // close ioureach class    
};
