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
  
  // number of elements in intersection
  int intrs;

  // number of elements in intersection
  int ivintrs;  

  // parallelotopic template and inverse
  IvMatrixNNd ptemp, dirtemp;

  // parallelotopic bounds
  IvVectorNd pbounds;
  IvVectorNd dirbounds;

  // discrete state action matrix at origin and its inverse
  IvMatrixNNd stmatorigin, invmatorigin;

  // iou of parallelotope bounds
  IvVectorNd piou[3*StateDim][MaxDivs];
    
  // array of optimal divisions vectors
  IntVectorNd DivVecs[100];
      
  /* iou array is 2-dimensional array of zonotopes.
     First dimension contains intersections. Second dimension
     contains unions. */
  vector<uzonNd> iou;

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

    ifstream dirfile( "src/pywrite/dirmat.txt" );
    for( int i=0; i<N; ++i ){
      for( int j=0; j<N; ++j ){
	double val;
	dirfile >> val;
	ptemp(i,j) = Interval(val, val);
	double valnext = EvRe[j](i);
	dirtemp(i,j) = Interval( valnext, valnext );
      }
    }
    dirfile.close();
    srand(0);
  }

  //======================================================================
  // Methods
  //======================================================================


  //----------------------------------------------------------------------
  /* Compute optimal division vectors along different directions 
     and the number of intersections */
  //----------------------------------------------------------------------
  void SetDivVecs( vector<uzonNd>&newiou, bool scalarobjective = false ){
    intrs = 0;
    bool unique;
    newiou.clear();

    IntVectorNd allDivVecs[100];
    VectorNd rev[100], imv[100];

    {
#pragma omp parallel for
      for( int i=0; i<N; ++i){
	if (scalarobjective){
	    OptErr E = OptDivision(bounds,LogDivs);
	    allDivVecs[i] = E.divs;
	  }
	else{
	  OptErr E = OptDivision(bounds,EvRe[i],EvIm[i],LogDivs);
	  allDivVecs[i] = E.divs;
	}
      }
    }

    for(int i=0; i<N; ++i){
      DivVecs[intrs] = allDivVecs[i];
      // retain unique vectors in the list of optimum division vectors
      // vectors should not be ones
      unique = true;
      unique = unique && (DivVecs[intrs].lpNorm<Eigen::Infinity>() > 1);
      for(int l=0; l<intrs; ++l){
	unique = unique && ((DivVecs[intrs]-DivVecs[l]).lpNorm<Eigen::Infinity>() > 0);
      }
      if(unique){
	uzonNd uz;
	//uz.flag(i) = true;
	newiou.push_back( uz );
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
  
  //----------------------------------------------------------------------
  // Method: compute iou representation
  //----------------------------------------------------------------------
  void SetIou(){
    // divide state into iou
    vector<uzonNd> newiou;
    SetDivVecs(newiou, false);
    int divs = pow(2,LogDivs);
        
# pragma omp parallel for
    for(int i=0; i<newiou.size(); i++){
      newiou[i].sets.resize(divs); // reset number of zonotopes in the union
#pragma omp parallel for
      for(int j=0; j<divs; j++){
	int q, d;
	q = 0;
	d = 0;
	Interval iv, delta;
	IvVectorNd addvect;
	// reset zonotope
	q = j;	
	for(int k=0; k<N; k++){
	  d = DivVecs[i](k);
	  iv = bounds(k);
	  delta = (Interval(iv.upper(),iv.upper())-Interval(iv.lower(),iv.lower()))/d;
	  addvect(k) = hull(iv.lower() + (q%d)*delta,iv.lower() + ((q%d) + 1)*delta);
	  q /= d;
	}
	// set zonotope value
	newiou[i].sets[j].MinSum(addvect);
	newiou[i].sets[j].reduceOrder(1);
	newiou[i].sets[j].setBounds();
      }
    }
    for(int i=0; i<newiou.size(); ++i)
      {
	iou.push_back( newiou[i] );
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
    for(int i=0; i<iou.size(); i++){
      U = iou[i].sets[0].bounds;
      V = piou[i][0];
      for(int j=1; j<divs; j++){
	U = join(U, iou[i].sets[j].bounds);
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
    pbounds = ptemp*L.region;
    MaxBounds = bounds;
    SimTime = delta;
  }  

  //----------------------------------------------------------------------
  // next IoU
  //----------------------------------------------------------------------
  void nextIou(){
    int divs = pow(2,LogDivs);
    int intersections = iou.size();
#pragma omp parallel for collapse(2)
    for(int j=0; ( j<intersections ); ++j){
      for(int k=0; k<divs; ++k){
	LinVals L;
	L.state = iou[j].sets[k].bounds;	
	L.state = meet(L.state,bounds);
	DisLin(L,true);
	iou[j].sets[k].prod(L.StMatDis);
	//IvVectorNd addvect = L.InpMatDis*Inp + L.ErrDis;
	IvVectorNd addvect = L.ErrDis;
	iou[j].sets[k].MinSum( addvect );
	iou[j].sets[k].setBounds();
	iou[j].sets[k].newReduceOrder(zonOrder);	
	//IvVectorNd rectbounds = L.StMatDis*L.state + addvect;
	//iou[j].sets[k].bounds = meet( iou[j].sets[k].bounds, rectbounds );		
	piou[j][k] = iou[j].sets[k].getBounds( ptemp );
	IvVectorNd dirbounds = ptemp*iou[j].sets[k].bounds;
	piou[j][k] = meet( piou[j][k], dirbounds );

	// update bounds
	if( k==0 )
	  {
	    iou[j].bounds = iou[j].sets[k].bounds;
	  }
	else
	  {
	    iou[j].bounds = join( iou[j].sets[k].bounds, iou[j].bounds );
	  }
      }
    }
    SetBounds();
  }

  //----------------------------------------------------------------------
  // refine iou of zonotopes
  //----------------------------------------------------------------------
  static bool sortuzon(uzonNd Z1, uzonNd Z2)
  {
    bool out;
    if(Z1.cost<=Z2.cost)
      {
	out = true;
      }
    else
      {
	out = false;
      }
    return out;
  }

  void refine()
  {
    int maxintrs = N;
    if(iou.size()<=maxintrs)
      {
	return;
      }
    
    vector<uzonNd> newiou = iou;
    iou.clear();
    IvVectorNd ioubounds;
    bool satbounds = false;

    for(int i=0; not satbounds && iou.size()<maxintrs; ++i)
      {
	int iousize = newiou.size();
#pragma omp parallel for
	for(int j=0; j<iousize; ++j)
	  {
	    IvVectorNd b;
	    // compute cost of union zonotope in newiou
	    newiou[j].cost = 0;
	    if(i==0)
	      {
		b = newiou[j].bounds;
	      }
	    else
	      {
		b = meet(ioubounds, newiou[j].bounds);
	      }
	    for(int k=0; k<N; ++k)
	      {
		newiou[j].cost += pow( width(b(k))/( width(bounds(k)) + 1e-10 ), 1 );
	      }
	  }
	// sort newiou based on cost
	sort(newiou.begin(), newiou.end(), sortuzon);	
	// update iou
	iou.push_back( newiou[0] );
	// update ioubounds
	if(i==0)
	  {
	    ioubounds = newiou[0].bounds;
	  }
	else
	  {
	    ioubounds = meet( ioubounds, newiou[0].bounds );
	  }
	// check whether ioubounds have reached saturation
	satbounds = true;
	for(int j=0; j<N; ++j)
	  {
	    satbounds = satbounds and (ioubounds(j).upper() <= bounds(j).upper()) and (ioubounds(j).lower() >= bounds(j).lower() );
	  }
	// update newiou (not iou) by removing first element
	newiou.erase(newiou.begin());
      }
    
//    // refine each zonotope in iou
//     int divs = pow(2,LogDivs);
//     int numintrs = iou.size();
// #pragma omp parallel for collapse(2)
//     for(int i=0; i<numintrs; ++i)
//       {
// 	for(int j=0; j<divs; ++j)
// 	  {
// 	    iou[i].sets[j].refine( bounds, 100 );
// 	  }
//       }    
  }

    
  
  //----------------------------------------------------------------------
  // flowpipe computation for a time period
  //----------------------------------------------------------------------
  void flow(double T, double MaxFlowTime){
    // set iou
    //SetIou();
    
    double FlowTime = 0;
    bool do_iter = true;
    
    while(do_iter){
      
       // append new iou
      SetIou();
      
      // compute next iou
      nextIou();

      // refine iou
      refine();
      cout << "iou size is " << iou.size() << "\n";     
      
      // update clocks
      SimTime += TimeStep;
      cout << SimTime.upper() << "\n";
      FlowTime += TimeStep;

      // update flowpipe and bounds
      ++FlowItr;
      (*(FlowItr)).bounds = pbounds;
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
    (*(FlowItr)).bounds = pbounds;

    /* perform recursive flowpipe computation with intermediate
       iou resets */
    while(SimTime.upper()<T){
      // compute flowpipe until IOU reset
      flow(T,IouResetTime);
    }
  }

  void simulate(double T){
    //omp_set_num_threads(100);
    //omp_set_nested(3);
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
