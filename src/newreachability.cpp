# include "linearization.cpp"


//****************************************************************************************************
// Zonotope class
//****************************************************************************************************
class ZonNd{
public:
  int dim, order;
  IvVectorNd bounds, refinebounds;
  
  vector<IvMatrixNNd> GenMat;
  IvVectorNd center, coeffBounds;

  //----------------------------------------------------------------------
  // constructor
  //----------------------------------------------------------------------
  ZonNd(){
    dim = StateDim;
    order = 2;
    center *= 0;
    GenMat.resize( order );
    for( int i = 0; i<order; ++i){
      GenMat[i] *= 0;
    }
    // create [-1,1]^n vector
    for(int i=0; i<dim; ++i){
      coeffBounds(i) = Interval(-1.0,1.0);
    }
  }

  //======================================================================
  // Methods
  //======================================================================

  // reset zonotope with a specified order
  void resetOrder( int neworder ){
    order = neworder;
    center *= 0;
    bounds *= 0;
    refinebounds *= 0;
    GenMat.resize( order );
    for(int i=0; i<order; ++i){
      GenMat[i] *= 0;
    }
  }
  
  //----------------------------------------------------------------------
  // Pre-multiply with a matrix
  //----------------------------------------------------------------------
  void prod(const IvMatrixNNd &M){
    center = M*center;
    for(vector<IvMatrixNNd>::iterator itr=GenMat.begin(); itr!= GenMat.end(); ++itr){      
      *itr = M*(*itr);
    }
  }

  //----------------------------------------------------------------------
  // Minkowski sum with an interval vector
  //----------------------------------------------------------------------  
  void MinSum(const IvVectorNd &x){
    IvVectorNd c = middle(x);
    IvVectorNd y = x-c;
    // update center
    center += c;    
    IvVectorNd b = *(GenMat.rbegin())*coeffBounds;
    // shift generators left
    vector<IvMatrixNNd>::reverse_iterator itr = GenMat.rbegin();
    vector<IvMatrixNNd>::reverse_iterator next_itr = GenMat.rbegin();
    ++next_itr;
    vector<IvMatrixNNd>::reverse_iterator end_itr = GenMat.rend();
    --end_itr;
    while(itr!=end_itr){      
      *itr = *(next_itr);
      ++itr;
      ++next_itr;
    }
    // compute first generator matrix
      *(GenMat.begin()) *= 0;
    for(int i=0; i<dim; ++i){
      Interval r = b(i)+y(i);
      (*(GenMat.begin()))(i,i) = (r.upper()-r.lower())/2;
    }
  }

  //----------------------------------------------------------------------
  // set bounds
  //----------------------------------------------------------------------
  void setBounds( IvMatrixNNd A){
    bounds = getBounds();
    refinebounds = getBounds( A );
  }

  // compute projection bounds on the zonotope
  IvVectorNd getBounds(){
    IvVectorNd out = center*1;
    for(vector<IvMatrixNNd>::iterator itr=GenMat.begin(); itr!=GenMat.end(); ++itr){
      out += (*itr)*coeffBounds;
    }
    return out;
  }

  IvVectorNd getBounds( const IvMatrixNNd &A ){
    IvVectorNd out = A*center;
    for( int i=0; i<order; ++i ){
      out += ( A*GenMat[i] )*coeffBounds;
    }
    return out;
  }

  //----------------------------------------------------------------------
  // refine with a threshold
  //----------------------------------------------------------------------
  void refine( double thres, IvMatrixNNd A, IvMatrixNNd invA ){
    // change basis of zonotope
    prod( A );

    // get bounds in transformed basis
    IvVectorNd compvect = getBounds();
    
    IvVectorNd addvect; // initialize sum vector for refinement
    IvMatrixNNd resetmat; // initialize reset matrix for refinement

    // compute above matrices
    for( int i=0; i<dim; ++i ){
      if( ( compvect(i).upper() - refinebounds(i).upper() ) > thres*width( refinebounds(i) ) || ( refinebounds(i).lower() - compvect(i).lower() ) > thres*width( refinebounds(i) ) ){
	addvect(i) += refinebounds(i);
	resetmat(i,i) *= 0;
      }
      else{
	resetmat(i,i) += 1;
      }
    }
    prod( resetmat );
    MinSum( addvect );

    // reset basis to original
    prod( invA );
  }


  // close ZonNd class
};


//****************************************************************************************************
// Intersection of Unions reach set class
//****************************************************************************************************
class ioureach: public nonlinear{
public:
  //======================================================================
  // attributes
  //======================================================================
  
  // initial set
  IvVectorNd InitState;
  
  // current bounds on the reachable set
  IvVectorNd bounds;
  
  // maximum bounds on the reachable set in the simulated time horizon
  IvVectorNd MaxBounds;

  // compute directional bounds used for refinement
  IvVectorNd refinebounds;
  
  // time interval;
  Interval SimTime;
    
  // log of number of elements in union
  int LogDivs;
  
  // number of elements in union
  int intrs;

  // order of zonotope
  int zonOrder;

  // threshold for refinement
  double refinethreshold;
    
  // array of optimal divisions vectors
  IntVectorNd DivVecs[StateDim];
    
  // union storage capacity of iou array
  static const int MaxDivs = 128;
  
  /* iou array is 2-dimensional array of zonotopes.
     First dimension contains intersections. Second dimension
     contains unions. */
  ZonNd iou[StateDim][MaxDivs];

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

  // matrices used for refinement of zonotope
  IvMatrixNNd refinedirs, invrefinedirs;

  //----------------------------------------------------------------------
  // constructor
  //----------------------------------------------------------------------
  ioureach(const IvVectorNd &State, const IvVectorMd &Input, const IvVectorKd &parvals, double tstep, int k)
    :nonlinear(Input, parvals){
    InitState = State;
    TimeStep = tstep;
    zonOrder = (iou[0][0]).order;
    LogDivs = k;
    refinethreshold = 100;
    SimNum = 0;
    doBloat = true;

    // set refdirs
    ifstream eigBasisFile( "src/pywrite/eigRealBasis.txt" );
    double readval; 
    for( int i=0; i<N; ++i ){
      for( int j=0; j<N; ++j ){
	eigBasisFile >> readval;
	refinedirs(i,j) += readval;
      }
    }
    eigBasisFile.close();
    invrefinedirs = refinedirs.inverse();
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
    for(int i=0; i<N; ++i){
      OptErr E = OptDivision(bounds,EvRe[i],EvIm[i],LogDivs);
      DivVecs[intrs] = E.divs;
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
  
  //----------------------------------------------------------------------
  // Method: compute iou representation
  //----------------------------------------------------------------------
  void SetIou(){
    SetDivVecs();
    // divide state into iou
    int divs = pow(2,LogDivs);
    int q, d;
    q = 0;
    d = 0;
    Interval iv, delta;
    for(int i=0; i<intrs; i++){
      for(int j=0; j<divs; j++){
	// reset zonotope
	iou[i][j] = ZonNd();
	iou[i][j].resetOrder( zonOrder );
	q = j;
	IvVectorNd addvect;
	for(int k=0; k<N; k++){
	  d = DivVecs[i](k);
	  iv = bounds(k);
	  delta = (Interval(iv.upper(),iv.upper())-Interval(iv.lower(),iv.lower()))/d;
	  addvect(k) = hull(iv.lower() + (q%d)*delta,iv.lower() + ((q%d) + 1)*delta);
	  q /= d;
	}
	// set zonotope value
	iou[i][j].MinSum(addvect);
	iou[i][j].setBounds( refinedirs );
      }
    }
  }


  //----------------------------------------------------------------------
  // Compute bounds from iou
  //----------------------------------------------------------------------
  void SetBounds(){
    IvVectorNd U;
    int divs = pow(2,LogDivs);
    // initialize bounds as union of elements in first intersection
    U = iou[0][0].bounds;
    for(int i=1; i<divs; i++){
      U = join(U,iou[0][i].bounds);
    }
    bounds = U*1; // bounds initialized
    // recursively get U and intersect with bounds
    for(int i=1; i<intrs; i++){
      U = iou[i][0].bounds;
      for(int j=1; j<divs; j++){
	U = join(U, iou[i][j].bounds);
      }
      bounds = meet(bounds,U);
    }
    MaxBounds = join( MaxBounds, bounds );
  }

    //----------------------------------------------------------------------
  // Compute bounds from iou
  //----------------------------------------------------------------------
  void SetRefineBounds(){
    IvVectorNd U;
    int divs = pow(2,LogDivs);
    // initialize bounds as union of elements in first intersection
    U = iou[0][0].refinebounds;
    for(int i=1; i<divs; i++){
      U = join(U,iou[0][i].refinebounds);
    }
    refinebounds = U*1; // bounds initialized
    // recursively get U and intersect with bounds
    for(int i=1; i<intrs; i++){
      U = iou[i][0].refinebounds;
      for(int j=1; j<divs; j++){
	U = join(U, iou[i][j].refinebounds);
      }
      refinebounds = meet(refinebounds,U);
    }
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
    MaxBounds = bounds;
    SimTime = delta;
  }

  //----------------------------------------------------------------------
  // Next IOU set
  //----------------------------------------------------------------------
  void NextIou(){
    int divs = pow(2,LogDivs);
    #pragma omp parallel for collapse(2)
    for(int j=0; j<intrs; ++j){
      for(int k=0; k<divs; ++k){
	iou[j][k].bounds = meet( bounds, iou[j][k].bounds );
	iou[j][k].refinebounds = meet( refinebounds, iou[j][k].refinebounds );
	iou[j][k].refine( refinethreshold, refinedirs, invrefinedirs );
	LinVals L;
	L.state = iou[j][k].bounds;
	L.state = meet(L.state,bounds);
	DisLin(L,true);
	iou[j][k].prod(L.StMatDis);
	IvVectorNd addvect = L.InpMatDis*Inp + L.ErrDis;
	iou[j][k].MinSum( addvect );
	iou[j][k].setBounds( refinedirs );
	// IvVectorNd boxBounds = L.StMatDis*L.state + addvect;
	// iou[j][k].bounds = meet( boxBounds, iou[j][k].bounds );
      }
    }
    SetBounds();
    SetRefineBounds();
  }


  //----------------------------------------------------------------------
  // flowpipe computation for a time period
  //----------------------------------------------------------------------
  void flow(double T, double MaxFlowTime){
    // set iou
    SetIou();

    SetBounds();
    SetRefineBounds();
    
    double FlowTime = 0;
    bool do_iter = true;
    
    while(do_iter){
      // compute next iou
      NextIou();
      
      // update clocks
      SimTime += TimeStep;
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
