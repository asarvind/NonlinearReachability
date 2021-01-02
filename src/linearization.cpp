# include "settings.cpp"

class nonlinear{
public:
  /* Dimensions of state and input */
  int N,M;
  
  /* Time step for discrete time approximation */
  double TimeStep;
  
  /* Input to dynamical system, its center */
  IvVectorMd Inp, InpCenter;
  
  /* Declare identity matrix */
  // this matrix is initialized in the contructor
  IvMatrixNNd eyeN;

  // Array of eigenvectors
  VectorNd EvRe[StateDim];
  VectorNd EvIm[StateDim];

  //----------------------------------------------------------------------
  // constructor
  //----------------------------------------------------------------------
  nonlinear(const IvVectorMd &Input){
    N = StateDim;
    M = InputDim;
    TimeStep = StepSize;
    // set identity matrix
    for(int i=0; i<StateDim; i++){
      for(int j=0; j<StateDim; j++){
	if(i==j){
	  eyeN(i,j) = Interval(1.0,1.0);
	}
	else{
	  eyeN(i,j) = Interval(0.0,0.0);
	}
      }
    }    
    Inp = Input*1;
    for(int i=0; i<M; i++){
      double c = (Inp(i).lower()+Inp(i).upper())/2.0;
      InpCenter(i) = Interval(c,c);
    }
    
    // set eigenvectors
    for(int i=0; i<N; ++i){
      for(int j=0; j<N; ++j){
	EvRe[i](j) = ReEV[j][i];
	EvIm[i](j) = ImEV[j][i];
      }
    }
  }
  
  // method to compute bounds on vector field
  IvVectorNd Field(const IvVectorNd &State){
    return VectorField(State, Inp);
  }
  // method to compute state action matrix in continuous time
  IvMatrixNNd ContStMat(const IvVectorNd &center){
    return StateMat(center, InpCenter);
  }
  // method to compute input action matrix in continuous time
  IvMatrixNMd ContInpMat(const IvVectorNd &center){
    return InputMat(center, InpCenter);
  }
  // method to compute continuous time taylor error (quadratic)
  IvVectorNd TaylorErr(const IvVectorNd &State, const IvVectorNd &center){
    IvVectorNd StError = State - center;
    return ContError(State,Inp,StError,InpCenter);
  }
  Interval TaylorErr(const IvVectorNd &State, const IvVectorNd &StError, int dim){
    return ContDimError(State,Inp,StError,InpCenter,dim);
  }
  Interval TaylorErr(const IvVectorNd &State, const IvVectorNd &StError, const VectorNd &ReDir, const VectorNd &ImDir){
    Interval out = Interval(0,0);
    IvVectorNd ErrVect = ContError(State,Inp,StError,InpCenter);
    for(int i=0; i<StateDim; ++i){
      out += pow(ErrVect(i)*ReDir(i), 2) + pow(ErrVect(i)*ImDir(i), 2);
    }
    out = sqrt(out);
    out = hull(-1*out,out);
    return out;
  }

  // method compute field at a point
  IvVectorNd PointField(const IvVectorNd &center){
    IvVectorNd out =  VectorField(center, InpCenter);
    for( int i=0; i<N; ++i ){
      if ( isnan( out(i).upper() ) || isnan( out(i).lower() ) ){
	cout<< out(i).lower() << " " << out(i).upper() << " field error\n";
	  exit(0);
	}      
    }
    return out;
  }
  
  // Structure that contains objects obtained from linearization
  struct LinVals{
    IvVectorNd region; // valid region of linearization
    IvVectorNd state; // bounds on set of states
    IvVectorNd center; // center at which to linearize
    IvVectorNd shift; // shift in continuous linearization error
    IvMatrixNNd StMatCont; // state action matrix
    IvMatrixNMd InpMatCont; // Input action matrix
    IvVectorNd ErrTaylor; // taylor expansion error
    IvVectorNd ErrCont; // error
    IvMatrixNNd StMatDis; // discrete time state action matrix
    IvMatrixNMd InpMatDis; // discrete time input action matrix
    IvVectorNd ErrDis; // discrete time error
    IvMatrixNNd InitStMat; // initial discrete state action matrix
    IvMatrixNMd InitInpMat; // initial discrete input action matrix
    IvVectorNd MatShiftErr; // error due to shift in state action matrix
  };

  
  // method to perform intantaneous (not regional) continuous time linearization
  void ContLin( LinVals &L, bool recompute ){
    if (recompute){
      L.center = middle(L.state);
      L.region = L.state;
      L.shift = PointField(L.center);
      L.StMatCont = ContStMat(L.center);
      L.InpMatCont = ContInpMat(L.center);
    }
    L.ErrTaylor = TaylorErr(L.region,L.center);
    L.ErrCont = L.shift + L.ErrTaylor;
  }

  // method to perform discrete time linearization
  void DisLin(LinVals &L, bool recompute){    
    // time interval
    Interval delta = Interval(0,TimeStep);
    /* perform continuous time linearization */
    LinRegion(L);
    IvMatrixNNd A = L.StMatCont;
    IvMatrixNNd SqA = A*A;
    // discrete time action matrices
    L.StMatDis = eyeN + A*TimeStep + A*A*pow(TimeStep,2)/2;
    L.InpMatDis = L.InpMatCont*TimeStep + A*L.InpMatCont*pow(TimeStep,2)/2;
    /* discrete time linearization error*/
    // square of delta
    Interval epsilon = pow(delta,3);
    // Discrete time error
    L.ErrDis =
      L.ErrCont*TimeStep + ( A*L.ErrCont*pow(TimeStep,2)/2 + SqA*L.ErrCont*epsilon/2 +
			     SqA*A*(L.region-L.center)*epsilon/6 +
			     SqA*L.InpMatCont*(Inp-InpCenter)*epsilon/2 );
    // add shift due to non-centered action
    L.ErrDis += L.center - L.StMatDis*L.center;
  }

  // method to compute valid region of linearization
  void LinRegion(LinVals &L){    
    // time interval
    Interval delta = Interval(0,TimeStep);
    /* perform continuous time linearization */
    ContLin(L,true);
    IvMatrixNNd A = L.StMatCont;
    IvMatrixNNd SqA = A*A;
    bool valid = false;
    for(int i=0; i<20; ++i){
      // discrete time action matrices
      L.StMatDis = eyeN + A*delta + A*A*pow(delta,2)/2;
      L.InpMatDis = L.InpMatCont*delta + A*L.InpMatCont*pow(delta,2)/2;
      /* discrete time linearization error*/
      // square of delta
      Interval epsilon = pow(delta,3);
      // intial state and input action matrices
      // Discrete time error
      L.ErrDis =
	L.ErrCont*delta + ( A*L.ErrCont*pow(delta,2)/2 + SqA*L.ErrCont*epsilon/2 +
			       SqA*A*(L.region-L.center)*epsilon/6 +
			       SqA*L.InpMatCont*(Inp-InpCenter)*epsilon/2 );
      // add shift due to non-centered action
      L.ErrDis += L.center - L.StMatDis*L.center;
      // update region
      IvVectorNd NextRegion = L.StMatDis*L.state + L.InpMatDis*(Inp-InpCenter) + L.ErrDis;
      if(is_subset(NextRegion,L.region)){
	valid = true;
	break;
      }
      else{
	L.region = join(L.region,NextRegion);
	ContLin(L,false);
	valid = false;
      }
    }
    if(not valid){
      cout << "could not find valid linearization region, try reducing time step\n";
      exit(0);
    }

    // check nan error
    for( int i=0; i<N; ++i ){
      Interval out = L.region(i);
      if ( isnan( out.upper() ) || isnan( out.lower() ) ){
	cout<< out.lower() << " " << out.upper() << " linearization error\n";
	  exit(0);
      }
    }
  }
  
  //----------------------------------------------------------------------
  // method to compute optimal division vector of state region
  //----------------------------------------------------------------------
  struct OptErr{
    double err;
    IntVectorNd divs;
  };
  OptErr OptDivision(const IvVectorNd &state, int dim, int K){
    OptErr out;
    out.err = numeric_limits<double>::infinity();
    double error;
    for(int i=0; i<N; i++){
      out.divs[i] = 1;
    }
    // optimum index
    int optind = 0; 
    // error variables
    IvVectorNd OptStateErr = state - middle(state);
    IvVectorNd StError, StartErr;
    // loop to compute optimal error
    for(int i=0; i<K; i++){
      optind = 0;
      // loop to compute error
      StartErr = OptStateErr*1;
      for(int j=0; j<N; j++){
	StError = StartErr*1;
  	StError[j] /= 2.0;
	error = TaylorErr(state,StError,dim).upper();
	if(error<out.err){
	  out.err = error;
	  optind = j;
	  OptStateErr = StError*1;
	}
      }
      if(out.err>0){
	out.divs[optind] *= 2;
      }
      else{
	break;
      }
    }
    return out;
  }

  //----------------------------------------------------------------------
  // optimal division vectors along a given dimension
  //----------------------------------------------------------------------
  OptErr OptDivision(const IvVectorNd &state, const VectorNd &ReDir, const VectorNd &ImDir, int K){
    OptErr out;
    out.err = numeric_limits<double>::infinity();
    double error;
    for(int i=0; i<N; i++){
      out.divs[i] = 1;
    }
    // optimum index
    int optind = 0; 
    // error variables
    IvVectorNd OptStateErr = state - middle(state);
    IvVectorNd StError, StartErr;
    // loop to compute optimal error
    for(int i=0; i<K; i++){
      optind = 0;
      // loop to compute error
      StartErr = OptStateErr*1;
      for(int j=0; j<N; j++){
	StError = StartErr*1;
  	StError[j] /= 2.0;
	error = TaylorErr(state,StError,ReDir,ImDir).upper();
	if(error<out.err){
	  out.err = error;
	  optind = j;
	  OptStateErr = StError*1;
	}
      }
      if(out.err>0){
	out.divs[optind] *= 2;
      }
      else{
	break;
      }
    }
    return out;
  }

  //----------------------------------------------------------------------
  // optimal division vectors for a measure along all directions
  //----------------------------------------------------------------------
  OptErr OptDivision(const IvVectorNd &state, int K){    
    OptErr out;
    out.err = numeric_limits<double>::infinity();
    double error;
    for(int i=0; i<N; i++){
      out.divs[i] = 1;
    }
    // optimum index
    int optind = 0; 
    // error variables
    IvVectorNd OptStateErr = state - middle(state);

    // compute error without splitting
    VectorNd baseErr = radius( ContError(state,Inp,OptStateErr,InpCenter) );

    IvVectorNd StError, StartErr;
    // loop to compute optimal error
    for(int i=0; i<K; i++){
      optind = 0;
      // loop to compute error
      StartErr = OptStateErr*1;
      for(int j=0; j<N; j++){
	StError = StartErr*1;
  	StError[j] /= 2.0;
	VectorNd splitErr = radius( ContError(state,Inp,StError,InpCenter) );

	// loop to compute infinity norm of relative error
	error = 0;
	for( int k = 0; k<N; ++k ){
	  error = max( splitErr(k)/max( baseErr(k), 1e-10 ), error );
	}
	if(error<out.err){
	  out.err = error;
	  optind = j;
	  OptStateErr = StError*1;
	}
      }
      if(out.err>0){
	out.divs[optind] *= 2;
      }
      else{
	break;
      }
    }
    return out;
  }
  
  
  // close class hybridize
  };


