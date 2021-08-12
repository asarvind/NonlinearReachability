//****************************************************************************************************
// Zonotope class
//****************************************************************************************************
class ZonNd{
public:
  int dim;
  IvVectorNd bounds;
  
  Eigen::Matrix< Interval, Eigen::Dynamic, Eigen::Dynamic > GenMat;
  IvVectorNd center;

  //----------------------------------------------------------------------
  // constructor
  //----------------------------------------------------------------------
  ZonNd(){
    center *= 0;
    dim = StateDim;
    GenMat.resize( dim, dim );
    for( int i = 0; i<dim; ++i){
      for( int j=0; j<dim; ++j ){
	GenMat(i,j) = 0;
      }
    }  
  }

  //======================================================================
  // Methods
  //======================================================================

  //----------------------------------------------------------------------
  // reduce order of zonotope below a threshold
  //----------------------------------------------------------------------

  double measuretoreduce( IvVectorNd x )
  {
    IvVectorNd z;
    for( int i=0; i<StateDim; ++i){
      z(i) = x(i)/( bounds(i)+1e-12 );       
    }
    IvVectorNd y = z/(z.norm() + 1e-12);
    for( int i=0; i<StateDim; ++i){
      y(i) = 1-y(i);       
    }    
    return y.norm().upper()*z.norm().upper();
  }
  
  void reduceOrder( int order ){
    if( GenMat.cols() > dim*order ){
      // Eigen::Matrix< Interval, Eigen::Dynamic, Eigen::Dynamic > newGenMat = GenMat;
      // int prevcols = newGenMat.cols();

      // // assign vector of magnitudes and indices for columns
      // vector< pair< double, int > > vals;
      // vals.resize( prevcols );
      // for( int i=0; i<prevcols; ++i ){
      // 	IvVectorNd colvect = newGenMat.block( 0,i,dim,1 ) ;
      // 	vals[i].first = measuretoreduce( colvect );
      // 	vals[i].second = i;
      // }

      // // sort the vals vector in descending order
      // //sort( vals.begin(), vals.end() );
      // //reverse( vals.begin(), vals.end() );

      // // reorder zonotope
      // for( int i=0; i<dim; ++i ){
      // 	for( int j=0; j<prevcols; ++j ){
      // 	  GenMat(i,j) = newGenMat(i, vals[j].second);
      // 	}
      // }

      // compute box approximation of last genmat
      Eigen::Matrix<Interval, Eigen::Dynamic, 1 > coeffVect;
      int cls = GenMat.cols() - dim*order + dim;
      coeffVect.resize( cls, 1 );
      for( int i=0; i<cls; ++i ){
	coeffVect(i,0) = Interval(-1,1);
      }
      IvVectorNd b = GenMat.block( 0, dim*order-dim, dim, cls )*coeffVect;

      IvMatrixNNd errmat;
      for( int i=0; i<dim; ++i ){
	for( int j=0; j<dim; ++j ){
	  if( i != j){
	    errmat(i,j) = Interval(0, 0);
	  }
	  else{
	    errmat(i,j) = ( width( b(i) )/2 )*Interval(1,1);
	  }
	}
      }

      Eigen::Matrix< Interval, Eigen::Dynamic, Eigen::Dynamic > newGenMat;
      if(order >1){	
	newGenMat.resize( dim, dim*order-dim );
	newGenMat << GenMat.block(0, 0, dim, dim*order - dim );
      }
      GenMat.resize( dim, order*dim );
      
      // reduce order
      if( order > 1 ){
	GenMat << errmat, newGenMat;
      }
      else{
	GenMat << errmat;
      }
    }
  }

  // void convertToParallelotope()
  // {
  //   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> newGenMat;
  //   newGenMat.resize(StateDim, GenMat.cols());
  //   for(int i=0; i<dim; ++i)
  //     {
  // 	for(int j=0; j<GenMat.cols(); ++j)
  // 	  {
  // 	    newGenMat(i,j) = GenMat(i,j).upper();
  // 	  }
  //     }
  //   MatrixNNd sqmat = newGenMat*( newGenMat.transpose() );
  //   Eigen::SelfAdjointEigenSolver<MatrixNNd> decop;
  //   decop.compute(sqmat);
  //   MatrixNNd eigv = decop.eigenvectors();
  //   Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic> newCoeffVect;
  //   newCoeffVect.resize(GenMat.cols(), 1);
  //   IvMatrixNNd ptopeGenMat;
  //   for(int i=0; i<GenMat.cols(); ++i)
  //     {
  // 	newCoeffVect(i) = Interval(-1,1);
  //     }
  //   for(int i=0; i<dim; ++i){
  //     for(int j=0; j<dim; ++j)
  // 	{
  // 	  ptopeGenMat(i,j) = eigv(i,j);
  // 	}
  //   }
  //   IvVectorNd scalevect = ( ( ptopeGenMat.inverse() )*GenMat )*newCoeffVect;
  //   for(int i=0; i<dim; ++i)
  //     {
  // 	ptopeGenMat.block(0,i,dim,1) *= scalevect(i).upper();
  //     }
  //   GenMat = ptopeGenMat;
  // }
  
  //----------------------------------------------------------------------
  // Pre-multiply with a matrix
  //----------------------------------------------------------------------
  void prod(const IvMatrixNNd &M){
    center = M*center;
    GenMat = M*GenMat;
  }

  //----------------------------------------------------------------------
  // Minkowski sum with an interval vector
  //----------------------------------------------------------------------  
  void MinSum(const IvVectorNd &x){
    IvVectorNd c = middle(x);
    // update center
    center += c;

    // convert span to box representation
    Eigen::Matrix< Interval, Eigen::Dynamic, Eigen::Dynamic > errmat;
    errmat.resize( dim, dim );

    for( int i=0; i<dim; ++i ){
      for( int j=0; j<dim; ++j ){
	if( i != j){
	  errmat(i,j) = Interval(0,0);
	}
	else{
	  errmat(i,j) = ( width( x(i) )/2 )*Interval(1,1);
	}
      }
    }
    
    // concatenate matrix at the end
    Eigen::Matrix< Interval, Eigen::Dynamic, Eigen::Dynamic > newGenMat = GenMat;
    GenMat.resize( dim, GenMat.cols() + dim );
    GenMat << errmat, newGenMat;
  }

  //----------------------------------------------------------------------
  // set bounds
  //----------------------------------------------------------------------
  void setBounds(){
    bounds = getBounds();
  }

  // compute interval bounds on the zonotope
  IvVectorNd getBounds(){
    IvVectorNd out = center*1;
    
    Eigen::Matrix<Interval, Eigen::Dynamic, 1 > coeffVect;
    coeffVect.resize( GenMat.cols(), 1 );
    for( int i=0; i<GenMat.cols(); ++i ){
      coeffVect(i) = Interval(-1,1);
    }
    out += GenMat*coeffVect;
    return out;
  }

  // compute directional bounds on the zonotope
  IvVectorNd getBounds( IvMatrixNNd ptemp ){
    IvVectorNd out = ptemp*center;
    
    Eigen::Matrix<Interval, Eigen::Dynamic, 1 > coeffVect;
    coeffVect.resize( GenMat.cols(), 1 );
    for( int i=0; i<GenMat.cols(); ++i ){
      coeffVect(i) = Interval(-1,1);
    }
    out += (ptemp*GenMat)*coeffVect;
    return out;
  }

  //----------------------------------------------------------------------
  // refine zonotope
  //----------------------------------------------------------------------
  void refine( IvVectorNd u, IvVectorNd thresbounds, IvMatrixNNd ptemp, IvMatrixNNd invptemp ){
    prod( ptemp );

    IvVectorNd addvect;
    addvect *= 0;
    int ncols = GenMat.cols();
    for( int i=0; i<dim; ++i ){
      if( ( u(i).upper() > thresbounds(i).upper() ) || ( u(i).lower() < thresbounds(i).lower() ) ){
	GenMat.block(i,0,1,ncols) *= 0;
	center(i) *= 0;
	addvect(i) = thresbounds(i);
      }
    }
    
    MinSum( addvect );
    prod( invptemp );
  }

  // close class
};

//------------------------------------------------------------------------------------------<<<<

//------------------------------------------------------------------------------------------>>>>
// Minkowski sum of two zonotopes
//------------------------------------------------------------------------------------------

ZonNd MinSum( ZonNd X, ZonNd Y )
{
  ZonNd out;
  out.center = X.center + Y.center;
  out.GenMat.resize( X.dim, X.GenMat.cols() + Y.GenMat.cols() );
  out.GenMat << Y.GenMat, X.GenMat;
  return out;
}

//------------------------------------------------------------------------------------------<<<<


