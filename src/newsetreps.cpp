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
  
  void reduceOrder( int order ){
    if( GenMat.cols() > dim*order ){
      Eigen::Matrix< Interval, Eigen::Dynamic, Eigen::Dynamic > newGenMat = GenMat;
      int prevcols = newGenMat.cols();

      // assign vector of magnitudes and indices for columns
      vector< pair< double, int > > vals;
      vals.resize( prevcols );
      for( int i=0; i<prevcols; ++i ){
	vals[i].first = ( newGenMat.block( 0,i,dim,1 ).norm() ).upper();
	vals[i].second = i;
      }

      // sort the vals vector in descending order
      sort( vals.begin(), vals.end() );
      reverse( vals.begin(), vals.end() );

      // reorder zonotope
      for( int i=0; i<dim; ++i ){
	for( int j=0; j<prevcols; ++j ){
	  GenMat(i,j) = newGenMat(i, vals[j].second);
	}
      }

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

      Eigen::Matrix< Interval, Eigen::Dynamic, Eigen::Dynamic > leftGenMat;
      if(order >1){	
	leftGenMat.resize( dim, dim*order-dim );
	leftGenMat << GenMat.block(0, 0, dim, dim*order - dim );
      }
      GenMat.resize( dim, order*dim );
      
      // reduce order
      if( order > 1 ){
	GenMat << leftGenMat, errmat;
      }
      else{
	GenMat << errmat;
      }
    }
  }
  
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
    GenMat << newGenMat, errmat;
  }

  //----------------------------------------------------------------------
  // set bounds
  //----------------------------------------------------------------------
  void setBounds(){
    bounds = getBounds();
  }

  // compute bounds on the zonotope
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

  // close class
};
