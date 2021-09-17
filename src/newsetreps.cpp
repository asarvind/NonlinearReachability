void removeRow(Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}



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
  // reduce order of zonotope
  //----------------------------------------------------------------------
  
  void reduceOrder( int order ){
    if( GenMat.cols() > dim*order ){

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

  void newReduceOrder(int order)
  {
    if( GenMat.cols() > dim*order )
      {
	Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic> G(dim, GenMat.cols());
	G = GenMat;
	// compute normalized matrix
	Eigen::MatrixXd M(G.rows(), G.cols());
	for(int i=0; i<G.rows(); ++i)
	  {
	    double r = (G.block(i,0,1,G.cols()).lpNorm<1>()).upper() + 1e-6;
	    for(int j=0; j<G.cols(); ++j)
	      {
		M(i,j) = G(i,j).upper()/r;
	      }
	  }
	
	// rearrange columns
	vector< pair<double, int> > normvals;
	normvals.resize(G.cols());
	for(int i=0; i<G.cols(); ++i)
	  {
	    VectorNd thiscol = M.block(0,i,dim,1);
	    normvals[i] = make_pair( thiscol.lpNorm<1>()-thiscol.lpNorm<Eigen::Infinity>(), i );
	  }
	sort(normvals.begin(), normvals.end(), greater< pair<double,int> >());
	for(int i=0; i<G.cols(); ++i)
	  {
	    GenMat.block(0,i,dim,1) = G.block(0,normvals[i].second,dim,1);
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
	
	// compute rest of generator matrix
	Eigen::Matrix< Interval, Eigen::Dynamic, Eigen::Dynamic > newGenMat;
	if(order >1){	
	  newGenMat.resize( dim, dim*order-dim );
	  newGenMat << GenMat.block(0, 0, dim, dim*order - dim );
	}
	      
	// reduce order
	GenMat.resize( dim, order*dim );
	if( order > 1 ){
	  GenMat << errmat, newGenMat;
	}
	else{
	  GenMat << errmat;
	}	
      }
  }

  void pruneRows(int num)
  {
    IvVectorNd addvect;
    addvect *= 0;
    int zonorder = GenMat.cols();
    for(int i=0; i<dim; ++i)
      {
	// compute cost pairs
	vector<pair<double,int>> vals;
	vals.resize(zonorder);
	for(int j=0; j<zonorder; ++j)
	  {
	    vals[j].first = std::abs(GenMat(i,j).upper());
	    vals[j].second = j;
	  }
	// rearrange cost pairs
	sort(vals.begin(), vals.end());
	// prune and assign remainder
	for(int j=0; j<zonorder-num*dim; ++j)
	  {
	    addvect(i) += GenMat(i,vals[j].second)*Interval(-1,1);
	    GenMat(i,vals[j].second) = 0;
	  }
      }
    MinSum(addvect);
  }
  
  //----------------------------------------------------------------------
  // Left-multiply with a matrix
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
    
    Eigen::Matrix<Interval, Eigen::Dynamic, 1> coeffVect;
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
  void refine( IvVectorNd ioubounds, double rho)
  {
    // refine bounds
    IvVectorNd zonbounds = getBounds();
    IvVectorNd addvect;
    addvect *= 0;
    int ncols = GenMat.cols();
    for( int i=0; i<dim; ++i ){
      if( width( zonbounds(i) ) > ( 1+rho )*width( ioubounds(i) ) ){
	GenMat.block(i,0,1,ncols) *= 0;
	center(i) *= 0;
	addvect(i) = ioubounds(i);
      }
    }
    
    MinSum( addvect );
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

//******************************************************************************************>>>>
//union zonotope along with validity vector
//******************************************************************************************
struct uzonNd
{
  //----------------------------------------------------------------------
  // members
  //----------------------------------------------------------------------

  std::vector<ZonNd> sets;

  IvVectorNd bounds;

  double cost;

  uzonNd()
  {
    cost = 0;
  }
  
};

//------------------------------------------------------------------------------------------<<<<





