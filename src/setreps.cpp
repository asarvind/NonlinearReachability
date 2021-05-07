//****************************************************************************************************
// Zonotope class
//****************************************************************************************************
class ZonNd{
public:
  int dim, order;
  IvVectorNd bounds;
  
  vector<IvMatrixNNd> GenMat;
  IvVectorNd center, CoeffVect;

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
      CoeffVect(i) = Interval(-1.0,1.0);
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
  void MinSum( const IvVectorNd &x ){
    IvVectorNd c = middle(x);
    // update center
    center += c;


    IvVectorNd U = (GenMat[order-1] - GenMat[order -2])*CoeffVect;
    IvMatrixNNd M;
    M *= 0;
    for(int i=0; i<dim; ++i){
      M(i,i) = width( x(i) )*Interval(1,1) + U(i);
    }
    
    GenMat.insert( GenMat.begin(), M );
    IvMatrixNNd K = ( GenMat[order-1] )*2;
            
    GenMat.insert( GenMat.begin(), K );
    GenMat.resize( order );
  }
  
  void newMinSum(const IvVectorNd &x){
    IvVectorNd c = middle(x);
    IvVectorNd y = x-c;
    // update center
    center += c;    
    IvVectorNd b = *(GenMat.rbegin())*CoeffVect;
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
  void setBounds(){
    bounds = getBounds();
  }

  // compute bounds on the zonotope
  IvVectorNd getBounds(){
    IvVectorNd out = center*1;
    for(vector<IvMatrixNNd>::iterator itr=GenMat.begin(); itr!=GenMat.end(); ++itr){
      out += (*itr)*CoeffVect;
    }
    return out;
  }

  // close class
};
