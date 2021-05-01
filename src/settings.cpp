# include <iostream>
# include <boost/numeric/interval.hpp>
# include "eigen/Eigen/Dense"
# include "vector"
# include <omp.h>
# include <fstream>
# include <algorithm> // for random shuffle
# include <bits/stdc++.h>
# include <ctime>

using namespace std;
using namespace boost::numeric;


/* Specialize is_convertible in Eigen::internal to allow matrix multiplication
   on interval matrices. */
template<typename X, typename S, typename P>
struct is_convertible<X,interval<S,P> > {
  enum { value = is_convertible<X,S>::value };
};
template<typename S, typename P1, typename P2>
struct is_convertible<interval<S,P1>,interval<S,P2> > {
  enum { value = true };
};

// Define interval datatype with double precision
typedef boost::numeric::interval<
    double,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_transc_std<double> >,
      boost::numeric::interval_lib::checking_base<double> > >
    Interval;

// Overload functions to add interval with integer
Interval operator+(const Interval &x, const int &y){
  return x+double(y);
}
Interval operator+(const int &y,const Interval &x){
  return x+double(y);
}

// Overload functions to subtract interval with integer
Interval operator-(const Interval &x, const int &y){
  return x-double(y);
}
Interval operator-(const int &y,const Interval &x){
  return x-double(y);
}


// Overload functions to multiply interval with integer
Interval operator*(const Interval &x, const int &y){
  return x*double(y);
}
Interval operator*(const int &y, const Interval &x){
  return x*double(y);
}

// Overload functions to divide interval with integer and vice-versa
Interval operator/(const Interval &x, const int &y){
  return x/double(y);
}
Interval operator/(const int &y,const Interval &x){
  return double(y)/x;
}
// Define functions to compute inverse trigonometric and trigonometric
// mycos
Interval mycos(const Interval &x){
  Interval y = cos(x + Interval(-1,1)*1e-10*x);
  if((isnan(y.lower())) || (isnan(y.upper()))){
    // cout << x.lower() << " " << x.upper() << " yes\n";
    Interval z(1.570,1.571);
    y = sin(z-x);
    if((isnan(y.lower())) || (isnan(y.upper()))){
      cout << x.lower() << " " << x.upper() << " cos error\n";
      exit(0);
    }
  }
  return y;
}

Interval mysin(const Interval &x){
  Interval y = sin(x);
  if((isnan(y.lower())) || (isnan(y.upper()))){
    // cout << x.lower() << " " << x.upper() << " yes\n";
    Interval z(1.570,1.571);
    y = cos(z-x);
    if((isnan(y.lower())) || (isnan(y.upper()))){
      cout << x.lower() << " " << x.upper() << " sin error\n";
      exit(0);
    }
  }
  return y;
}

Interval mytan(const Interval &x){
  Interval y = mysin(x)/mycos(x);
  return y;
}

// secant
Interval sec(const Interval &x){
  Interval y = 1/mycos(x);
  if((isnan(y.lower())) || (isnan(y.upper()))){
    cout << x.lower() << " " << x.upper() << " sec error\n";
    exit(0);
  }
  return y;
}

// cosecant
Interval csc(const Interval &x){
  Interval y =  1/mysin(x);
  if((isnan(y.lower())) || (isnan(y.upper()))){
    cout << x.lower() << " " << x.upper() << " csc error\n";
    exit(0);
  }
  return y;
}

// cotangent
Interval cot(const Interval &x){
  Interval y =  mycos(x)/mysin(x);
  if((isnan(y.lower())) || (isnan(y.upper()))){
    cout << x.lower() << " " << x.upper() << " cot error\n";
    exit(0);
  }
  return y;
}

// Functions to compute powers
Interval pow(const Interval &x, const double &y){
  Interval out =  exp(y*log(x));
  if((isnan(out.lower())) || (isnan(out.upper()))){
    cout << x.lower() << " " << x.upper() << " log error\n";
    exit(0);
  }
  return out;
}
Interval pow(const Interval &x, const Interval &y){
  Interval out = exp(y*log(x));
  if((isnan(out.lower())) || (isnan(out.upper()))){
    cout << x.lower() << " " << x.upper() << " log error\n";
    exit(0);
  }
  return out;
}




// Include dynamics

// Constants StateDim, InputDim and pardim are extracted from dynamics by running python scripts
#include "pywrite/Constants.cpp"

// Define types of vectors and Matrices
typedef Eigen::Matrix<Interval,StateDim,StateDim> IvMatrixNNd;
typedef Eigen::Matrix<Interval,StateDim,InputDim> IvMatrixNMd;
typedef Eigen::Matrix<Interval,StateDim,1> IvVectorNd;
typedef Eigen::Matrix<Interval,InputDim,1> IvVectorMd;
typedef Eigen::Matrix<int,StateDim,1> IntVectorNd;
typedef Eigen::Matrix<double,StateDim,1> VectorNd;
typedef Eigen::Matrix<Interval,pardim,1> IvVectorKd;
typedef Eigen::Matrix<Interval,tempRows,1> IvVectorLd; // interval vector of bounds on polytope type
typedef Eigen::Matrix<Interval,Eigen::Dynamic,StateDim> IvMatrixLNd; // template matrix type

// Function to compute vector field, continuous state and input matrices
// and continuous error from linearization.
#include "pywrite/VectorField.cpp"
#include "pywrite/StateMat.cpp"
#include "pywrite/InputMat.cpp"
#include "pywrite/ContError.cpp"
#include "pywrite/DimError.cpp"

//======================================================================
// Operations on interval state vectors
//======================================================================
  
// join of interval vectors
IvVectorNd join(IvVectorNd x, IvVectorNd y){
  IvVectorNd out;
  for(int i=0; i<StateDim; i++){
    out(i) = hull(x(i),y(i));
  }
  return out;
}

IvVectorLd join(IvVectorLd x, IvVectorLd y){
  IvVectorLd out;
  for(int i=0; i<tempRows; i++){
    out(i) = hull(x(i),y(i));
  }
  return out;
}


// meet of interval vectors
IvVectorNd meet(IvVectorNd x, IvVectorNd y){
  IvVectorNd out;
  for(int i=0; i<StateDim; i++){
    out(i) = intersect(x(i),y(i));
    if ( isnan( out(i).upper() ) || isnan( out(i).lower() ) ){
      cout << "blunder";
      cout<< x(i).lower() << " " << x(i).upper() << "\n";
      cout<< y(i).lower() << " " << y(i).upper() << "\n intersection error\n";
      exit(0);
    }
    if (not overlap(x(i),y(i)) ){
	cout << x(i).lower() << " " << x(i).upper() << "\n";
	cout << y(i).lower() << " " << y(i).upper() << "\n";
	cout << "meet error\n";
	exit(0);
      }
  }
  return out;
}

IvVectorLd meet(IvVectorLd x, IvVectorLd y){
  IvVectorLd out;
  for(int i=0; i<tempRows; i++){
    out(i) = intersect(x(i),y(i));
    if ( isnan( out(i).upper() ) || isnan( out(i).lower() ) ){
      cout << "blunder";
      cout<< x(i).lower() << " " << x(i).upper() << "\n";
      cout<< y(i).lower() << " " << y(i).upper() << "\n intersection error\n";
      exit(0);
    }
    if ( not overlap(x(i),y(i)) ){
	cout << x(i).lower() << " " << x(i).upper() << "\n";
	cout << y(i).lower() << " " << y(i).upper() << "\n";
	cout << "meet error\n";
	exit(0);
      }   
  }
  return out;
}

// method to compute center of a state interval vector
IvVectorNd middle(const IvVectorNd &x){
  IvVectorNd out;
  for(int i=0; i<StateDim; i++){
    double l = x( i ).lower();
    double r = x( i ).upper();
    out(i) = ( Interval( r, r ) + Interval( l, l ) )/2.0;
    if ( isnan( out(i).upper() ) || isnan( out(i).lower() ) ){
      cout<< out(i).lower() << " " << out(i).upper() << " center error\n";
      exit(0);
    }      
  }
  return out;
}


// radius of interval vector
VectorNd radius(const IvVectorNd &x){
  VectorNd out;
  for(int i=0; i<StateDim; ++i){
    out(i) = (x(i).upper()-x(i).lower())/2.0;
  }
  return out;
}

// check containment between interval vectors
bool is_subset(const IvVectorNd &x, const IvVectorNd &y){
  bool out = true;
  for(int i=0; i<StateDim; ++i){
    out = out && (subset(x(i),y(i)));
  }
  return out;
}

// check non-empty intersection between interval vectors
bool is_overlap(const IvVectorNd &x, const IvVectorNd &y){
  bool out = true;
  for(int i=0; i<StateDim; ++i){
    out = out && (overlap(x(i),y(i)));
  }
  return out;
}

//======================================================================
// Operations on double vectors
//======================================================================

// compare two double vectors
bool not_larger(const VectorNd &x, const VectorNd &y){
  bool out = true;
  for(int i=0; i<StateDim; ++i){
    out = out && x(i)<=y(i);
  }
  return out;
}


