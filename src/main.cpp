# include "reachability.cpp"

int main(){
  ifstream readlower, readupper, simhypar;
  double lb, ub, num, den;
  string line;
  
  // read state initial bounds
  readlower.open( "src/pywrite/stlb.txt" );
  readupper.open( "src/pywrite/stub.txt" );
  int ind = 0;
  IvVectorNd initstate;
  
  for( int i = 0; i<StateDim; ++i){
    readlower >> lb;
    readupper >> ub;
    initstate( i ) = Interval( lb, ub );
  }
  readlower.close();
  readupper.close();
    
  // read input bounds
  readlower.open( "src/pywrite/inplb.txt" );
  readupper.open( "src/pywrite/inpub.txt" );
  ind = 0;
  IvVectorMd inpbounds;
  for( int i = 0; i<InputDim; ++i ){
    readlower >> lb;
    readupper >> ub;
    inpbounds( i ) = Interval( lb, ub );
  }
  readlower.close();
  readupper.close();

  // read parameter values
  readlower.open( "src/pywrite/parden.txt" );
  readupper.open( "src/pywrite/parnum.txt" );
  ind = 0;
  IvVectorKd parvals;
  for( int i = 0; i<pardim; ++i ){
    readlower >> den;
    readupper >> num;
    parvals( i ) = Interval( num, num )/Interval(den, den);
  }
  readlower.close();
  readupper.close();

  // read simulation hyperparameters
  simhypar.open( "src/pywrite/simpars.txt" );
  
  double T; // simulation time horizon length
  simhypar >> T;

  int logDivs;
  simhypar >> logDivs; // log_e( no. divisions )

  double tStep; // time step size
  simhypar >> tStep;

  int zonOrder;
  simhypar >> zonOrder;  // order of zonotope

  simhypar.close();

  // bloating
  ioureach bloatobj( initstate, inpbounds, parvals, tStep/10, logDivs );
  bloatobj.zonOrder = zonOrder;
  bloatobj.simulate( tStep );

  // create reachset object
  ioureach reachobj( initstate, inpbounds, parvals, tStep, logDivs );
  reachobj.zonOrder = zonOrder;
  reachobj.doBloat = false;
  reachobj.bounds = bloatobj.MaxBounds;
  reachobj.MaxBounds = bloatobj.MaxBounds;
  reachobj.SimTime = Interval( 0, tStep );

  // simulate
  reachobj.simulate( T );

  // save traces
  reachobj.SaveTraces();
}
