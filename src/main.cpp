# include "reachability.cpp"

int main(){
  ifstream readlower, readupper, simhypar;
  double lb, ub;
  string line;
  
  // read state initial lower bounds
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
    
  // read input initial lower bounds
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

  // read simulation hyperparameters
  simhypar.open( "src/pywrite/simpars.txt" );
  
  double T; // simulation time horizon length
  simhypar >> T;

  int logDivs;
  simhypar >> logDivs; // log_e( no. divisions )

  double tStep; // time step size
  simhypar >> tStep;

  simhypar.close();

  // bloating
  ioureach bloatobj( initstate, inpbounds, tStep/10, logDivs );
  bloatobj.simulate( tStep );

  // create reachset object
  ioureach reachobj( initstate, inpbounds, tStep, logDivs );
  reachobj.doBloat = false;
  reachobj.bounds = bloatobj.MaxBounds;
  reachobj.MaxBounds = bloatobj.MaxBounds;
  reachobj.SimTime = Interval( 0, tStep );

  // simulate
  reachobj.simulate( T );

  // save traces
  reachobj.SaveTraces();
}
