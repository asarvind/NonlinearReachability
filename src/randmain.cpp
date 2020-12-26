# include "unionreach.cpp"

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
  
  double T; // simulation time horizon length
  cout << "enter length of time horizon\n";
  cin >> T;
  cout << "\n";

  int logDivs;
  cout << "enter log_2( no. of divisions )\n";
  cin >> logDivs; // log_2( no. divisions )
  cout << "\n";

  double tStep; // time step size
  cout << "enter time step size\n";
  cin >> tStep;
  cout << "\n";

  // get seed number for pseudorandomnumber generation
  int seedNum;
  cout << "enter pseudorandom generator initialization seed (integer):\n";
  cin >> seedNum;
  
  // bloating
  ioureach bloatobj( initstate, inpbounds, tStep/10, logDivs );
  bloatobj.simulate( tStep );

  // create reachset object
  ioureach reachobj( initstate, inpbounds, tStep, logDivs );
  reachobj.doBloat = false;
  reachobj.bounds = bloatobj.MaxBounds;
  reachobj.MaxBounds = bloatobj.MaxBounds;
  reachobj.SimTime = Interval( 0, tStep );

  // set seed
  reachobj.seedNum = seedNum;
  
  // simulate
  reachobj.RandSim( T );

  // save traces
  reachobj.SaveTraces();
}
