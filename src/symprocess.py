from sympy import *
from sympy.codegen.ast import Assignment
import time # for compute time calculation
import numpy as np
import scipy.linalg as scl
import math

import os

#****************************************************************************************************
# class to process dynamics and write c++ functions
#****************************************************************************************************

class writetocpp:
    def __init__(self,dynamics, tempOrder = 10):
        # specify all symbolic variables
        stVars = list( var( list( dynamics["field"].keys() ) ) )
        inpVars = list( var( list( dynamics["inp"].keys() ) ) )
        allVars = stVars + inpVars
        self.paramVars = var( list( dynamics["param"].keys() ) )
        self.stVars = stVars
        self.inpVars = inpVars
        self.allVars = allVars

        self.dirmat = dynamics["directions"]

        # define initial state bounds
        self.initState = dynamics[ "initState" ]
        
        # define symbolic vector field
        self.field = {}
        for ky in dynamics["field"]:
            self.field[ ky ] = nsimplify( eval( dynamics[ "field" ][ ky ] ) )

        # define input
        self.inp = dynamics["inp"]

        # define parameters
        self.param = dynamics["param"]

        # dimensions of state, input and parameters
        self.N = len( self.field )
        self.M = len( dynamics["inp"] )
        self.K = len( dynamics["param"] )
        self.tempOrder = tempOrder 
        self.L = 1 # will be reset later
        
        # write vector field as a matrix
        self.matField = list( self.field.values() )

        # declare Taylor error variables
        self.errVars = {}
        for myvar in allVars:
            self.errVars[ str(myvar) ] =  var( "e_"+str(myvar) )        # convert to vector
        self.errVarsMat = Matrix( self.N+self.M, 1, list( self.errVars.values() ) )

        # state action matrix
        self.stJacobian = Matrix(self.matField).jacobian( Matrix( self.N, 1, stVars ))

        # input action matrix
        self.inpJacobian = Matrix(self.matField).jacobian( Matrix( self.M, 1, inpVars ))

        # Taylor error
        self.linerr = {}
                                                          
        for ky in self.field:
            self.linerr[ ky ] = 0.5*self.errVarsMat.transpose()*hessian( self.field[ ky ], tuple( allVars ) )*self.errVarsMat

        # remove old pywrite and results folder
        os.system( "rm -rf src/pywrite" )
        os.system( "rm -rf results" )

        # create pywrite and results folder
        os.system( "mkdir -p src/pywrite" )
        os.system( "mkdir -p results" )
               
        # filenames for writing
        self.ConstantsFile = "src/pywrite/Constants.cpp"
        self.VectorFieldFile = "src/pywrite/VectorField.cpp"
        self.StateMatFile = "src/pywrite/StateMat.cpp"
        self.InputMatFile = "src/pywrite/InputMat.cpp"
        self.ContErrorFile = "src/pywrite/ContError.cpp"
        self.DimErrorFile = "src/pywrite/DimError.cpp"
        self.stLbFile = "src/pywrite/stlb.txt"
        self.stUbFile = "src/pywrite/stub.txt"
        self.inpLbFile = "src/pywrite/inplb.txt"
        self.inpUbFile = "src/pywrite/inpub.txt"

        # write all files
        self.WriteTemp()
        self.WriteConstants()
        self.WriteVectorField()
        self.WriteStateMat()
        self.WriteInputMat()
        self.WriteContError()
        self.WriteDimError()        

    def WriteConstants(self):
        f = open(self.ConstantsFile,"w")
        # write dimension of vector field
        f.write("const int StateDim = "+str(self.N)+";\n")
        # write dimension of input
        f.write("const int InputDim = "+str(self.M)+";\n")
        # write number of parameters
        f.write("const int pardim = "+str(self.K)+";\n")
        # write dimension of template
        f.write("const int tempRows = "+str(self.L)+";\n")

    def WriteVectorField(self):
        # clear previous content on file
        f = open(self.VectorFieldFile, "w")
        f.write("// Function to compute bounds on vector field\n")
        f.close()
        
        # Initialize Function
        f = open(self.VectorFieldFile,"a")
        
        # open function 
        f.write("IvVectorNd VectorField(const IvVectorNd &state, const IvVectorMd &inp, const IvVectorKd &parvals)"+"{\n")
        
        # define output of function
        f.write("// define output\n")
        f.write("IvVectorNd out;\n")
        
        # assign values to state symbols
        f.write("//assign values to state symbols\n")
        ind = 0;
        for ky in self.field:
            f.write("Interval "+ ky +"= "+"state("+str(ind)+");\n")
            ind += 1
            
        # assign values to input symbols
        f.write("// assign values to input symbols \n")
        ind = 0
        for ky in self.inp:
            f.write("Interval "+ ky +"= "+"inp("+str(ind)+");\n")
            ind += 1
            
        # assign values to parameters
        f.write("// assign values to input symbols \n")
        ind = 0
        for ky in self.param:
            f.write("Interval "+ ky +"= "+"parvals("+str(ind)+");\n")
            ind += 1
        
        # compute vector field bounds
        f.write("// compute bounds on vector field\n")
        ind = 0
        for ky in self.field:
            mystr = "out("+str(ind)+")"+"= "+cxxcode(self.field[ky]).replace("std::","")+";\n"
            for repl in ["cos","sin","tan"]:
                mystr = mystr.replace(repl,"my"+repl)
            f.write(mystr)
            ind += 1
        f.write("\n")
        # return output
        f.write("return out;\n")
        # close function
        f.write("}")
        f.close()

    def WriteStateMat(self):
        # clear previous content on file
        f = open(self.StateMatFile, "w")
        f.write("// Function to compute continuous state change matrix\n")
        f.close()
        # Initialize Function
        f = open(self.StateMatFile,"a")
        # open function 
        f.write("IvMatrixNNd StateMat(const IvVectorNd &state, const IvVectorMd &inp, const IvVectorKd &parvals)"+"{\n")
        # define output of function
        f.write("// define output\n")
        f.write("IvMatrixNNd out;\n")
        
        # assign values to state symbols
        f.write("//assign values to state symbols\n")
        ind = 0;
        for ky in self.field:
            f.write("Interval "+ ky +"= "+"state("+str(ind)+");\n")
            ind += 1
            
        f.write("// assign values to input symbols \n")
        ind = 0
        # assign values to input symbols
        for ky in self.inp:
            f.write("Interval "+ ky + "= "+"inp("+str(ind)+");\n")
            ind += 1

        # assign values to parameters
        f.write("// assign values to parameters \n")
        ind = 0
        for ky in self.param:
            f.write("Interval "+ ky +"= "+"parvals("+str(ind)+");\n")
            ind += 1

        # compute matrix
        f.write("// compute matrix\n")
        for i in range(self.N):
            for j in range(self.N):
                mystr = "out("+str(i)+","+str(j)+")"+"= "+cxxcode(self.stJacobian[i,j]).replace("std::","")+";\n"
                for repl in ["cos","sin","tan"]:
                    mystr = mystr.replace(repl,"my"+repl)
                f.write(mystr)
        f.write("\n")
        # return output
        f.write("return out;\n")
        # close function
        f.write("}")
        f.close()

    def WriteInputMat(self):
        # clear previous content on file
        f = open(self.InputMatFile, "w")
        f.write("// Function to compute continuous input action matrix\n")
        f.close()
        # Initialize Function
        f = open(self.InputMatFile,"a")
        # open function 
        f.write("IvMatrixNMd InputMat(const IvVectorNd &state, const IvVectorMd &inp, const IvVectorKd &parvals)"+"{\n")
        # define output of function
        f.write("// define output\n")
        f.write("IvMatrixNMd out;\n")
        
        # assign values to state symbols
        f.write("//assign values to state symbols\n")
        ind = 0;
        for ky in self.field:
            f.write("Interval "+ ky +"= "+"state("+str(ind)+");\n")
            ind += 1
        f.write("// assign values to input symbols \n")
        ind = 0
        
        # assign values to input symbols
        for ky in self.inp:
            f.write("Interval "+ ky +"= "+"inp("+str(ind)+");\n")
            ind += 1

        # assign values to parameters
        f.write("// assign values to parameters \n")
        ind = 0
        for ky in self.param:
            f.write("Interval "+ ky +"= "+"parvals("+str(ind)+");\n")
            ind += 1

        # compute matrix
        f.write("// compute matrix\n")
        for i in range(self.N):
            for j in range(self.M):
                mystr = "out("+str(i)+","+str(j)+")"+"= "+cxxcode(self.inpJacobian[i,j]).replace("std::","")+";\n"
                for repl in ["cos","sin","tan"]:
                    mystr = mystr.replace(repl,"my"+repl)
                f.write(mystr)
        f.write("\n")
        # return output
        f.write("return out;\n")
        # close function
        f.write("}")
        f.close()

    def WriteContError(self):
        # clear previous content on file
        f = open(self.ContErrorFile, "w")
        f.write("// Function to compute continuous time linearization error\n")
        f.close()
        # Initialize Function
        f = open(self.ContErrorFile,"a")
        # open function 
        f.write("IvVectorNd ContError(const IvVectorNd &state, const IvVectorMd &inp,"+
        "const IvVectorNd &StError, const IvVectorMd &InpCenter, const IvVectorKd &parvals)"+"{\n")
        # define output of function
        f.write("// define output\n")
        f.write("IvVectorNd out;\n")
        
        # assign values to state symbols
        f.write("//assign values to state symbols\n")
        ind = 0;
        for ky in self.field:
            f.write("Interval "+ ky +"= "+"state("+str(ind)+");\n")
            ind += 1
            
        # assign values to input symbols
        f.write("// assign values to input symbols \n")
        ind = 0
        for ky in self.inp:
            f.write("Interval "+ ky +"= "+"inp("+str(ind)+");\n")
            ind += 1
            
        # assign values to parameters
        f.write("// assign values to parameters \n")
        ind = 0
        for ky in self.param:
            f.write("Interval "+ ky +"= "+"parvals("+str(ind)+");\n")
            ind += 1

        # assign values to error symbols
        f.write("// assign values to error symbols \n")
        ind = 0
        for ky in self.field:
            f.write("Interval "+str(self.errVars[ky])+"= "+"StError("+str(ind)+");\n")
            ind += 1
        ind = 0
        for ky in self.inp:
            f.write("Interval "+str(self.errVars[ky])+"= "+"inp("+str(ind)+")-"+"InpCenter("+str(ind)+");\n")
            ind += 1
            
        # compute linearization error
        f.write("// compute continuous time linearization error\n")
        ind = 0
        for ky in self.linerr:
            mystr = "out("+str(ind)+")"+"= "+cxxcode(self.linerr[ky][0]).replace("std::","")+";\n";
            for repl in ["cos","sin","tan"]:
                mystr = mystr.replace(repl,"my"+repl)
            f.write(mystr)
            ind += 1        
        f.write("\n")
        # return output
        f.write("return out;\n")
        # close function
        f.write("}")
        f.close()

    def WriteDimError(self):
        # clear previous content on file
        f = open(self.DimErrorFile, "w")
        f.write("// Function to compute continuous time linearization error\n")
        f.close()
        # Initialize Function
        f = open(self.DimErrorFile,"a")
        # open function 
        f.write("Interval ContDimError(const IvVectorNd &state, const IvVectorMd &inp,"+
        "const IvVectorNd &StError, const IvVectorMd &InpCenter, const IvVectorKd &parvals, int dim)"+"{\n")
        # define output of function
        f.write("// define output\n")
        f.write("Interval out;\n")
        # assign values to state symbols
        f.write("//assign values to state symbols\n")
        ind = 0;
        for ky in self.field:
            f.write("Interval "+ ky +"= "+"state("+str(ind)+");\n")
            ind += 1
        # assign values to input symbols
        f.write("// assign values to input symbols \n")
        ind = 0
        for v in self.inp:
            f.write("Interval "+ v +"= "+"inp("+str(ind)+");\n")
            ind += 1
            
        # assign values to parameters
        f.write("// assign values to parameters \n")
        ind = 0
        for ky in self.param:
            f.write("Interval "+ ky +"= "+"parvals("+str(ind)+");\n")
            ind += 1
            
        # assign values to error symbols
        f.write("// assign values to error symbols \n")
        ind = 0
        for ky in self.field:
            f.write("Interval "+str(self.errVars[ky])+"= "+"StError("+str(ind)+");\n")
            ind += 1
        ind = 0
        for ky in self.inp:
            f.write("Interval "+str(self.errVars[ky])+"= "+"inp("+str(ind)+")-"+"InpCenter("+str(ind)+");\n")
            ind += 1
            
        # compute linearization error
        f.write("// compute continuous time linearization error\n")
        ind = 0
        for ky in self.linerr:
            mystr = "out = "+cxxcode(self.linerr[ky][0]).replace("std::","")+";\n";
            for repl in ["cos","sin","tan"]:
                mystr = mystr.replace(repl,"my"+repl)
            # open if statement
            f.write("\n if(dim=="+str(ind)+"){\n")
            f.write(mystr)
            # close if statement
            f.write("}\n")
            ind += 1        
        f.write("\n")
        # return output
        f.write("return out;\n")
        # close function
        f.write("}")
        f.close()

    def WriteTemp(self):
        # compute value replacements in state action matrix
        repl = []
        for ky in self.field:
            repl.append( (self.initState[ ky ][0]+self.initState[ ky ][1])/2 )
        for ky in self.inp:
            repl.append( (self.inp[ ky ][0]+self.inp[ ky ][0])/2 )
        for v in self.paramVars:
            repl.append( self.param[ str(v) ] )

        # compute state action matrix at origin
        stmatfn = lambdify( self.allVars + self.paramVars, self.stJacobian, modules="numpy" )
        stmatorigin = stmatfn( *repl )

        # polar decomposition
        Q, _ = scl.polar(stmatorigin)
        
        # compute eigenvectors
        _,V = np.linalg.eig(stmatorigin);
        ReV = np.real(V);
        ImV = np.imag(V);

        # write eigenvectors to txt file (1st column is real part, 2nd column is imag part)
        np.savetxt( "src/pywrite/eigRe.txt", np.transpose(stmatorigin), delimiter = " " )
        np.savetxt( "src/pywrite/eigIm.txt", ImV*0, delimiter = " " )

        # write directions matrix at origin and its pseudo inverse
        np.savetxt( "src/pywrite/dirmat.txt", self.dirmat, delimiter = " " )
        
        
#end-class====================================================================================================

#====================================================================================================
# function to process dynamics
#====================================================================================================
def preprocess( model, server = None, compiler = "g++-11" ):
    symstart = time.time() # start timer for symbolic processing
    # perform symbolic processing and write c++ files
    writetocpp( model )
    symend = time.time()
    
    # compile c++ file
    if server is None:
        compstart = time.time()
        if compiler == "g++-9":
            os.system( "make compile9" )
        else:
            os.system( "make compile11" )
        compend = time.time()
        print( "time for symbolic processing and compiling is ", compend+symend-compstart-symstart  )
        
    else:
        pth = server[ "toolpath" ]
        # remove pywrite and results directory in host
        execstr = "ssh " + server[ "name" ] + " \"" + "rm -rf "
        execstr += pth + "/src/pywrite "
        execstr += pth + "/results; "

        # make fresh pywrite directory
        execstr += "mkdir "
        execstr += pth + "/src/pywrite; "

        # make fresh results directory
        execstr += "mkdir "
        execstr += pth + "/results\" "
        os.system( execstr )
        
        # copy files to pywrite directory
        execstr = "scp " + "src/pywrite/* " + server[ "name" ]
        execstr += ":" + pth + "/src/pywrite/"
        os.system( execstr )

        # compile inside server
        execstr = "ssh " + server[ "name" ] + " \"cd " + pth 
        execstr += "/;" + " make compile9; cd\""
        compstart = time.time()
        os.system( execstr )
        compend = time.time()

        print( "time for symbolic processing and compiling is ", compend+symend-compstart-symstart, " secs" )

#====================================================================================================

#====================================================================================================
# function to simulate
#====================================================================================================
def setvals( model, hypar, server = None ):
    # write initial state bounds
    lbobj = open( "src/pywrite/stlb.txt", "w" )
    ubobj = open( "src/pywrite/stub.txt", "w" )
    for mykey in model[ "initState" ]:
        lbobj.write( str( model[ "initState" ][mykey][0] ) + " ")
        ubobj.write( str( model[ "initState" ][mykey][1] ) + " ")
    lbobj.close()
    ubobj.close()

    # write input bounds
    lbobj = open( "src/pywrite/inplb.txt", "w" )
    ubobj = open( "src/pywrite/inpub.txt", "w" )
    for mykey in model[ "inp" ]:
        lbobj.write( str( model[ "inp" ][mykey][0] ) + " ")
        ubobj.write( str( model[ "inp" ][mykey][1] ) + " ")
    lbobj.close()
    ubobj.close()

    # write values of numerator and denominator of parameters
    numobj = open( "src/pywrite/parnum.txt", "w" )
    denobj = open( "src/pywrite/parden.txt", "w" )
    for mykey in model[ "param" ]:
        numpar, denpar = fraction( model[ "param" ][ mykey ] )
        numobj.write( str( numpar ) + " " )
        denobj.write( str( denpar ) + " " )
    # close file objects
    numobj.close()
    denobj.close()

    # write simulation hyperparameters
    simwrite = open( "src/pywrite/simpars.txt", "w" )
    simstr = str( hypar[ "maxTime" ] ) + " " + str( hypar[ "logDivs" ] ) + " " + str( hypar[ "timeStep" ] ) + " " \
        + str( hypar[ "zonOrder" ] )
    simwrite.write( simstr )
    simwrite.close()

    # run executable
    if server is None:
        os.system( "make execute" )        
    else:
        pth = server[ "toolpath" ]
        # copy files to server
        execstr = "scp " + "src/pywrite/*.txt " + server[ "name" ]
        execstr += ":" + pth + "/src/pywrite/"
        os.system( execstr )
        
        # execute inside server
        #execstr =  "ssh " + server[ "name" ] + " \"cd "
        #execstr += pth + "/; " + "make execute\"; cd;"
        #os.system( execstr )


    



