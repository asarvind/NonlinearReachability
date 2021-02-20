from sympy import *
from sympy.codegen.ast import Assignment
import time # for compute time calculation
import numpy as np
import scipy.linalg as scl
import math

#****************************************************************************************************
# class to process dynamics and write c++ functions
#****************************************************************************************************

class writetocpp:
    def __init__(self,dynamics,state_vars,inp_vars,params=[],vals=[]):
        self.state_vars = state_vars
        self.inp_vars = inp_vars
        self.params = params
        self.vals  = vals
        self.err_var = {}
        new_vars = []
        for v in self.state_vars+self.inp_vars:
            self.err_var[str(v)] = var("e_"+str(v))
            new_vars.append(self.err_var[str(v)])
        # indices of state variables
        self.state_ind = {}
        ind = 0
        for v in self.state_vars:
            self.state_ind[str(v)] = ind
            ind += 1
        # dimensions of state, input and parameters
        self.N = len(self.state_vars)
        self.M = len(self.inp_vars)
        self.K = len( self.params )
        
        # time constants
        self.TimeStep = 0.01 # step size
        # store dynamics along each dimension and their corresponding lambda functions
        self.dynamics = {}
        # std string for replacement
        for v in self.state_vars:
            self.dynamics[str(v)] = dynamics[str(v)]
        # symbolic state coefficient matrix  after linearization
        # also compute concrete linearization matrix at origin 
        self.StMat = np.zeros((self.N,self.N))
        varis = self.inp_vars+self.state_vars
        
        # create list of values of state+input+parameter
        initStUb = np.loadtxt("src/pywrite/stub.txt")
        initStLb = np.loadtxt("src/pywrite/stlb.txt")
        initInpLb = np.loadtxt("src/pywrite/inplb.txt")
        initInpUb = np.loadtxt("src/pywrite/inpub.txt")
        initlist = []
        for i in range(initStUb.size):
            initlist.append( ( initStUb[i] + initStLb[i] )/2.0 )
        if initInpUb.size>1:
            for i in range(initInpUb.size):
                initlist.append( ( initInpUb[i] + initInpLb[i] )/2.0 )
        else:
              initlist.append( (initInpUb+initInpLb)/2 )  
        initlist += self.vals
        
        self.statemat = [[None for j in range(self.N)] for i in range(self.N)]
        for i in range(self.N):
            for j in range(self.N):
                self.statemat[i][j] = simplify(dynamics[str(self.state_vars[i])].diff(self.state_vars[j]))
                StFun = lambdify(self.state_vars+self.inp_vars+self.params,self.statemat[i][j])
                self.StMat[i,j] = StFun(*(initlist))
        # inp coefficient matrix after expressions
        self.inpmat = [[None for j in range(self.N)] for i in range(self.N)]
        self.InpMat = np.zeros((self.N,self.M))
        for j in range(self.M):
            for i in range(self.N):
                self.inpmat[i][j] = simplify(dynamics[str(self.state_vars[i])].diff(self.inp_vars[j]))
                InpFun = lambdify(self.state_vars+self.inp_vars+self.params,self.inpmat[i][j])
                self.InpMat[i,j] = InpFun(*(initlist))
        # taylor error
        self.linerr = {}
        for i in range(self.N):
            ky = str(self.state_vars[i])
            out_expr = 0
            for r in self.state_vars+self.inp_vars:
                for s in self.state_vars+self.inp_vars:
                    expr = (dynamics[ky].diff(r,s))/2
                    out_expr += expr*self.err_var[str(r)]*self.err_var[str(s)]
            self.linerr[ky] = simplify(out_expr)
            
        # filenames for writing
        self.ConstantsFile = "src/pywrite/Constants.cpp"
        self.VectorFieldFile = "src/pywrite/VectorField.cpp"
        self.StateMatFile = "src/pywrite/StateMat.cpp"
        self.InputMatFile = "src/pywrite/InputMat.cpp"
        self.ContErrorFile = "src/pywrite/ContError.cpp"
        self.DimErrorFile = "src/pywrite/DimError.cpp"
        # write constants
        f = open(self.ConstantsFile,"w")
        # write dimension of vector field
        f.write("const int StateDim = "+str(self.N)+";\n")
        # write dimension of input
        f.write("const int InputDim = "+str(self.M)+";\n")
        # write number of parameters
        f.write("const int pardim = "+str(self.K)+";\n")

        # write time constants
        self.TimeConstants()
        f.write("const double StepSize = "+str(self.TimeStep)+";\n")
        f.close()
        self.WriteEig()
        self.WriteVectorField()
        self.WriteStateMat()
        self.WriteInputMat()
        self.WriteContError()
        self.WriteDimError()
        
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
        for v in self.state_vars:
            f.write("Interval "+str(v)+"= "+"state("+str(ind)+");\n")
            ind += 1
            
        # assign values to input symbols
        f.write("// assign values to input symbols \n")
        ind = 0
        for v in self.inp_vars:
            f.write("Interval "+str(v)+"= "+"inp("+str(ind)+");\n")
            ind += 1
            
        # assign values to parameters
        f.write("// assign values to input symbols \n")
        ind = 0
        for v in self.params:
            f.write("Interval "+str(v)+"= "+"parvals("+str(ind)+");\n")
            ind += 1
        
        # compute vector field bounds
        f.write("// compute bounds on vector field\n")
        ind = 0
        for v in self.state_vars:
            ky = str(v)
            mystr = "out("+str(ind)+")"+"= "+cxxcode(self.dynamics[ky]).replace("std::","")+";\n"
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
        for v in self.state_vars:
            f.write("Interval "+str(v)+"= "+"state("+str(ind)+");\n")
            ind += 1
            
        f.write("// assign values to input symbols \n")
        ind = 0
        # assign values to input symbols
        for v in self.inp_vars:
            f.write("Interval "+str(v)+"= "+"inp("+str(ind)+");\n")
            ind += 1

        # assign values to parameters
        f.write("// assign values to parameters \n")
        ind = 0
        for v in self.params:
            f.write("Interval "+str(v)+"= "+"parvals("+str(ind)+");\n")
            ind += 1

        # compute matrix
        f.write("// compute matrix\n")
        for i in range(self.N):
            for j in range(self.N):
                mystr = "out("+str(i)+","+str(j)+")"+"= "+cxxcode(self.statemat[i][j]).replace("std::","")+";\n"
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
        for v in self.state_vars:
            f.write("Interval "+str(v)+"= "+"state("+str(ind)+");\n")
            ind += 1
        f.write("// assign values to input symbols \n")
        ind = 0
        
        # assign values to input symbols
        for v in self.inp_vars:
            f.write("Interval "+str(v)+"= "+"inp("+str(ind)+");\n")
            ind += 1

        # assign values to parameters
        f.write("// assign values to parameters \n")
        ind = 0
        for v in self.params:
            f.write("Interval "+str(v)+"= "+"parvals("+str(ind)+");\n")
            ind += 1

        # compute matrix
        f.write("// compute matrix\n")
        for i in range(self.N):
            for j in range(self.M):
                mystr = "out("+str(i)+","+str(j)+")"+"= "+cxxcode(self.inpmat[i][j]).replace("std::","")+";\n"
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
        for v in self.state_vars:
            f.write("Interval "+str(v)+"= "+"state("+str(ind)+");\n")
            ind += 1
            
        # assign values to input symbols
        f.write("// assign values to input symbols \n")
        ind = 0
        for v in self.inp_vars:
            f.write("Interval "+str(v)+"= "+"inp("+str(ind)+");\n")
            ind += 1
            
        # assign values to parameters
        f.write("// assign values to parameters \n")
        ind = 0
        for v in self.params:
            f.write("Interval "+str(v)+"= "+"parvals("+str(ind)+");\n")
            ind += 1

        # assign values to error symbols
        f.write("// assign values to error symbols \n")
        ind = 0
        for v in self.state_vars:
            ky = str(v)
            f.write("Interval "+str(self.err_var[ky])+"= "+"StError("+str(ind)+");\n")
            ind += 1
        ind = 0
        for v in self.inp_vars:
            ky = str(v)
            f.write("Interval "+str(self.err_var[ky])+"= "+"inp("+str(ind)+")-"+"InpCenter("+str(ind)+");\n")
            ind += 1                      
        # compute linearization error
        f.write("// compute continuous time linearization error\n")
        ind = 0
        for v in self.state_vars:
            ky = str(v)
            mystr = "out("+str(ind)+")"+"= "+cxxcode(self.linerr[ky]).replace("std::","")+";\n";
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
        for v in self.state_vars:
            f.write("Interval "+str(v)+"= "+"state("+str(ind)+");\n")
            ind += 1
        # assign values to input symbols
        f.write("// assign values to input symbols \n")
        ind = 0
        for v in self.inp_vars:
            f.write("Interval "+str(v)+"= "+"inp("+str(ind)+");\n")
            ind += 1
            
        # assign values to parameters
        f.write("// assign values to parameters \n")
        ind = 0
        for v in self.params:
            f.write("Interval "+str(v)+"= "+"parvals("+str(ind)+");\n")
            ind += 1
            
        # assign values to error symbols
        f.write("// assign values to error symbols \n")
        ind = 0
        for v in self.state_vars:
            ky = str(v)
            f.write("Interval "+str(self.err_var[ky])+"= "+"StError("+str(ind)+");\n")
            ind += 1
        ind = 0
        for v in self.inp_vars:
            ky = str(v)
            f.write("Interval "+str(self.err_var[ky])+"= "+"inp("+str(ind)+")-"+"InpCenter("+str(ind)+");\n")
            ind += 1
            
        # compute linearization error
        f.write("// compute continuous time linearization error\n")
        ind = 0
        for v in self.state_vars:
            ky = str(v)
            mystr = "out = "+cxxcode(self.linerr[ky]).replace("std::","")+";\n";
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
        

    def TimeConstants(self):
        # step size
        d = (np.absolute(np.linalg.matrix_power(self.StMat,3))).sum(1)
        x = (6*np.amin((np.reciprocal(d))*0.01))**(1/3)
        d = (np.absolute(np.matmul(np.matmul(self.StMat,self.StMat),self.InpMat))).sum(1)
        y = (2*np.amin(np.reciprocal(d))*0.01)**(1/3)
        # self.TimeStep = np.minimum(x,y)
        # self.TimeStep = np.minimum(self.TimeStep,1)
        # n = math.ceil(-1*math.log10(self.TimeStep))
        # self.TimeStep = round(self.TimeStep,n)
        self.TimeStep = 0.1

    def WriteEig(self):
        _,V = np.linalg.eig(self.StMat);
        ReV = np.real(V);
        ImV = np.imag(V);
        f = open(self.ConstantsFile,"a")
        re_arr_str = "\n" + "const double ReEV[StateDim][StateDim] = {"
        im_arr_str = "\n" + "const double ImEV[StateDim][StateDim] = {"
        for i in range(self.N):
            re_col_str = " {"
            im_col_str = " {"
            for j in range(self.N):
                re_col_str += " " + str(ReV[i,j]) + ", ";
                im_col_str += " " + str(ImV[i,j]) + ", ";
            re_col_str += "},"
            im_col_str += "},"
            re_arr_str += re_col_str
            im_arr_str += im_col_str
        re_arr_str += " };\n"
        im_arr_str += " };\n"
        f.write(re_arr_str)
        f.write(im_arr_str)
        f.close();
        
#end-class====================================================================================================

