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
    def __init__(self,dynamics,state_vars,inp_vars,params=[],assignments=[],vals=[]):
        self.state_vars = state_vars
        self.inp_vars = inp_vars
        mydict = parDynamics(dynamics, state_vars, inp_vars,params,assignments,vals)
        self.assignments = mydict["assignments"]
        self.params = mydict["params"]
        self.vals  = mydict["vals"]
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
        # dimensions of state and input space
        self.N = len(self.state_vars)
        self.M = len(self.inp_vars)
        # time constants
        self.TimeStep = 0.01 # step size
        # store dynamics along each dimension and their corresponding lambda functions
        self.dynamics = {}
        # std string for replacement
        for v in self.state_vars:
            self.dynamics[str(v)] = mydict["dynamics"][str(v)]
        # state coefficient matrix  after linearization
        # also compute linearization matrix at origin 
        self.StMat = np.zeros((self.N,self.N))
        varis = self.inp_vars+self.state_vars
        self.statemat = [[None for j in range(self.N)] for i in range(self.N)]
        for i in range(self.N):
            for j in range(self.N):
                self.statemat[i][j] = simplify(dynamics[str(self.state_vars[i])].diff(self.state_vars[j]))
                StFun = lambdify(self.state_vars+self.inp_vars+self.params,self.statemat[i][j])
                self.StMat[i,j] = StFun(*(np.zeros(self.N+self.M).tolist() + self.vals))
        # inp coefficient matrix after expressions
        self.inpmat = [[None for j in range(self.N)] for i in range(self.N)]
        self.InpMat = np.zeros((self.N,self.M))
        for j in range(self.M):
            for i in range(self.N):
                self.inpmat[i][j] = simplify(dynamics[str(self.state_vars[i])].diff(self.inp_vars[j]))
                InpFun = lambdify(self.state_vars+self.inp_vars+self.params,self.inpmat[i][j])
                self.InpMat[i,j] = InpFun(*(np.zeros(self.N+self.M).tolist() + self.vals))
        self.linerr = {}
        for i in range(self.N):
            ky = str(self.state_vars[i])
            out_expr = 0
            for r in self.state_vars+self.inp_vars:
                for s in self.state_vars+self.inp_vars:
                    if simplify(r-s) == 0:
                        expr = (dynamics[ky].diff(r,s))/4
                    else:
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
        # write time constants
        self.TimeConstants()
        f.write("const double StepSize = "+str(self.TimeStep)+";\n")
        # write parametric assignments
        for asng in self.assignments:
            f.write("const Interval " + asng + "\n")
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
        f.write("IvVectorNd VectorField(const IvVectorNd &state, const IvVectorMd &inp)"+"{\n")
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
        f.write("IvMatrixNNd StateMat(const IvVectorNd &state, const IvVectorMd &inp)"+"{\n")
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
        f.write("IvMatrixNMd InputMat(const IvVectorNd &state, const IvVectorMd &inp)"+"{\n")
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
        "const IvVectorNd &StError, const IvVectorMd &InpCenter)"+"{\n")
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
        "const IvVectorNd &StError, const IvVectorMd &InpCenter, int dim)"+"{\n")
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
        self.TimeStep = np.minimum(x,y)
        self.TimeStep = np.minimum(self.TimeStep,1)
        n = math.ceil(-1*math.log10(self.TimeStep))
        self.TimeStep = round(self.TimeStep,n)

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

#****************************************************************************************************
# function to convert dynamics with floats to parametric dynamics
#****************************************************************************************************

def parDynamics(dynamics, state_vars, inp_vars, params, assignments, vals):
    out = {}
    out["dynamics"] = dynamics;
    out["assignments"] = assignments;
    out["params"] = params;
    out["vals"] = vals;
    s = {}
    ide_list = []
    token = 0;

    for ky in dynamics:
        for exp in preorder_traversal( dynamics[ky] ):
            if(isinstance(exp,(Float,Integer))):
                ide = "var_"+str(exp)
                if ide not in ide_list:
                    ide_list.append(ide)
                    s[ide] = var("par"+"_"+str(token))
                    out["params"].append(s[ide])
                    out["assignments"].append( str( s[ide] ) + " = " + str(exp) + ";" )
                    out["vals"].append( float(exp) )
                    out["dynamics"][ky] = out["dynamics"][ky].subs( exp, s[ide] )
                    token += 1
    return out;
                    

#end-function====================================================================================================
