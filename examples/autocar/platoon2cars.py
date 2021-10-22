import  numpy as np

#====================================================================================================
# define autocar model
#====================================================================================================
# initialize autocar
autocar = {}

# define field of autocar
field = {}
field["x1"] = "x4*cos(x5+x7)"
#field["x2"] = "x4*sin(x5+x7)"
field["x3"] = "-K0*(x5+x7+x3)"
field["x4"] = "K2*(5.5-x4) + u"
field["x5"] = "x6"
field["x6"] = "mu*m/(Iz*(lr+lf))*(lf*CSf*(g*lr-u*hcg)*x3+(lr*CSr*(g*lf+u*hcg)-lf*CSf*(g*lr-u*hcg))*x7-(lf*lf*CSf*(g*lr-u*hcg) + lr*lr*CSr*(g*lf+u*hcg))*x6/x4)"
field["x7"]  = "mu/(x4*(lr+lf))*(CSf*(g*lr-u*hcg)*x3+(CSr*(g*lf+u*hcg)-CSf*(g*lr-u*hcg))*x7-(lf*CSf*(g*lr-u*hcg) + lr*CSr*(g*lf+u*hcg))*x6/x4)-x6"

# define field of autocar
field["y1"] = "y4*cos(y5+y7)"
#field["y2"] = "y4*sin(y5+y7)"
field["y3"] = "-K0*(y5+y7+y3)"
field["y4"] = "( K2*(x4-y4) + K1*(x1-y1 - 20) )"
field["y5"] = "y6"
field["y6"] = "mu*m/(Iz*(lr+lf))*(lf*CSf*(g*lr-u*hcg)*y3+(lr*CSr*(g*lf+u*hcg)-lf*CSf*(g*lr-u*hcg))*y7-(lf*lf*CSf*(g*lr-u*hcg) + lr*lr*CSr*(g*lf+u*hcg))*y6/y4)"
field["y7"]  = "mu/(y4*(lr+lf))*(CSf*(g*lr-u*hcg)*y3+(CSr*(g*lf+u*hcg)-CSf*(g*lr-u*hcg))*y7-(lf*CSf*(g*lr-u*hcg) + lr*CSr*(g*lf+u*hcg))*y6/y4)-y6"

autocar["field"] = field

# define parameters of autocar
autocar["param"] =  {
    "g" : 9.81,
    "m" : 1093.3,
    "mu" : 1.0489,
    "lf" : 1.156,
    "lr" : 1.422,
    "hcg" : 0.574,
    "Iz" : 1791.6,
    "CSf" : 20.89,
    "CSr" : 20.89,
    "K0" : 4,
    "K1" : 0,
    "K2" : 3,
}

# define input of autocar
autocar["inp"] = {
    "u" : [-0.01, 0.01]
}

# define initial state of autocar
initState = {}
initState[ "x1" ] = [0, 0]
#initState[ "x2" ] = [0, 0]
initState[ "x3" ] = [-0.5, 0.5]
initState[ "x4" ] = [8, 9]
initState[ "x5" ] = [-0.2, 0.2]
initState[ "x6" ] = [-0.25, 0.25]
initState[ "x7" ] = [-0.2, 0.2]

initState[ "y1" ] = [-25 , -25]
#initState[ "y2" ] = [-0.5, 0.5]
initState[ "y3" ] = [-0.1, 0.1]
initState[ "y4" ] = [5, 9]
initState[ "y5" ] = [-0.05, 0.05]
initState[ "y6" ] = [-0.1, 0.1]
initState[ "y7" ] = [-0.05, 0.05]

autocar[ "initState" ] = initState

# directions of flowpipe bounds
dirmat = np.identity(10);
dirmat[6,0] = 1;
dirmat[6,6] = -1;
autocar["directions"] = dirmat;

#====================================================================================================
# define hyperparameters
#====================================================================================================
hypar = {}
hypar[ "timeStep" ] = 0.005
hypar[ "maxTime" ] = 5
hypar[ "zonOrder" ] = 100
hypar[ "logDivs" ] = 0

#====================================================================================================
# set server
#====================================================================================================
server = {}
server["name"] = "amazon"
server["toolpath"] = "~/main/NonlinearReachability"
