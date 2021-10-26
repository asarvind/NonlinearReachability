import  numpy as np

#====================================================================================================
# define autocar model
#====================================================================================================
# initialize autocar
model = {}

# define field of autocar
field = {}
field["x1"] = "x4*cos(x5+x7)"
#field["x2"] = "x4*sin(x5+x7)"
field["x3"] = "-K0*(x5+x7+x3)"
field["x4"] = "(K2*(5.5-x4) + u1)"
field["x5"] = "x6"
field["x6"] = "mu*m/(Iz*(lr+lf))*(lf*CSf*(g*lr-(K2*(5.5-x4) + u1)*hcg)*x3+(lr*CSr*(g*lf+(K2*(5.5-x4) + u1)*hcg)-lf*CSf*(g*lr-(K2*(5.5-x4) + u1)*hcg))*x7-(lf*lf*CSf*(g*lr-(K2*(5.5-x4) + u1)*hcg) + lr*lr*CSr*(g*lf+(K2*(5.5-x4) + u1)*hcg))*x6/x4)"
field["x7"]  = "mu/(x4*(lr+lf))*(CSf*(g*lr-(K2*(5.5-x4) + u1)*hcg)*x3+(CSr*(g*lf+(K2*(5.5-x4) + u1)*hcg)-CSf*(g*lr-(K2*(5.5-x4) + u1)*hcg))*x7-(lf*CSf*(g*lr-(K2*(5.5-x4) + u1)*hcg) + lr*CSr*(g*lf+(K2*(5.5-x4) + u1)*hcg))*x6/x4)-x6"

# define field of autocar
field["y1"] = "y4*cos(y5+y7)"
#field["y2"] = "y4*sin(y5+y7)"
field["y3"] = "-K0*(y5+y7+y3)"
field["y4"] = "( K2*(x4-y4) )"
field["y5"] = "y6"
field["y6"] = "mu*m/(Iz*(lr+lf))*(lf*CSf*(g*lr-(K2*(x4-y4)+u2)*hcg)*y3+(lr*CSr*(g*lf+(K2*(x4-y4)+u2)*hcg)-lf*CSf*(g*lr-(K2*(x4-y4)+u2)*hcg))*y7-(lf*lf*CSf*(g*lr-(K2*(x4-y4)+u2)*hcg) + lr*lr*CSr*(g*lf+(K2*(x4-y4)+u2)*hcg))*y6/y4)"
field["y7"]  = "mu/(y4*(lr+lf))*(CSf*(g*lr-(K2*(x4-y4)+u2)*hcg)*y3+(CSr*(g*lf+(K2*(x4-y4)+u2)*hcg)-CSf*(g*lr-(K2*(x4-y4)+u2)*hcg))*y7-(lf*CSf*(g*lr-(K2*(x4-y4)+u2)*hcg) + lr*CSr*(g*lf+(K2*(x4-y4)+u2)*hcg))*y6/y4)-y6"

model["field"] = field

# define parameters of autocar
model["param"] =  {
    "g" : 9.81,
    "m" : 1093.3,
    "mu" : 1.0489,
    "lf" : 1.156,
    "lr" : 1.422,
    "hcg" : 0.1,
    #"hcg" : 0.574,
    "Iz" : 1791.6,
    "CSf" : 20.89,
    "CSr" : 20.89,
    "K0" : 4,
    "K1" : 0,
    "K2" : 3,
}

# define input of autocar
model["inp"] = {
    "u1" : [-0.01, 0.01],
    "u2" : [-0.01, 0.01]
}

# define initial state of autocar
initState = {}
s = 1
initState[ "x1" ] = [0, 0]
#initState[ "x2" ] = [0, 0]
initState[ "x3" ] = [-0.5*s, 0.5*s]
initState[ "x4" ] = [8, 9]
initState[ "x5" ] = [-0.3*s, 0.3*s]
initState[ "x6" ] = [-0.25*s, 0.25*s]
initState[ "x7" ] = [-0.3*s, 0.3*s]

initState[ "y1" ] = [-25 , -25]
#initState[ "y2" ] = [-0.5, 0.5]
initState[ "y3" ] = [-0.1*s, 0.1*s]
initState[ "y4" ] = [5, 9]
initState[ "y5" ] = [-0.05*s, 0.05*s]
initState[ "y6" ] = [-0.1*s, 0.1*s]
initState[ "y7" ] = [-0.05*s, 0.05*s]

model[ "initState" ] = initState

# directions of flowpipe bounds
dirmat = np.identity(12);
dirmat[0,6] = -1;
model["directions"] = dirmat;

#====================================================================================================
# define hyperparameters
#====================================================================================================
hypar = {}
hypar[ "timeStep" ] = 0.005
hypar[ "maxTime" ] = 5
hypar[ "zonOrder" ] = 50
hypar[ "logDivs" ] = 1

#====================================================================================================
# set server
#====================================================================================================
server = {}
server["name"] = "amazon"
server["toolpath"] = "~/main/NonlinearReachability"
