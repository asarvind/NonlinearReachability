from src.symprocess import *
from examples.autocar.platoon2cars import autocar, hypar, server

server = None
compiler = "g++-9"

preprocess( autocar, server )
setvals( autocar, hypar, server )
