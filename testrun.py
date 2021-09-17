from src.symprocess import *
from examples.autocar.platoon2cars import autocar, hypar, server

server = None

preprocess( autocar, server )
setvals( autocar, hypar, server )
