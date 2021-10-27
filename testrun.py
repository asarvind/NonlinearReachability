from src.symprocess import *
#from examples.autocar.platoon2cars import autocar, hypar, server
from examples.autocar.platoon2cars import model, hypar, server

server = None

# preprocess( autocar, server )
# setvals( autocar, hypar, server )
preprocess( model, server )
setvals( model, hypar, server )
