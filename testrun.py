from src.symprocess import *
from examples.unicycle.model import model, hypar, server

server = None

preprocess( model, server )
setvals( model, hypar, server )
