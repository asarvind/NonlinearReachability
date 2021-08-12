from src.symprocess import *
from examples.autocar.model import autocar, hypar, server

server = None

preprocess( autocar, server )
setvals( autocar, hypar, server )
