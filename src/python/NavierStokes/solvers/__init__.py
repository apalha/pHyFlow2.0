"""
    Import the solver algorithm  the user wants. All the available solver 
    are saved in this folder. The available solvers are:
        
        - Chorin projection scheme      :: ['chorin'] in chorin.py
        
"""

#__author__      = "Lento Manickathan <l.manickathan@student.tudelft.nl>"
#__date__        = "2013-07-10"
#__copyright__   = "Copyright (C) 2013 " + __author__
#__license__     = "GNU GPL version 3 or any later version"


# Wrapper for solver classes
def Solver(solverType):
    "Return solver instance for given solver name"
    exec("from %s import Solver as NamedSolver" % solverType)
    return NamedSolver(solverType)
	
