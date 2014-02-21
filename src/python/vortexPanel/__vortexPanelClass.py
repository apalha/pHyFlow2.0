#-*- coding: utf-8 -*-
__doc__ = """

VORTEXPANELCLASS
================

The main class for the vortex-panel coupled.

* Note: No Turbulence scheme implemented (only laminar flow)

Description
-----------
This class wraps the vortex blobs and panel together.

Methodology
-----------


Implemented NS Algorithms
-------------------------


:First Added:   2014-02-21
:Last Modified: 2014-02-21                         
:Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
:License:       GNU GPL version 3 or any later version
"""

#   Copyright (C) 2014 Lento Manickathan                                                                         
#   
#   This file is part of pHyFlow.                                                                                                      
#   
#   pHyFlow is free software: you can redistribute it and/or modify                                                                    
#   it under the terms of the GNU Lesser General Public License as published by                                                       
#   the Free Software Foundation, either version 3 of the License, or                                                                 
#   (at your option) any later version.                                                                                               
#   
#   pHyFlow is distributed in the hope that it will be useful,                                                                         
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                    
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                                      
#   GNU Lesser General Public License for more details.                                                                               
#      
#   You should have received a copy of the GNU Lesser General Public License                                                          
#   along with pHyFlow. If not, see <http://www.gnu.org/licenses/>.    

__all__ = ['VortexPanel']

# External packages
import numpy as _numpy

# Import pHyFlow
from pHyFlow.vortexMethod import base as _base
from pHyFlow.vortexMethod import vmOptions as _vmOptions

from pHyFlow.vortex import VortexBlobs as _VortexBlobs
from pHyFlow.panel import Panels as _Panels

class VortexPanel(object):
    r"""
    
    Usage
    -----
    .. code-block:: python
        
    Parameters
    ----------

    Attribute
    ---------
   
    Methods
    -------
   

    :First Added:   2014-02-21
    :Last Modified: 2014-02-21                         
    :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
    :License:       GNU GPL version 3 or any later version   
   
    """
    """
    Revisions
    ---------


    """
    
    def __init__(self,blobs,panels):
        
        
        #---------------------------------------------------------------------
        # Check/Set input parameters
            
        # Initialize the parameters
        self.__set('blobs',blobs)
        self.__set('panels', panels)
            
            
    def __set(self,varName,var):
        """
        Function to set parameters.
        """
        
        if varName == 'blobs':
            if type(var) != _VortexBlobs:
                raise ValueError("'blobs' should be of type pHyFlow.vortex.VortexBlobs. It is %s" % type(var))
                
            self.blobs = var
            
        elif varName == 'panels':
            if type(var) != _Panels:
                raise ValueError("'blobs' should be of type pHyFlow.vortex.VortexBlobs. It is %s" % type(var))
                
            self.panels = var
            
          
    