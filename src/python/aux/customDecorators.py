"""
This module contains a list of custom decorators
"""
# Copyright (C) 2014 Lento Manickathan
#
# This file is part of pHyFlow.
#
# pHyFlow is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pHyFlow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with pHyFlow. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2014-02-28
# Last changed: 2014-03-04


# Standard property (only has get function)
def simpleGetProperty(func):
    # Get the function name
    name = func.__name__
    
    # Define set error function
    def setError(instance, value):
        raise AttributeError("Not possible to manually set '%s' !" % name)    
        
    # Define delete error function
    def deleteError(instance):
        raise AttributeError("Not possible to manually delete '%s' !" % name)           
    
    # Make the property with error set and delete
    prop = property(fget=func, fset=setError, fdel=deleteError, doc=func.__doc__)
    return prop
