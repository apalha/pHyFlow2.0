"""
This module contains a list of custom decorators


"""



# Standard property (only has get function)
def simpleGetProperty(func):
    #    locals = sys._getframe(1).f_locals
    name = func.__name__

    def setError(instance, value):
        raise AttributeError("Not possible to manually set '%s' !" % name)    
        
    def deleteError(instance):
        raise AttributeError("Not possible to manually delete '%s' !" % name)           
    
    #    prop = locals.get(name)
    prop = property(fget=func, fset=setError, fdel=deleteError, doc=func.__doc__)
    return prop
