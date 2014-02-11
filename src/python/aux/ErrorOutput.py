__doc__ = """

ErrorOutput master class.

Description
-----------
Implements standard error message format output across all modules of pHyFlow.


Implemented Algorithms
----------------------


References
----------



:First added:   2014-02-06
:Last updated:  2013-02-06
:Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
:License:       GNU GPL version 3 or any later version

"""

"""
Reviews:
-------


"""


# define all text colors
__TEXT_BLACK = '\033[0;30m'
__TEXT_RED = '\033[0;31m'
__TEXT_GREEN = '\033[0;32m'
__TEXT_YELLOW = '\033[0;33m'
__TEXT_BLUE = '\033[0;34m'
__TEXT_PURPLE = '\033[0;35m'
__TEXT_CYAN = '\033[0;36m'
__TEXT_LIGHTGRAY = '\033[0;37m'
__TEXT_DARKGRAY = '\033[1;30m'
__TEXT_BOLDRED = '\033[1;31m'
__TEXT_BOLDGREEN = '\033[1;32m'
__TEXT_BOLDYELLOW = '\033[1;33m'
__TEXT_BOLDBLUE = '\033[1;34m'
__TEXT_BOLDPURPLE = '\033[1;35m'
__TEXT_BOLDCYAN = '\033[1;36m'
__TEXT_WHITE = '\033[1;37m'

# define all background colors
__BACKGROUND_BLACK = '\033[0;40m'
__BACKGROUND_RED = '\033[0;41m'
__BACKGROUND_GREEN = '\033[0;42m'
__BACKGROUND_YELLOW = '\033[0;43m'
__BACKGROUND_BLUE = '\033[0;44m'
__BACKGROUND_PURPLE = '\033[0;45m'
__BACKGROUND_CYAN = '\033[0;46m'
__BACKGROUND_LIGHTGRAY = '\033[0;47m'
__BACKGROUND_DARKGRAY = '\033[1;40m'
__BACKGROUND_BOLDRED = '\033[1;41m'
__BACKGROUND_BOLDGREEN = '\033[1;42m'
__BACKGROUND_BOLDYELLOW = '\033[1;43m'
__BACKGROUND_BOLDBLUE = '\033[1;44m'
__BACKGROUND_BOLDPURPLE = '\033[1;45m'
__BACKGROUND_BOLDCYAN = '\033[1;46m'
__BACKGROUND_WHITE = '\033[1;47m'

# define the termination code
__ENDC = '\033[0m'

# define the functions for specific output format

def errorMessage(message,errorType,errorOrigin):
    r"""
    Generates the string of an error message, using a specific formatting for error messages. Displays the error message,
    the type of error and the origin of the error.

    Usage
    -----
    .. code-block errorMessage(message,errorType,errorOrigin)

    Parameters
    ----------
        message : string
                  the string containing the text that explains the error.
                  shape: (1,)
        errorType : string
                    the string containing the text that identifies the errorType.
                    It can be a standard python error like ValueError, TypeError, etc.,
                    or a user defined error.
                    shape: (1,)
        errorOrigin : string
                      the string containing the text that identifies the origin of the error.
                      That is, the function where the error was generated.
                      shape: (1,)

    Returns
    -------
        formattedMessage : string
                           the formatted string containing the error message ready to be printed.
                           shape: (1,)


    :First Added:   2014-02-06
    :Last Modified: 2014-02-06
    :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
    :License:       GNU GPL version 3 or any later version

    """

    """
    Reviews:

    """

    # format the message with the errorType in bold red, the message in red and the origin of the error in black
    # errorType : message
    # Origin : error Origin
    formattedMessage = __TEXT_BOLDRED + errorType + ' : ' + __ENDC + __TEXT_RED + message + __ENDC + '\n' + __TEXT_BOLDBLUE + 'Origin : ' + __ENDC + errorOrigin

    return formattedMessage