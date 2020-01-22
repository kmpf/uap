from logging import getLogger
logger = getLogger('uap_logger')

class UAPError(StandardError):
    pass
    """Exception raised for ayn error in the UAP.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, *args):
        logger.error(*args)
        message = 'No info.'
        if len(args)>0 and isinstance(args[0], str):
            message = args[0]
        super(UAPError, self).__init__(message)
