from logging import getLogger
logger = getLogger('uap_logger')

class UAPError(StandardError):
    pass
    """Exception raised for ayn error in the UAP.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        logger.error(message)
        super(UAPError, self).__init__(message)
