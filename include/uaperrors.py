from logging import getLogger
logger = getLogger('uap_logger')


class UAPError(Exception):
    """Exception raised for any error in the UAP.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        super(UAPError, self).__init__(message)


class StepError(UAPError):
    """Exception raised for any error in a UAP step.

    Attributes:
        step -- the AbstractStep instance it is called in
        message -- explanation of the error
    """

    def __init__(self, step, message):
        step_message = '[%s] %s' % (step.get_step_name(), message)
        super(StepError, self).__init__(step_message)
