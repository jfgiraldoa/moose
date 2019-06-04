#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

from collections import namedtuple

def initStatus():
    status = namedtuple('status', 'status color code')
    return status

class StatusSystemError(Exception):
    pass

class StatusSystem(object):
    """
    A Class for supplying statuses, with status text color and corresponding exit codes.

    Syntax:
    status = StatusSystem()

    status.getStatus()
      returns a named tuple:
      status(status='NA', color='GREY', code=0x0)

    status.setStatus(status.fail)
      returns a named tuple:
      status(status='FAIL', color='RED', code=0x80)


    Available statuses:

      no_status   the default status when instanced
      success     exit code 0, passing
      skip        exit code 0, skipped
      silent      exit code 0, skipped and also silent
      fail        exit code 0x80, an general error
      diff        exit code 0x81, an error due to exodiff, csvdiff
      deleted     exit code 0x83, some tests are instructed not to run (never fixable), but when instructed to do so they fail
      error       exit code 0x84, an error caused by the TestHarness itself
      timeout     exit code 0x1, an error caused by a test exceeding it's max_time
      hold        exit code 0, used to identify the state of a test as it moves around in the TestHarness
      queued      exit code 0, used to identify the state of a test as it moves around in the TestHarness
      running     exit code 0, used to identify the state of a test as it moves around in the TestHarness
      finished    exit code 0, used to identify the state of a test as it moves around in the TestHarness
    """
    status = initStatus()

    # Default statuses
    no_status = status(status='NA', color='GREY', code=0x0)

    # exit-zero statuses
    success = status(status='OK', color='GREEN', code=0x0)
    skip = status(status='SKIP', color='GREY', code=0x0)
    silent = status(status='SILENT', color='GREY', code=0x0)

    # non-zero statuses
    fail = status(status='FAIL', color='RED', code=0x80)
    diff = status(status='DIFF', color='YELLOW', code=0x81)
    deleted = status(status='DELETED', color='RED', code=0x83)
    error  = status(status='ERROR', color='RED', code=0x80)
    timeout  = status(status='TIMEOUT', color='RED', code=0x1)

    # Pending statuses
    hold  = status(status='HOLD', color='CYAN', code=0x0)
    queued  = status(status='QUEUED', color='CYAN', code=0x0)
    running  = status(status='RUNNING', color='CYAN', code=0x0)

    # all-encompassing finished status
    finished  = status(status='FINISHED', color='GREY', code=0x0)

    __all_statuses = [no_status,
                      success,
                      skip,
                      silent,
                      fail,
                      diff,
                      deleted,
                      error,
                      timeout,
                      hold,
                      queued,
                      running,
                      finished]

    __exit_nonzero_statuses = [fail,
                               diff,
                               deleted,
                               error,
                               timeout]

    __exit_zero_statuses = [success,
                            skip,
                            silent]

    __pending_statuses = [hold,
                          queued,
                          running]

    def __init__(self):
        self.__status = self.no_status

    def createStatus(self, status_key='NA'):
        """ return a specific status object based on supplied status name """
        for status in self.__all_statuses:
            if status_key == status.status:
                return status

    def getStatus(self):
        """
        Return the status object.
        """
        return self.__status

    def getAllStatuses(self):
        """ return list of named tuples containing all status types """
        return self.__all_statuses

    def getFailingStatuses(self):
        """ return list of named tuples containing failing status types """
        return self.__exit_nonzero_statuses

    def getSuccessStatuses(self):
        """ return list of named tuples containing exit code zero status types """
        return self.__exit_zero_statuses

    def getPendingStatuses(self):
        """ return list of named tuples containing pending status types """
        return self.__pending_statuses

    def setStatus(self, status=no_status):
        """
        Set the current status to status. If status is not supplied, 'no_status' is implied.
        There is a validation check during this process to ensure the named tuple adheres to
        this class's set statuses.
        """
        if self.isValid(status):
            self.__status = status
        else:
            raise StatusSystemError('Invalid status! %s' % (str(status)))
        return self.__status

    def isValid(self, status):
        original = set(self.no_status._asdict().keys())
        altered = set(status._asdict().keys())
        if not original.difference(altered) or status in self.__all_statuses:
            return True
