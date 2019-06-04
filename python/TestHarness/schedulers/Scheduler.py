#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

from TestHarness.JobDAG import JobDAG
from FactorySystem.MooseObject import MooseObject
import os, traceback
from time import sleep
from timeit import default_timer as clock
from multiprocessing.pool import ThreadPool
import threading # for thread locking and thread timers

class SchedulerError(Exception):
    pass

class Scheduler(MooseObject):
    """
    Base class for handling jobs asynchronously. To use this class, call .schedule()
    and supply a list of testers to schedule. Each group of testers supplied will begin
    running immediately.

    Syntax:
       .schedule([list of tester objects])

    A list of testers will be added to a queue and begin calling their derived run method.
    You can continue to add more testers to the queue in this fashion.

    Once all jobs have been scheduled, call .waitFinish() to wait until all jobs have
    finished.
    """

    @staticmethod
    def validParams():
        params = MooseObject.validParams()
        params.addRequiredParam('average_load',  64.0, "Average load to allow")
        params.addRequiredParam('max_processes', None, "Hard limit of maxium processes to use")
        params.addParam('min_reported_time', 10, "The minimum time elapsed before a job is reported as taking to long to run.")

        return params

    # This is what will be checked for when we look for valid schedulers
    IS_SCHEDULER = True

    def __init__(self, harness, params):
        MooseObject.__init__(self, harness, params)

        ## The test harness to run callbacks on
        self.harness = harness

        # Retrieve and store the TestHarness options for use in this object
        self.options = harness.getOptions()

        # The Scheduler class can be initialized with no "max_processes" argument and it'll default
        # to a soft limit. If however a max_processes is passed we'll treat it as a hard limit.
        # The difference is whether or not we allow single jobs to exceed the number of slots.
        if params['max_processes'] == None:
            self.available_slots = 1
            self.soft_limit = True
        else:
            self.available_slots = params['max_processes'] # hard limit
            self.soft_limit = False

        self.average_load = params['average_load']

        self.min_report_time = params['min_reported_time']

        # Initialize run_pool based on available slots
        self.run_pool = ThreadPool(processes=self.available_slots)

        # Initialize status_pool to only use 1 process (to prevent status messages from getting clobbered)
        self.status_pool = ThreadPool(processes=1)

        # Slot lock when processing resource allocations and modifying slots_in_use
        self.slot_lock = threading.Lock()

        # Job lock when modifying a jobs status
        self.activity_lock = threading.Lock()

        # A combination of processors + threads (-j/-n) currently in use, that a job requires
        self.slots_in_use = 0

        # Set containing all scheduled jobs
        self.__scheduled_jobs = set([])

        # Set containing jobs entering the run_pool
        self.__job_bank = set([])

        # Total running Job and Test failures encountered
        self.__failures = 0

        # Allow threads to set a global exception
        self.__error_state = False

        # Private set of jobs currently running
        self.__active_jobs = set([])

        # Jobs that are taking longer to finish than the alloted time are reported back early to inform
        # the user 'stuff' is still running. Jobs entering this set will not be reported again.
        self.jobs_reported = set([])

        # The last time the scheduler reported something
        self.last_reported_time = clock()

        # Sets of threading objects created by jobs entering and exiting the queues. When scheduler.waitFinish()
        # is called, and both thread pools are empty, the pools shut down, and the call to waitFinish() returns.
        self.__status_pool_lock = threading.Lock()
        self.__runner_pool_lock = threading.Lock()
        self.__status_pool_jobs = set([])
        self.__runner_pool_jobs = set([])

        # True when scheduler.waitFinish() is called. This alerts the scheduler, no more jobs are
        # to be scheduled. KeyboardInterrupts are then handled by the thread pools.
        self.__waiting = False

    def triggerErrorState(self):
        self.__error_state = True
        self.run_pool.close()
        self.status_pool.close()

    def killRemaining(self, keyboard=False):
        """ Method to kill running jobs """
        with self.activity_lock:
            for job in self.__active_jobs:
                job.killProcess()
        if keyboard:
            self.triggerErrorState()
            self.harness.keyboard_interrupt()
        else:
            self.triggerErrorState()

    def retrieveJobs(self):
        """ return all the jobs the scheduler was tasked to perform work for """
        return self.__scheduled_jobs

    def schedulerError(self):
        """ boolean if the scheduler prematurely exited """
        return self.__error_state and not self.maxFailures()

    def maxFailures(self):
        """ Boolean for hitting max failures """
        return ((self.options.valgrind_mode and self.__failures >= self.options.valgrind_max_fails)
                or self.__failures >= self.options.max_fails)

    def run(self, job):
        """ Call derived run method """
        return

    def notifyFinishedSchedulers(self):
        """ Notify derived schedulers we are finished """
        return

    def augmentJobs(self, Jobs):
        """
        Allow derived schedulers to augment Jobs before they perform work.
        Note: This occurs before we perform a job count sanity check. So
        any additions or subtractions to the number of jobs will result in
        an exception.
        """
        return

    def waitFinish(self):
        """
        Inform the Scheduler there are no further jobs to schedule.
        Return once all jobs have completed.
        """
        self.__waiting = True
        try:
            # wait until there is an error, or if all the queus are empty
            waiting_on_status_pool = True
            waiting_on_runner_pool = True

            while (waiting_on_status_pool or waiting_on_runner_pool) and self.__job_bank:

                if self.__error_state:
                    break

                with self.__status_pool_lock:
                    waiting_on_status_pool = sum(1 for x in self.__status_pool_jobs if not x.ready())
                with self.__runner_pool_lock:
                    waiting_on_runner_pool = sum(1 for x in self.__runner_pool_jobs if not x.ready())

                sleep(0.1)

            # Completed all jobs sanity check
            if not self.__error_state and self.__job_bank:
                raise SchedulerError('Scheduler exiting with different amount of work than what was tasked!')

            if not self.__error_state:
                self.run_pool.close()
                self.run_pool.join()
                self.status_pool.close()
                self.status_pool.join()

            # allow derived schedulers to perform any exit routines
            self.notifyFinishedSchedulers()

        except KeyboardInterrupt:
            self.killRemaining(keyboard=True)

    def schedule(self, testers):
        """
        Generate and submit a group of testers to a thread pool queue for execution.
        """
        # If we are not to schedule any more jobs for some reason, return now
        if self.__error_state:
            return

        # Instance our job DAG, create jobs, and a private lock for this group of jobs (testers)
        Jobs = JobDAG(self.options)
        j_dag = Jobs.createJobs(testers)
        j_lock = threading.Lock()

        # Allow derived schedulers access to the jobs before they launch
        self.augmentJobs(Jobs)

        # job-count to tester-count sanity check
        if j_dag.size() != len(testers):
            raise SchedulerError('Scheduler was going to run a different amount of testers than what was received (something bad happened)!')

        # Store all jobs in the global job bank. As jobs finish, they will be removed from
        # this set. This will function as our final sanity check on 100% job completion
        with j_lock:
            self.__job_bank.update(j_dag.topological_sort())

        # Store all scheduled jobs
        self.__scheduled_jobs.update(j_dag.topological_sort())

        # Launch these jobs to perform work
        self.queueJobs(Jobs, j_lock)

    def queueJobs(self, Jobs, j_lock):
        """
        Determine which queue jobs should enter. Finished jobs are placed in the status
        pool to be printed while all others are placed in the runner pool to perform work.

        A finished job will trigger a change to the Job DAG, which will allow additional
        jobs to become available and ready to enter the runner pool (dependency jobs).
        """
        with j_lock:
            concurrent_jobs = Jobs.getJobsAndAdvance()
            for job in concurrent_jobs:
                if job.isFinished():
                    if not self.status_pool._state:
                        with self.__status_pool_lock:
                            self.__status_pool_jobs.add(self.status_pool.apply_async(self.jobStatus, (job, Jobs, j_lock)))

                elif job.isHold():
                    if not self.run_pool._state:
                        with self.__runner_pool_lock:
                            job.setStatus(job.queued)
                            self.__runner_pool_jobs.add(self.run_pool.apply_async(self.runJob, (job, Jobs, j_lock)))

    def getLoad(self):
        """ Method to return current load average """
        loadAverage = 0.0
        try:
            loadAverage = os.getloadavg()[0]
        except AttributeError:
            pass      # getloadavg() not available in this implementation of os
        return loadAverage

    def satisfyLoad(self):
        """ Method for controlling load average """
        while self.slots_in_use > 1 and self.getLoad() >= self.average_load:
            sleep(1.0)

    def reserveSlots(self, job, j_lock):
        """
        Method which allocates resources to perform the job. Returns bool if job
        should be allowed to run based on available resources.
        """
        # comply with load average
        if self.options.load:
            self.satisfyLoad()

        with self.slot_lock:
            can_run = False
            if self.slots_in_use + job.getSlots() <= self.available_slots:
                can_run = True

            # Check for insufficient slots -soft limit
            elif job.getSlots() > self.available_slots and self.soft_limit:
                job.addCaveats('OVERSIZED')
                can_run = True

            # Check for insufficient slots -hard limit (skip this job)
            elif job.getSlots() > self.available_slots and not self.soft_limit:
                job.addCaveats('insufficient slots')
                with j_lock:
                    job.setStatus(job.skip)

            if can_run:
                self.slots_in_use += job.getSlots()
        return can_run

    def handleTimeoutJob(self, job, j_lock):
        """ Handle jobs that have timed out """
        with j_lock:
            if job.isRunning():
                job.setStatus(job.timeout, 'TIMEOUT')
                job.killProcess()

    def handleLongRunningJob(self, job, Jobs, j_lock):
        """ Handle jobs that have not reported in the alotted time """
        with self.__status_pool_lock:
            self.__status_pool_jobs.add(self.status_pool.apply_async(self.jobStatus, (job, Jobs, j_lock)))

    def jobStatus(self, job, Jobs, j_lock):
        """
        Instruct the TestHarness to print the status of job. This is a serial
        threaded operation, so as to prevent clobbering of text being printed
        to stdout.
        """
        # The pool is closing down due to a failure, or this job has previously been handled:
        #
        # A job which triggers the long_running timer, has a chance to finish before this
        # slower serialized status pool, has a chance to process it. Meaning two of the same
        # jobs now exist in this queue, with a finished status. This method can only work on
        # a finished job object once (a set removal operation occurs to signify scheduled job
        # completion as a sanity check).
        if self.status_pool._state or job not in self.__job_bank:
            return

        # Peform within a try, to allow keyboard ctrl-c
        try:
            with j_lock:
                if job.isRunning():
                    # already reported this job once before
                    if job in self.jobs_reported:
                        return

                    # this job will be reported as 'RUNNING'
                    elif clock() - self.last_reported_time >= self.min_report_time:
                        job.addCaveats('FINISHED')

                        with self.activity_lock:
                            self.jobs_reported.add(job)

                    # TestHarness has not yet been inactive long enough to warrant a report
                    else:
                        # adjust the next report time based on delta of last report time
                        adjusted_interval = max(1, self.min_report_time - max(1, clock() - self.last_reported_time))
                        job.report_timer = threading.Timer(adjusted_interval,
                                                           self.handleLongRunningJob,
                                                           (job, Jobs, j_lock,))
                        job.report_timer.start()
                        return

                # Inform the TestHarness of job status
                self.harness.handleJobStatus(job)

                # Reset activity clock
                if not job.isSilent():
                    self.last_reported_time = clock()

                if job.isFail():
                    self.__failures += 1

                if job.isFinished():
                    if job in self.__job_bank:
                        self.__job_bank.remove(job)
                    else:
                        raise SchedulerError('job accountability failure while working with: %s' % (job.getTestName()))

            # Max failure threshold reached, begin shutdown
            if self.maxFailures():
                self.killRemaining()

        except Exception:
            print('statusWorker Exception: %s' % (traceback.format_exc()))
            self.killRemaining()

        except KeyboardInterrupt:
            self.killRemaining(keyboard=True)

    def runJob(self, job, Jobs, j_lock):
        """ Method the run_pool calls when an available thread becomes ready """
        # Its possible, the queue is just trying to empty. Allow it to do so
        # with out generating overhead
        if self.__error_state:
            return

        try:
            # see if we have enough slots to start this job
            if self.reserveSlots(job, j_lock):
                with j_lock:
                    job.setStatus(job.running)

                with self.activity_lock:
                    self.__active_jobs.add(job)

                timeout_timer = threading.Timer(float(job.getMaxTime()),
                                                self.handleTimeoutJob,
                                                (job, j_lock,))

                job.report_timer = threading.Timer(self.min_report_time,
                                                   self.handleLongRunningJob,
                                                   (job, Jobs, j_lock,))

                job.report_timer.start()
                timeout_timer.start()
                self.run(job) # Hand execution over to derived scheduler
                timeout_timer.cancel()

                # Recover worker count before attempting to queue more jobs
                with self.slot_lock:
                    self.slots_in_use = max(0, self.slots_in_use - job.getSlots())

                # Stop the long running timer
                job.report_timer.cancel()

                # All done
                with j_lock:
                    job.setStatus(job.finished)

                with self.activity_lock:
                    self.__active_jobs.remove(job)

            # Not enough slots to run the job...
            else:
                # ...currently, place back on hold before placing it back into the queue
                if not job.isFinished():
                    with j_lock:
                        job.setStatus(job.hold)
                    sleep(.1)

            # Job is done (or needs to re-enter the queue)
            self.queueJobs(Jobs, j_lock)

        except Exception:
            print('runWorker Exception: %s' % (traceback.format_exc()))
            self.killRemaining()

        except KeyboardInterrupt:
            self.killRemaining(keyboard=True)
