from subprocess import check_output
from collections import Counter
from time import sleep
import gc
class LSF:
    '''Abstraction layer for the LSF utility, used to submit jobs from python scripts. \n
    Each LSF class instance has two class fields:
    cmd: is the full shell command to be submitted
    Example: cmd = python hello_world.py

    job_id: Initially set to the empty string,  six figure number used to identify the job for checking it's status.'''

    def submit_command (self, cpu_cores = 1, memory = '', queue = "new-all.q", output = "/dev/null"):
        '''Submits the command command string found in self.cmd
        Parameters:
        cpu_cores - (Integer) number of cores to use
        memory - (Integer) amount of memory in MB
        queue - (String) the name of the queue'''

        queue_string = " -q "+ queue+" "
        nodes_span_string = " span[hosts=1]"
        direct_output = " -u " +output +" " if not output == "" else " "
        mem_string = " rusage[mem="+str(memory)+"] " if not memory=='' else ' '
        cpu_usage_string = " -n " + str(cpu_cores) if cpu_cores > 1 else " "
        bsub_prefix = "bsub "+queue_string + cpu_usage_string
        bsub_prefix += direct_output
        bsub_prefix += " -R '"+mem_string + nodes_span_string+"' "
        self.job_id = str (check_output ([bsub_prefix +" "+ self.cmd], shell=True)).split ('<') [1].split ('>') [0]

    def job_status_running (self):
        '''Checks to see if this specific job is finished.\n Returns False if job is PENDING or RUNNING, else True'''
        if not self.job_id:
            return True
        try:
            gc.collect()
            return  [w for w in str (check_output (['bjobs ' + self.job_id], shell=True)).split (' ') if w != ''] [9] in ['RUN', 'PEND','SSUSP']
        except IndexError:
            return False
        except OSError:
            return True


def wait_for_jobs(job_list):
    '''Takes a list of jobs submitted to the queue and blocks the process untill all jobs in the list are done (or failed)'''
    while any([j.job_status_running() for j in job_list]):
        sleep(60)
        gc.collect()
    sleep(60)


def count_jobs(job_list):
    while Counter([j.job_status_running() for j in job_list])[True] > 100:
        sleep(300)
        gc.collect()
