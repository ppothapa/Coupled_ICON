import subprocess
import re
import sys

from batch_job import BatchJob

class PBSJob(BatchJob):
    def __init__(self, cmd, cwd):
        super().__init__(cmd, cwd)
        self.system = "PBS"

    def submit(self, script):
        submit_cmd = self.cmd.split()

        if len(self.parents) > 0:
            parent_ids = [p.jobid for p in self.parents]
            submit_cmd += ["--after", ",".join(parent_ids)]

        submit_cmd.append(script)
        print("submitting PBS job: '{}'".format(" ".join(submit_cmd)))
        qsub = subprocess.Popen(submit_cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.cwd, encoding="UTF-8")
        
        try:
            self.jobid = re.findall(r"\d{5,100}.\w{3,100}", qsub.stdout.readline())[0]
        except:
            for line in qsub.stderr.readlines():
                print(line)

        if not self.jobid:
            print("Parsing jobid from pbs job failed, got {}".format(self.jobid))
            sys.exit(1)

        # qsub will return immediately. That's why we pass along qwait
        self.job = subprocess.Popen(["qwait", self.jobid])
