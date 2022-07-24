import subprocess
import re
import sys

from batch_job import BatchJob

class SlurmJob(BatchJob):
    def __init__(self, cmd, cwd):
        super().__init__(cmd, cwd)
        self.system = "Slurm"

    def submit(self, script):
        # --wait so that the subprocess waits for slurm completion
        submit_cmd = self.cmd.split() + ["--wait"]

        if len(self.parents) > 0:
            parent_ids = [p.jobid for p in self.parents]
            submit_cmd.append("--dependency=afterany:{}".format(",".join(parent_ids)))

        submit_cmd.append(script)
        print("submitting slurm job: '{}'".format(" ".join(submit_cmd)))
        sp = subprocess.Popen(submit_cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.cwd, encoding="UTF-8")
        try:
            self.jobid = re.findall(r"\d{5,100}", sp.stdout.readline())[0]
        except:
            for line in sp.stderr.readlines():
                print(line)

        if not self.jobid or not self.jobid.isnumeric():
            print("Parsing jobid from slurm job failed, got {}".format(self.jobid))
            sys.exit(1)
        self.job = sp
