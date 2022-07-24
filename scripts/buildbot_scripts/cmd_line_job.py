import subprocess
import sys

from batch_job import BatchJob

class CmdLineJob(BatchJob):
    def __init__(self, cmd, cwd):
        super().__init__(cmd, cwd)
        self.system = "Commandline"

    # Command line jobs cannot be submitted in parallel - wait at submission
    def wait(self):
        print("Command line jobs run sequentially. Continue.")

    def submit(self, script):
        if len(self.parents) > 0:
            print("Dependencies are not supported for {}-jobs".format(self.system))
            sys.exit(1)

        # store output as LOG-file
        full_cmd = "./{}".format(script)
        with open("{}/LOG.{}.o".format(self.cwd, script), "wb") as out:
            self.job = subprocess.Popen(full_cmd, shell=False, stdout=out, stderr=out, cwd=self.cwd, encoding="UTF-8")
            # wait for job to finish before starting next one
            self.returncode = self.job.wait()

