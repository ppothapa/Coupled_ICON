import subprocess
import re
import sys

from batch_job import BatchJob

debugOutput = False

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
        print(f"submitting PBS job: '{submit_cmd}'")
        qsub = subprocess.Popen(submit_cmd,
                                shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                cwd=self.cwd,
                                encoding="UTF-8")

        try:
            stdout = qsub.stdout.readlines()[0]
            if debugOutput:
                print(f"|qsub: stdout = {stdout}|")
            self.jobid = re.findall(r"(\d+)\.\w+", stdout)[0]
            if debugOutput:
                print(f"|qsub: jobId = {self.jobid}|")
        except:
            for line in qsub.stderr.readlines():
                print(f"|qsub STDERR: {line}|")

        if not self.jobid:
            print(f"Parsing jobid from pbs job failed, got {self.jobid}")
            sys.exit(1)

        # qsub will return immediately. That's why we pass along qwait
        qwaitCmd = 'qwait {}'.format(self.jobid)
        if debugOutput:
            print('|qwait call: "{}"|'.format(qwaitCmd))
        self.job = subprocess.Popen(qwaitCmd,
                                    shell=True,
                                    stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    encoding="utf-8")

    def wait(self):
        # PBS does not return the exitcode of the submitted script. instead it prints it to stdout
        ret = self.job.communicate()
        stdout = ret[0]
        stderr = ret[1]

        if debugOutput:
            print(f"|qwait:stdout: {stdout}|")
            print(f"|qwait:stderr: {stderr}|")

        try:
            exit_code = int(stderr.split()[-1])
            if debugOutput:
                print(f"qwait: exit_code = |{exit_code}|")
        except:
            print('pbs_job.py: Could not get a proper return value from "qwait"')
            print(f"qwait: stderr = |{stderr}|")
            exit_code = 1

        self.returncode = exit_code

        return exit_code
