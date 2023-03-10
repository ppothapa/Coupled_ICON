class BatchJob(object):
    def __init__(self, cmd, cwd):
        self.system = "undefined"
        self.cmd = cmd
        self.cwd = cwd
        self.jobid = None
        self.job = None
        self.returncode = None
        self.parents = []

    def add_parent(self, parent):
        self.parents.append(parent)

    def wait(self):
        self.returncode = self.job.wait()

    def cancel(self):
        self.job.cancel()
