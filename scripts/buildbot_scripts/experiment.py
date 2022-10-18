from util import config_dict_to_string, config_dict_to_list
from icon_paths import run_path, base_path

from pathlib import Path
import subprocess
import sys

class Experiment(object):
    def __init__(self, name, builder, run_flags=None):
        self.name = name
        self.run_flags = run_flags

        self.builder = builder

        # The runscript will be generated with the filename. This means
        # self.name needs to be transformed: "path/to/filename" -> "filename"
        # by pathlibs "name" functionality.
        self.run_name = Path(self.name).name if Path(self.name).suffix == ".run" else "{}.run".format(Path(self.name).name)

        self.parents = []
        self.children = []

        self.batch_job = None

        return

    def to_string(self):
        out = self.name
        if(self.run_flags): out += " run_flags: {}".format(config_dict_to_string(self.run_flags))
        if(len(self.parents) > 0): out += " depends: {}".format(" ".join(["{} ({})".format(p.name, p.builder.name) for p in self.parents]))
        return out

    def add_child(self, child):
        print("adding child {} ({}) to experiment {} ({})".format(child.name, child.builder.name, self.name, self.builder.name))
        self.children.append(child)

    def add_parent(self, parent):
        print("adding parent {} ({}) to experiment {} ({})".format(parent.name, parent.builder.name, self.name, self.builder.name))
        self.parents.append(parent)

    def get_run_name(self, relative=True):
        return self.run_name if relative else run_path / self.run_name

    def submit(self, from_builder, parent_job=None):
        # if no job is submitted, return empty process list
        processes = []

        # the type of batch system must be specified by a helper class externally
        if not self.batch_job:
            print("{}: no information on the batch system given, please set 'batch_job' before submitting a job".format(self.run_name))
            sys.exit(1)

        # add the submitting parent job as dependency
        if parent_job:
            print("setting {} as parent job for {}".format(parent_job.jobid, self.run_name))
            self.batch_job.add_parent(parent_job)

        # make sure the job is not already running
        if self.batch_job.jobid:
            print("{}: already submitted with jobid {}".format(self.run_name, self.batch_job.jobid))
            return processes

        # TODO: implement cross-builder dependencies
        if from_builder and (self.builder.name != from_builder.name):
            print("cross builder dependency from {} to {} for experiment {}".format(self.builder.name, from_builder.name, self.name))
            print("cross builder dependencies are not yet allowed")
            sys.exit(1)

        # if this experiment has no dependencies or all parent jobs are collected, submit
        if len(self.parents) == 0 or (len(self.batch_job.parents) == len(self.parents)):
            p = self.batch_job.submit(self.run_name)
            processes += [p]
            for c in self.children:
                p = c.submit(self.builder, self.batch_job)
                if p: processes += p

        # waiting for all parent jobs to be submitted
        else:
            print("{}: waiting for {} parent jobs to be submitted...".format(self.name, len(self.parents) - len(self.batch_job.parents)))
            print("dependencies: {}".format(" ".join([exp.name for exp in self.parents])))

        return processes

    def make_runscript(self):
        # full path including experiment name
        exp_path = run_path / Path(self.name)

        # check if path exists. Note: this is already checked when the experiment is added
        if not exp_path.exists():
            print("could not find experiment {}".format(self.name))
            return 1

        # check if this is a hardcoded experiment
        if exp_path.suffix == ".run":
            # if experiment is located in the run directory, do nothing. Else add symlink
            if not exp_path.parent.absolute() == run_path:
                new_path = run_path / exp_path.name
                # remove symlink if already exists
                if new_path.exists(): new_path.unlink()
                # new (absolute) path in the run directory
                new_path.symlink_to(exp_path)

                print("linking {} to {}".format(new_path.relative_to(base_path), exp_path.relative_to(base_path)))

            # return relative path to runscript
            status = 0
        else:
            # get filename -> get last element of filename (i.e. check.atm_amip -> atm_amip)
            exp_name = Path(self.name).name.split(".")[-1]
            # check if this script will run icon
            cmd = ["./run/make_target_runscript"]
            cmd.append("in_script={}".format(self.name))
            if Path(self.name).name.split(".")[0] in ["check", "exp"]:
                cmd.append("in_script=exec.iconrun")
            cmd.append("out_script={}".format(self.run_name))
            cmd.append("EXPNAME={}".format(exp_name))
            cmd += config_dict_to_list(self.run_flags) # append the list of run flags
            sp = subprocess.run(cmd, shell=False, cwd=base_path, stderr=subprocess.STDOUT, encoding="UTF-8")
            status = sp.returncode
        
        return status

