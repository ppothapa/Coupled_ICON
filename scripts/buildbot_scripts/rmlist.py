#!/usr/bin/python3.2
# -*- coding: utf-8 -*-
#==============================================================================
# delete an experiment lists
#==============================================================================
from model_paths import *
import argparse

parser = argparse.ArgumentParser(description='Delete an experiment list.')
parser.add_argument('name', type=str, help='the name of the list')
args = parser.parse_args()

#print(args.name)
name=args.name
paths = model_paths()

paths.deleteThisExperimentList(name)
print("Experiment list "+name+" is deleted.")

