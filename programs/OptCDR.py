#!/usr/bin/env python

# Where the IPRO Suite is installed
InstallFolder = "/storage/work/mfa5147/IPRO_Suite/"

# Name of this file
__name__ = "IPRO Suite Program to De Novo Design the CDRs for an Antigen"

# Documentation string
__doc__ = """
Written in 2014 by Matthew Grisewood of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This program de novo designs the binding portions (CDRs) of antibodies to have 
high specificity and affinity against any targeted epitope for an antigen."""

# Load the needed PYTHON Modules
import os
import sys
sys.path.append(InstallFolder + "modules/")
import EXPERIMENT
import IPRO_FUNCTIONS
import OPTCDR

# Start by loading an Experiment
experiment = EXPERIMENT.Experiment()
# Initialize the libraries for the OptCDR experiment
dummy = IPRO_FUNCTIONS.INITIALIZE(experiment)

# Go through each of the libraries
for ln in range(1, experiment["Optcdr Libraries"] + 1):
    # Make the initial library for the random antigen position
    OPTCDR.DO(experiment, ln)
