#/usr/bin/env python

# The name of this file
__name__ = "General OptGraft Functions"
# Documentation
__doc__ = """
Written in 2015 by Matthew Grisewood  of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions for running OptGraft. Each of the main functions
requires an Experiment class object as its input."""

# Include the PYTHON modules
import os
# Include the IPRO Suite STANDARDS module
from STANDARDS import *
# And include other IPRO Suite modules


# An error class for problems specific to this module
class OptgraftError(IPRO_Error):
    """An error for the problems that come up in the OPTGRAFT module"""
    def __init__(self, error = ''):
        """The initialization of the OptgraftError class."""
        IPRO_Error.__init__(self, error)

def create_binding_locations():
    """mjg5185 edit"""
    pass

def make_IPRO_experiment(experiment, folder):
    """Create a new Experiment class object to handle additional mutations
    nearby the binding site."""
    # Change the folder
    os.chdir(folder)
    # Copy the "input_files" folder
    os.system("cp -r " + experiment["Folder"] + "input_files/ .")
    # Make the "structures" folder- mjg5185 edit
    # Make a "results" folder
    os.mkdir("results")
    # Export the 'Experiment_Details.txt" file to this folder
    experiment.output(local = True)
    # Read through the file and make the required adjustments
    text = ""
    # Keep flags to accomplish various formatting tasks
    # mjg5185 Edits (refer to OptCDR.py) once IO_ASK & IO functions are
    # completed
