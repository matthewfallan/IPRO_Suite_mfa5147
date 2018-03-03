import os
import shutil
import sys
sys.path.append("/storage/home/mfa5147/work/IPRO_Suite/modules/")

import REFINEMENT
import EXPERIMENT

directory = sys.argv[1]
os.chdir(directory)
if not os.path.isdir("refinement"):
    os.mkdir("refinement")

experiment = EXPERIMENT.Experiment()

# Remove directories left by any previous refinements aborted before proper completion.
# These directories would otherwise cause the new refinement to fail.
refinement = os.path.join(experiment["Folder"], "refinement")
if os.path.isdir(refinement):
    shutil.rmtree(refinement)
sharing = os.path.join(experiment["Folder"], "sharing")
if os.path.isdir(sharing):
    shutil.rmtree(sharing)

iteration = 0
REFINEMENT.Start(experiment, iteration)
REFINEMENT.DO(experiment)
