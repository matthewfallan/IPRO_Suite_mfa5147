import sys
sys.path.append("/storage/home/mfa5147/work/IPRO_Suite/modules/")

from CHARMM import Relaxation
from MOLECULES import *

pdb_file = sys.argv[1]

os.chdir("/storage/home/mfa5147/work/IPRO_Suite/test")
mols = [MoleculeFile(pdb_file)[0]]
bindingAssembly = DesignGroup(1, mols)

experiment = {}
experiment["User"] = "mfa5147"
experiment["CHARMM Topology Files"] = ["expanded_topology.rtf"]
experiment["CHARMM Parameter Files"] = ["expanded_parameter.prm"]
experiment["Use Solvation"] = False
experiment["CHARMM Iterations"] = 100

Relaxation(bindingAssembly, experiment)
