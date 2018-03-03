#/usr/bin/env python

# The name of this file
__name__ = "General OptCDR Functions"
# Documentation
__doc__ = """
Written in 2015 by Matthew Grisewood of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions for running OptCDR. Each of the main functions
requires an Experiment class object as its input."""

# Include the PYTHON modules
import os
import random
import sys
import time
# Include the IPRO Suite STANDARDS module
from STANDARDS import *
# And include other IPRO Suite modules
import CHARMM
import EXPERIMENT
import GAMS
#import CPLEX
import IO_OUTPUT
import IPRO_FUNCTIONS
import MOLECULES
import REFINEMENT
import ROTAMERS
import SHARING

# An error class for problems specific to this module
class OptcdrError(IPRO_Error):
    """An error for problems in the OPTCDR module"""
    def __init__(self, error = ''):
        """The initialization of the OptcdrError class."""
        IPRO_Error.__init__(self, error)

def load_clusters(experiment, cdrs = ["H1", "H2", "H3", "L1", "L2", "L3"]):
    """This function loads the information pertaining to the OptCDR canonical
    clusters."""
    # Store the information in this dictionary
    clusters = {}
    # Access the files with all of the clustering information
    files = os.listdir(experiment['Optcdr Cluster Folder'])
    # Go through each of the files
    for file in files:
        # Skip any swap files
        if file[-4:] == ".swp":
            continue
        # Read through each file
        f = open(experiment['Optcdr Cluster Folder'] + file)
        for line in f:
            # Split the line on white space
            items = line.split()
            # If it is basic information about the cluster, store it
            if len(items) > 0 and items[0] in cdrs:
                # CDR Cluster INT   Length: INT   Model: NAME   Members: INT
                # 0   1       2     3       4     5      6      7        8
                cdr = items[0]
                clusterNo = int(items[2])
                clusterLength = int(items[4])
                modelName = items[6]
                members = int(items[8])
                if cdr not in clusters.keys():
                    clusters[cdr] = {}
                # Store the information for this cluster number
                clusters[cdr][clusterNo] = {"Length": clusterLength, "Model": \
                                            modelName, "Members": members}
        f.close()
    # Store the clustering information
    return clusters

def molecule_name_association(experiment, design = []):
    """This function associates a Molecule name with its corresponding CDR
    name."""
    # Create a dictionary for the associations
    associations = {}
    # Gather the chains being designed in this experiment
    chains = experiment["Optcdr Chains"]
    # Sort them
    chains.sort()
    # Create a list of CDRs
    cdrs = []
    # Go through the chains
    for chain in chains:
        # Each chain has 3 CDRs
        for i in range(1, 4):
            cdrs.append(chain[0].upper() + str(i))
    # If this is an affinity maturation, Design Molecules are handled a little
    # differently
    if design != [] and len(design) != len(cdrs):
        for mol in experiment[0]:
            if mol.design:
                design.append(mol.name)
    else:
        for mol in experiment["Molecules"]:
            if mol[0] == None:
                design.append(mol[1])
    design.sort()
    # Create a temporary list to manage situations where either (or both) of the
    # frameworks are provided
    temp = []
    # If a heavy framework is provided
    if ("heavy" in experiment["Optcdr Chains"]) and \
                                ("heavy" in experiment["Optcdr Frameworks"]):
        temp.extend([design[0]] * 3)                              
    # If a light framework is provided
    if ("light" in experiment["Optcdr Chains"]) and \
                                ("light" in experiment["Optcdr Frameworks"]):
        temp.extend([design[-1]] * 3)
    # If necessary, overwrite the design list
    if temp != []:
        design = temp
    # The length of the two lists must match
    if len(design) != len(cdrs):
        text = "There is an unexpected number of Design Molecules in the "
        text += "experiment."
        raise OptcdrError(text)
    # Associate the two pieces of information
    for i, j in zip(design, cdrs):
        associations[j] = i
    return associations

def optcdr_usage_limits(rotamers, experiment):
    """Obtain the amino acid usage information for use during rotamer
    optimization."""
    # Load the CDR-Molecule name associations to determine which CDR is free
    associations = molecule_name_association(experiment, [])
    # Store the residue numbers that are free to move
    numbers = {}
    # Go through the Binding Assemblies
    for group in experiment:
        for mol in group:
            # Go though the residues
            for res in mol:
                # If the residue is free to move, store its information
                if res.freedom == "FREE":
                    # Get the residue's number
                    name = res.name
                    # If the name is composed of only digits
                    if name.isdigit():
                        pass
                    # If the last character is a letter
                    elif name[:-1].isdigit() and name[-1].isalpha():
                        name = name[:-1] 
                    # If this is the first instance of the molecule, create a
                    # dictionary key
                    if mol.name not in numbers.keys():
                        numbers[mol.name] = []
                    # Store the residue name
                    numbers[mol.name].append(int(name))
    # There should only be one Molecule that is at least partially free
    number_keys = list(numbers.keys())
    if len(number_keys) != 1:
        text = "There is more than 1 Molecule that is at least partially free "
        text += "within an initial antibody refinement."
        raise OptcdrError(text)    
    # Determine the CDR that is free
    ranges = {1: range(27, 39), 2: range(56, 66), 3: range(105, 118)}
    # Go through the CDR ranges
    for cdrNo in ranges.keys():
        # Break the loop if the number falls within the range
        if numbers[number_keys[0]][0] in ranges[cdrNo]:
            break
    # Ensure that all free residues are in this range
    allInOne = True
    # Go throguh the residue numbers
    for number in numbers[number_keys[0]]:
        # If the number does not fall within the span, there is a problem
        if number not in ranges[cdrNo]:
            allInOne = False
            break
    # If all free residues do not fall within the same CDR, raise an error
    if not allInOne:
        text = str(number) + " is not in the range for CDR" + str(cdrNo)
        raise OptcdrError(text)
    # Determine if this is the light or heavy chain
    for chain in experiment["Optcdr Chains"]:
        # Create the hypothetical CDR 
        cdr = chain[0].upper() + str(cdrNo)
        # Use the "associations" dictionary to see if the Molecule in "numbers"
        # is used for this CDR. The Molecule names are not redundant so it can
        # only be used for one chain (but can span multiple CDRs)
        if associations[cdr] == number_keys[0]:
            break
    # Get the current folder name to determine which library this is
    ln = int(os.getcwd().split("library")[1])
    # Get the CDRs for this experiment
    cdrs = list(associations.keys())
    cdrs.sort()
    # Get the index for the CDR
    index = cdrs.index(cdr) + 1
    # Get the canonical structure from the solutions
    solution = experiment["Scores"][ln-1][1]
    canonical = solution[index]
    # Store the data in this dictionary
    information = {}
    # Open the file with the composition rules
    name = cdr + "_" + str(canonical) + "_Composition_Rules.txt"
    f = open(experiment['Optcdr Usage Pattern'] + name, "r")
    # Go through the file
    for line in f:
        # Split the line on white space
        items = line.split()
        # Store the data in the dictionary
        # The data is stored as information[category] = [min, max]
        information[items[0][:-1]] = [items[2], items[4]]    
    # Close the file
    f.close()
    # Gather the permission information
    name = cdr + "_" + str(canonical) + "_permittedKinds.txt"
    permissions = {}
    f = open(experiment['Optcdr Usage Pattern'] + name, "r")
    # Go thorugh the permissions file
    for line in f:
        # Split the line on white space
        items = line.split()
        # Store the permissions for the CDR
        # Format is permissions[residueName] = [AA1, AA2, ..., AAN]
        permissions[items[0]] = items[1:]
    # Close the file
    f.close()
    # Go through the Binding Assemblies
    for group in experiment:
        # Go through the Molecules
        for mol in group:
            # Go though the Residues
            for res in mol:
                # Set the amino acid permissions
                # If this is the molecule with the CDR and the residueName is in
                # the CDR, set the permissions using the permissions dictionary
                if (mol.name == associations[cdr]) and (res.name in \
                                                        permissions.keys()):
                    # Allow this residue to be mutated
                    res.design = True
                    # Specify the amino acids allowed
                    res.permittedKinds = permissions[res.name]
                # Otherwise, it should be an empty list    
                else:       
                    res.design = False
                    res.permittedKinds = []
    return information    

def free_cdr(experiment, associations, cdr, maturation = None):
    """Edit the freedom attribute of each residue so that only the specified CDR
    is permitted to move. Also edit the permissions attribute so the amino acid
    can only change within the CDR being examined to a specified list of amino
    acids"""
    # Get the current folder name to determine which library this is
    ln = int(os.getcwd().split("library")[1].split("/")[0])
    # Get the CDRs for this experiment
    cdrs = list(associations.keys())
    cdrs.sort()
    # Get the index for the CDR
    index = cdrs.index(cdr) + 1
    # Get the canonical structure from the solutions
    solution = experiment["Scores"][ln-1][1]
    canonical = solution[index]
    # Gather the permission information
    name = cdr + "_" + str(canonical) + "_permittedKinds.txt"
    permissions = {}
    f = open(experiment['Optcdr Usage Pattern'] + name, "r")
    # Go thorugh the permissions file
    for line in f:
        # Split the line on white space
        items = line.split()
        # Store the permissions for the CDR
        # Format is permissions[residueName] = [AA1, AA2, ..., AAN]
        permissions[items[0]] = items[1:]
    # Close the file
    f.close()
    # Store the perturbation information in this dictionary 
    perturbed = []
    # If this is an affinity maturation, use the IPRO Experiment class object
    if maturation != None:
        exp = maturation
    # Otherwise, use the OptCDR Experiment class object      
    else:
        exp = experiment
    # Go through all of the Binding Assemblies
    for group in exp:
        # Go through the Molecules within the Binding Assembly
        for mol in group:
            # Go through each Residue of the Molecule
            for res in mol:
                # Get the residue's number
                name = res.name
                # If the name is composed of only digits
                if name.isdigit():
                    pass
                # If the last character is a letter
                elif name[:-1].isdigit() and name[-1].isalpha():
                    name = name[:-1]
                # Otherwise, the residue is not properly named so raise an error
                else:
                    text = name + " is not a proper residue name"
                    raise OptcdrError(text)
                # Convert the number to an integer
                name = int(name)
                # If the Residue is in a Target Molecule, fix it in place
                if not mol.design:
                    res.freedom = "FIXED"
                    res.permission = "FIXED"
                # If this CDR is not in the Molecule, fix everything
                elif associations[cdr] != mol.name:
                    res.freedom = "FIXED"
                    res.permission = "FIXED"
                # If the residue number is within the CDR, then allow it to move
                # and receive mutations
                elif ((cdr[1] == "1") and (name in range(27, 39))) or \
                     ((cdr[1] == "2") and (name in range(56, 66))) or \
                     ((cdr[1] == "3") and (name in range(105, 118))):
                    res.freedom = "FREE"
                    res.permission = "ROTAMER"
                    perturbed.append([mol.name, res.name])
                # Otherwise, fix it in place and do not allow amino acid changes
                else:
                    res.freedom = "FIXED"
                    res.permission = "FIXED"
                # Set the amino acid permissions
                # If this is the molecule with the CDR and the residueName is in
                # the CDR, set the permissions using the permissions dictionary
                if (mol.name == associations[cdr]) and (res.name in \
                                                        permissions.keys()):
                    # Allow this residue to be mutated
                    res.design = True
                    # Specify the amino acids allowed
                    res.permittedKinds = permissions[res.name]
                # Otherwise, it should be an empty list    
                else:
                    res.design = False
                    res.permittedKinds = []
    return perturbed

def initialize_antibodies(experiment, canonicals):
    """Refine the construction of the antibody and the position of the 
    antigen in the antibody's binding pocket"""
    # Store the order in which the CDRs should be refined
    order = ['H3', 'L3', 'H2', 'H1', 'L1', 'L2']
    # Get the molecule name-CDR associations
    associations = molecule_name_association(experiment, [])
    # Do the refinement protocol 2x for each CDR
    for i in range(2):
        # Go through the CDRs in order
        for cdr in order:
            # If this CDR is not included in the experiment, skip it
            if cdr not in canonicals.keys():
                continue
            # Only allow this CDR to move in the energy minimization
            dummy = free_cdr(experiment, associations, cdr)
            # Find the optimal set of rotamers with the additional
            # OptCDR-specific constraints
            IPRO_FUNCTIONS.Optimal_Rotamers(experiment)
        # Run a rigid-body docking protocol
        IPRO_FUNCTIONS.Docking(experiment, experiment["Docking Frequency"])

def orient_framework(molecule, reference):
    """Properly orient framework molecules so that CDRs can be attached to them
    within an OptCDR experiment."""
    # If the dimensions of the framework molecule and reference molecule are
    # different, raise an error
    if molecule.dimensions != reference.dimensions:
        text = "The dimensions of the framework reference and the framework "
        text += "Molecule class object must match."
        raise OptcdrError(text)
    # Identify the backbone atoms for the attachment points
    atoms1 = []
    atoms2 = []
    # Go through the residues in the framework Molecule
    for res in molecule:
        # If the residue name is one of the attachment point positions, store
        # the atoms
        if res.name.isdigit() and int(res.name) in [26, 39, 55, 66, 104, 118]:
            # Go through the backbone atoms
            for name in ['N', 'CA', 'C']:
                # Store the attachment point atom in the reference structure
                atoms1.append(reference[res.name][name])
                # Store the attachment point atom in the framework molecule
                atoms2.append(res[name])
    # Make sure there are corresponding lengths of atom lists
    if len(atoms1) != len(atoms2):
        text = "Orienting the framework molecule requires the same number of "
        text += "attachment point atoms for both the framework Molecule as "
        text += "well as the framework reference."
        raise OptcdrError(text)
    # Center the reference attachment point atoms
    ref_coors = [0.0] * reference.dimensions
    # Center the framework attachment point atoms
    frame_coors = [0.0] * molecule.dimensions
    # Go through the lists of atoms
    for ref_atom, frame_atom in zip(atoms1, atoms2):
        # Go through the dimensions of the atoms
        for dim in range(reference.dimensions):
            # Add the dimension to the list controlling the centering of the
            # molecules
            ref_coors[dim] += ref_atom[dim]
            frame_coors[dim] += frame_atom[dim]
    # Go back through the lists and divide by the total number of atoms (to get
    # the geometric center of the atoms)
    for dim in range(reference.dimensions):
        ref_coors[dim] /= len(atoms1)
        frame_coors[dim] /= len(atoms2)
    # Calculate the rotation matrix
    rmatrix = MOLECULES.calculate_rmatrix(atoms1, atoms2)
    # Center the framework molecule
    MOLECULES.move(molecule, frame_coors)
    # Rotate the framework molecule
    MOLECULES.rotate(molecule, rmatrix)
    # Move the framework molecule back to the reference's coordinates
    MOLECULES.move(molecule, ref_coors, "+")

def include_framework(frameMol):
    """Create the antibody structure including the framework regions"""
    # Store the framework regions in a list of residues
    residues = []
    # Go through the framework PDB residues
    for res in frameMol:
        # Extract the residue number
        name = res.name
        # If the name is composed of only digits
        if name.isdigit():
            n = int(name)
        # If the last character is a letter
        elif name[:-1].isdigit() and name[-1].isalpha():
            n = int(name[:-1])
        # If this position is within a CDR or not between 1 and 129, skip it
        if (n >= 27 and n <= 38) or (n >= 56 and n <= 65) or \
        (n >= 105 and n <= 117) or n < 1 or n > 129:
            continue
        # If this is the first residue of a framework region, initialize a
        # string to store the atoms
        if n in [1, 39, 66, 118]:
            # If this is not FR1, append the previous framework to the list
            if n != 1:
                residues.append(atoms)
            atoms = ""
        # Add this residue to the string of residues      
        atoms += str(res)        
    # Add FR4 to the list 
    residues.append(atoms)
    return residues

def build_antibodies(experiment, canonicals, ln):
    """This function constructs the antibody structures for a given OptCDR
    library."""
    # Store all of the antibodies in this list
    antibodies = []
    # Get the list of CDRs to be considered
    cdrs = list(canonicals.keys())
    cdrs.sort()
    # Get the optimal set of canonical structures for this library
    solution = experiment["Scores"][ln-1][1]
    # Go through the antibody chains being designed
    chains = experiment["Optcdr Chains"]
    chains.sort()
    # Find the reference framework molecule
    if experiment["Optcdr Frameworks"] != {}:
        file = experiment["Optcdr Framework Reference"].split("/")[-1]
        path = experiment["Optcdr Framework Reference"].replace(file, '')
        reference = MOLECULES.MoleculeFile(file, path) 
    for chain in chains:
        # Store the molecules in a list
        molecules = []
        # Go through the CDR numbers
        for i in range(1, 4):
            # Extract the CDR name
            cdr = chain[0].upper() + str(i)  
            # Get the index for the CDR in the solution dictionary
            index = cdrs.index(cdr) + 1
            # Append the canonical structure molecule to the list of molecules
            molecules.append(canonicals[cdr][solution[index]])
        # If a framework has been specified, add it
        name = chain.lower()
        if name in experiment["Optcdr Frameworks"]:
            # Extract the Molecule class object
            molecule = experiment["Optcdr Frameworks"][name][0][2]
            # Properly orient the framework
            orient_framework(molecule, reference[chain[0].upper()])
            # Obtain the list of framework regions
            frameworks = include_framework(molecule)
            # Go through the molecules and convert them to text
            texts = []
            for mol in molecules:
                # Skip the first and last residues of each CDR since they are
                # attach points
                text = ""
                for rn in range(len(mol)):
                    if rn not in [0, len(mol) - 1]:
                        text += str(mol[rn])
                # Append this CDR to the list of CDRs           
                texts.append(text)
            # Make sure there are 4 framework regions and 3 CDRs
            if len(frameworks) != 4 or len(texts) != 3:
                text = "The framework molecule does not include all of the "
                text += "necessary regions"
                raise OptcdrError(text)
            # Concatenate the text into a single string
            atoms = ''
            # Add FR1, CDR1, FR2, CDR2, FR3, CDR3 in that order 
            for i in range(0, 3):
                atoms += frameworks[i] + texts[i] 
            # Add FR4      
            atoms += frameworks[3]
            # Overwrite the molecules list with this single molecule
            molecules =[MOLECULES.Molecule(atoms, 1, 1, chain[0].upper(), True,\
                         experiment["Force Field"], experiment["File Format"])]
        # Add the list of molecules to the list of antibodies
        antibodies.extend(molecules)
    # Generate the expected number of Design Molecules
    dms = 0
    for chain in chains:
        if chain in experiment["Optcdr Frameworks"]:
            dms += 1
        else:
            dms += 3
    # Create a formatted list
    formatted_antibodies = []
    # If a single molecule is present for each chain, use the framework name
    if len(antibodies) == len(chains):
        for antibody in antibodies:
            formatted_antibodies.append([None, antibody.name, antibody])
    # If there are 3 CDRs for each chain, generate a unique name for each CDR
    elif len(antibodies) == dms:  
        # Create a formatted list
        formatted_antibodies = []
        # Create a string of possible names
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        # Use a counter to go through the alphabet
        an = 0
        # Go through the antibodies
        for antibody in antibodies:
            # Use a while loop to find an appropriate name for the molecule
            goOn = True
            while goOn:
                # Increment the counter
                an += 1
                # If the name has not been used, store it
                if alphabet[an] not in experiment[0]:
                    goOn = False
                    # Store the formatted molecule list
                    antibody.name = alphabet[an]
                    formatted_antibodies.append([None, alphabet[an], antibody])
    # Otherwise, something is wrong so raise an error
    else:
        text = "There is an unexpected number of antibodies generated for "
        text += "library " + str(ln)
        raise OptcdrError(text)
    # Update the experiment to include the antibody details
    experiment["Molecules"].extend(formatted_antibodies)
    experiment.make_DesignGroups()
    existing = experiment['Summary']
    experiment.finish_creation()
    experiment['Summary'] += existing
    # Update the OptCDR Scores
    SHARING.load_scores(experiment, experiment["Folder"] + \
                               "results/initial/")
    # Output these molecules to the 'Current' folder
    SHARING.output_Current(experiment, './Current/')
    # Refine the initial antibody structure
    initialize_antibodies(experiment, canonicals)

def optcdr_restraints(experiment):
    """Automatically create restraints that fix framework regions in place and 
    restrict movement in the CDRs"""
    # Extract any previously-created restraints
    if "Restraints" in experiment:
        restraints = experiment["Restraints"]
    # If there are not any previously-created restraints, begin with an empty
    # dictionary
    else:
        restraints = {}
    # Avoid having multiple position constraints on the Design Molecules, which
    # were created in the OptCDR experiment. Atoms in Design Molecules could not
    # have been previously fixed, had distance restraints, or dihedral
    # restraints since the Design Molecules did not previosuly exist
    positions = []
    # Go through the position restraints
    if "Position" in restraints:
        # Go through the restraints
        for restraint in restraints["Position"]:
            # If the restraint is on all Molecules, edit it to just be Target 
            # Molecules
            if restraint[1] == "all":
                # Go through the Target Molecules
                for molecule in experiment[0]:
                    if not molecule.design:
                        # Add a restraint
                        positions.append([restraint[0], molecule.name,
                        restraint[2], restraint[3], restraint[4]])
            # Otherwise, add the restraint
            else:
                positions.append(restraint)
    # Create a list to keep track of the CDR positions
    cdr_positions = range(27, 39)
    cdr_positions.extend(range(56, 66))
    cdr_positions.extend(range(105, 118))  
    # If 'Fixed Atoms' is not already in the restraints dictionary, then add it
    if "Fixed Atoms" not in restraints.keys():
        restraints["Fixed Atoms"] = {"all": {}}
    # If there is not already a key within the 'Fixed Atoms' dictionary for
    # 'all' Binding Assemblies, add it
    if "all" not in restraints["Fixed Atoms"].keys():
        restraints["Fixed Atoms"]["all"] = {}
    # Go through the Design Molecules
    for molecule in experiment[0]:
        if molecule.design:
            # Go through each residue in the Molecule
            for res in molecule:
                # Get the residue's position
                name = res.name
                # If the name is composed of only digits
                if name.isdigit():
                    pass
                # If the last character is a letter
                elif name[:-1].isdigit() and name[-1].isalpha():
                    name = name[:-1] 
                # If this position is within a CDR, add a position restraint
                if int(name) in cdr_positions:
                    # The position restraint should occur on all atoms with a
                    # force constant of 0.05
                    positions.append(["all", molecule.name, res.name, "all", 
                    0.05])
                # Otherwise, fix the atoms of this residue in place
                else:
                    # Store the atoms in a list
                    atoms = []
                    for atom in res:
                        atoms.append(atom.name)
                    # If this Molecule is not in the "Fixed Atoms" dictionary,
                    # add it
                    if molecule.name not in \
                    restraints["Fixed Atoms"]["all"].keys():
                        restraints["Fixed Atoms"]["all"][molecule.name] = {}
                    restraints["Fixed Atoms"]["all"][molecule.name][res.name]=\
                    atoms                                                
    # Add position restraints              
    restraints["Position"] = positions 
    # Create a temporary dictionary to output the restraints
    exp = {"Restraints": restraints, "Design Positions": []}
    for name in ["Molecules", "Design Groups", "Force Field", "File Format"]:
        exp[name] = experiment[name]
    # Generate the text for the Restraints     
    text = IO_OUTPUT.Restraints(exp)
    return text

def make_IPRO_experiment(experiment, folder):
    """Create a new Experiment class object to handle the affinity maturation"""
    # Change the folder
    os.chdir(folder)
    # Make an "input_files" folder
    os.system("cp -r " + experiment["Folder"] + "input_files/ ." )
    # Make a "structures" folder
    os.system("cp -r " + experiment["Folder"] + "structures/ ." )
    # Make a "results" folder
    os.mkdir("results")
    # Copy the "Current" folder from the library initialization to this newly
    # created folder
    command = "cp -r " + experiment["Folder"] + "initial_" + folder  
    command += "/Current/ ./temp_Current"
    os.system(command)
    # Rename the Design Molecules
    for molecule in experiment["Molecules"]:
        if molecule[0] == None:
            new = "Molecule" + molecule[2].name + ".pdb"
            molecule[0] = new
    # Export the 'Experiment_Details.txt' file to this folder
    experiment.output(local = True)
    # Read through the file and make required adjustments
    text = ""
    # Keep flags to accomplish various formatting tasks
    # Skip the blank line after the OptCDR section
    skipLine = False
    # Keep track if the Design Molecules were specified
    dms = False
    # Keep track if the Design Positions were specified
    dps = False
    # Keep track if the restraints were specified
    rest = False
    with open('Experiment_Details.txt', 'r') as file:
        for line in file:
            # Get the index of the first ':'. 
            try:
                i = line.index(':')
            # If there isn't a ":", keep the line (except OptCDR information)
            # for formatting purposes       
            except ValueError:
                # Skip the OptCDR title
                if line.strip() == "How to run OptCDR":
                    skipLine = True
                # Skip antigen lines
                elif line.strip() == "Antigen epitope residues":
                    skipLine = True
                # Skip any Restraints lines
                elif line.strip() in ["Atoms that may never move", \
                                      "Restraints on Atom Positions", \
                                      "Restraints on Atom-Atom distances", \
                                      "Restraints on Dihedral Angles"]:
                    skipLine = True                           
                # Also skip the next blank space following the OptCDR
                # information
                elif line.strip() == "" and skipLine:
                    skipLine = False
                # If the Design Positions should be specified, list them
                elif dps:
                    # Only add the Design Position info once
                    dps = False
                    # Add this line
                    text += line
                    # Add the proper title
                    text += "Residues that are Permitted to Mutate\n"
                    # Store the left-hand side text
                    lhs = "Design Position: ".ljust(30)
                    # Get the CDR positions
                    positions = range(27, 39)
                    positions.extend(range(56, 66))
                    positions.extend(range(105, 118))
                    # Go through the Molecules
                    for mol in experiment[0]:
                        # Only add text for Design Molecules
                        if mol.design:
                            # Go through the residues in the Design Molecule
                            for res in mol:
                                # Get the residue's position
                                name = res.name
                                # If the name is composed of only digits
                                if name.isdigit():
                                    pass
                                # If the last character is a letter
                                elif name[:-1].isdigit() and name[-1].isalpha():
                                    name = name[:-1]
                                if int(name) in positions:
                                    text += lhs + "Residue " + res.name 
                                    text += " in Molecule " + mol.name + "\n"
                    # Add an extra line
                    text += "\n"
                # Otherwise, add the line      
                else:
                    text += line
                continue
            # Split the line on that value
            attribute = line[:i].strip()
            info = line[i+1:].strip() 
            lhs = attribute + ": "
            lhs = lhs.ljust(30)
            # If this is the type of experiment
            if attribute == "Type":
                text += lhs + "IPRO\n"
            # If this is the name of the experiment
            elif attribute == "Name":
                text += lhs + folder + "\n"
            # If this is the path to the experiment's directory
            elif attribute == "Folder":
                text += lhs + info + folder + "/\n"
            # Add Design Molecules
            elif attribute == "Molecule":
                # If the text for the Design Molecules has not been added
                if not dms:
                    # Move to the "structures" folder
                    os.chdir("structures")
                    # Go through the Molecules
                    for molecule in experiment["Molecules"]:
                        # Extract the Molecule
                        mol = molecule[2]
                        # Only format Design Molecules
                        if mol.design:
                            # Add the text for the Design Molecule
                            text += lhs + "Molecule " + mol.name + " from file"
                            text += " Molecule" + mol.name + ".pdb is "
                            text += "Design Molecule " + mol.name + "\n"
                            # Output the Molecule
                            mol.output(None, experiment["File Format"], \
                                                            experiment["User"])
                    # Move back to the proper folder
                    os.chdir("../")
                    # Make it known that the Design Molecules have been
                    # specified
                    dms = True
                # Add the Target Molecule(s)
                if info.split()[6] == "Target":
                    text += line
            # Add Design Positions
            elif attribute == "Design Group":
                # Edit the flag so that Design Positions are added next
                dps = True
                # Add the Binding Assembly line
                text += line
            # If the line contains restraint information, skip it
            elif attribute in ["Fixed Atoms", "Position Restraint", \
                                "Distance Restraint", "Dihedral Restraint"]:
                pass                      
            # Skip any OptCDR-specific information
            elif attribute in ["Home Directory", "Canonical Folder", 
                               "Clash File", "Position File", "Cluster Folder",
                               "Framework Reference", "Optcdr Chain", 
                               "Optcdr Positions", "Optcdr Libraries",
                               "Antigen Rotation", "Heavy Framework",
                               "Light Framework", "Epitope Position",
                               "Usage Pattern"]:
                pass                        
            # Otherwise, just keep the same
            else:
                text += line
    # Add the text from the restraints at the end of the file
    text += optcdr_restraints(experiment)
    # Use the stored text to write a new 'Experiment_Details.txt' file
    with open('Experiment_Details.txt', 'w') as file:
        file.write(text) 
    # Return to the OptCDR home directory
    os.chdir("../")

def initialize_libraries(experiment, ln):
    """Create copies of the canonical structures and properly oriented antigens
    in each library."""
    # Move into the folder to do the intial calculations in
    folder = "initial_library" + str(ln)
    os.chdir(folder)    
    # Create a time stamp for beginning the calculations
    experiment["Summary"] = "Library " + str(ln)  + " Initialization\n"
    experiment["Summary"] += "Started" + SHARING.time_stamp()
    # Find the proper number of coordinates to consider
    N = len(experiment["Movements"][ln])/2
    # Go through each antigen
    for mol in experiment[0]:
        # Apply the proper rotation
        for cn in range(N):
            # Create a generic vector of zeros of the appropriate length
            vector = [0.0] * N
            # Place a value of 1.0 in the correct location in the vector
            vector[cn] = 1.0
            # Find the angle to rotate the antigens by
            angle = experiment["Movements"][ln][N+cn]
            # Rotate each of the antigens by the appropriate angle
            rmatrix = MOLECULES.calculate_rmatrix(angle, vector)
            MOLECULES.rotate(mol, rmatrix)
        # Translate each antigen by the appropriate amount
        MOLECULES.move(mol, experiment["Movements"][ln][:N], '+')
    # Update the reference folder with these updated coordinates
    SHARING.output_Current(experiment, "./Current/") 
    # Load the canonical structures
    canonicals = IPRO_FUNCTIONS.load_canonicals(experiment)
    cdrs = list(canonicals.keys())
    cdrs.sort()
    # Load the clashes
    clashes = IPRO_FUNCTIONS.load_clashes(experiment, cdrs) 
    # Load the C++ scores
    raw_scores = IPRO_FUNCTIONS.load_scores(experiment["Folder"])
    # Look for alternate solutions using integer cuts
    goOn = True
    # Store the solutions in a list
    solutions = [experiment["Scores"][ln-1]]
    # Keep searching for alternate solutions until the quality of the result is
    # worse
    while goOn:
        # Resolve the MILP using integer cuts
        if useCPLEX:
            #solution = CPLEX.optcdr_canonicals(canonicals, clashes, \
            #                                 raw_scores[ln], solutions)
            pass
        else:
            solution = GAMS.optcdr_canonicals(canonicals, clashes, \
                                              raw_scores[ln], solutions)
        # If the solution found has an equal objective value to the first, store
        # it and re-run the MILP
        if solution["Score"] == experiment["Scores"][ln-1][1]["Score"]:
            solutions.append([experiment["Scores"][ln-1][0], solution])
        # Otherwise, break out of the loop and analyze the results
        else:
            goOn = False
    # Update the library based on the most members for the cluster
    best = 0
    # Skip this if there is only one solution after applying the integer cuts
    if len(solutions) > 1:
        # Load the clusters
        cdrs = list(canonicals.keys())
        cdrs.sort()
        clusters = load_clusters(experiment, cdrs)
        # Initialize the variables to store the solution with the most cluster
        # members
        best = None
        amount = 0
        # Go through the solutions
        for i, solution in enumerate(solutions):
            # Store the total number of members throughout the CDRs
            total = 0
            # Go through the CDRs
            for j, cdr in enumerate(cdrs):
                # Extract the number of members from the "clusters" dictionary 
                members = clusters[cdr][solution[1][j+1]]["Members"]
                # 30 is the number where the permitted amino acids change from
                # "of the same type" to "only those observed" at each position
                if members > 30:
                    members = 30
                # Add the number of members to the total for this solution
                total += members
            # If applicable, update the "best" solution found and its
            # corresponding total number of members
            if total > amount:
                best = i
                amount = total
    # Update the library based on the most structures
    experiment["Scores"][ln-1] = solutions[best]
    # If the set of canonical structures has changed, update the referenced
    # values
    if best != 0:
        SHARING.output_scores(experiment, experiment["Folder"] + "Current/", ln)
    # Copy the necessary files
    SHARING.copy_standard_files(experiment, solv = True) 
    # Generate the antibody structures
    build_antibodies(experiment, canonicals, ln) 
    # Go back to the home directory
    os.chdir("../")
    # Try to create a new folder to handle the IPRO affinity maturation
    folder = "library" + str(ln)
    try:
        os.mkdir(folder)
    # If the folder already exists, delete it and make a new one. This is the
    # proper procedure since the library should only be there if the
    # initialization has already finished
    except OSError:
        os.system("rm -rf " + folder)
        os.mkdir(folder)
    # Create a new Experiment class object to handle the IPRO affinity maturation
    make_IPRO_experiment(experiment, folder)
    # Delete the initialization folder
    os.system("rm -rf initial_" + folder) 
    # Update the summary file
    # Create a summary file
    experiment["Summary"] += "Ended" + SHARING.time_stamp()
    name = SHARING.summary_name(SHARING.get_current())
    f = open(name, "a")
    f.write(experiment["Summary"])
    f.close()
 
def initial_check(experiment, ln):
    """Do an initial check about doing calculations for a library"""
    # Start sharing
    SHARING.Start(experiment)
    # This is the name of the folder that we are searching for
    folder = "library" + str(ln)
    # If the results are already completed
    if folder in os.listdir(experiment["Folder"] + "results/"):
        SHARING.End(experiment)
        return False, False
    # Or if the affinity maturation is ongoing
    elif folder in os.listdir(experiment["Folder"]):
        SHARING.End(experiment)
        return False, True
    # If the initial folder for that library does not yet exist, construct it
    folder = "initial_"  + folder
    if folder not in os.listdir(experiment["Folder"]):
        os.mkdir(folder)
        os.mkdir(folder + "/Current")
    # Otherwise, another processor is performing the initialization
    else:
        SHARING.End(experiment)
        return False, False
    # Remove any existing design molecule information
    molecules = []
    for mol in experiment["Molecules"]:
        if mol[0] != None:
            molecules.append(mol)
    # Recreate the initial Molecules           
    experiment["Molecules"] = molecules
    experiment.make_DesignGroups()
    experiment.finish_creation()
    # Load the unpositioned antigen information and scoring information
    SHARING.update_Current(experiment, experiment["Folder"] + \
                           "results/initial/")
    SHARING.load_scores(experiment, experiment["Folder"] + \
                               "results/initial/")
    # End sharing
    SHARING.End(experiment)
    return True, True

def Backbone_Perturbation(experiment, gn, optcdr):
    """Perform a backbone iteration for an OptCDR Affinity Maturation."""
    # Make sure the experiment is an Experiment
    if not isinstance(experiment, EXPERIMENT.Experiment):
        text = "The Backbone Perturbation function requires an Experiment class"
        text += " object as its input"
        raise IPRO_Error(text)
    # Find out if a refinement is happening
    refinement = IPRO_FUNCTIONS.refinement_check(experiment, gn)
    if refinement:
        return refinement
    # Assign the closest rotamers to all Residues
    # IPRO_FUNCTIONS.Closest_Rotamers(experiment, gn)
    # Start timing this
    startTime = time.time()
    # Determine which Group to make the perturbation selections for
    if gn == None:
        N = 1
    else:
        N = gn
    # Generate the Molecule name - CDR relationships
    # Get the design Molecules
    design = []
    for mol in experiment[0]:
        if mol.design:
            design.append(mol.name) 
    associations = molecule_name_association(optcdr, design)
    # Create a list of the possible CDRs
    cdrs = list(associations.keys())
    cdrs.sort()
    # Randomly select a CDR
    cdr = cdrs[random.randint(0, len(cdrs) - 1)]
    # Change the permissions and freedoms for the CDR
    perturbed = free_cdr(optcdr, associations, cdr, experiment) 
    # Add a summary for the perturbation
    experiment["Summary"] += "Perturbing the " + cdr + " CDR\n"
    # Generate the random perturbation angles
    angles = IPRO_FUNCTIONS.generate_angles(experiment[N], perturbed)
    # Loop through the Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Time the group
        groupTime = time.time()
        # Mutate to Glycines
        ROTAMERS.Glycine(group, experiment)
        # Perturb the backbone
        if experiment["Force Field"] == "CHARMM":
            CHARMM.Perturbation(group, angles, experiment)
        else:
            text = "The Backbone Perturbation function does not support the "
            text += str(experiment["Force Field"]) + " force field."
            raise IPRO_Error(text)
        # Update the summary
        if gn == None:
            experiment["Summary"] += SHARING.summary_update(groupTime, \
                        "Perturbing Design Group " + str(group.number))
    experiment["Summary"] += SHARING.summary_update(startTime, \
                             "The Backbone Perturbation function")
    return refinement

def IPRO_ITERATION(experiment, optcdr_experiment, gn = None):
    """Do an iteration of IPRO for an OptCDR Affinity Maturation."""
    # Start the iteration
    iteration, refinement = IPRO_FUNCTIONS.Start_Iteration(experiment, gn)
    if not refinement and iteration > IPRO_FUNCTIONS.max_iterations(experiment):
        return iteration
    # Do the steps of an iteration
    if not refinement:
        refinement = Backbone_Perturbation(experiment, gn, optcdr_experiment)
    if not refinement:
        refinement = IPRO_FUNCTIONS.Optimal_Rotamers(experiment, gn)
    if not refinement:
        docked, refinement = IPRO_FUNCTIONS.Docking(experiment, iteration, gn)
    if not refinement:
        refinement = IPRO_FUNCTIONS.Relaxation(experiment, gn, True)
    if not refinement:
        energies, refinement = IPRO_FUNCTIONS.Calculate_Energy(experiment, gn)
    if not refinement:
        refinement = IPRO_FUNCTIONS.End_Iteration(experiment, energies, iteration, gn)
    # If a refinement has been started and an iteration folder was created,
    # delete the iteration folder
    if refinement and iteration > 0:
        os.chdir("../")
        os.system("rm -rf iteration" + str(iteration))
    return iteration

def Affinity_Maturation(experiment, ln):
    """Run the IPRO calculations analogous to in vivo affinity maturation."""
    # Change the foler to the appropriate library folder
    os.chdir("library" + str(ln))
    # Create a new Experiment class object
    ipro_experiment = EXPERIMENT.Experiment()
    # Update the Current information
    for folder in ["./temp_Current/", "./initialize/Current/", "./Current/"]:
        try:
            SHARING.update_Current(ipro_experiment, folder)
            break
        except (IOError, OSError):
            pass
    # Delete the temporary Current folder, if it exists
    try:         
        os.system("rm -rf temp_Current")
    except OSError:
        pass
    # Initialize the IPRO experiment
    dummy = IPRO_FUNCTIONS.INITIALIZE(ipro_experiment)
    # Declare an iteration counter variable and set it to 0
    iteration = 0
    # Do the appropriate number of iterations
    while iteration <= ipro_experiment["IPRO Iterations"]:
        # If appropriate, do a refinement
        REFINEMENT.DO(ipro_experiment)
        # Do an IPRO Iteration
        iteration = IPRO_ITERATION(ipro_experiment, experiment)
    # Wait for IPRO to be finished, assisting with any Refinements that get
    # started while waiting
    IPRO_FUNCTIONS.Wait(ipro_experiment)
    # Go back to the OptCDR experiment home folder
    os.chdir("../")

def check_finish(experiment, ln):
    """If the library results are completed but have not been formatted, create
    a system flag to format the final results."""
    # Generate the system flag
    flag = False
    # If the library is in the results folder, it has already been formatted
    name = "library" + str(ln) 
    if name in os.listdir(experiment["Folder"] + "results/"):
        return flag
    # Create a shortcut for the folder name
    folder = experiment["Folder"] + name + "/"
    # Get the total number of iterations completed
    iteration = SHARING.iteration_counter(folder, False)
    # If the iterations completed is identical to the amount that should have
    # been completed, the library is finished
    if iteration == experiment['IPRO Iterations']:
        flag = True    
    return flag

def Finish(experiment, ln):
    """Perform the final OptCDR calculations to properly format the results."""
    # Move to the "results" folder within the experiment's home directory
    os.chdir(experiment["Folder"] + "results/")
    # Make a folder of the best structures in each library
    list = os.listdir("./")
    # If a "best" folder is not already in the "results" folder, make it
    if "best" not in list:
        os.mkdir("best")
    # Move to the "best" folder
    os.chdir("best")
    # Make a folder for the library
    os.mkdir("library" + str(ln))
    os.chdir("library" + str(ln))
    # Find the best iteration in the library's results folder
    folder = experiment["Folder"] + "library" + str(ln) + "/results/"
    list = os.listdir(folder)
    best = 0
    # Go through the information in the "results" folder
    for name in list:
        if name.startswith("iteration"):
            # Get the iteration number
            iteration = int(name[9:])
            # If it is higher than "best", then store its value
            if iteration > best:
                best = iteration
    # Copy the information from the "best" in that folder into the experiment's
    # home results folder
    folder += "iteration" + str(best) + "/"
    # List the files within this folder
    files = os.listdir(folder)
    # Copy each file to the experiment's results "best" folder
    for file in files:
        os.system("cp " + folder + file + " ./") 
    # List the sequence information and energy information in the summary file
    text = "LIBRARY " + str(ln) + " RESULTS\n"
    # Gather the total number of groups to have their information output
    groups = len(experiment)
    # Create a list of all Target Molecules in the experiment
    target_molecules = []
    # Go through all of the Molecules in the experiment
    for molecule in experiment[0]:
        # If it is a Target Molecule
        if not molecule.design:
            # Then store it
            target_molecules.append(molecule.name)
    # Now gather all of the Design Molecules
    molecules = []
    # Go through the files
    for file in files:
        # If it is a Molecule File, get the name of the Molecule
        name = file.split(".")[0][-1]
        # If it is in the 1st Binding Assembly (to avoid redundancy), store it
        # if it is not in the list of Target Molecules, meaning it is a Design
        # Molecule
        if file.startswith("Group1_Molecule") and name not in target_molecules:
            molecules.append(name)
    molecules.sort()
    # Create a Summary of the amino acids used within each CDR, as well as the
    # canonical structures used to make the CDRs
    # List the canonical structure information
    # Get the optimal set of canonical structures
    solution = experiment["Scores"][ln-1][1]
    # Output the score
    canonical = "The score for the set of canonical structures used is "
    canonical += str(solution["Score"]) + "\n"
    # Store the position information for each of the CDRs
    ranges = {1: range(27, 39), 2: range(56, 66), 3: range(105, 118)}
    # Go thorugh each of the CDRs and output the canonical structure used
    associations = molecule_name_association(experiment, molecules)
    cdrs = list(associations.keys())
    cdrs.sort()
    # Store the sequence information in this string
    sequence = ""
    for num, cdr in enumerate(cdrs):
        # Add the canonical structure information
        canonical += "The " + cdr + " CDR used canonical structure #"
        canonical += str(solution[num+1]) + "\n" 
        # Get the appropriate Molecule for the CDR
        name = "Group1_Molecule" + associations[cdr] + ".pdb"
        mol = MOLECULES.MoleculeFile(name)[0]
        # Go through all of the residues
        for res in mol:
            # Get its name so that its position may be extracted
            rName = res.name
            # If the name is composed of only digits
            if rName.isdigit():
                pass
            # If the last character is a letter
            elif rName[:-1].isdigit() and rName[-1].isalpha():
                rName = rName[:-1] 
            # Convert the name to an integer
            rName = int(rName)
            # If this position lies within the CDR position, add its sequence
            # information
            if rName in ranges[int(cdr[-1])]:
                sequence += cdr + " Residue " + str(rName) + " in Molecule "
                sequence += mol.name + ": " + res.kind + "\n"
    # Store the Energy information
    energy = ""
    # Go through the Binding Assemblies
    for gn in range(1, groups + 1):
        # Open the Energy file
        name = "Group" + str(gn) + "_Energies.txt"
        f = open(name, "r")
        # Go through the file
        for line in f:
            # Split the line on white space
            items = line.split()
            # Add the text to the energy string
            energy += "The " + items[0] + " " + items[1][:-1] + " of Design "
            energy += "Group " + str(gn) + " is " + items[2] + " kcal / mol\n"    
        # Close the file
        f.close()
    # Change back to the Experiment's home directory
    os.chdir(experiment["Folder"])
    # Add all of this information to the Summary file
    experiment["Summary"] += canonical + sequence + energy + "\n\n"
    name = SHARING.summary_name(SHARING.get_current())
    f = open(name, "a")
    f.write(experiment["Summary"])
    f.close()    
    # Move the library to the results folder
    command = "mv library" + str(ln) + " results/"    
    os.system(command)
    # If this is the final library, delete the SCORES.txt file
    if ln == experiment['Optcdr Libraries']:
        os.system("rm SCORES.txt")

def DO(experiment, ln):
    """Make a particular library based on the antigen positioning information"""
    # Determine if the library needs to be constructed
    do_init, do_mature = initial_check(experiment, ln)
    # If it should be, then do so
    if do_init:
        # Initialize the library
        initialize_libraries(experiment, ln)
    # If an affinity maturation should be done
    if do_mature:
        # Affinity maturation of the antibody       
        Affinity_Maturation(experiment, ln)
    # Check to see if everything is finished
    finished = check_finish(experiment, ln)
    # Finish the experiment
    if finished:
        Finish(experiment, ln)
        sys.exit(0)
