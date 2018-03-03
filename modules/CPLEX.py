#!/usr/bin/env python

# Name of the file
__name__ = "IPRO Suite CPLEX Interface"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions for accessing the CPLEX Modeling environment"""

# Include normal PYTHON modules
import os
import sys
# Include all contents from the STANDARDS module
from STANDARDS import *
import OPTCDR
import EXPERIMENT
# Allow access to the CPLEX solver
sys.path.append(CPLEXFolder)
import cplex
import shutil
import DEIMMUNIZATION

class CPLEX_Error(IPRO_Error):
    """An Error for problems in the IPRO Suite CPLEX Module"""
    def __init__(self, error = ''):
        """The initialization of the CPLEX_Error"""
        IPRO_Error.__init__(self, error)

def cut_previous_selection(model, rotamers, previous):
    """Create an integer cut to eliminate previous solutions"""
    # Use a lot of error checking to make sure the generated equations are valid
    if isinstance(previous, list) and len(previous) > 0:
        # Loop through the options in the list
        for I, solution in enumerate(previous):
            # If it isn't a dictionary of answers, skip it
            if not isinstance(solution, dict) or len(solution) != len(rotamers):
                continue
            problem = False
            for rn in solution:
                if not isinstance(solution[rn], int) or rn not in rotamers or \
                not 0 <= solution[rn] <= len(rotamers[rn]):
                    problem = True
                    break
            if problem:
                continue
            # We've validated that there is a choice for every Residue and that
            # each choice refers to a specific rotamer, so it is possible to
            # make the restraint. Create a list of variables and coefficients
            vars = []
            coefs = []
            # Go through the Residues in the solution
            for rn in solution:
                # Get the kind of the relevant rotamer
                kind = rotamers[rn][solution[rn]].kind
                # Go through the rotamers that match this kind and include them
                # in the sum
                for i, rot in enumerate(rotamers[rn]):
                    if rot.kind == kind:
                        vars.append("X_" + str(rn) + "_" + str(i+1))
                        coefs.append(1)
            # Add an integer cut to the model to eliminate the possibility of
            # this solution being selected again
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, \
                                coefs)], senses = ["L"], rhs = [len(solution)-1])

def make_immunization_restraints(model, rotamers, residues):
    """Create an integer cut to eliminate non hummanizaiton design"""
    spots = []
    molName = ""
    indexDict = {}
    for molecule in residues:
        for i, residue in enumerate(molecule):
            if residue.permission == "ROTAMER":
                spots.append(residue.name)
                molName = residue.moleculeName
                indexDict[(residue.moleculeName, residue.name)] = residue.number
    if len(spots) == 0:
        return
    if not len(spots) in [1, 2, 3]:
        text = "The number of perturbation residues for OptMAVEn is not correct"
        raise CPLEX_Error(text)
    molecule = residues[molName]
    # Need to copy the database file to current folder
    path = "/gpfs/group/cdm/IPRO_Suite/databases/OptMAVEn/"
    fileName = "Human_9mer_Sequences.txt"
    shutil.copyfile(path + fileName, fileName)
    #humanSeqs = DEIMMUNIZATION.load_human_sequences()
    cuts = DEIMMUNIZATION.make_humanization_cuts(molecule, spots)
    #print "cuts: ", cuts
    if isinstance(cuts, list) and len(cuts) > 0:
        for constraint in cuts:
            #print constraint
            # We've validated that there is a choice for every Residue and that
            # each choice refers to a specific rotamer, so it is possible to
            # make the restraint. Create a list of variables and coefficients
            vars = []
            coefs = []
                # Go through the rotamers that match this kind and include them
                # in the sum
            residueName = constraint.keys()
            length = len(constraint)
            for name in residueName:
                residueNumber = indexDict[(molName, name)]
                for i, rotamer in enumerate(rotamers[residueNumber]):
                    if rotamer.kind == constraint[name]:
                        vars.append("X_" + str(residueNumber) + "_" + str(i+1))
                        coefs.append(1)
            # Add an integer cut to the model to eliminate the possibility of
            # this solution being selected again
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, \
                                coefs)], senses = ["L"], rhs = [length-1])

    else:
        pass
        #text = "The humanization cuts are not correct"
        #raise CPLEX_Error(text)

def match_sequences(model, rotamers, experiment):
    """Make sure Dimer sequences match"""
    # Make sure the Experiment contains dimer information
    if isinstance(experiment, (EXPERIMENT.Experiment, dict)) and "Dimers" in \
    experiment:
        # Separate the rotamers by Molecule name, Residue name, and amino acid
        # kind
        sorted = {}
        for N in rotamers:
            for i, rot in enumerate(rotamers[N]):
                if rot.moleculeName not in sorted:
                    sorted[rot.moleculeName] = {}
                if rot.name not in sorted[rot.moleculeName]:
                    sorted[rot.moleculeName][rot.name] = {}
                if rot.kind not in sorted[rot.moleculeName][rot.name]:
                    sorted[rot.moleculeName][rot.name][rot.kind] = []
                # Store the rotamer's position number and rotamer number
                sorted[rot.moleculeName][rot.name][rot.kind].append([N, i+1])
        # Now go through the Dimers
        for pair in experiment["Dimers"]:
            # If there are no rotamers for one of the Molecules, it can be
            # skipped
            if pair[0] not in sorted or pair[1] not in sorted:
                continue
            # Only create equations when there are more than one kind of amino
            # acid allowed for each of the Residues
            for rn in sorted[pair[0]]:
                if rn not in sorted[pair[1]] or len(sorted[pair[0]][rn]) < 2:
                    continue
                # go through each kind of amino acid. We know they match up, but
                # check that anyway
                for kind in sorted[pair[0]][rn]:
                    if kind not in sorted[pair[1]][rn]:
                        continue
                    # Store the created variables and coefficients
                    vars = []
                    coefs = []
                    # Include all of the rotamers from the first Molecule
                    for data in sorted[pair[0]][rn][kind]:
                        vars.append("X_" + str(data[0]) + "_" + str(data[1]))
                        coefs.append(1)
                    # Include all of the rotamers from the second Molecule, but
                    # subtract them instead
                    for data in sorted[pair[1]][rn][kind]:
                        vars.append("X_" + str(data[0]) + "_" + str(data[1]))
                        coefs.append(-1)
                    # Store the linear constraint for this kind of amino acid at
                    # this pair of dimer positions in the model
                    model.linear_constraints.add(lin_expr = \
                    [cplex.SparsePair(vars, coefs)], senses = ["E"], rhs = [0])

def optcdr_usage(model, rotamers, experiment):
    """Amino acids incorporated and total number of charged amino acids is
    restricted to be below one standard deviation from their mean in a database
    of natural antibodies."""
    # Use a try statement to determine if this is an OptCDR experiment & the
    # usage constraints should be implemented. The try statement is to to ensure
    # the experiment is given
    try:
        if experiment["Type"] == "OptCDR":
            # Store the constraints in these lists and strings
            constr_lhs = []
            constr_rhs = []
            constr_senses = []
            # Store the positions that have rotamers
            spots = list(rotamers.keys())
            spots.sort()
            # Gather the data pertaining to the amino acid usage in OptCDR
            data = OPTCDR.optcdr_usage_limits(rotamers, experiment)
            # Go through the data
            info = list(data.keys())
            info.sort()
            # Go through each composition constraint
            for comp in info:
                # Get the list of amino acids meeting the composition
                # If it is an amino acid, only that type of amino acid meets the
                # criteria
                if comp in aminoAcids[experiment["File Format"]]:
                    aminos = [comp]
                # Otherwise, multiple amino acids can meet the criteria
                elif comp in propertiesAA[experiment["File Format"]]:
                    aminos = propertiesAA[experiment["File Format"]][comp]
                # If "comp" is neither, raise an error
                else:
                    text = comp + " is an unrecognized composition type"
                    raise CPLEX_Error(text)
                # If this is a histidine, add the protonation states
                if comp == "HIS":
                    aminos.append("HSD")
                # Create a flag to determine if any rotamers meet the
                # composition criteria
                use = False
                # Store the binary variables in this list
                bns = []
                # Go through the positions with rotamers
                for spot in spots:
                    # Go through the list of rotamers
                    for i in range(len(rotamers[spot])):
                        # Reference the rotamer from the dictionary
                        rotamer = rotamers[spot][i]
                        # If this rotamer is an amino acid involved in the
                        # composition limitation, add it to the constraint
                        if rotamer.kind in aminos:
                            # State that a rotamer of the proper composition has
                            # been found
                            use = True
                            # Add the binary variable to the list
                            bn = "X_" + str(spot) + "_" + str(i+1)
                            bns.append(bn)
                # Do not add the constraints if no rotamers met the composition
                # conditions
                if not use:
                    continue
                # Otherwise, add the constraints      
                constr_lhs.append([bns, [1]*len(bns)])
                constr_rhs.append(float(data[comp][0]))
                constr_senses.append("G")
                constr_lhs.append([bns, [1]*len(bns)])
                constr_rhs.append(float(data[comp][1]))
                constr_senses.append("L")
            # If constraints exist, add them to the CPLEX model
            if constr_lhs != []:
                model.linear_constraints.add(lin_expr = constr_lhs, senses = \
                                             constr_senses, rhs = constr_rhs)
    except KeyError:
        pass

def special_rotamer_restraints(model, rotamers, experiment = None, previous = \
                               None):
    """Create special integer cuts for use in the CPLEX model"""
    # Add any previous solutions to the model
    cut_previous_selection(model, rotamers, previous)
    # And match the sequences between dimers
    match_sequences(model, rotamers, experiment)
    # Add the OptCDR amino acid usage pattern constraints
    optcdr_usage(model, rotamers, experiment)

def make_rotamer_selector(residues, rotamers, RCE, RRE, experiment = None, previous = \
                          None):
    """Make a CPLEX model to select an optimal combination of rotamers"""
    # First, make a CPLEX variable
    model = cplex.Cplex()
    # Set it up as an MILP that minimizes a value
    model.set_problem_type(cplex.Cplex.problem_type.MILP)
    model.objective.set_sense(model.objective.sense.minimize)
    # Create the objective function, which minimizes the sum of the selected
    # rotamers with the constant portions of the system and with each other
    objV = []
    objC = []
    # Go through the positions receiving rotamers
    spots = list(rotamers.keys())
    spots.sort()
    for i in spots:
        for r in range(1, len(rotamers[i]) + 1):
            # Store the binary variable for this position
            objV.append("X_" + str(i) + "_" + str(r))
            # And its energy with the constant portions of the system. Have a
            # large maximum penalty
            e = RCE[i][r - 1]
            if e > 10000:
                e = 10000.0
            objC.append(e)
    # Now create the 'double' binary variable term (i.e the ones that contain
    # the rotamer rotamer energy information)
    for data in RRE:
        objV.append("Z_"+str(data[0]) + "_" + str(data[1]) + "_" + str(data[2])\
                    + "_" + str(data[3]))
        e = data[4]
        # EDIT by mfa5147:
        # IPRO is now programmed to ignore rotamer-rotamer interaction energies
        # and optimize rotamer-constant interaction energies alone, but because
        # ignoring rotamer-rotamer interaction energies could lead to clashes
        # between rotamers, if the interaction energy is above a threshhold at
        # which the rotamers are considered to clash, then the large positive
        # interaction energy is retained to prevent both rotamers from being
        # chosen simulataneously.
        # Upper threshold (to prevent the problem from blowing up)
        threshold_upper = 10000.0
        # Clash threshold (above which rotamers are considered to clash and
        # their interaction energies are not ignored.
        threshold_clash = 100.0
        if e > threshold_upper:
            e = threshold_upper
        elif e < threshold_clash:
            e = 0.0 
        objC.append(e)
    # Put the objective function in the CPLEX model
    model.variables.add(names = objV, obj = objC, lb = [0] * len(objV), ub = \
                        [1] * len(objV), types = ['B'] * len(objV))
    # Include the restraint to select exactly one rotamer at each position
    for i in spots:
        vars = []
        for r in range(1, len(rotamers[i]) + 1):
            vars.append("X_" + str(i) + "_" + str(r))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])
    # Create the two linearization restraints that relate Z to X
    # For each i, r, and j, the sum over s of Z(i,r,j,s) = x(i,r) ->
    # (sum over s of Z(i,r,j,s)) - x(i,r) = 0
    for i in spots:
        for r in range(1, len(rotamers[i]) + 1):
            for j in spots:
                if j <= i:
                    continue
                vars = ["X_" + str(i) + "_" + str(r)]
                coefs = [-1]
                for s in range(1, len(rotamers[j]) + 1):
                    vars.append("Z_" + str(i) + "_" + str(r) + "_" + str(j) + \
                                "_" + str(s))
                    coefs.append(1)
                model.linear_constraints.add(lin_expr=[cplex.SparsePair(vars, \
                                             coefs)], senses = ["E"], rhs = [0])
    # for each j, s, and i, the sum over r of Z(i,r,j,s) = x(j,s)
    for i in spots:
        for j in spots:
            if j <= i:
                continue
            for s in range(1, len(rotamers[j]) + 1):
                vars = ["X_" + str(j) + "_" + str(s)]
                coefs = [-1]
                for r in range(1, len(rotamers[i]) + 1):
                    vars.append("Z_" + str(i) + "_" + str(r) + "_" + str(j) + \
                                "_" + str(s))
                    coefs.append(1)
                model.linear_constraints.add(lin_expr=[cplex.SparsePair(vars, \
                                             coefs)], senses = ["E"], rhs = [0])
    # Add in any special restraints
    special_rotamer_restraints(model, rotamers, experiment, previous)
    # Add in humanization restraints
    if experiment["Type"] == "OptMAVEN":
        make_immunization_restraints(model, rotamers, residues)
    # Return the CPLEX model
    return model

def optimal_rotamer_selector(residues, rotamers, RCE, RRE, experiment = None, previous = \
                             None):
    """Use CPLEX to select an optimal combination of rotamers"""
    # Make the CPLEX model
    model = make_rotamer_selector(residues, rotamers, RCE, RRE, experiment, previous)
    # Suppress the printed output to the screen
    # model.set_results_stream(None)
    # Set the time limit (seconds) for cplex optimizer
    if experiment != None and experiment["Type"] == "OptCDR":
        model.parameters.timelimit.set(7200)
    else:     
        model.parameters.timelimit.set(900)
    # Solve the model
    model.solve()
    # Get the status of the solution
    status = model.solution.get_status()
    if status not in [101, 108]:
        text = "The optimal rotamer selector function did not find an optimal "
        text += "solution. Good luck determining why."
        raise CPLEX_Error(text)
    objective = model.solution.get_objective_value()
    # Determine what rotamers were selected and calculate the energy of
    solution = {}
    spots = list(rotamers.keys())
    spots.sort()
    # Go through the positions and rotamers for that position
    for i in spots:
        for r in range(1, len(rotamers[i]) + 1):
            # Get the value of the binary variable
            x = model.solution.get_values("X_" + str(i) + "_" + str(r))
            if x > 0.01:
                if i not in solution:
                    solution[i] = r - 1
                else:
                    text = "The CPLEX solution has multiple selections for "
                    text += "Residue Number " + str(i)
                    raise CPLEX_Error(text)
    return objective, solution

def make_canonical_selector(canonicals, clashes, scores, results = []):
    """Use CPLEX to select an optimal set of canonical structures."""
    # First, make a CPLEX variable
    model = cplex.Cplex()
    # Set it up as an MILP that maximizes the objective function value
    model.set_problem_type(cplex.Cplex.problem_type.MILP)
    model.objective.set_sense(model.objective.sense.maximize)
    # Identify the cdrs that are being used in the experiment
    cdrs = list(canonicals.keys())
    cdrs.sort()
    # Store the parameters and their associated values in lists
    parameter_names = []
    parameter_values = []
    # Go through each CDR
    for cdr in cdrs:
        # Find the number associated with the CDR
        index = cdrs.index(cdr) + 1
        # Go through each canonical structure within the CDR
        for i in range(1, len(canonicals[cdr].keys()) + 1):
            # Store the parameter name
            par_name = "Y_" + str(index) + "_" + str(i)
            parameter_names.append(par_name)
            # Store the parameter value
            parameter_values.append(scores[(index, i)])
    # Add the variables (binary times parameter value) to the CPLEX model
    lbs = [0] * len(parameter_names)
    ubs = [1] * len(parameter_names)
    Types = ['B'] * len(parameter_names)
    model.variables.add(names = parameter_names, obj = parameter_values, \
                        lb = lbs, ub = ubs, types = Types)
    # Store the constraints in these lists
    constr_lhs = []
    constr_rhs = []
    # Store the inequalities/equalities as a string
    constr_senses = []
    # Add the clash constraints for the MILP formulation
    for set in clashes:
        # Get the indices for the CDRs involved in the clash
        cdr1 = cdrs.index(set[0]) + 1
        cdr2 = cdrs.index(set[2]) + 1
        # Get the numbers of the canonical structures
        cn1 = set[1]
        cn2 = set[3]
        # Get the two binary variables' names
        bn1 = "Y_" + str(cdr1) + "_" + cn1
        bn2 = "Y_" + str(cdr2) + "_" + cn2
        # Add the clash constraint to the constraint list
        constr_lhs.append([[bn1, bn2],[1,1]])
        constr_rhs.append(1)
        constr_senses.append("L")
    # Now generate the integer cuts to eliminate previously chosen solutions    
    for solution in results:
        # Store the binary variables for the solution
        bns = []
        # Go through each CDR
        for number in range(1, len(cdrs) + 1):
            # Get the binary variable of the previous solution
            bn = "Y_" + str(number) + "_" + str(solution[1][number]) 
            # Add the binary to the list
            bns.append(bn)
        # Add the integer cut to the constraint list      
        constr_lhs.append([bns, [1]*len(bns)])
        value = len(cdrs) - 1
        constr_rhs.append(value)
        constr_senses.append("L")
    # Ensure that only one canonical structure is selected for each CDR
    # Go through each CDR
    for cdr in cdrs:
        # Find the number associated with the CDR
        index = cdrs.index(cdr) + 1
        # Store the binary variables in this list
        bns = []
        # Go through each canonical structure within the CDR
        for i in range(1, len(canonicals[cdr].keys()) + 1): 
            # Get the associated binary variable
            bn = "Y_" + str(index) + "_" + str(i)
            # Add this binary variable to the list of binary variables
            bns.append(bn)
        # Add the constraint to the constraint list
        constr_lhs.append([bns, [1]*len(bns)])
        constr_rhs.append(1)
        constr_senses.append("E")
    # Add the restraints for the MILP        
    model.linear_constraints.add(lin_expr = constr_lhs, senses = constr_senses,
                                 rhs = constr_rhs)
    # Return the CPLEX model
    return model

def optcdr_canonicals(canonicals, clashes, scores, cuts = []):
    """Use CPLEX to select the optimal set of canonical structures."""
    # Create the CPLEX model
    model = make_canonical_selector(canonicals, clashes, scores, cuts)
    # Suppress the printed output to the screen
    model.set_results_stream(None)
    # Set the time limit for the CPLEX optimizer
    model.parameters.timelimit.set(900)
    # Solve the model
    model.solve()
    # Get the status of the solution
    status = model.solution.get_status()
    if status not in [101, 108]:
        text = "The optcdr_canonicals function in CPLEX was unable to converge"
        text += " to a solution. Remove the supression of the model output "
        text += "in the CPLEX.py file and re-run for more information."
        raise CPLEX_Error(text)
    # Get the objective function value       
    objective = model.solution.get_objective_value()
    # Store the solution
    solution = {}
    # Identify the cdrs that are being used in the experiment
    cdrs = len(canonicals.keys())
    cdrs.sort()
    # Store the CDRs' corresponding canonical structures 
    for cdr in cdrs:
        # Find the number associated with the CDR
        index = cdrs.index(cdr) + 1
        # Go through each canonical structure within the CDR
        for i in range(1, len(canonicals[cdr].keys()) + 1):    
            # Get the binary variable
            bn = "Y_" + str(index) + "_" + str(i)
            # Get the binary variable value
            y = model.solution.get_values(bn)
            # Store the selected canonical structures
            if y == 1.0:
                # Raise an error if there was more than 1 canonical structure
                # selected for this CDR
                if index in solution.keys():
                    text = "There was more than one canonical structure "
                    text += "selected per CDR. This indicates a problem with "
                    text += "the MILP formulation within CPLEX."
                    raise CPLEX_Error(text)
                else:
                    solution[index] = i
    # Raise an error if there were not enough canonical structures selected
    if len(solution.keys()) != len(cdrs):
        text = "There were not canonical structures selected for each CDR in "
        text += "the OptCDR canonical selector. This indicates a problem with "
        text += "the MILP formulation within CPLEX."
        raise CPLEX_Error(text)
    # Add the total score to the solution
    solution["Score"] = int(objective)
    # Return the identified solution
    return solution

def make_OptMAVEn_selector(energies, struCuts, solutionCuts):
    """Make a CPLEX model to select an optimal combination of MAPs parts"""
    # Generate a CPLEX model
    model = cplex.Cplex()
    # Set it up as an MILP minimizing the energy
    model.set_problem_type(cplex.Cplex.problem_type.MILP)
    model.objective.set_sense(model.objective.sense.minimize)
    # Make a list of the parts, divided by region
    parts = {}
    for data in energies:
        if data[0] not in parts:
            parts[data[0]] = []
        parts[data[0]].append(data[1])
    # Determine the types of domains that are being created
    domains = []
    for part in parts:
        if part.startswith("H") and "H" not in domains:
            domains.append("H")
        elif part.startswith(("K", "L")) and "L" not in domains:
            domains.append("L")
    # Create the objective function
    # Store the variables here
    objV = []
    # and their coefficients (i.e. energies) here
    objC = []
    # go through the energies
    for data in energies:
        objV.append("X_" + data[0] + "_" + str(data[1]))
        objC.append(data[2])
    model.variables.add(names = objV, obj = objC, lb = [0] * len(objV), ub = \
                        [1] * len(objV), types = ['B'] * len(objV))
    # Now create the restraints. First, make sure exactly one heavy V* structure
    # is selected
    if "H" in domains:
        vars = []
        for i in parts["HV"]:
            vars.append("X_HV_" + str(i))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])
        vars = []
        for i in parts["HCDR3"]:
            vars.append("X_HCDR3_" + str(i))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])
        vars = []
        for i in parts["HJ"]:
            vars.append("X_HJ_" + str(i))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])

    # Now make sure only one light V* is selected
    if "L" in domains:
        vars = []
        for p in ["KV", "LV"]:
            if p not in parts:
                continue
            for i in parts[p]:
                vars.append("X_" + p + "_" + str(i))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])

        # Now make sure exactly one CDR3 and J* structure is selected for each type
        # of variable domain
        for L in ["K", "L"]:
            if L + "V" not in parts:
                continue
            # Store the V* info here
            Vvars = []
            Vcoefs = []
            for i in parts[L + "V"]:
                Vvars.append("X_" + L + "V_" + str(i))
                Vcoefs.append(-1)
            # Do this for each of the other two regions of structure
            for R in ['CDR3', 'J']:
                if L + R not in parts:
                    continue
                vars = []
                coefs = []
                for i in parts[L+R]:
                    vars.append("X_" + L + R + "_" + str(i))
                    coefs.append(1)
                # Include the V* information in the restraint
                vars.extend(Vvars)
                coefs.extend(Vcoefs)
                # Make the restraint
                model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, \
                                            coefs)], senses = ["E"], rhs = [0])
    # Include the integer cuts on structure combinations
    for cut in struCuts:
        if cut[0] in parts and cut[1] in parts[cut[0]] and \
           cut[2] in parts and cut[3] in parts[cut[2]]:
            # Say that those two cannot be selected together
            vars = ["X_" + cut[0] + "_" + str(cut[1]), \
                    "X_" + cut[2] + "_" + str(cut[3])]
            coefs = [1, 1]
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, \
                                         coefs)], senses = ["L"], rhs = [1])
    # Include the integer cuts from previous optimized solutions for finding the next optimized one
    # First delete the HJ or LJ or KJ from the solution cuts because they almost contribute nothing to the energy
    # Otherwise you can get almost the same parts solution except the 'J' parts are different
    newSolutionCuts = []
    for cut in solutionCuts:
        index = 0
        newCut = []
        length = len(cut)
        while(index < length-1):
            if not cut[index][1] == 'J':
                newCut.append(cut[index])
                newCut.append(cut[index + 1])
            index += 2
        newSolutionCuts.append(newCut)
    for cut in newSolutionCuts:
        if cut[0] in parts and cut[1] in parts[cut[0]] and \
           cut[2] in parts and cut[3] in parts[cut[2]] and \
           cut[4] in parts and cut[5] in parts[cut[4]] and \
           cut[6] in parts and cut[7] in parts[cut[6]]:
            vars = ["X_" + cut[0] + "_" + str(cut[1]), \
                    "X_" + cut[2] + "_" + str(cut[3]), \
                    "X_" + cut[4] + "_" + str(cut[5]), \
                    "X_" + cut[6] + "_" + str(cut[7])]
            coefs = [1, 1, 1, 1]
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, \
                                         coefs)], senses = ["L"], rhs = [3])

    # Return the model and the dictionary of parts
    return model, parts

def OptMAVEn_selector(energies, struCuts, solutionCuts):
    """Use CPLEX to select an optimal combination of MAPs parts"""
    # Make the model and get a dictionary containing the different kinds of
    # parts
    model, parts = make_OptMAVEn_selector(energies, struCuts, solutionCuts)
    # Suppress the output that gets printed to the screen
    #model.set_results_stream(None)
    # Solve the model
    model.solve()
    # Get the status of the solution
    status = model.solution.get_status()
    if status != 101 and status != 102:
        text = "The OptMAVEn selector of structures did not work. Good luck "
        text += "determining why."
        raise CPLEX_Error(text)
    # Get the objective value
    objective = model.solution.get_objective_value()
    # Store the solution here
    solution = {}
    for part in parts:
        for i in parts[part]:
            # The name of this variable
            name = "X_" + part + "_" + str(i)
            # The value of the variable
            x = model.solution.get_values(name)
            if x > 0.01:
                if part not in solution:
                    solution[part] = i
                else:
                    text = "The CPLEX solution has multiple parts selected for "
                    text += part
                    raise CPLEX_Error(text)
    return solution, objective
