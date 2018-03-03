#!/usr/bin/env python

# The name of this file
__name__ = "GAMS IPRO Suite Module"
# The documentation string
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions related to the use of GAMS for optimization
purposes in the IPRO suite of programs."""

# Include PYTHON modules
import os
import sys
# Include the STANDARDS module, which contains all "default" variables, the
# "supported" lists, the backboneAtoms and aminoAcids dictionaries, and the
# IPRO_Error class
from STANDARDS import *
import OPTCDR

# Create an error class for problems in the module
class GAMS_Error(IPRO_Error):
    """An error for problems in the GAMS module of the IPRO suite."""
    def __init__(self, error = ''):
        """The initialization of the GAMS_Error class."""
        IPRO_Error.__init__(self, error)

# Have a function that runs a GAMS script
def execute_GAMS_script(script, inputFile, outputFile):
    """Run a GAMS script."""
    # Wrap everything in a try statement to catch any errors
    try:
        # Make the input file
        f = open(inputFile, "w")
        f.write(script)
        f.close()
        # Run GAMS without printed output
        os.system(GAMSCommand + " " + inputFile + " lo=0 o /dev/null")
        # Open the output file
        f = open(outputFile, "r")
        lines = f.readlines()
        f.close()
        return lines
    # If there was an error, complain about it
    except (IOError, TypeError, OSError):
        text = "Running a GAMS script has failed. Here are the expected input "
        text += "and output file names:\nInput: " + str(inputFile)+"\nOutput: "
        text += str(outputFile)
        raise GAMS_Error(text)

def cut_previous(rotamers, previous, names, equations):
    """Use integer cuts to remove previous solutions from being possible."""
    # If the previous entry isn't a list, just return the names and equations
    # without modification
    if not isinstance(previous, list):
        return names, equations
    # Loop through the entries
    for I, group in enumerate(previous):
        # If it is not a dictionary, skip it
        if not isinstance(group, dict):
            continue
        # Make sure each entry in the dictionary corresponds to a residue
        # and a rotamer
        problem = False
        for spot in group:
            if spot not in rotamers:
                problem = True
                break
            if not isinstance(group[spot], int) or not 0 <= group[spot] < \
            len(rotamers[spot]):
                problem = True
                break
        if problem:
            continue
        # It has been confirmed that each Residue has a specific rotamer
        # designated in this solution, so use it to create integer cuts.
        # Create the name of the cut
        name = "CUT" + str(I+1)
        names += "\n\t" + name
        equation = name + " ..\n"
        # Loop through the different spots
        for spot in group:
            # identify the kind of Residue that was selected
            kind = rotamers[spot][group[spot]].kind
            # Loop through the rotamers for this Residue and include all of
            # them that are of this kind
            for i, rot in enumerate(rotamers[spot]):
                if rot.kind == kind:
                    equation += 'x("' + str(spot) + '","' + str(i+1)
                    equation += '") + '
        # Cut the last two characters off of the equation string to get rid
        # of the last +_
        equations += equation[:-2] + '=l= ' + str(len(group) - 1) + ';\n\n'
    # Return the names and equations
    return names, equations

def match_sequences(rotamers, experiment, names, equations):
    """Make sure the sequences of Dimers stay the same."""
    # Use a try statement to check the experiment input for whether or not it
    # contains a Dimers specification. The try statement is to make sure the
    # experiment is something that the 'in' function works on
    try:
        if "Dimers" not in experiment:
            return names, equations
    except (KeyError, TypeError, AttributeError):
        return names, equations
    # We know that the experiment contains information about pairs of molecules
    # that should have matching sequences, so this function may continue
    # The number of matching equations made so far
    I = 0
    # Separate the Rotamers by Molecule name, Residue name, and amino acid kind
    sorted = {}
    for i in rotamers:
        for j, r in enumerate(rotamers[i]):
            if r.moleculeName not in sorted:
                sorted[r.moleculeName] = {}
            if r.name not in sorted[r.moleculeName]:
                sorted[r.moleculeName][r.name] = {}
            if r.kind not in sorted[r.moleculeName][r.name]:
                sorted[r.moleculeName][r.name][r.kind] = []
            # Store the rotamer's position number and rotamer number
            sorted[r.moleculeName][r.name][r.kind].append([r.number, j+1])
    # Loop through the pairs of dimers
    for pair in experiment["Dimers"]:
        # If either doesn't have rotamers, skip this pair
        if pair[0] not in sorted or pair[1] not in sorted:
            continue
        # Go through all of the positions in the first Molecule
        for spot in sorted[pair[0]]:
            # Only do equations when both Residues in both Molecules allow
            # more than one kind of amino acid. Refer to the dimer_match
            # function in IPRO_FUNCTIONS for more information.
            if spot not in sorted[pair[1]] or len(sorted[pair[0]][spot]) <= 1 \
            or len(sorted[pair[1]][spot]) <= 1:
                continue
            # Go through each kind of amino acid. They should be guaranteed to
            # match
            for kind in aminoAcids[experiment["File Format"]]:
                # If this residue kind does not exist within the rotamers at
                # either position, skip the residue kind
                if kind not in sorted[pair[0]][spot].keys() and kind not in \
                sorted[pair[1]][spot].keys():
                    continue
                # Increment the equation counter
                I += 1
                # Generate a name for this equation
                name = "MATCH" + str(I)
                # Include this name in the list of equation names
                names += "\n\t" + name
                # Start creating the equation
                equation = name + " ..\n"
                # If this residue kind is absent in one set of rotamers, then it
                # should not be selected
                if kind not in sorted[pair[0]][spot].keys():
                    equation += '0 =e= '
                else:    
                    # Include the rotamers of this kind for the first Molecule
                    # on the LHS of the equation
                    for data in sorted[pair[0]][spot][kind]:
                        equation += 'x("' + str(data[0]) + '","' + \
                        str(data[1]) + '") + '
                    # Cut off the last two characters and put an equals
                    # statement
                    equation = equation[:-2] + "=e= "
                # If this residue kind is absent in the other set of rotamers,
                # then it should also not be selected
                if kind not in sorted[pair[1]][spot].keys():
                    equations += equation + '0 ;\n\n'
                else:     
                    # Put the RHS of the equation
                    for data in sorted[pair[1]][spot][kind]:
                        equation += 'x("' + str(data[0]) + '","' + \
                        str(data[1]) + '") + '
                    # Store the equation, cutting off the last two characters
                    equations += equation[:-2] + ";\n\n"
    # Return the names of the equations as well as the equations themselves
    return names, equations

def optcdr_usage(rotamers, experiment, names, equations):
    """Amino acids incorporated and total number of charged amino acids is
    restricted to be below one standard deviation from their mean in a database
    of natural antibodies."""
    # Use a try statement to determine if this is an OptCDR experiment & the
    # usage constraints should be implemented. The try statement is to to ensure
    # the experiment is given
    try:
        if experiment["Type"] != "OptCDR":
            return names, equations
    except KeyError:
        return names, equations
    # Store the positions that have rotamers
    spots = list(rotamers.keys())
    spots.sort()
    # Gather the data pertaining to amino acid usage in OptCDR
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
        # If "comp" is in neither, raise an error
        else:
            text = comp + " is an unrecognized composition type"
            raise GAMS_Error(text)
        # If this is a histidine, add the protonation states
        if comp == "HIS":
            aminos.append("HSD")
        # Initialize a mininum usage constraint equation
        minEqn = comp.upper() + "_MIN..\n"
        # Intialize a maximum usage constraint equation
        maxEqn = comp.upper() + "_MAX..\n"
        # Create a flag to determine if any rotamers meet the composition
        # criteria
        use = False
        # Go through the positions with rotamers
        for spot in spots:
            # Go through the list of rotamers
            for i in range(len(rotamers[spot])):
                # Reference the rotamer from the dictionary
                rotamer = rotamers[spot][i]
                # If this rotamer is an amino acid involved in the composition,
                # add it to the constraints
                if rotamer.kind in aminos:
                    # State that a rotamer of the proper composititon has been
                    # found
                    use = True
                    # Sum the binary variables
                    minEqn += "x('" + str(spot) + "', '" + str(i+1) + "') + "
                    maxEqn += "x('" + str(spot) + "', '" + str(i+1) + "') + "
        # Do not add the equations if there are no rotamers specified
        if not use:
            continue
        # Add the equation name
        names += "\n\t" + comp.upper() + "_MIN"
        names += "\n\t" + comp.upper() + "_MAX"
        # Add the minimum usage value to the equation            
        minEqn = minEqn[:-3] + " =g= " + str(data[comp][0]) + ";\n\n"
        # Add the maximum usage value to the equation
        maxEqn = maxEqn[:-3] + " =l= " + str(data[comp][1]) + ";\n\n"
        # Add these equations to the string of all equations
        equations += minEqn + maxEqn
    # Return the names of the equations and the actual equations
    return names, equations

def glycine_limitation(rotamers, experiment, names, equations):
    """Restrict the number of Glycines to the wild-type usage."""
    # Only use this restriction on specified experiment names
    if (not experiment["Name"].startswith("Optzyme_mcr_r")) or (not \
    experiment["Name"].startswith("MetH_R")):
        return names, equations    
    else:
        from MOLECULES import MoleculeFile
    # Store the positions that have rotamers
    spots = list(rotamers.keys())
    spots.sort()
    # Store the new equation information in this dictionary
    glyEqns = {}
    # Store the wild-type usage for each Molecule in this dictionary
    WT_glys = {}
    for spot in spots:
        # Go through the list of rotamers
        for i in range(len(rotamers[spot])):
            # Reference the rotamer from this dictionary
            rotamer = rotamers[spot][i]
            mN = rotamer.moleculeName
            # If this is a new Molecule
            if rotamer.moleculeName not in glyEqns.keys():
                glyEqns[mN] = ""
                WT_glys[mN] = 0
                WT = MoleculeFile("Molecule" + rotamer.moleculeName + ".pdb", \
                experiment["Folder"] + "structures/")
            # If this is the first rotamer at this spot, check the WT residue
            if i == 0:
                if WT[mN][rotamer.name].kind == "GLY":
                    WT_glys[mN] += 1
            # Store the binary variable information
            if rotamer.kind == "GLY":
                glyEqns[mN] += "x('" + str(spot) + "', '" + str(i+1) + "') + "
    # Write the proper restraint names
    mNs = list(glyEqns.keys())
    mNs.sort()
    for mN in mNs:
        if WT_glys[mN] == 0:
            continue
        names += "\n\tGLY_" + mN
        equations += "GLY_" + mN + "..\n" + glyEqns[mN][:-3] + " =l= "
        equations += str(WT_glys[mN]) + ";\n\n"
    # Return the names of the equations and the equations themselves
    return names, equations

# Have a function that creates special restraints for use in GAMS based on the
# type of IPRO suite program that is being run.
def special_rotamer_restraints(residues, rotamers, experiment = None, \
                               previous = None):
    """Create special restraints for use in selecting optimal rotamers."""
    # Store the names of the equations in this string
    names = ''
    # And the equation declarations themselves in this one
    equations = ''
    # Call the different special restraint functions
    names, equations = cut_previous(rotamers, previous, names, equations)
    names, equations = match_sequences(rotamers, experiment, names, equations)
    names, equations = optcdr_usage(rotamers, experiment, names, equations)
    names, equations = glycine_limitation(rotamers, experiment, names,
    equations)
    return names, equations

def make_rotamer_selector(residues, rotamers, RCE, RRE, experiment = None, \
                          previous = None):
    """Make a GAMS script to select an optimal combination of rotamers."""
    # Wrap everything in a try statement so that if there is any problem a
    # particular error message can be given
    try:
        # Retrieve any special equations
        eq_names, eqs = special_rotamer_restraints(residues, rotamers, \
                                                   experiment, previous)
        # Declare names for the input and output files for gams
        inputFile = "select_rotamers.gms"
        outputFile = "selected_rotamers.txt"
        # Identify the different positions that are receiving rotamers
        spots = list(rotamers.keys())
        spots.sort()
        # Identify the maximum number of rotamers at any single position
        max = 0
        for spot in spots:
            if len(rotamers[spot]) > max:
                max = len(rotamers[spot])
        # Start creating the script
        script = """* GAMS script for selecting an optimal rotamer arrangement

$INLINECOM /* */

option subsystems;

OPTIONS
\tsysout = off
\tsolprint = off
\treslim = 900
\titerlim = 1000000
\tdomlim = 0
\tlimcol = 0
\tlimrow = 0
\toptca = 0.0
\toptcr = 0.0
\twork = 50000000;

SETS
\ti\tthe positions with rotamers\t/
"""
        for spot in spots:
            script += "\t\t\t" + str(spot) + "\n"
        script += "/\n\tr\trotamers per residue /1*" + str(max) + """/;

ALIAS(i,j);
ALIAS(r,s);

PARAMETERS
\tE_rc(i,r)
\tE_rr(i,r,j,s)
\texist(i,r);

"""
        # Now its time to actually extract the Rotamer-Constant and
        # Rotamer-Rotamer energies from the RCE and RRE dictionaries,
        # respectively. Have a large maximum energy value to include, which is
        # necessary to keep GAMS from treating the problem as infeasible
        cutoff = 10000.0
        # First declare the E_rc (rotamer-constant portions) values
        for spot in spots:
            # Loop through the possible rotamers
            for i in range(max):
                # If this rotamer exists
                if i < len(RCE[spot]):
                    # If the value is greater than a cutoff (large, but needed
                    # to keep GAMS in the feasible range)
                    if RCE[spot][i] > cutoff:
                        RCE[spot][i] = cutoff
                    # Include the E_rc value
                    script += 'E_rc("' + str(spot) + '","' + str(i+1) + '") = '
                    script += format(RCE[spot][i], '.3f') + ";\n"
                    # And declare that this rotamer exists
                    script += 'exist("' + str(spot) +'","'+str(i+1)+'") = 1;\n'
                # If the rotamer doesn't exist, declare that
                else:
                    script += 'exist("'+str(spot)+'","'+str(i+1)+'") = 0;\n'
        # Now include the rotamer-rotamer energies
        for group in RRE:
            # If the energy is too large, use the cutoff value
            if group[4] > cutoff:
                group[4] = cutoff
            script += 'E_rr("' + str(group[0]) + '","' + str(group[1]) + '","'
            script += str(group[2]) + '","' + str(group[3]) + '") = '
            script += format(group[4], '.3f') + ';\n'
        # Keep declaring things in the script. I apologize that some of these
        # lines go longer than 80 characters, but that's allowed in GAMS and its
        # the easiest way to do this.
        script += """
VARIABLES
\tobj
\tz(i,r,j,s);

BINARY VARIABLES
\tx(i,r);

EQUATIONS
\tOBJECTIVE
\tROTAMER_CHOICE
\tZ_RESTRAINT1
\tZ_RESTRAINT2""" + eq_names + """;

OBJECTIVE ..
obj =e= sum((i,r)$(exist(i,r)), x(i,r) * E_rc(i,r)) + sum((i,r,j,s)$(ord(i) lt ord(j) and (exist(i,r) and exist(j,s))), z(i,r,j,s) * E_rr(i,r,j,s));

ROTAMER_CHOICE(i) ..
sum(r$(exist(i,r)), x(i,r)) =e= 1;

Z_RESTRAINT1(i,r,j)$(ord(i) lt ord(j) and exist(i,r)) ..
x(i,r) =e= sum(s$(exist(j,s)), z(i,r,j,s));

Z_RESTRAINT2(i,j,s)$(ord(i) lt ord(j) and exist(j,s)) ..
x(j,s) =e= sum(r$(exist(i,r)), z(i,r,j,s));

""" + eqs + """MODEL best_rotamers /
\tOBJECTIVE
\tROTAMER_CHOICE
\tZ_RESTRAINT1
\tZ_RESTRAINT2""" + eq_names + """/;
best_rotamers.workspace = 1500;

file results /""" + outputFile + """/;
put results;

z.lo(i,r,j,s) = 0;
z.up(i,r,j,s) = 1;

SOLVE best_rotamers USING mip MINIMIZING obj;

put best_rotamers.modelstat /;
put obj.l /;

LOOP(i,
        LOOP(r$(x.l(i,r)),
                put i.tl, r.tl /;
        );
);
"""
        # Return the script, input file, and output file
        return script, inputFile, outputFile
    # If there was any error
    except KeyError:
        pass
    """
    except (IOError, OSError, KeyError, TypeError, ValueError, AttributeError, \
            IPRO_Error):
        text = "There was an error in the make_rotamer_selector function in "
        text += "the " + __name__ + ". Good luck determining why."
        raise GAMS_Error(text)
    """

def optimal_rotamer_selector(residues, rotamers, RCE, RRE, experiment = None, \
                             previous = None):
    """Use GAMS to select an optimal arrangement of rotamers."""
    # Create the GAMS script
    script, inputFile, outputFile = make_rotamer_selector(residues, rotamers, \
                                    RCE, RRE, experiment, previous)
    # Run GAMS
    lines = execute_GAMS_script(script, inputFile, outputFile)
    # Store the solution here
    solution = {}
    for i, line in enumerate(lines):
        if i == 0:
            status = int(float(line) + 0.01)
        elif i == 1:
            obj = float(line)
        else:
            solution[int(line.split()[0])] = int(line.split()[1]) - 1
    # If the status is not acceptable (refer to pages 88-9 in "GAMS-A User's
    # Guide" by Richard E. Rosenthal for meanings), raise an error. 1 is an
    # optimal solution and 8 is an integer solution
    if status not in [1, 8]:
        text = "The optimal_rotamer_selection function in GAMS did not converge"
        text += " to an appropriate solution. Good luck determining why."
        raise GAMS_Error(text)
    # With the solution collected, delete the input and output files
    #os.remove(inputFile)
    #os.remove(outputFile)
    # Return the status, objective function value, and identified solution
    return obj, solution

def make_canonical_selector(canonicals, clashes, scores, results = []):
    """Make a GAMS script to select an optimal set of canonical structures
    within OptCDR."""
    # Wrap the script generation within a "try" statement so an error can be
    # provided if the construction fails
    try:
        # Declare names for the GAMS input and output files
        inputFile = "select_canonicals.gms"
        outputFile = "selected_canonicals.txt"
        # Identify the cdrs that are being used in the experiment
        cdrs = list(canonicals.keys())
        cdrs.sort()
        # Determine the maximum number of canonical structures possible in any
        # CDR
        max = 0
        for cdr in cdrs:
            if len(canonicals[cdr].keys()) > max:
                max = len(canonicals[cdr].keys())
        # Create the parameters for the formulation
        t1 = 'S("'
        t2 = 'exist("'
        t3 = '","'
        t4 = '") = '
        t5 =';\n'
        # Store the parameters in a string
        parameters = ''
        # Go through each CDR
        for cdr in cdrs:
            # Go from canonical structure #1 to the maximum number of canonical
            # structures possible
            for i in range(1, max + 1):
                # Find the number associated with the CDR
                index = cdrs.index(cdr) + 1
                # If the canonical structure exists, write the appropriate
                # parameters
                if (index, i) in scores:
                    # Write the score for the canonical structure
                    parameters += t1 + str(index) + t3 + str(i) + t4
                    parameters += str(scores[(index, i)]) + t5;
                    # Write that the canonical structure exists
                    parameters += t2 + str(index) + t3 + str(i) + t4 + '1' + t5
                # Otherwise, state that it does not exist     
                else:
                    parameters += t2 + str(index) + t3 + str(i) + t4 + '0' + t5
        # Create the integer cut equations and clash constraints for the MILP 
        # formulation              
        t1 = 'y("'
        t2 = '") + '
        t6 = '") =l= 1;\n'
        t7 = '") =l= ' + str(len(cdrs) - 1) + ';\n'
        # Store the equations in a string             
        equations = ''
        # Store the equation names in a string
        equation_names = ''
        # Keep track of the number of generated equations
        count = 0
        # Go through all of the clashes
        for set in clashes:
            # Increment the counter
            count += 1
            # Get the indices for the CDRs involved in the clash
            cdr1 = cdrs.index(set[0]) + 1
            cdr2 = cdrs.index(set[2]) + 1
            # Get the numbers of the canonical structures
            cn1 = set[1]
            cn2 = set[3]
            # Create a name for the clash constraint
            equation_names += "\n\tCLASH" + str(count)
            # Create the actual constraint
            equations += "\nCLASH" + str(count) + " ..\n"
            equations += t1 + str(cdr1) + t3 + cn1 + t2 + t1 + str(cdr2) + t3
            equations += cn2 + t6
        # Generate integer cuts to eliminate previously chosen solutions
        # Reset the counter
        count = 0
        # Now generate the integer cuts to eliminate previously chosen solutions
        for solution in results:
            # Increment the counter
            count += 1
            # Create a name for the integer cut
            equation_names += "\n\tIC" + str(count)
            # Create the integer cut
            equations += "\nIC" + str(count) + " ..\n"
            # Go through each CDR
            for number in range(1, len(cdrs) + 1):
                # Write the previously-selected canonical structures
                equations += t1 + str(number) + t3 + str(solution[1][number]) 
                if number != len(cdrs):
                    equations += t2
                else:
                    equations += t7
        # Create a script to be executed within GAMS
        script = """* Code for selecting the optimal set of canonical structures in OptCDR

$INLINECOM /* */

option subsystems;

OPTIONS
\tsysout = off
\tsolprint = off
\treslim = 900
\titerlim = 1000000
\tdomlim = 0
\tlimcol = 0
\tlimrow = 0
\toptca = 0.0
\toptcr = 0.0
\twork = 50000000;

SETS
c   The antibody CDR                /1*""" + str(len(cdrs)) + """/
r   The canonical structure         /1*""" + str(max) + """/

PARAMETERS
\tS(c,r)
\texist(c,r);

""" + parameters + """

VARIABLES
\tz;

BINARY VARIABLES
\ty(c,r);

EQUATIONS
\tOBJ
\tCANONICAL_SELECTION""" + equation_names + """;

OBJ ..
z =e= sum(c, sum(r$(exist(c,r)), y(c,r)*S(c,r)));

CANONICAL_SELECTION(c) ..
sum(r$(exist(c,r)), y(c,r)) =e= 1;
""" + equations + """
MODEL best_canonicals /
\tOBJ
\tCANONICAL_SELECTION""" + equation_names + """
/;

file results /""" + outputFile + """/;
put results;

SOLVE best_canonicals USING mip MAXIMIZING z;

put best_canonicals.modelstat /;
put z.l /;

LOOP(c,
        LOOP(r$(y.l(c,r)),
                put c.tl, r.tl /;
        );
);
"""
    # Return the script, input file name, and output file name
        return script, inputFile, outputFile
    # If there was any kind of error
    except (IOError, OSError, KeyError, TypeError, ValueError, AttributeError, \
        IPRO_Error):
        text = "There was an error in the make_canonical_selector function "
        text += "in " + __name__ + ". Good luck determining why."
        raise GAMS_Error(text)

def optcdr_canonicals(canonicals, clashes, scores, cuts = []):
    """Use GAMS to select the optimal set of canonical structures."""
    # Create the GAMS script
    script, inputFile, outputFile = make_canonical_selector(canonicals, \
                                    clashes, scores, cuts)
    # Run GAMS
    lines = execute_GAMS_script(script, inputFile, outputFile)
    # Store the solution
    solution = {}
    # Go through the output file
    for i, line in enumerate(lines):
        # If it is the first line, it is the model status
        if i == 0:
            status = int(float(line) + 0.01)
        # If this is the second line, it is the objective function value
        elif i == 1:
            obj = int(float(line))
        # If it is any other line, it is the CDR and the appropriate canonical
        # structure. Store this information.
        else:
            items = line.split()
            solution[int(items[0])] = int(items[1])
    # If the status is not acceptable (refer to pages 88 & 89 in "GAMS- A User's
    # Guide" by Richard E. Rosenthal for meanings), raise an error. 1 is an
    # optimal solution and 8 is an integer solution
    if status not in [1, 8]:
        text = "The optcdr_canonicals function in GAMS was unable to converge "
        text += "to a solution. Please examine the GAMS output for more "
        text += "information."
        raise GAMS_Error(text)
    # With the MILP solved, delete the GAMS input and output files
    #os.remove(inputFile)
    #os.remove(outputFile)
    # Add the total score to the solution
    solution["Score"] = obj
    # Return the identified solution
    return solution     
