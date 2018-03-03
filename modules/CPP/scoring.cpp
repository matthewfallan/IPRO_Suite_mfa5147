/* Written in 2014 by Matthew Grisewood of the Costas Maranas Lab in the
 * Chemical Engineering Department of the Pennsylvania State University
 
 * This file generates random movements for the antigen(s) and scores their
 * interactions with the canonical structures */

// Include the needed C++ files
# include <iostream>
# include <fstream>
# include <sstream>
# include <vector>
# include <string>
# include <time.h>
# include <stdlib.h>
// Include the MOLECULES header file
# include "MOLECULES.h"
# include "DOCKING_FUNCTIONS.h"

// Use the standard namespace
using namespace std;

// Include the Molecule, Residue, and Atom classes
using MOLECULES::Atom;
using MOLECULES::Residue;
using MOLECULES::Molecule;

using DOCKING_FUNCTIONS::calculate_rmatrix;

// Gather the random antigen positions generated in Python
vector<vector<float> > gather_movements(){
    string line;
    vector<vector<float> > movements;
    vector<float> x_move, y_move, z_move, x_rotate, y_rotate, z_rotate;
    float value;
    ifstream in1("MOVEMENTS.txt");
    in1 >> value;
    while (!in1.eof()){
        x_move.push_back(value);
        in1 >> value;
        y_move.push_back(value);
        in1 >> value;
        z_move.push_back(value);
        in1 >> value;
        x_rotate.push_back(value);
        in1 >> value;
        y_rotate.push_back(value);
        in1 >> value;
        z_rotate.push_back(value);
        in1 >> value;
    }
    in1.close();
    movements.push_back(x_move);
    movements.push_back(y_move);
    movements.push_back(z_move);
    movements.push_back(x_rotate);
    movements.push_back(y_rotate);
    movements.push_back(z_rotate);
    return movements;
}

// Create the scoring function
int OptCDR_scoring(Molecule &cdr, vector<Residue> &antigen){
    // This function calculates the VDW clash score between the CDR and the
    // antigen
    int score = 0;
    for(int i=0; i<cdr.residues.size(); i++){
        for(int j=0; j<cdr.residues[i].atoms.size(); j++){
            float closest = 500;
            float radius = 0;
            for(int k=0; k<antigen.size(); k++){
                for(int l=0; l<antigen[k].atoms.size(); l++){
                    float dis = pow(cdr.residues[i].atoms[j].x -
                                antigen[k].atoms[l].x, 2);  
                    dis += pow(cdr.residues[i].atoms[j].y - 
                           antigen[k].atoms[l].y, 2);
                    dis += pow(cdr.residues[i].atoms[j].z - 
                           antigen[k].atoms[l].z, 2);
                    dis = sqrt(dis);
                    if (dis < closest){
                        closest = dis;
                        radius = antigen[k].atoms[l].vdwRadius;
                    }
                }
            }
            radius += cdr.residues[i].atoms[j].vdwRadius;
            if (closest < radius){
                score -= 5;
            }
            else if(closest < 8.0){
                score += 1;
            }
        }
    }
    return score;    
}

// Start the main program
int main(){
    // Load the random antigen movements and rotations
    vector<vector<float> > movements;
    vector<float> x_move, y_move, z_move, x_rotate, y_rotate, z_rotate;
    movements = gather_movements();
    x_move = movements[0];
    y_move = movements[1];
    z_move = movements[2];
    x_rotate = movements[3];
    y_rotate = movements[4];
    z_rotate = movements[5];
    int positions = x_move.size();
    // Open the file that contains the antigens' atoms
    ifstream in2("ANTIGENS.txt"); 
    // Create the variables needed to read in the contents of the file
    vector<Atom> atoms; 
    // Create a vector to store the antigen structures
    vector<Residue> antigens;
    // Read in the file
    string line;
    getline(in2, line);
    // Create a dummy molecule name
    string molName = "NONE";
    // Create a dummy residue number
    int resNum = -1000000;
    while(!in2.eof()){
        // Store the line as an Atom
        Atom atom; atom.load(line);
        // If this molecule name is different from the last or the residue
        // number is different from the last, store the residue
        if ((atom.moleculeName != molName) || (atom.residueNumber != resNum)){
            // Don't store the molecule if there are no stored atoms
            if ((molName != "NONE") && (resNum != -1000000)){
                // Store the Residue
                Residue res;
                res.load(atoms);
                // Append the residue to the list of antigen residues
                antigens.push_back(res);
                // Clear the list of atoms
                atoms.clear();
            }
            // Store the new Molecule name
            molName = atom.moleculeName;
            // Store the new Residue number
            resNum = atom.residueNumber;
        }
        // Append the Atom to the list of Atoms
        atoms.push_back(atom);
        // Read in the next line
        getline(in2, line);}
    // Close the file
    in2.close();
    // Store the Residue
    Residue res;
    res.load(atoms);
    // Append the Residue to the list of antigen residues
    antigens.push_back(res);

    // Open the file that contains the canonical structures' atoms
    ifstream in3("CANONICALS.txt");
    // Create the variables needed to read in the contents of the file
    atoms.clear();
    // Create a vector to store the CDR structures
    vector<Molecule> CDR;
    vector<vector<Molecule> > canonicals;
    vector<string> names;
    // Read in the file
    getline(in3, line);
    // Create a dummy Molecule name and dummy Rotamer number
    molName = "NONE";
    int rotNumber = -1;
    while(!in3.eof()){
        // Store the line as an Atom
        Atom atom; atom.load(line); 
        if ((atom.moleculeName != molName) || (atom.rotamerNumber !=
        rotNumber)){
            if ((molName != "NONE") && (rotNumber != -1)){
                Molecule mol;
                mol.load(atoms);
                CDR.push_back(mol);
                if (atom.moleculeName != molName){
                    canonicals.push_back(CDR);
                    names.push_back(molName);
                    CDR.clear();
                }    
                atoms.clear();
            }        
        }
        // Append the Atom to the list of Atoms
        atoms.push_back(atom);
        molName = atom.moleculeName;
        rotNumber = atom.rotamerNumber;
        getline(in3, line);}
    // Close the file      
    in3.close();
    Molecule mol;
    mol.load(atoms);
    CDR.push_back(mol);
    canonicals.push_back(CDR);
    names.push_back(molName);
    vector<vector<Molecule> > temp;
    string alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" ;
    int i = 0;
    while(temp.size() != canonicals.size()){
        for(int j=0; j<canonicals.size(); j++){
            if(canonicals[j][0].name[0] == alpha[i]){
                temp.push_back(canonicals[j]);
                i++;
                break;
            }
        }
    }
    canonicals = temp;

    // Calculate the score for each canonical structure for each random antigen
    // movement
    vector<vector<vector<int> > > scores;
    if(positions != x_move.size()){
        string temp = "The number of random positions generated does not ";
        temp += "match the desired number.";
        cout << temp << endl;
        exit(EXIT_FAILURE);
    } 
    // Go through each movement
    for(int i=0; i<positions; i++){
        // Make a copy of the antigen
        vector<Residue> current = antigens; 
        // Rotate
        float rmatrix [3][3];
        calculate_rmatrix(x_rotate[i], 'x', rmatrix);
        for(int j=0; j<current.size(); j++){
            current[j].rotate(rmatrix);
        }
        calculate_rmatrix(y_rotate[i], 'y', rmatrix);
        for(int j=0; j<current.size(); j++){
            current[j].rotate(rmatrix);
        }
        calculate_rmatrix(z_rotate[i], 'z', rmatrix);
        for(int j=0; j<current.size(); j++){
            current[j].rotate(rmatrix);
        }
        // Move
        float coords[3] = {x_move[i], y_move[i], z_move[i]};
        for(int j=0; j<current.size(); j++){
            current[j].move(coords, '+');
        }
        vector<vector<int> > position_scores;
        for (int j=0; j<canonicals.size(); j++){
            vector<int> cdr_scores;
            for (int k=0; k<canonicals[j].size(); k++){
                int score =  OptCDR_scoring(canonicals[j][k], current); 
                cdr_scores.push_back(score);
            }
            position_scores.push_back(cdr_scores);     
        }
        scores.push_back(position_scores);    
    }      
    ofstream out2("SCORES.txt");
    for(int i=0; i<scores.size(); i++){
        for(int j=0; j<scores[i].size(); j++){
            for(int k=0; k<scores[i][j].size(); k++){
                out2 << i+1 << " " << j+1 << " " << k+1 << " " <<
                scores[i][j][k] << endl;
            }
        }
    }
    out2.close();
    
    // End the program
    return 0;
}
