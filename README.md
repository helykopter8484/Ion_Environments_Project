# Ion Environments Project
This repository contains code for the Computational Geometry based analysis of ion environments in three dimensional protein structures, created for the thesis project.

### Prerequisites
The scripts are written in Python and use a few other standard libraries and biopython.
* Python
* Biopython

#### Driver.py
This is the main module which is used to call other modules and perform the analysis.
The sequence identity cutoff value is the main variable parameter

#### PDB.py
This module contains the PDB Class which has methods to:
* Download PDB Files
* Describe PDB Files
* Get unique chains from PDB files, based on ATOM information (not sequence based)
* Extract specific chains from a PDB file
* Create unique PDB entries based on a specific HETATM records for a PDB chain

#### Subprocess.py
This module uses the subprocess module to call upon a few Java programs.

#### Calculations.py
This module performs calculations relevent to the computational geometry based analysis.
