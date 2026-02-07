# Dependencies (ase)
from ase.io import read, write
from ase.visualize import view

# Add a MOF into the simulation through reading the .cif file

# Find the center of the MOF's pores in order to add the TFSi anion
# Use 2 atoms as a reference for center coordinates and then center inside based on another elements position

# Add the TFSi anion inside the MOF pores (2 anions)
# Rotate accoridngly

# Place Li cations near the Mg sites in the MOF
# Li should be lined up near the Nitrogen on the anion for proper bonding
# Base position based on the Mg site and Nitrogen

# Sourround the MOF with  solvent (CH3OCH3)
# Define a region that solvent sould not appear
# Randomly generate a certain amount of solvent around the simulation box
# If the molecule is touching the MOF, generate a new positon
# Loop until all molecules of solvent are added

# Perform an energy minimization to relax the solution
# Output the relaxed system to an xyz file for analysis
