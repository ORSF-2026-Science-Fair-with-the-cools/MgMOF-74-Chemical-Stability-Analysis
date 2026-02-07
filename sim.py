# Dependencies (ase)
from ase.io import read, write
from ase.build import molecule
from ase.visualize import view

# Add a MOF into the simulation through reading the .cif file
mof = read('MgMOF-74.cif')
view(mof)

# Find the center of the MOF's pores in order to add the TFSI anion
for atom in mof:
    print(atom.index, atom.symbol, atom.position)

# Use 2 atoms (Mg #1 and Mg #9) as a reference for center coordinates and then center inside based on another elements position (Mg #0)

# Read TFSI anion from the ntf2_pack.xyz file
# Add the TFSI anion inside the MOF pores (2 anions)
# Rotate accoridngly

# Place Li cations near the Mg sites in the MOF
# Li should be lined up near the Nitrogen on the anion for proper bonding
# Base position based on the Mg site and Nitrogen

# Sourround the MOF with DME solvent (CH3OCH3)
# Define a region that solvent sould not appear
# Randomly generate a certain amount of solvent around the simulation box
# If the molecule is touching the MOF, generate a new positon
# Loop until all molecules of solvent are added

# Perform an energy minimization to relax the solution
# Output the relaxed system to an xyz file for analysis
