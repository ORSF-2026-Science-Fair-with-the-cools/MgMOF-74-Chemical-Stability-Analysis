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
pos_mg1 = mof.get_positions()[1]
pos_mg9 = mof.get_positions()[9]
center_point = (pos_mg1 + pos_mg9) / 2

# Read TFSI anion from the ntf2_pack.xyz file
# Add the TFSI anion inside the MOF pores (2 anions)
# Rotate accoridngly
tfsi = read('ntf2_pack.xyz')
tfsi.translate(center_point) # Move to the calculated center
mof.extend(tfsi)

# Place Li cations near the Mg sites in the MOF
# Li should be lined up near the Nitrogen on the anion for proper bonding
# Base position based on the Mg site and Nitrogen

# Sourround the MOF with DME solvent (CH3OCH3)
dme = molecule('CH3OCH3')
n_solvent = 20  # Number of molecules to add
min_dist = 2.0  # Minimum distance in Angstroms to avoid overlap
# Define a region that solvent sould not appear
# Randomly generate a certain amount of solvent around the simulation box
# If the molecule is touching the MOF, generate a new positon
# Loop until all molecules of solvent are added
added = 0
while added < n_solvent:
    # Generate random position within cell boundaries
    cell_dims = mof.get_cell().lengths()
    rand_pos = np.random.rand(3) * cell_dims

# Perform an energy minimization to relax the solution
# Output the relaxed system to an xyz file for analysis
write('system_initial.xyz', mof)
print(f"Successfully added {added} solvent molecules.")

