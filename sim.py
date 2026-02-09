# Dependencies (ase)
import numpy as np
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
pos_mg0 = mof[0].position
pos_mg1 = mof[1].position
pos_mg9 = mof[9].position
center_point = 0.5 * (pos_mg1 + pos_mg9)

# Read TFSI anion from the ntf2_pack.xyz file
# Add the TFSI anion inside the MOF pores (2 anions)
# Rotate accoridngly
tfsi = read('ntf2_pack.xyz')

offsets = [
    np.array([0.0, 0.0, 0.0]),
    np.array([0.0, 0.0, 4.0])
]

for off in offsets:
    tfsi_copy = tfsi.copy()
    tfsi_copy.rotate(45, 'z')
    tfsi_copy.translate(center_point + off)
    mof.extend(tfsi_copy)

# Place Li cations near the Mg sites in the MOF
# Li should be lined up near the Nitrogen on the anion for proper bonding
# Base position based on the Mg site and Nitrogen

# Sourround the MOF with DME solvent (CH3OCH3)
dme = molecule('CH3OCH3')
dme.translate(-dme.get_center_of_mass())

n_solvent = 20
min_dist = 2.0
cell = mof.get_cell().lengths()

def too_close(mol, system):
    for a in mol:
        for b in system:
            if np.linalg.norm(a.position - b.position) < min_dist:
                return True
    return False

added = 0
attempts = 0

while added < n_solvent and attempts < 5000:
    attempts += 1

    dme_try = dme.copy()
    dme_try.rotate(np.random.uniform(0, 360), 'x')
    dme_try.rotate(np.random.uniform(0, 360), 'y')
    dme_try.rotate(np.random.uniform(0, 360), 'z')

    rand_pos = np.random.rand(3) * cell
    dme_try.translate(rand_pos)

    if not too_close(dme_try, mof):
        mof.extend(dme_try)
        added += 1
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
view(mof)


