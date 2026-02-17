# Dependencies (ase)
import numpy as np
from ase import Atom
from ase.io import read, write
from ase.build import molecule
from ase.visualize import view
from ase.geometry import get_distances

# Add a MOF into the simulation through reading the .cif file
mof = read('MgMOF-74.cif')
mof.rotate(29, 'z')

# Find the center of the MOF's pores in order to add the TFSI anion
for atom in mof:
    print(atom.index, atom.symbol, atom.position)

# Use 2 atoms (Mg #1 and Mg #10) as a reference for center coordinates
pos_mg1 = mof[1].position
pos_mg10 = mof[10].position
pos_mg14 = mof[14].position
pos_mg5 = mof[5].position
center_points = [0.5 * (pos_mg1 + pos_mg10),
                 0.5 * (pos_mg14 + pos_mg5)]

# Read TFSI anion from the ntf2_pack.xyz file
tfsi = read('ntf2_pack.xyz')

tfsi.translate(-tfsi.get_center_of_mass())

offsets = [np.array([0.0, 0.0, -1.0]),
           np.array([0.0, 0.0, 1.0])] # Offsets in order to center the anion

rotations = [np.array([120, 0, -25]), 
             np.array([0, -20, 65])]

# Add the TFSI anion inside the MOF pores (2 anions)
for i in range(2):
    tfsi_copy = tfsi.copy()
    tfsi_copy.rotate(rotations[i][2], 'z') # Rotate accordingly
    tfsi_copy.rotate(rotations[i][1], 'y')
    tfsi_copy.rotate(rotations[i][0], 'x')
    tfsi_copy.translate(center_points[i] + offsets[i])
    mof.extend(tfsi_copy)

# Place Li cations near the Mg sites in the MOF
# Li should be lined up near the Nitrogen on the anion for proper bonding
# Base position based on the Mg site and Nitrogen (Mg #1, Mg #4 and N #167, N #182)
pos_n167 = mof[167].position
pos_n182 = mof[182].position
pos_mg4 = mof[4].position

ion_offsets = [np.array([0.5, -1.0, 0.0]),
               np.array([-0.25, 1.0, 0.0])]

ion_positions = [0.5 * (pos_mg1 + pos_n167),
                 0.5 * (pos_mg4 + pos_n182)]

Li1 = Atom("Li", position=(ion_positions[0] + ion_offsets[0]), charge=1)
Li2 = Atom("Li", position=(ion_positions[1] + ion_offsets[1]), charge=1)
mof.extend(Li1)
mof.extend(Li2)

# Sourround the MOF with DME solvent (CH3OCH3)
dme = molecule('CH3OCH3')
dme.translate(-dme.get_center_of_mass())

#Define region at which solvent can be placed
rangeX = [-22, 22]
rangeY = [-1.5, 27]
rangeZ = [-5.5, 12.5]

#define minimum distance (angstroms) and # of molecules as well as number of attempts before giving up
n_solvent = 20
min_dist = 1.2
attempts = 5

#For each molecule
#for i in range(n_solvent):
#    for x in range(attempts):
#        #Place at random position and rotation within range
#        #If it is too close, retry (up to 5 times)

# Perform an energy minimization to relax the solution
# Output the relaxed system to an xyz file for analysis
#write('system_initial.xyz', mof)
#print(f"Successfully added {added} solvent molecules.")
view(mof)
