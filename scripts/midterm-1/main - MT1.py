import copy
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.colors as mplcolors
import math
from functions import *

# Protein ID and data folder in which save data
pdb_id = '5lxe'
# pdb_id = '3lye'
data_folder = '../../data/midterm-1'

# Retrieve the protein from PDB
pdbl = PDBList()
pdbl.retrieve_pdb_file(pdb_id, pdir=data_folder, file_format='pdb')
structure = PDBParser(QUIET=True).get_structure(pdb_id, '{}/pdb{}.ent'.format(data_folder, pdb_id))

# Save a PDB file containing only a list of residues with id " "
io = PDBIO()
io.set_structure(structure[0])
io.save("{}/pdb{}_cut.ent".format(data_folder, pdb_id), select=Select())

pdb_id = '{}_cut'.format(pdb_id)
structure = PDBParser(QUIET=True).get_structure(pdb_id, '{}/pdb{}.ent'.format(data_folder,pdb_id))

# Generation of distance matrix and plot of heatmap
sequence_separation = [0, np.inf]
dist_matrix = get_distance_matrix(structure[0]['A'], sequence_separation)

plot_heatmap(dist_matrix, pdb_id, data_folder)

# Count the number of residues for different sequence separation ranges
sequence_separation = [0, 6]
dist_matrix = get_distance_matrix(structure[0]['A'], sequence_separation)
contact_map = (dist_matrix[:] < 5).astype(float)
print("Number of residues for range [0, 6]: ", int(np.sum(contact_map[:]) / 2))

sequence_separation = [7, 12]
dist_matrix = get_distance_matrix(structure[0]['A'], sequence_separation)
contact_map = (dist_matrix[:] < 5).astype(float)
count = np.sum(contact_map[:]) / 2
print("Number of residues for range [7, 12]: ", int(np.sum(contact_map[:]) / 2))

sequence_separation = [13, 24]
dist_matrix = get_distance_matrix(structure[0]['A'], sequence_separation)
contact_map = (dist_matrix[:] < 5).astype(float)
print("Number of residues for range [13, 24]: ", int(np.sum(contact_map[:]) / 2))

sequence_separation = [25, np.inf]
dist_matrix = get_distance_matrix(structure[0]['A'], sequence_separation)
contact_map = (dist_matrix[:] < 5).astype(float)
print("Number of residues for range [25, inf]: ", int(np.sum(contact_map[:]) / 2))

# Build the peptides (reveal structure holes) and Calculate PSI and PHI
ppb = PPBuilder()
rama = {}  # { chain : [[residue_1, ...], [phi_residue_1, ...], [psi_residue_2, ...] ] }
chain = structure[0]['A']
for pp in ppb.build_peptides(chain):
    phi_psi = pp.get_phi_psi_list()  # [(phi_residue_1, psi_residue_1), ...]

    for i, residue in enumerate(pp):
        if phi_psi[i][0] is not None and phi_psi[i][1] is not None:
            rama.setdefault(chain.id, [[], [], []])
            rama[chain.id][0].append(residue)
            rama[chain.id][1].append(math.degrees(phi_psi[i][0]))
            rama[chain.id][2].append(math.degrees(phi_psi[i][1]))

# Ramachandran regions: 2 = 90%, 1 = 60%
regions_matrix = []
with open("{}/ramachandran.dat".format(data_folder)) as f:
    for line in f:
        if line:
            regions_matrix.append([int(ele) for ele in line.strip().split()])

# Plot Ramachandran regions
# fig, ax = plt.subplots(figsize=(12, 12))
cmap = mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF'])

f, axes = plt.subplots(1, len(rama))
axes = np.array(axes).reshape(-1)
print("\nOutliers in the Ramachandran plot:")
for ax, chain_id in zip(axes, rama):

    # Plot regions (60 percentile & 90 percentile)
    im = ax.imshow(regions_matrix, cmap=cmap, extent=(-180, 180, -180, 180))

    # Plot points
    ax.scatter(rama[chain_id][1], rama[chain_id][2], s=3, alpha=0.5)

    ax.set_xlabel('phi')
    ax.set_ylabel('psi')

    # Count and print outliers
    for residue, phi, psi in zip(*rama[chain_id]):
        phi_col = int(phi) + 179  # column of the Ramachandran range matrix
        psi_row = -1 * int(psi) + 179  # row of the Ramachandran range matrix

        if regions_matrix[psi_row][phi_col] == 0:
            print(residue, phi, psi, phi_col, psi_row, regions_matrix[psi_row][phi_col])

plt.tight_layout()  # Remove figure padding
plt.savefig('data/ramachandran_regions.png', bbox_inches='tight')