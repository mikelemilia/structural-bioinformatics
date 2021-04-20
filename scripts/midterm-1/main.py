import matplotlib.colors as mplcolors
import math

import numpy as np

from functions import *

# Protein ID and data folder in which save data
pdb_id = '3lye'
# pdb_id = '5lxe'
data_folder = 'data/midterm-1'
output_folder = 'output/midterm-1'

# Retrieve the protein from PDB
pdbl = PDBList()
pdbl.retrieve_pdb_file(pdb_id, pdir=data_folder, file_format='pdb')
structure = PDBParser(QUIET=True).get_structure(pdb_id, '{}/pdb{}.ent'.format(data_folder, pdb_id))

# Save a PDB file containing only a list of residues with id " "
io = PDBIO()
io.set_structure(structure[0])
io.save("{}/pdb{}_cut.ent".format(data_folder, pdb_id), select=Select())

pdb_id = '{}_cut'.format(pdb_id)
structure = PDBParser(QUIET=True).get_structure(pdb_id, '{}/pdb{}.ent'.format(data_folder, pdb_id))

# Generation of distance matrix and plot of heatmap
sequence_separation = [0, np.inf]
dist_matrix = get_distance_matrix(structure[0]['A'], sequence_separation)

plot_heatmap(dist_matrix, pdb_id, output_folder)



# DISTANCE MATRIX COMPARISON

dist_matrix_prof = get_distance_prof(structure[0]['A'])

dist_matrix_prof = np.asmatrix(dist_matrix_prof)
dist_matrix = np.asmatrix(dist_matrix)

# with open('distance_outfile.txt', 'wb') as f:
#     for i in range(dist_matrix.shape[0]):
#         np.savetxt(f, dist_matrix_prof[i, :] - dist_matrix[i, :], fmt='%.2f')

# CONTACT MAP COMPARISON

contact_map_prof = np.asmatrix((dist_matrix_prof[:] < 5).astype(float))
contact_map = np.asmatrix((dist_matrix[:] < 5).astype(float))

# with open('contact_outfile.txt', 'wb') as f:
#     for i in range(contact_map.shape[0]):
#         np.savetxt(f, contact_map_prof[i, :] - contact_map[i, :], fmt='%.2f')

# CONTACT MAP [7, 12] COMPARISON

dist_matrix_7 = get_distance_matrix(structure[0]['A'], [7, 12])
contact_map_7 = np.asmatrix((dist_matrix_7[:] < 5).astype(float))
contact_map_prof_7 = np.asmatrix((np.triu(contact_map_prof, 7) - np.triu(contact_map_prof, 13)))

# with open('contact_7_outfile.txt', 'wb') as f:
#     for i in range(contact_map_7.shape[0]):
#         np.savetxt(f, contact_map_prof_7[i, :] - contact_map_7[i, :], fmt='%.2f')




# Count the number of residues for different sequence separation ranges
sequence_separations = [[0, 6], [7, 12], [13, 24], [25, np.inf]]

for sequence_separation in sequence_separations:
    dist_matrix = get_distance_matrix(structure[0]['A'], sequence_separation)
    contact_map = (dist_matrix[:] < 5).astype(float)
    print("Number of residues for range {}: {}".format(sequence_separation, int(np.sum(np.triu(contact_map[:])))))
    if sequence_separation[0] == 0:
        print(dist_matrix)
        print("of them, {} are on the diagonal".format(contact_map.shape[0]))

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
fig, ax = plt.subplots(figsize=(12, 12))
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
plt.savefig('{}/ramachandran_regions_{}.png'.format(output_folder, pdb_id), bbox_inches='tight')
