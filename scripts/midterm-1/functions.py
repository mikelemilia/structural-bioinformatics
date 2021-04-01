import copy
from Bio.PDB import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


# Save a PDB chain
class Select(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "

    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"


def get_distance_matrix(residues, seq_sep):
    distances = []
    atoms = []

    for residue in residues:
        found = False
        for atom in residue.get_atoms():
            if atom.get_id() == 'CB':
                found = True
                atoms.append(atom)
                break
        if found == False: atoms.append(residue["CA"])

    i = 0
    for residue1 in residues:
        row = []
        j = 0
        for residue2 in residues:
            if seq_sep[0] <= abs(residue1.id[1] - residue2.id[1]) <= seq_sep[1]:
                row.append(atoms[i] - atoms[j])
            else:
                row.append(None)
            j += 1
        distances.append(row)
        i += 1

    return np.array(distances, dtype=float)


def plot_heatmap(dist_matrix, pdb_id, output_folder):
    current_cmap = copy.copy(matplotlib.cm.get_cmap())
    current_cmap.set_bad(color='white')

    fig, ax = plt.subplots(figsize=(12, 12))
    im = ax.imshow(dist_matrix)
    fig.colorbar(im, fraction=0.03, pad=0.05)
    plt.savefig('{}/ca_distances_{}.png'.format(output_folder, pdb_id), bbox_inches='tight')

    # Plot contact map
    contact_map = (dist_matrix[:] < 5).astype(
        float)  # Calculate the contact map based on a distance threshold 8 Angstrom
    fig, ax = plt.subplots(figsize=(24, 24))
    im = ax.imshow(contact_map)

    # Set ticks
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))

    plt.savefig('{}/ca_contacts_{}.png'.format(output_folder, pdb_id), bbox_inches='tight')
