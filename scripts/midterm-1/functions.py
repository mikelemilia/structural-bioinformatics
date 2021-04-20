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


def get_distance_prof(residues):
    # Calculate the distance matrix
    distances = []
    for residue1 in residues:
        if residue1.id[0] == " " and residue1.has_id("CA"):  # Exclude hetero/water residues
            row = []
            for residue2 in residues:
                if residue2.id[0] == " " and residue2.has_id("CA"):  # Exclude hetero/water residues
                    atom1 = residue1["CB"] if residue1.has_id('CB') else residue1["CA"]
                    atom2 = residue2["CB"] if residue2.has_id('CB') else residue2["CA"]
                    row.append(atom1 - atom2)
            distances.append(row)
    return np.array(distances, dtype=float)


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
        if found is False:
            atoms.append(residue["CA"])

    for i in range(len(list(residues))):
        row = []
        for j in range(len(list(residues))):
            if seq_sep[0] <= abs(list(residues)[i].id[1] - list(residues)[j].id[1]) <= seq_sep[1]:
                row.append(atoms[i] - atoms[j])
                if seq_sep[0] == 7:
                    print(abs(list(residues)[i].id[1] - list(residues)[j].id[1]))
            else:
                row.append(None)
        distances.append(row)

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
