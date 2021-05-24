import numpy as np
from Bio.PDB import NeighborSearch

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

from iupred_data import aa_list, p_matrix, m_matrix



def get_distance_matrix(residues):
    # Calculate the distance matrix
    distances = []
    for residue1 in residues:
        if residue1.id[0] == " " and residue1.has_id("CA"):  # Exclude hetero/water residues
            row = []
            for residue2 in residues:
                if residue2.id[0] == " " and residue2.has_id("CA"):  # Exclude hetero/water residues
                    row.append(residue1["CA"] - residue2["CA"])
            distances.append(row)

    np.array(distances).reshape(len(residues), -1)

    return distances


def get_contact_map(residues):
    contacts = [[0 for _ in residues] for _ in residues]
    ns = NeighborSearch([atom for residue in residues for atom in residue.get_atoms()])
    for residue1, residue2 in ns.search_all(3.5, level="R"):  # level="R" returns pairs of residues in contact considering all atoms
        if abs(residue1.id[1] - residue2.id[1]) >= 2:
            contacts[residues.index(residue1)][residues.index(residue2)] = 1
            contacts[residues.index(residue2)][residues.index(residue1)] = 1

    return np.array(contacts)

def iupred(seq, sequence_separation=2, window_size=100, window_size_smooth=10):
    '''
    Calculate residue IUPRED energy considering neighbouring residues (windows_size) and
    smoothing by window_size_smooth
    :param seq: a string of aminoacids
    :param sequence_separation: neighbours min distance
    :param window_size: neighbours max distance
    :param window_size_smooth: sliding average window size
    :return: row prediction, smoothed prediction
    '''

    pred = []
    pred_smooth = []

    p_mat = np.array(p_matrix)

    indices = [aa_list.index(aa) for aa in list(seq)]  # Transform sequence into indexes as in the P matrix
    for i, aa_index in enumerate(indices):

        # Get the slice i-100/i+100 excluding adjacent positions (+/-1)
        start_before = max(0, i - window_size)
        end_before = max(0, i - sequence_separation)
        start_after = min(len(indices) - 1, i + sequence_separation)
        end_after = min(len(indices) - 1, i + window_size)
        indices_local = indices[start_before: end_before] + indices[start_after: end_after]
        # print(i, aa_index, aa_list[aa_index], len(indices), len(indices_local), i, start_before, end_before, start_after, end_after)

        # Count amino acids in the window
        row = np.full((20,), 0)
        for index in indices_local:
            row[index] += 1
        # print(row)

        row = row / len(indices_local)  # calculate AA frequency
        # print(row)

        row = row * p_mat[aa_index, ]  # calculate energy
        # print(row)

        aa_energy = np.sum(row)
        # print(i, seq[i], aa_energy)

        pred.append(aa_energy)

    # Smooth the prediction (moving average)
    for i in range(len(pred)):
        frag = pred[max(0, i - window_size_smooth): min(i + window_size_smooth, len(pred))]
        pred_smooth.append(sum(frag) / len(frag))

    return pred, pred_smooth


def plot_contact_map(contact_map, output_folder, pdb_id, chain):
    fig, ax = plt.subplots(figsize=(24, 24))
    im = ax.imshow(contact_map)

    # Set ticks
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))

    plt.savefig('{}/ca_contacts_{}_{}.png'.format(output_folder, pdb_id, chain), bbox_inches='tight')
