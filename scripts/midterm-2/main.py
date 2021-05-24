from Bio.PDB import *
from Bio.SeqUtils import seq1

from functions import *

if __name__ == "__main__":

    # IUPRED prediction
    pdb_ids = [('2loj', 'A'), ('2m1s', 'A'), ('2kkw', 'A')]

    pdbl = PDBList()
    for pdb_id, chain_id in pdb_ids:
        pdbl.retrieve_pdb_file(pdb_id, pdir='data/midterm-2', file_format='pdb')

        structure = PDBParser(QUIET=True).get_structure(pdb_id, "data/midterm-2/pdb{}.ent".format(pdb_id))
        residues = [residue for residue in structure[0][chain_id] if residue.id[0] == " "]

        dist_matrix = get_distance_matrix(residues)
        contact_map = get_contact_map(residues)

        plot_contact_map(contact_map, 'output/midterm-2', pdb_id, chain_id)

        seq = "".join([seq1(residue.get_resname()) for residue in residues])

        # Punto 2

        exact_energy = []
        exact_energy_smooth = []
        m_mat = np.matrix(m_matrix)
        for i in range(len(seq)):

            energy_temp = 0

            contact = contact_map[i, :]
            idx = (np.array(contact) == 1)

            interaction_i = [seq[j] for j in range(len(idx)) if idx[j] == True]

            k = aa_list.index(seq[i])
            for elem in interaction_i:
                h = aa_list.index(elem)

                energy_temp += m_mat[k, h]

            exact_energy.append(energy_temp)


        # Smooth the prediction (moving average)
        for i in range(len(exact_energy)):
            frag = exact_energy[max(0, i - 10): min(i + 10, len(exact_energy))]
            exact_energy_smooth.append(sum(frag) / len(frag))


        fig, ax = plt.subplots(figsize=(12, 6))
        ax.set_title("{}_{}".format(pdb_id, chain_id))
        ax.axhline()
        ax.plot(np.arange(len(seq)), exact_energy, ls='--')
        ax.plot(np.arange(len(seq)), exact_energy_smooth, ls='-')

        plt.tight_layout()  # Remove figure padding
        plt.savefig('output/midterm-2/mypred_{}_{}.png'.format(pdb_id, chain_id), bbox_inches='tight')

        # Punto 3
        fig, ax = plt.subplots(figsize=(12, 6))
        pred, pred_smooth = iupred(seq)  # prediction
        ax.set_title("{}_{}".format(pdb_id, chain_id))
        ax.axhline()
        ax.plot(np.arange(len(seq)), pred, ls='--')
        ax.plot(np.arange(len(seq)), pred_smooth, ls='-')

        plt.tight_layout()  # Remove figure padding
        plt.savefig('output/midterm-2/iupred_{}_{}.png'.format(pdb_id, chain_id), bbox_inches='tight')

        # Punto 4
        count = np.sum(np.array(pred_smooth) >= 0, axis=0)
        total = len(pred_smooth)
        ratio = count / total
        print(count, total, ratio)
        print(pred_smooth)


        count2 = np.sum(np.array(exact_energy_smooth) >= 0, axis=0)
        total2 = len(exact_energy_smooth)
        ratio2 = count2 / total2
        print(count2, total2, ratio2)
        print(exact_energy_smooth)

