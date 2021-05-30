from Bio.PDB import *
from Bio.SeqUtils import seq1

from functions import *

if __name__ == "__main__":

    # IUPRED prediction
    pdb_ids = [('2loj', 'A'), ('2m1s', 'A'), ('2kkw', 'A'), ('1tac', 'A')]

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

            # Index of the residue in position i
            k = aa_list.index(seq[i])
            energy_temp = 0

            # List of the contacts of i
            contact = contact_map[i, :]
            # DA DECOMMENTARE PER FAR VENIRE IL RISULTATO DEL PROF
            # for j in range(i):
            #     contact[j] = 0

            idx = np.array(contact) == 1
            interaction_i = [seq[j] for j in range(len(idx)) if idx[j] == True]

            for elem in interaction_i:
                # Index of the resiude in contact
                h = aa_list.index(elem)
                # Add the current contribution
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
        pred_energy, pred_energy_smooth = iupred(seq)  # prediction

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.set_title("{}_{}".format(pdb_id, chain_id))
        ax.axhline()
        ax.plot(np.arange(len(seq)), pred_energy, ls='--')
        ax.plot(np.arange(len(seq)), pred_energy_smooth, ls='-')

        plt.tight_layout()  # Remove figure padding
        plt.savefig('output/midterm-2/iupred_{}_{}.png'.format(pdb_id, chain_id), bbox_inches='tight')

        # Punto 4

        count_exact = np.sum(np.array(exact_energy_smooth) >= 0, axis=0)
        total_exact = len(exact_energy_smooth)
        ratio_exact = count_exact / total_exact
        print("Disorder content for exact energy: ", count_exact, total_exact, ratio_exact)
        # print(exact_energy_smooth)

        count_pred = np.sum(np.array(pred_energy_smooth) >= 0, axis=0)
        total_pred = len(pred_energy_smooth)
        ratio_pred = count_pred / total_pred
        print("Disorder content for predicted energy: ", count_pred, total_pred, ratio_pred)
        # print(pred_smooth)

        print("\n")
