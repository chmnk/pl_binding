import sys, getpass, os
import urllib, json
import urllib.request

username = getpass.getuser()
sys.path.append('/home/{0}/pylig'.format(username))


from PDB.rcsb_downloader import pdbe_assembly_download, validation_download
from databases.db_comparison.data_set_combined_paths import *


class UTILS:
    @staticmethod
    def try_mkdir(dir_path):
        try:
            os.mkdir(dir_path)
        except FileExistsError:
            pass

    @staticmethod
    def read_ids(fname_from):
        data = []
        with open(fname_from) as F:
            for ln in F:
                data.append(ln.replace('\n', ''))
        return data

    @staticmethod
    def write_ids(data, fname_to):
        F = open(fname_to, 'w')
        for pdbcode in data:
            F.write(pdbcode + '\n')
        F.close()

    @staticmethod
    def read_counter(fname):
        cntr = {}
        with open(fname) as F:
            for ln_ in F:
                ln = ln_.replace('\n', '').split()
                cntr[ln[0]] = int(ln[1])
        return cntr

    @staticmethod
    def write_sorted_dict_simple(data, fname_to):
        F = open(fname_to, 'w')
        for k in sorted(list(data.keys())):
            F.write('{0} {1}\n'.format(k, data[k]))
        F.close()

    @staticmethod
    def read_sorted_dict_simple(fname_from):
        data = {}
        with open(fname_from) as F:
            for ln_ in F:
                ln = ln_.replace('\n', '').split()
                data[ln[0]] = ln[1]
        return data

    @staticmethod
    def json_from_url(url, verbose=True):
        try:
            response = urllib.request.urlopen(url).read().decode()
            return json.loads(response)
        except urllib.error.HTTPError:
            if verbose:
                print('could not open', url)
            return None

    @staticmethod
    def get_rscc_from_validation_info(pdbcode, ligcode, chain_id, sa_id, res_id,
                                      no_chain_id=False, no_sa_id=False, no_res_id=False):
        if not os.path.exists('{0}/pdb-evaluations/{1}.xml'.format(PATH_main, pdbcode)):
            return None
        with open('{0}/pdb-evaluations/{1}.xml'.format(PATH_main, pdbcode), 'r') as F:
            for ln_ in F:
                ln = ln_.lstrip('\t')
                if ln.startswith('<ModelledSubgroup '):
                    l_id = ln.split(' resname="')[1].split('"')[0]
                    c_id = ln.split(' chain="')[1].split('"')[0]
                    r_id = int(ln.split(' resnum="')[1].split('"')[0])
                    s_id = ln.split(' said="')[1].split('"')[0]
                    if l_id == ligcode and (chain_id == c_id or no_chain_id) and (r_id == res_id or no_res_id) and (s_id == sa_id or no_sa_id):
                        if 'rscc' not in ln:
                            return 'no density'
                        rscc_id = float(ln.split(' rscc="')[1].split('"')[0])
                        return rscc_id
        return None

    @staticmethod
    def get_all_outliers_validation_info(pdbcode, threshold=10.):
        if not os.path.exists('{0}/pdb-evaluations/{1}.xml'.format(PATH_main, pdbcode)):
            return None
        outliers = {}
        with open('{0}/pdb-evaluations/{1}.xml'.format(PATH_main, pdbcode), 'r') as F:
            mg_open = False
            mg_line = ''
            for ln_ in F:
                ln = ln_.lstrip('\t')
                if ln.startswith('<ModelledSubgroup '):
                    mg_open = True
                    mg_line = ln
                elif ln.startswith('</ModelledSubgroup'):
                    mg_open = False
                    mg_line = ''
                elif mg_open and ln.startswith('<mog-angle-outlier'):
                    l_id = mg_line.split(' resname="')[1].split('"')[0]
                    c_id = mg_line.split(' chain="')[1].split('"')[0]
                    r_id = int(mg_line.split(' resnum="')[1].split('"')[0])
                    s_id = mg_line.split(' said="')[1].split('"')[0]
                    k = tuple([l_id, c_id, r_id, s_id])
                    zscore = float(mg_line.split(' Zscore="')[1].split('"')[0])
                    if abs(zscore) > threshold:
                        if k not in outliers:
                            outliers[k] = []
                        outliers[k].append([zscore, ln])
        return outliers

    @staticmethod
    def get_outliers_validation_info(pdbcode, ligcode, chain_id, sa_id, res_id, threshold=10.,
                                      no_chain_id=False, no_sa_id=False, no_res_id=False):
        if not os.path.exists('{0}/pdb-evaluations/{1}.xml'.format(PATH_main, pdbcode)):
            return None
        outliers = {}
        with open('{0}/pdb-evaluations/{1}.xml'.format(PATH_main, pdbcode), 'r') as F:
            mg_open = False
            mg_line = ''
            for ln_ in F:
                ln = ln_.lstrip('\t')
                if ln.startswith('<ModelledSubgroup '):
                    mg_open = True
                    mg_line = ln
                elif ln.startswith('</ModelledSubgroup'):
                    mg_open = False
                    mg_line = ''
                elif mg_open and ln.startswith('<mog-bond-outlier'):
                    l_id = mg_line.split(' resname="')[1].split('"')[0]
                    c_id = mg_line.split(' chain="')[1].split('"')[0]
                    r_id = int(mg_line.split(' resnum="')[1].split('"')[0])
                    s_id = mg_line.split(' said="')[1].split('"')[0]
                    if l_id == ligcode and (chain_id == c_id or no_chain_id) and (r_id == res_id or no_res_id) and (s_id == sa_id or no_sa_id):
                        k = tuple([l_id, c_id, r_id, s_id])
                        zscore = float(ln.split('Zscore="')[1].split('"')[0])
                        if abs(zscore) > threshold:
                            if k not in outliers:
                                outliers[k] = []
                            outliers[k].append([zscore, ln])
        return outliers

    @staticmethod
    def get_sa_id_validation_info(pdbcode, ligcode, chain_id, res_id):
        if not os.path.exists('{0}/pdb-evaluations/{1}.xml'.format(PATH_main, pdbcode)):
            return None
        with open('{0}/pdb-evaluations/{1}.xml'.format(PATH_main, pdbcode), 'r') as F:
            for ln_ in F:
                ln = ln_.lstrip('\t')
                if ln.startswith('<ModelledSubgroup '):
                    l_id = ln.split(' resname="')[1].split('"')[0]
                    c_id = ln.split(' chain="')[1].split('"')[0]
                    r_id = int(ln.split(' resnum="')[1].split('"')[0])
                    s_id = ln.split(' said="')[1].split('"')[0]
                    if l_id == ligcode and chain_id == c_id and r_id == res_id:
                        print(l_id, c_id, r_id, s_id)
                        return s_id

    @staticmethod
    def get_sa_id_lig_pdbe(pdbcode, ligcode, chain_id, res_id):
        ligs_data = UTILS.json_from_url('https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{0}'.format(pdbcode))[pdbcode]
        for lig_entry in ligs_data:
            l_id = lig_entry['chem_comp_id']
            c_id = lig_entry['chain_id']
            r_id = lig_entry['author_residue_number']
            sa_id = lig_entry['struct_asym_id']

            if l_id == ligcode and chain_id == c_id and r_id == res_id:
                return sa_id


    @staticmethod
    def get_sa_id_resi_pdbe(pdbcode, ligcode, chain_id, res_id):
        resi_data = UTILS.json_from_url('https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/{0}'.format(pdbcode))[pdbcode]
        for mol_data in resi_data['molecules']:
            for ch in mol_data['chains']:
                c_id = ch['chain_id']
                sa_id = ch['struct_asym_id']
                if c_id != chain_id:
                    continue
                for res_entry in ch['residues']:
                    l_id = res_entry['residue_name']
                    r_id = res_entry['author_residue_number']

                    if l_id == ligcode and chain_id == c_id and r_id == res_id:
                        return sa_id