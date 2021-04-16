import sys, getpass, os

username = getpass.getuser()
sys.path.append('/home/{0}/pylig'.format(username))


from PDB.rcsb_downloader import pdbe_assembly_download, validation_download
from databases.db_comparison.data_set_combined_paths import *
from databases.db_comparison.data_utils import *


def download_with_cofactors():
    F = open('/home/maria/data/combined_dataset/cofactor_proteins_cur.txt', 'r')
    pdbcodes = F.readlines()
    F.close()
    pdbcodes = [p[:4] for p in pdbcodes]
    pdbe_assembly_download(pdbcodes, PATH_assemblies)


def download_cf_validation():
    existing_vr = [f[:4] for f in os.listdir('/home/maria/data/combined_dataset/pdb-evaluations/')]
    F = open('/home/maria/data/combined_dataset/cofactor_proteins_cur.txt', 'r')
    pdbcodes = F.readlines()
    F.close()
    pdbcodes = [p[:4] for p in pdbcodes]
    pdbcodes = [p[:4] for p in pdbcodes if p not in existing_vr]
    # print(len(pdbcodes))
    validation_download(ids=pdbcodes, path_to='/home/maria/data/combined_dataset/pdb-evaluations/')


def download_bm_validation():
    existing_vr = [f[:4] for f in os.listdir('/home/maria/data/combined_dataset/pdb-evaluations/')]
    F = open('{0}pdbcodes_bindingmoad.txt'.format(PATH_main), 'r')
    pdbcodes = F.readlines()
    F.close()
    pdbcodes = [p[:4] for p in pdbcodes if p not in existing_vr]
    # print(len(pdbcodes))
    validation_download(ids=pdbcodes, path_to='/home/maria/data/combined_dataset/pdb-evaluations/')


def download_bm_assembly():
    F = open('{0}pdbcodes_bindingmoad.txt'.format(PATH_main), 'r')
    pdbcodes = F.readlines()
    F.close()
    pdbcodes = [p[:4] for p in pdbcodes]
    pdbe_assembly_download(pdbcodes, PATH_assemblies)


def download_mr_validation():
    existing_vr = [f[:4] for f in os.listdir('/home/maria/data/combined_dataset/pdb-evaluations/')]
    F = open('{0}has_mod_res_lig.txt'.format(PATH_main), 'r')
    pdbcodes = F.readlines()
    F.close()
    pdbcodes = [p[:4] for p in pdbcodes]
    pdbcodes = [p[:4] for p in pdbcodes if p not in existing_vr]
    # print(len(pdbcodes))
    validation_download(ids=pdbcodes, path_to='/home/maria/data/combined_dataset/pdb-evaluations/')


def download_mr_assembly():
    F = open('{0}has_mod_res_lig.txt'.format(PATH_main), 'r')
    pdbcodes = F.readlines()
    F.close()
    pdbcodes = [p[:4] for p in pdbcodes]
    pdbe_assembly_download(pdbcodes, PATH_assemblies)


def download_mod_res_info():
    i_ln = 0

    ok_structures = []
    with open('{0}all_pdb_ok_resolution.txt'.format(PATH_main), 'r') as F:
        for ln in F:
            ok_structures.append(ln.replace('\n', ''))

    # check alternate_conformers!
    F = open('{0}has_mod_res.txt'.format(PATH_main), 'w')
    has_mod_res = {}
    for i_p, p in enumerate(ok_structures):
        if i_p % 1000 == 0:
            print(i_p)
        mod_res_data = UTILS.json_from_url('https://www.ebi.ac.uk/pdbe/api/pdb/entry/modified_AA_or_NA/{0}'.format(p),
                                     verbose=False)
        if mod_res_data is None:
            continue
        #     print(p, mod_res_data)
        has_mod_res[p] = []
        for lig_entry in mod_res_data[p]:
            ligcode = lig_entry['chem_comp_id']
            c_id = lig_entry['chain_id']
            res_id = lig_entry['author_residue_number']
            sa_id = lig_entry['struct_asym_id']
            F.write('{0} {1} {2} {3} {4}\n'.format(p, ligcode, c_id, res_id, sa_id))
            has_mod_res[p].append([ligcode, c_id, res_id, sa_id])
            print(p, has_mod_res[p][-1])
    F.close()


download_bm_validation()
download_bm_assembly()