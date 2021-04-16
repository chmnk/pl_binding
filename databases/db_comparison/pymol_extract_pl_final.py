import sys, getpass, os
import os, getpass, socket, sys
import urllib, json
import re
from collections import Counter
from scipy.spatial.distance import cdist
import numpy as np
import scipy
import re
from matplotlib.pyplot import *
import pickle
import copy
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

username = getpass.getuser()
sys.path.append('/home/{0}/pylig'.format(username))

from pymol import cmd

from databases.db_comparison.data_set_combined_paths import *
from databases.db_comparison.data_utils import *


def copy_pd_simple():
    pdbcodes = UTILS.read_ids('{0}pdbcodes_pb_good.txt'.format(PATH_main))
    casf_2013 = UTILS.read_ids('{0}casf_2013_codes.dat'.format(PATH_data))
    casf_2016 = UTILS.read_ids('{0}casf_2016_codes.dat'.format(PATH_data))
    for p in pdbcodes:
        if p not in casf_2013 and p not in casf_2016:
            os.system('cp -r {0}{1} {2}/'.format(PATH_pdbbind_general, p, PATH_pb19_bm20))


def copy_pd_list(pdbcodes):
    for p in pdbcodes:
        os.system('cp -r {0}{1} {2}/'.format(PATH_pb19_bm20, p, PATH_pb19_bm20_mr_cf_v1))


def extract_bindingmoad():
    pdbs_to_save = []  # each pdb file is unique (i.e. has 1 ligand) in this dataset
    casf_2013 = UTILS.read_ids('{0}casf_2013_codes.dat'.format(PATH_data))
    casf_2016 = UTILS.read_ids('{0}casf_2016_codes.dat'.format(PATH_data))
    with open('{0}pdbcodes_ligcodes_bm_not_pb_good.txt'.format(PATH_main)) as F:
        for ln_ in F:
            ln = ln_.replace('\n', '').split()
            n_resi = int(ln[1])
            r_id0 = int(ln[3 + n_resi])
            pdbs_to_save.append([ln[0], ln[2: 2 + n_resi], ln[2 + n_resi], [r_id0 + i for i in range(n_resi)]])
    for pdbcode, ligcodes, chain_id, res_ids in pdbs_to_save:
        if pdbcode in casf_2013:
            continue
        if pdbcode in casf_2016:
            continue
        sa_id = UTILS.get_sa_id_validation_info(pdbcode, ligcodes[0], chain_id, res_ids[0])
        if sa_id is None:
            sa_id = UTILS.get_sa_id_lig_pdbe(pdbcode, ligcodes[0], chain_id, res_ids[0])
            if sa_id is None:
                sa_id = UTILS.get_sa_id_resi_pdbe(pdbcode, ligcodes[0], chain_id, res_ids[0])
                if sa_id is None:
                    skip_resi = {'2q5a': -1, '1oby': -1, '1w9o': -1, '5t5j': -1, '5yd5': -1, '6rt6': -1, '6rt7': -1,
                                 '1mfl': -4,
                                 '1rsu': -2, '1w9e': -2, '1y7l': -2, '2j6o': -3, '3db3': -2, '3p36': -2, '3qzs': -2,
                                 '5f5k': -3, '5yd3': -2, '6hhp': -3, '6hmt': -3}
                    if pdbcode in skip_resi:
                        res_ids = [r_id + skip_resi[pdbcode] for r_id in res_ids]
                    sa_id = UTILS.get_sa_id_resi_pdbe(pdbcode, ligcodes[0], chain_id, res_ids[0])
                    if sa_id is None:
                        print(pdbcode, ligcodes, chain_id, res_ids)
                        print(sa_id)
        path_dir_out = '{0}{1}/'.format(PATH_pb19_bm20, pdbcode)
        try:
            os.mkdir(path_dir_out)
        except FileExistsError:
            pass

        cmd.load('{0}{1}-assembly-1.cif'.format(PATH_assemblies, pdbcode), 'my_com')
        for i, ligcode_res_id in enumerate(zip(ligcodes, res_ids)):
            ligcode, res_id = ligcode_res_id
            cmd.create('lig_{0}'.format(i), 'my_com and resn {0} and resi {1} and segi {2} and (alt A or alt \'\')'.format(ligcode, res_id, sa_id))
            cmd.remove('my_com and resn {0} and resi {1} and segi {2}'.format(ligcode, res_id, sa_id))
        cmd.create('all_li', 'lig_*')
        cmd.delete('lig_*')
        cmd.save('{0}{1}_protein.pdb'.format(path_dir_out, pdbcode), 'my_com')
        cmd.save('{0}{1}_ligand.mol2'.format(path_dir_out, pdbcode), 'all_li')
        cmd.delete('all')


load_pb = lambda pdbcode: [cmd.load('../pdbbind/2019/general/{0}/{0}_ligand.mol2'.format(pdbcode)), cmd.load('../pdbbind/2019/general/{0}/{0}_protein.pdb'.format(pdbcode)), cmd.origin(pdbcode + '_ligand')]
load_pa = lambda pdbcode: cmd.load('./pdb-assemblies/{0}-assembly-1.cif'.format(pdbcode))


def get_occupancies_pb19_bm20_mr_cf():
    mr_cf_lig = pickle.load(open("{0}pb19_bm20_mr_cf_v1_mr_cf_dct".format(PATH_main), "rb"))
    mod_res_final = pickle.load(open("{0}pb19_bm20_mr_cf_v1_mr_dct".format(PATH_main), "rb"))
    cf_res_final = pickle.load(open("{0}pb19_bm20_mr_cf_v1_cf_dct".format(PATH_main), "rb"))
    all_mr_cf = copy.deepcopy(mr_cf_lig)
    for p, e_p in mod_res_final.items():
        all_mr_cf[p] = {}
        for l, e_l in e_p.items():
            all_mr_cf[p][l] = {'mr': e_l}
    for p, e_p in cf_res_final.items():
        all_mr_cf[p] = {}
        for l, e_l in e_p.items():
            all_mr_cf[p][l] = {'cf': e_l}

    F = open('pb19_bm20_mr_cf_v1_occupancies.txt', 'w')

    for p, e_p in all_mr_cf.items():
        F.write(p + '\n')
        cmd.load(PATH_assembly_patt.format(p), 'prot')
        cmd.remove('resname HOH')
        for l, e_l in e_p.items():
            lig_name = '_'.join([l[0], l[2], str(l[3])])
            myspace = {'occ': []}
            cmd.select('lig', 'segi {0}'.format(l[2]))
            cmd.iterate('lig', 'occ.append(q)', space=myspace)
            cmd.delete('lig')
            occ_lig = np.mean(myspace['occ'])
            F.write('{0} {1} {2} {3} {4}\n'.format(*l, occ_lig))
            if occ_lig < 0.8:
                print(p, l, occ_lig)
            e_ls = []
            if 'mr' in e_l:
                e_ls += e_l['mr']
            if 'cf' in e_l:
                e_ls += e_l['cf']

            for lig1 in e_ls:
                myspace = {'occ': []}
                resi = str(lig1[3]) if type(lig1[3]) == int else ','.join(map(str, lig1[3]))
                cmd.select('lig1', 'segi {0} and resn {1} and resi {2}'.format(lig1[2], lig1[0].replace(',', '+'), resi))
                cmd.iterate('lig1', 'occ.append(q)', space=myspace)
                cmd.delete('lig1')
                occ_lig = np.mean(myspace['occ'])
                F.write('{0} {1} {2} {3} {4}\n'.format(*lig1, occ_lig))
                if occ_lig < 0.8:
                    print(p, l, lig1, occ_lig)
        cmd.delete('all')
    F.close()


def write_pb19_bm20_mr_cf():
    UTILS.try_mkdir(PATH_pb19_bm20_mr_cf_v1)
    allowed_pdbbind = pickle.load(open("{0}pb19_bm20_mr_cf_v1_pb_pdbcodes".format(PATH_main), "rb"))
    allowed_bm = pickle.load(open("{0}pb19_bm20_mr_cf_v1_bm_pdbcodes".format(PATH_main), "rb"))
    mr_cf_lig = pickle.load(open("{0}pb19_bm20_mr_cf_v1_mr_cf_dct".format(PATH_main), "rb"))
    mod_res_final = pickle.load(open("{0}pb19_bm20_mr_cf_v1_mr_dct".format(PATH_main), "rb"))
    cf_res_final = pickle.load(open("{0}pb19_bm20_mr_cf_v1_cf_dct".format(PATH_main), "rb"))
    all_mr_cf = copy.deepcopy(mr_cf_lig)
    for p, e_p in mod_res_final.items():
        all_mr_cf[p] = {}
        for l, e_l in e_p.items():
            all_mr_cf[p][l] = {'mr': e_l}
    for p, e_p in cf_res_final.items():
        all_mr_cf[p] = {}
        for l, e_l in e_p.items():
            all_mr_cf[p][l] = {'cf': e_l}

    # copy_pd_list(allowed_pdbbind)
    # copy_pd_list(allowed_bm)
    # cofactors_list = '' resname:resid:segid
    F = open('pb19_bm20_mr_cf_v1_mr_cf_info.txt', 'w')

    for p, e_p in all_mr_cf.items():
        print(p)
        cmd.load(PATH_assembly_patt.format(p), 'prot')
        cmd.remove('resname HOH')
        lig_names_cnt = Counter()
        for l, e_l in e_p.items():
            lig_name = '_'.join([l[0], l[2], str(l[3])]).replace(', ', ',')
            p_lig_name = '{0}_{1}'.format(p, lig_name)
            path_to = '{0}{1}/'.format(PATH_pb19_bm20_mr_cf_v1, p_lig_name)
            UTILS.try_mkdir(path_to)

            lig_names_cnt[lig_name] += 1
            if 'mr' in e_l:
                mr_tail = ','.join([':'.join([lig1[0], lig1[2], str(lig1[3])]) for lig1 in e_l['mr']])
                F.write('{0} mr {1}\n'.format(p_lig_name, mr_tail))
            if 'cf' in e_l:
                cf_tail = ','.join([':'.join([lig1[0], lig1[2], str(lig1[3])]) for lig1 in e_l['cf']])
                F.write('{0} cf {1}\n'.format(p_lig_name, cf_tail))

            myspace = {'c_ids': set([])}
            cmd.select('lig', 'segi {0}'.format(l[2]))
            print('lig atoms', cmd.count_atoms('lig'))
            cmd.select('around_lig', 'prot near_to 10 of lig')
            cmd.iterate('around_lig', 'c_ids.add(chain)', space=myspace)
            cmd.select('around_lig_all', 'chain ' + '+'.join(myspace['c_ids']))
            cmd.save('{0}{1}_protein.pdb'.format(path_to, p_lig_name), 'around_lig_all')
            cmd.save('{0}{1}_ligand.mol2'.format(path_to, p_lig_name), 'lig')
            cmd.delete('lig')
            cmd.delete('around_lig')
            cmd.delete('around_lig_all')

        if max(lig_names_cnt.values()) > 1:
            print(p, e_p)
        cmd.delete('all')
    F.close()




#TODO occupancy..... i.e. 3dna...

# extract_bindingmoad()
# copy_pd_simple()

write_pb19_bm20_mr_cf()
# get_occupancies_pb19_bm20_mr_cf()