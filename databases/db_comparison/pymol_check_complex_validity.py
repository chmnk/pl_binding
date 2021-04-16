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
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr

username = getpass.getuser()
sys.path.append('/home/{0}/pylig'.format(username))

from pymol import cmd

from databases.db_comparison.data_set_combined_paths import *


cofactor_ids = ['01A', '01K', '0AF', '0ET', '0HG', '0HH', '0UM', '0WD', '0XU', '0Y0', '0Y1', '0Y2', '18W', '1C4', '1CV', '1CZ', '1DG', '1HA', '1JO', '1JP', '1R4', '1TP', '1TY', '1U0', '1VU', '1XE', '1YJ', '29P', '2CP', '2MD', '2NE', '2TP', '2TY', '36A', '37H', '3AA', '3CD', '3CP', '3GC', '3H9', '3HC', '48T', '4AB', '4CA', '4CO', '4IK', '4LS', '4LU', '4YP', '5AU', '5GY', '62X', '6FA', '6HE', '6NR', '6V0', '76H', '76J', '76K', '76L', '76M', '7AP', '7HE', '8EF', '8EL', '8EO', '8FL', '8ID', '8JD', '8PA', '8Z2', 'A3D', 'ABY', 'ACO', 'AGQ', 'AHE', 'AMX', 'AP0', 'ASC', 'AT5', 'ATA', 'B12', 'BCA', 'BCO', 'BHS', 'BIO', 'BOB', 'BSJ', 'BTI', 'BTN', 'BYC', 'BYG', 'BYT', 'C2F', 'CA3', 'CA5', 'CA6', 'CA8', 'CAA', 'CAJ', 'CAO', 'CCH', 'CIC', 'CMC', 'CMX', 'CNC', 'CND', 'CO6', 'CO8', 'COA', 'COB', 'COD', 'COF', 'COH', 'COM', 'COO', 'COT', 'COW', 'COY', 'COZ', 'D7K', 'DCA', 'DCC', 'DDH', 'DG1', 'DHE', 'DN4', 'DPM', 'DTB', 'EAD', 'EEM', 'EN0', 'ENA', 'EPY', 'ESG', 'F43', 'FA8', 'FAA', 'FAB', 'FAD', 'FAE', 'FAM', 'FAO', 'FAS', 'FCG', 'FCX', 'FDA', 'FDE', 'FED', 'FFO', 'FMI', 'FMN', 'FNR', 'FNS', 'FON', 'FOZ', 'FRE', 'FSH', 'FYN', 'G27', 'GBI', 'GBP', 'GBX', 'GDN', 'GDS', 'GF5', 'GGC', 'GIP', 'GNB', 'GPR', 'GPS', 'GRA', 'GS8', 'GSB', 'GSF', 'GSH', 'GSM', 'GSN', 'GSO', 'GTB', 'GTD', 'GTS', 'GTX', 'GTY', 'GVX', 'H2B', 'H4B', 'HAG', 'HAS', 'HAX', 'HBI', 'HCC', 'HDD', 'HDE', 'HEA', 'HEB', 'HEC', 'HEM', 'HIF', 'HMG', 'HSC', 'HTL', 'HXC', 'IBG', 'ICY', 'IRF', 'ISW', 'JM2', 'JM5', 'JM7', 'K15', 'L9X', 'LEE', 'LNC', 'LPA', 'LPB', 'LZ6', 'M43', 'M6T', 'MCA', 'MCD', 'MCN', 'MDE', 'MDO', 'MGD', 'MH0', 'MLC', 'MNH', 'MNR', 'MPL', 'MQ7', 'MSS', 'MTE', 'MTQ', 'MTV', 'MYA', 'N01', 'N1T', 'N3T', 'NA0', 'NAD', 'NAE', 'NAI', 'NAJ', 'NAP', 'NAQ', 'NAX', 'NBD', 'NBP', 'NDC', 'NDE', 'NDO', 'NDP', 'NHD', 'NHM', 'NHQ', 'NHW', 'NMX', 'NOP', 'NPL', 'NPW', 'ODP', 'OXK', 'P1H', 'P2Q', 'P3Q', 'P5F', 'PAD', 'PCD', 'PDP', 'PLP', 'PLR', 'PMP', 'PNS', 'PP9', 'PQQ', 'PXP', 'PZP', 'R1T', 'RBF', 'RFL', 'RGE', 'S0N', 'S1T', 'SA8', 'SAD', 'SAE', 'SAH', 'SAM', 'SCA', 'SCD', 'SCO', 'SDX', 'SFD', 'SFG', 'SH0', 'SHT', 'SMM', 'SND', 'SOP', 'SRM', 'SX0', 'T1G', 'T5X', 'T6F', 'TAD', 'TAP', 'TC6', 'TD6', 'TD7', 'TD8', 'TD9', 'TDK', 'TDL', 'TDM', 'TDP', 'TDT', 'TDW', 'TGG', 'THD', 'THF', 'THG', 'THH', 'THV', 'THW', 'THY', 'TOQ', 'TP7', 'TP8', 'TPP', 'TPQ', 'TPU', 'TPW', 'TPZ', 'TQQ', 'TRQ', 'TS5', 'TT8', 'TXD', 'TXE', 'TXP', 'TXZ', 'TYQ', 'TYY', 'TZD', 'UAH', 'UQ1', 'UQ2', 'UQ5', 'UQ6', 'VWW', 'WCA', 'WSD', 'WWF', 'XAX', 'XP8', 'XP9', 'Y7Y', 'YNC', 'ZBF', 'ZEM', 'ZID', 'ZNH', 'ZOZ']

ignore_list_mr = ['6yfy',  # huge NMR
                  '1kqe', '2xdc', '5isw', '5isx',  # gramicidin
                  '1c0q', '1c0r', '1qd8', '1ghg', '1sho', '1ehi', '1aa5', '2dln', '1rrv', '5m2h', '5m2k', '5hnm', '4oak', '4mut', '4mur', '4mus', '4ecl', '1fvm'  # VANCOMYCIN
                  '1cw8',  # histone?..
                  '2dsr', '1gac',  # peptide,
                  '1go6', '1hhu',  # Balhimycin
                  '6kgx',  # GIANT structure
                 ]


def read_counter(fname):
    cntr = {}
    with open(fname) as F:
        for ln_ in F:
            ln = ln_.replace('\n', '').split()
            cntr[ln[0]] = int(ln[1])
    return cntr


def get_complex_info():
    """
    Parse the precomputed information on complex' cofactors, ligands and their quality.
    :return:
    """
    lig_info = {}

    with open('/home/maria/data/combined_dataset/lig_info.txt') as F:
        for ln_ in F:
            ln = ln_.replace('\n', '').replace('no density', 'no_density').split()
            pdbcode, ligcode, c_id, sa_id, res_id, is_carb, rscc_cur = ln
            if pdbcode not in lig_info.keys():
                lig_info[pdbcode] = {}
            if rscc_cur == 'None':
                rscc_cur = None
            elif rscc_cur != 'no_density':
                rscc_cur = float(rscc_cur)

            lig_info[pdbcode][(ligcode, c_id, sa_id, int(res_id))] = [ligcode, c_id, sa_id, int(res_id), bool(is_carb), rscc_cur]

    no_density = read_counter('/home/maria/data/combined_dataset/no_density.txt')
    ok_rscc = read_counter('/home/maria/data/combined_dataset/ok_rscc.txt')

    return sorted(lig_info.keys()), lig_info, ok_rscc, no_density


def get_complex_info_mr():
    """
    Parse the precomputed information on complex' modified residues, ligands and their quality.
    :return:
    """
    lig_info = {}
    mr_info = {}

    with open('/home/maria/data/combined_dataset/lig_info_mr.txt') as F:
        for ln_ in F:
            ln = ln_.replace('\n', '').replace('no density', 'no_density').split()
            pdbcode, ligcode, c_id, sa_id, res_id, is_carb, rscc_cur = ln
            if pdbcode not in lig_info.keys():
                lig_info[pdbcode] = {}
            if rscc_cur == 'None':
                rscc_cur = None
            elif rscc_cur != 'no_density':
                rscc_cur = float(rscc_cur)

            lig_info[pdbcode][(ligcode, c_id, sa_id, int(res_id))] = [ligcode, c_id, sa_id, int(res_id), bool(is_carb), rscc_cur]

    no_density = read_counter('/home/maria/data/combined_dataset/no_density_mr.txt')
    ok_rscc = read_counter('/home/maria/data/combined_dataset/ok_rscc_mr.txt')

    with open('/home/maria/data/combined_dataset/info_residues.txt') as F:
        for ln_ in F:
            ln = ln_.replace('\n', '').replace('no density', 'no_density').split()
            pdbcode, ligcode, c_id, sa_id, res_id, rscc_cur = ln
            if pdbcode not in mr_info.keys():
                mr_info[pdbcode] = {}
            if rscc_cur == 'None':
                rscc_cur = None
            elif rscc_cur != 'no_density':
                rscc_cur = float(rscc_cur)

            mr_info[pdbcode][(ligcode, c_id, sa_id, int(res_id))] = [ligcode, c_id, sa_id, int(res_id), rscc_cur]

    no_density_mr = read_counter('/home/maria/data/combined_dataset/no_density_residues.txt')
    ok_rscc_mr = read_counter('/home/maria/data/combined_dataset/ok_rscc_residues.txt')

    valid_keys = set(lig_info.keys()).intersection(set(mr_info.keys()))

    return sorted(list(valid_keys)), lig_info, ok_rscc, no_density, mr_info, ok_rscc_mr, no_density_mr


def check_dists_lig_cofactor():
    """
    Computes the minimum distance for each cofactor-ligand pair, saves it to a file.
    Cofactor is also considered a ligand if there are several cofactors.
    :return:
    """
    pdbcodes, lig_info, ok_rscc, no_density = get_complex_info()
    has_nice_pair = 0
    fname_out = '{0}pdb_lig_co_pairs.csv'.format(PATH_main)

    F_pairs = open(fname_out, 'w')
    F_pairs.write('pdbcode,ligcode1,c_id1,sa_id1,res_id1,ligcode2,c_id2,sa_id2,res_id2,dist_min\n')
    F_pairs.close()
    for i_pdb, pdbcode in enumerate(pdbcodes):
        print(i_pdb, pdbcode)
        cur_lig_crds = {}
        cmd.load('{0}{1}-assembly-1.cif'.format(PATH_assemblies, pdbcode), 'my_asm')
        cmd.remove('hydrogens')
        lig_info_cur = lig_info[pdbcode]
        for k_ in sorted(lig_info_cur.keys()):
            ligcode, c_id, sa_id, res_id, is_carb, rscc_cur = lig_info_cur[k_]
            if rscc_cur is not None and rscc_cur != 'no_density':
                if rscc_cur < 0.8:
                    continue
            # print('my_asm and resname {0} and segi {1} and resi {2}'.format(ligcode, sa_id, res_id))
            crd = cmd.get_coords('my_asm and resname {0} and segi {1} and resi {2}'.format(ligcode, sa_id, res_id))
            if crd is not None:
                cur_lig_crds[k_] = crd
            # print(cur_lig_crds[k_])
        cmd.delete('all')
        cur_ok_keys = sorted(list(cur_lig_crds.keys()))
        nice_pairs = []
        for i1, k1 in enumerate(cur_ok_keys):
            for i2, k2 in enumerate(cur_ok_keys):
                if i1 >= i2:
                    continue
                if k1[0] not in cofactor_ids and k2[0] not in cofactor_ids:
                    continue
                dst_min = np.min(cdist(cur_lig_crds[k1], cur_lig_crds[k2]))
                if dst_min <= 5.:
                    nice_pairs.append((k1, k2))
                F_pairs = open(fname_out, 'a')
                F_pairs.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}\n'.format(pdbcode,
                                                                                 k1[0], k1[1], k1[2], k1[3],
                                                                                 k2[0], k2[1], k2[2], k2[3],
                                                                                 dst_min))
                F_pairs.close()

        if len(nice_pairs) > 0:
            print(nice_pairs)
            has_nice_pair += 1

    print(len(pdbcodes))
    print(has_nice_pair)


def check_dists_lig_modres():
    """
    Computes the minimum distance for each modified residue-ligand pair, saves it to a file.
    :return:
    """
    pdbcodes, lig_info, ok_rscc, no_density, mr_info, ok_rscc_mr, no_density_mr = get_complex_info_mr()
    has_nice_pair = 0

    fname_out = '{0}pdb_lig_mr_pairs.csv'.format(PATH_main)

    F_pairs = open(fname_out, 'w')
    F_pairs.write('pdbcode,modres1,c_id1,sa_id1,res_id1,ligcode2,c_id2,sa_id2,res_id2,dist_min\n')
    F_pairs.close()

    for i_pdb, pdbcode in enumerate(pdbcodes):
        print(i_pdb, pdbcode)
        cur_lig_crds = {}
        cur_mr_crds = {}
        cmd.load('{0}{1}-assembly-1.cif'.format(PATH_assemblies, pdbcode), 'my_asm')
        cmd.remove('hydrogens')
        lig_info_cur = lig_info[pdbcode]
        for k_ in sorted(lig_info_cur.keys()):
            ligcode, c_id, sa_id, res_id, is_carb, rscc_cur = lig_info_cur[k_]
            if rscc_cur is not None and rscc_cur != 'no_density':
                if rscc_cur < 0.8:
                    continue
            # print('my_asm and resname {0} and segi {1} and resi {2}'.format(ligcode, sa_id, res_id))
            crd = cmd.get_coords('my_asm and resname {0} and segi {1} and resi {2}'.format(ligcode, sa_id, res_id))
            if crd is not None:
                cur_lig_crds[k_] = crd

        res_info_cur = mr_info[pdbcode]
        for k_ in sorted(res_info_cur.keys()):
            ligcode, c_id, sa_id, res_id, rscc_cur = res_info_cur[k_]
            if rscc_cur is not None and rscc_cur != 'no_density':
                if rscc_cur < 0.8:
                    continue
            crd = cmd.get_coords('my_asm and resname {0} and segi {1} and resi {2}'.format(ligcode, sa_id, res_id))
            if crd is not None:
                cur_mr_crds[k_] = crd

        cmd.delete('all')
        cur_ok_keys = sorted(list(cur_lig_crds.keys()))
        cur_ok_keys_mr = sorted(list(cur_mr_crds.keys()))
        nice_pairs = []
        for i1, k1 in enumerate(cur_ok_keys):
            for i2, k2 in enumerate(cur_ok_keys_mr):
                dst_min = np.min(cdist(cur_mr_crds[k2], cur_lig_crds[k1]))
                if dst_min <= 5.:
                    nice_pairs.append((k1, k2))
                F_pairs = open(fname_out, 'a')
                F_pairs.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}\n'.format(pdbcode,
                                                                                 k1[0], k1[1], k1[2], k1[3],
                                                                                 k2[0], k2[1], k2[2], k2[3],
                                                                                 dst_min))
                F_pairs.close()

        if len(nice_pairs) > 0:
            print(nice_pairs)
            has_nice_pair += 1

    print(len(pdbcodes))
    print(has_nice_pair)


def check_if_peptide_not_modres():
    mod_res_info = {}
    mod_res_good = {}
    # with open('{0}pdb_lig_mr_pairs_noncov_dmin5.2_30resmin.csv'.format(PATH_main)) as F:
    with open('{0}pdb_lig_mr_pairs.csv'.format(PATH_main)) as F:
        next(F)
        for ln_ in F:
            ln = ln_.replace('\n', '').split(',')
            pdbcode, ligcode, c_id_lig, sa_id_lig, resi_lig, resn, c_id_mr, sa_id_mr, resi_mr, dist = ln
            if pdbcode in ignore_list_mr:
                continue
            resi_lig, resi_mr, dist = int(resi_lig), int(resi_mr), float(dist)
            if pdbcode not in mod_res_info:
                mod_res_info[pdbcode] = {}
            if sa_id_lig not in mod_res_info[pdbcode]:
                mod_res_info[pdbcode][sa_id_lig] = []
            mod_res_info[pdbcode][sa_id_lig].append((ligcode, c_id_lig, sa_id_lig, resi_lig, resn, c_id_mr, sa_id_mr, resi_mr, dist))

    for pdbcode in sorted(list(mod_res_info.keys())):
        cmd.load('{0}{1}-assembly-1.cif'.format(PATH_assemblies, pdbcode), 'my_asm')
        cmd.remove('hydrogens')
        for sa_id_lig in sorted(list(mod_res_info[pdbcode].keys())):
            for entry in mod_res_info[pdbcode][sa_id_lig]:
                sa_id_mr = entry[6]
                dsts = [ele[-1] for ele in mod_res_info[pdbcode][sa_id_lig]]
                cmd.select('segi_cur', 'segi {0}'.format(sa_id_mr))
                cmd.select('has_ca', 'segi_cur and name CA')
                n_ca_segi = cmd.count_atoms('has_ca')
                if min(dsts) > 5.2:
                    continue
                if min(dsts) <= 1.8:  # covalent
                    continue
                    # print('probably covalent bonding:', min(dsts))
                    # print(pdbcode, n_ca_segi, "  ", sa_id, mod_res_info[pdbcode][sa_id])
                if n_ca_segi <= 30:
                    continue
                # if n_ca_segi <= 25:
                #     print(pdbcode, n_ca_segi, "  ", sa_id, mod_res_info[pdbcode][sa_id])
                cmd.delete('segi_cur')
                cmd.delete('has_ca')

                if entry[-1] > 5.2:
                    continue
                if pdbcode not in mod_res_good:
                    mod_res_good[pdbcode] = {}
                if sa_id_lig not in mod_res_good[pdbcode]:
                    mod_res_good[pdbcode][sa_id_lig] = []
                mod_res_good[pdbcode][sa_id_lig].append(entry)


        # for sa_id in sorted(list(mod_res_info[pdbcode].keys())):
        #     dsts = [ele[-1] for ele in mod_res_info[pdbcode][sa_id]]
        #     cmd.select('segi_cur', 'segi {0}'.format(sa_id))
        #     cmd.select('has_ca', 'segi_cur and name CA')
        #     n_ca_segi = cmd.count_atoms('has_ca')
        #     if min(dsts) > 5.2:
        #         continue
        #     if min(dsts) <= 1.8:  # covalent
        #         continue
        #         # print('probably covalent bonding:', min(dsts))
        #         # print(pdbcode, n_ca_segi, "  ", sa_id, mod_res_info[pdbcode][sa_id])
        #     if n_ca_segi <= 30:
        #         continue
        #     # if n_ca_segi <= 25:
        #     #     print(pdbcode, n_ca_segi, "  ", sa_id, mod_res_info[pdbcode][sa_id])
        #     cmd.delete('segi_cur')
        #     cmd.delete('has_ca')
        #     for entry in mod_res_info[pdbcode][sa_id]:
        #         if entry[-1] > 5.2:
        #             continue
        #         if pdbcode not in mod_res_good:
        #             mod_res_good[pdbcode] = {}
        #         if sa_id not in mod_res_good[pdbcode]:
        #             mod_res_good[pdbcode][sa_id] = []
        #         mod_res_good[pdbcode][sa_id].append(entry)
        cmd.delete('all')
    print(len(mod_res_info), len(mod_res_good))

    F = open('{0}pdb_lig_mr_pairs_noncov_dmin5.2_30resmin.csv'.format(PATH_main), 'w')
    for pdbcode in sorted(list(mod_res_good.keys())):
        for sa_id in sorted(list(mod_res_good[pdbcode].keys())):
            for entry in mod_res_good[pdbcode][sa_id]:
                F.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}\n'.format(pdbcode, *entry))
    F.close()



# check_dists_lig_modres()
check_if_peptide_not_modres()