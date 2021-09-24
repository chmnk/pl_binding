#!/usr/bin/python3
import sys
import getpass

username = getpass.getuser()
sys.path.append("/home/{0}/pylig".format(username))  # my source to get the relative path, add your source

from databases.pdbbind.pdbbind_paths import *

PATH_val = '/home/maria/data/pdbbind/pdb-evaluations/'
PATH_index, PATH_index_refined, PATH_index_detailed, PATH_pdbbind_general = get_pdbbind_paths('2019')


pdbcodes_refined = []
pdbcodes_general = []
with open(PATH_index_refined) as F:
    for ln_ in F:
        if ln_.startswith('#'):
            continue
        ln = ln_.split()
        if len(ln[0]) == 4:
            pdbcodes_refined.append(ln[0])
print(PATH_index_detailed)
with open(PATH_index_detailed) as F:
    for ln_ in F:
        # print(ln_)
        if ln_.startswith('#'):
            continue
        ln = ln_.split()
        if len(ln[0]) == 4:
            pdbcodes_general.append(ln[0])

outlier_data = {}
for pdbcode in pdbcodes_general:
    fname = '{0}{1}.xml'.format(PATH_val, pdbcode)
    if not os.path.exists(fname):
        print('not exists!')
        continue
    with open(fname) as F:
        for ln_ in F:
            if '\t\t<mog-bond-outlier' in ln_:
                ln = ln_.split('<mog-bond-outlier')[1].split()
                if float(ln[0].split('"')[1]) > 10.:
                    if fname not in outlier_data:
                        outlier_data[fname] = []
                    print(fname, ln)
    # bond-outlier Zscore="5.30"
print(len(outlier_data), len(pdbcodes_general))