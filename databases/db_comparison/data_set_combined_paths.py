import os, getpass, socket, sys

username = getpass.getuser()
hostname = socket.gethostname()

sys.path.append("/home/{0}/pylig".format(username))  # my source to get the relative path, add your source
from databases.bindingmoad.bindingmoad_paths import *
from databases.pdbbind.pdbbind_paths import *

PATH_index, PATH_index_refined, PATH_index_detailed, PATH_pdbbind_general = get_pdbbind_paths('2019')

PATH_CCD = '/home/{0}/data/pdbe_dict/components.cif'.format(username)
PATH_main = '/home/{0}/data/combined_dataset/'.format(username)
PATH_assemblies = '{0}/pdb-assemblies/'.format(PATH_main)
PATH_assembly_patt = '{0}/pdb-assemblies/{1}-assembly-1.cif'.format(PATH_main, '{0}')

PATH_co_example = '{0}co_examples/'.format(PATH_main)
PATH_pb19_bm20 = '{0}pb19_bm20/'.format(PATH_main)
PATH_pb19_bm20_mr_cf_v1 = '{0}pb19_bm20_mr_cf_v1/'.format(PATH_main)