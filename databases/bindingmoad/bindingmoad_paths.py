import os, getpass, socket

username = getpass.getuser()

hostname = socket.gethostname()

PATH_data = '/home/{0}/data/'.format(username)
if 'inria' in hostname or 'mistis' in hostname or 'node' in hostname:
    PATH_data = '/home/mkaducov/chmnk/data/'
PATH_bindingmoad = '{0}bindingmoad/'.format(PATH_data)
PATH_bindingmoad_cur = '{0}bindingmoad/2020/'.format(PATH_data)
# PATH_index_template = '{0}{1}/index/INDEX_general_PL_data.{1}'.format(PATH_pdbbind, '{0}')
# PATH_index_refined_template = '{0}{1}/index/INDEX_refined_data.{1}'.format(PATH_pdbbind, '{0}')
# PATH_index_detailed_template = '{0}{1}/index/INDEX_general_PL.{1}'.format(PATH_pdbbind, '{0}')
# PATH_pdbbind_general_template = '{0}{1}/general/'.format(PATH_pdbbind, '{0}')
# PATH_pdbbind_general_nocasf_template = '{0}{1}/general-no2013-no2016/'.format(PATH_pdbbind, '{0}')
# PATH_pdbbind_general_withcasf_template = '{0}{1}/general-with-coreset/'.format(PATH_pdbbind, '{0}')
# PATH_protlig_stat_template = '{0}{1}/protlig_stat/'.format(PATH_pdbbind, '{0}')


# def get_pdbbind_general(pdbbind_version='2016'):
#     PATH_general = PATH_pdbbind_general_template.format(pdbbind_version)
#     return PATH_general
#
#
# def get_pdbbind_paths(pdbbind_version='2016'):
#     global PATH_index, PATH_index_refined, PATH_index_detailed_template
#
#     PATH_index = PATH_index_template.format(pdbbind_version)
#     PATH_index_refined = PATH_index_refined_template.format(pdbbind_version)
#     PATH_index_detailed = PATH_index_detailed_template.format(pdbbind_version)
#     PATH_pdbbind_general = PATH_pdbbind_general_template.format(pdbbind_version)
#     return PATH_index, PATH_index_refined, PATH_index_detailed, PATH_pdbbind_general
#
#
# def get_pdbbind_nocasf_path(pdbbind_version='2016'):
#     PATH_pdbbind_general = PATH_pdbbind_general_nocasf_template.format(pdbbind_version)
#     return PATH_pdbbind_general
#
#
# def get_pdbbind_withcasf_path(pdbbind_version='2016'):
#     PATH_pdbbind_general = PATH_pdbbind_general_withcasf_template.format(pdbbind_version)
#     return PATH_pdbbind_general
#
#
# def get_protlig_stat_path(pdbbind_version='2016'):
#     PATH_protlig_stat = PATH_protlig_stat_template.format(pdbbind_version)
#     return PATH_protlig_stat
#
#
# def get_pdbbind_root(pdbbind_version='2016'):
#     return PATH_pdbbind + pdbbind_version + '/'