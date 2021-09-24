import os


# def fasta_download():
#     with open(PL_DOC_PATH + "INDEX_general_PL_name.2015", 'r') as Fgen:
#         for line in Fgen:
#             if line.startswith('#'):
#                 continue
#             pdbcode = line.split()[0]
#             # downloading with curl
#             os.system("curl -o " + FASTA_PATH + "{0}.fasta http://www.rcsb.org/pdb/files/fasta.txt?structureIdList={1}".format(pdbcode, pdbcode.upper()))

# A function for downloading
#   pdb, pdb.gz,
#   PDBx/mmCIF,
#   PDBML/XML,
#   Biological Assemblies,
#   Structure Factors,
#   NMR Restraints


def struct_download(ids, path_to, tail='pdb'):
    if path_to[-1] == '/':
        path_to = path_to[:-1]
    for cur_id in ids:
        if not os.path.isfile('{0}/{1}.{2}'.format(path_to, cur_id.lower(), tail)):
            print(cur_id.lower())
            print("curl -o {0}/{1}.{2} https://files.rcsb.org/download/{1}.{2}".format(path_to, cur_id.lower(), tail))
            os.system(
                "curl -o {0}/{1}.{2} https://files.rcsb.org/download/{1}.{2}".format(path_to, cur_id.lower(), tail))


def fasta_download(ids, path_to):
    if path_to[-1] == '/':
        path_to = path_to[:-1]
    for cur_id in ids:
        os.system(
            "curl -o {0}/{1}.fasta 'http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=FASTA&compression=NO&structureId={1}'".format(
                path_to, cur_id.lower()))


def validation_download(ids, path_to):
    if not os.path.exists(path_to):
        os.mkdir(path_to)
    else:
        pass
    if path_to[-1] == '/':
        path_to = path_to[:-1]
    for iid, cur_id in enumerate(ids):
        if iid % 100 == 0:
            print(iid, cur_id)
        if os.path.exists("{0}/{1}.xml".format(path_to, cur_id.lower())):
            continue
        os.system(
            "curl -o {0}/{1}.xml.gz https://files.rcsb.org/pub/pdb/validation_reports/{2}/{1}/{1}_validation.xml.gz".format(
                path_to, cur_id.lower(),  cur_id.lower()[1:3]))
        os.chdir(path_to)
        os.system('gzip -d {0}.xml.gz'.format(cur_id.lower()))


def pdbe_assembly_download(ids, path_to):
    if not os.path.exists(path_to):
        os.mkdir(path_to)
    else:
        pass
    if path_to[-1] == '/':
        path_to = path_to[:-1]
    for iid, cur_id in enumerate(ids):
        if iid % 100 == 0:
            print(iid, cur_id)
        if os.path.exists("{0}/{1}-assembly-1.cif".format(path_to, cur_id.lower())):
            continue
        os.system(
            "curl -L -o {0}/{1}-assembly-1.cif.gz 'http://www.ebi.ac.uk/pdbe/static/entry/download/{1}-assembly-1.cif.gz'".format(
                path_to, cur_id.lower()))
        os.chdir(path_to)
        os.system('gzip -d {0}-assembly-1.cif.gz'.format(cur_id.lower()))


# p_ids = []
# out_path = '/media/Disk2/lab/d3r/FXR_30_similarity/rcsb_pdb/'
# struct_download(p_ids, out_path)
# out_path = '/media/Disk2/lab/d3r/FXR_30_similarity/rcsb_fasta/'
# fasta_download(p_ids, out_path)

# F = open('/media/hdd1/data/pdbbind/2016/index/INDEX_general_PL_data.2016')
# struct_download([''], '/media/hdd1/data/pdbbind/2016/cif/', tail='cif')

# pdbcodes = [p.lower() for p in os.listdir('/home/maria/data/pdbbind/2019/general/') if len(p) == 4]
# validation_download(ids=pdbcodes, path_to='/home/maria/data/pdbbind/pdb-evaluations/')

