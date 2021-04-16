import getpass, socket, sys

username = getpass.getuser()

hostname = socket.gethostname()

sys.path.append("/home/{0}/pylig".format(username))  # my source to get the relative path, add your source
print(sys.path)
from databases.bindingmoad.bindingmoad_paths import *

print(len(os.listdir(PATH_bindingmoad_cur)))
fls = [f for f in os.listdir(PATH_bindingmoad_cur) if '.bio' in f]
for f in fls:
    try:
        os.mkdir(PATH_bindingmoad_cur + f.split('.bio')[0])
    except FileExistsError:
        pass
    os.system('mv {0}{1} {0}{2}/{3}'.format(PATH_bindingmoad_cur, f, f.split('.bio')[0], f + '.pdb'))