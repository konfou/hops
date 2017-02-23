import os
import shutil
import platform

app = 'hops'
system = platform.system()

excecutable = {'Darwin': 'command', 'Linux': 'sh', 'Windows': 'cmd'}

current_app_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), app)
app_dir = os.path.join(os.path.expanduser('~'), app)

shortcut = os.path.join(os.path.expanduser('~'), 'Desktop', app + '.' + excecutable[system])

# install to home

if os.path.isdir(app_dir):
    shutil.rmtree(app_dir)
    shutil.copytree(current_app_dir, app_dir)
else:
    shutil.copytree(current_app_dir, app_dir)

# create shortcut

w = open(shortcut, 'w')
w.write('python ' + app_dir)
w.close()

if system == 'Darwin':
    os.system('chmod 755 ' + shortcut)
elif system == 'Linux':
    os.system('chmod +x' + shortcut)

try:
    import yaml
    print '\nPackage yaml already installed.'
except ImportError:
    os.system("pip install pyyaml")
try:
    import numpy
    print '\nPackage numpy already installed.'
except ImportError:
    os.system("pip install numpy")
try:
    import scipy
    print '\nPackage scipy already installed.'
except ImportError:
    os.system("pip install scipy")
try:
    import matplotlib
    print '\nPackage matplotlib already installed.'
except ImportError:
    os.system("pip install matplotlib")
try:
    import quantities
    print '\nPackage quantities already installed.'
except ImportError:
    os.system("pip install quantities")
try:
    import pyfits
    print '\nPackage pyfits already installed.'
except ImportError:
    os.system("pip install pyfits")
try:
    import ephem
    print '\nPackage ephem already installed.'
except ImportError:
    os.system("pip install ephem")

print '\n'
raw_input('Installation completed. Press enter to exit.')