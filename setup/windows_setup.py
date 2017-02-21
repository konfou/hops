import os
import shutil

app = 'hops'

current_app_dir = os.path.join(os.path.split(os.path.abspath(os.path.dirname(__file__)))[0], app)
app_dir = os.path.join(os.path.expanduser('~'), app)

test_shortcut = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'run.command')
current_shortcut = os.path.join(os.path.abspath(os.path.dirname(__file__)), app + '.cmd')
shortcut = os.path.join(os.path.expanduser('~'), 'Desktop', app + '.cmd')

# install to home

if os.path.isdir(app_dir):
    shutil.rmtree(app_dir)
    shutil.copytree(current_app_dir, app_dir)
else:
    shutil.copytree(current_app_dir, app_dir)

# create shortcut

shutil.copy(test_shortcut, current_shortcut)

w = open(current_shortcut, 'a')
w.write('python ' + app_dir)
w.close()

shutil.move(current_shortcut, shortcut)

print '\n'
raw_input('Installation completed. Press enter to exit.')