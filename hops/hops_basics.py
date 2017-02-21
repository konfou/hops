import warnings
warnings.filterwarnings("ignore",
                        message='Matplotlib is building the font cache using fc-list. This may take a moment.')
warnings.filterwarnings("ignore",
                        message='The installed version of numexpr 2.4.4 is not supported in pandas and will be not be used')
from Tkinter import *
import tkFileDialog
from tkMessageBox import *
import pyfits as pf
import os
import glob
import shutil
import yaml
import numpy as np
import time
import ephem
import tools
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasBase, FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler, MouseEvent
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.offsetbox import AnchoredText
import gzip
import socket
import urllib
import ttk
import sys


__location__ = os.path.dirname(__file__)

holomon_logo = glob.glob(__location__ + '/holomon.gif')[0]
main_logfile = __location__ + '/log.yaml'
logfile = 'log.yaml'

def read_main_log(keyword, keyword2=None):
    x = yaml.load(open(main_logfile, 'r'))
    if keyword2:
        return x[keyword][keyword2]
    else:
        return x[keyword]


def write_main_log(keyword, value, keyword2=None):

    x = yaml.load(open(main_logfile, 'r'))
    if keyword2:
        x[keyword][keyword2] = value
    else:
        x[keyword] = value
    yaml.dump(x, open(main_logfile, 'w'), default_flow_style=False)


def copy_main_log():
    if not os.path.isfile(logfile):
        shutil.copy(main_logfile, "{0}/{1}".format(os.path.abspath('.'), logfile))


def read_log(keyword, keyword2=None):

    try:
        x = yaml.load(open(logfile, 'r'))
        if keyword2:
            return x[keyword][keyword2]
        else:
            return x[keyword]

    except KeyError:

        try:
            value = read_main_log(keyword, keyword2)
            write_log(keyword, value, keyword2)

            return value

        except KeyError:

            return False


def write_log(keyword, value, keyword2=None):
    x = yaml.load(open(logfile, 'r'))
    if keyword2:
        x[keyword][keyword2] = value
    else:
        x[keyword] = value

    yaml.dump(x, open(logfile, 'w'), default_flow_style=False)


def test_fits_keyword(fits_file, keyword):

    if len(fits_file) == 0:
        return [False, 'No keyword found']

    else:
        try:
            fits_file = glob.glob('*' + fits_file + '*.f*t*')[0]

            if pf.open(fits_file)[0].header[str(keyword)]:
                return [True, 'Keyword found']

            else:
                return [False, 'No keyword found']

        except (KeyError, IndexError):
            return [False, 'No keyword found']


def test_file_number(fits_file):

    if len(fits_file) == 0:
        test = 0
    else:
        test = len(glob.glob('*' + fits_file + '*.f*t*'))

    if test > 0:
        return [True, '{0} files found'.format(test)]
    else:
        return [False, 'No files found']


def test_coordinates(radec_string):

    if len(radec_string.split()) != 2:
        return [False, 'Wrong\ncoordinates']

    if len(radec_string.split()[0].split(':')) != 3:
        return [False, 'Wrong\ncoordinates']

    if len(radec_string.split()[0].split(':')[0]) != 2:
        return [False, 'Wrong\ncoordinates']
    else:
        try:
            if int(radec_string.split()[0].split(':')[0]) >= 24:
                return [False, 'Wrong\ncoordinates']
        except ValueError:
            return [False, 'Wrong\ncoordinates']

    if len(radec_string.split()[0].split(':')[1]) != 2:
        return [False, 'Wrong\ncoordinates']
    else:
        try:
            if int(radec_string.split()[0].split(':')[1]) >= 60:
                return [False, 'Wrong\ncoordinates']
        except ValueError:
            return [False, 'Wrong\ncoordinates']

    if len(radec_string.split()[0].split(':')[2]) < 2:
        return [False, 'Wrong\ncoordinates']
    else:
        try:
            if float(radec_string.split()[0].split(':')[2]) >= 60:
                return [False, 'Wrong\ncoordinates']
        except ValueError:
            return [False, 'Wrong\ncoordinates']

    if len(radec_string.split()[1][1:].split(':')) != 3:
        return [False, 'Wrong\ncoordinates']

    if len(radec_string.split()[1].split(':')[0]) != 3:
        return [False, 'Wrong\ncoordinates']
    else:
        try:
            if abs(int(radec_string.split()[1].split(':')[0])) >= 90:
                return [False, 'Wrong\ncoordinates']
        except ValueError:
            return [False, 'Wrong\ncoordinates']

    if len(radec_string.split()[1].split(':')[1]) != 2:
        return [False, 'Wrong\ncoordinates']
    else:
        try:
            if int(radec_string.split()[1].split(':')[1]) >= 60:
                return [False, 'Wrong\ncoordinates']
        except ValueError:
            return [False, 'Wrong\ncoordinates']

    if len(radec_string.split()[1].split(':')[2]) < 2:
        return [False, 'Wrong\ncoordinates']
    else:
        try:
            if float(radec_string.split()[1].split(':')[2]) >= 60:
                return [False, 'Wrong\ncoordinates']
        except ValueError:
            return [False, 'Wrong\ncoordinates']

    try:
        if ephem.Equatorial(radec_string.split()[0], radec_string.split()[1]):
            return [True, 'Coordinates\naccepted']

    except (ValueError, TypeError):
        return [False, 'Wrong\ncoordinates']


def initialise_window(window, window_name=None, exit_command=None):

    if not window_name:
        window_name = read_main_log('windows', 'software_window')

    if not exit_command:
        def exit_command():
            os._exit(-1)

    window.wm_title(window_name)
    window.protocol('WM_DELETE_WINDOW', exit_command)

    window.withdraw()


def setup_window(window, objects):

    main_font = tuple(read_main_log('windows', 'main_font'))
    button_font = tuple(read_main_log('windows', 'button_font'))
    entries_bd = read_main_log('windows', 'entries_bd')

    for row in range(len(objects)):
        if len(objects[row]) == 0:
            label_empty = Label(window, text='')
            label_empty.grid(row=row, column=100)
        else:
            for obj in objects[row]:

                if obj[0].winfo_class() == 'Button':
                    obj[0].configure(font=button_font)
                elif obj[0].winfo_class() == 'Entry':
                    obj[0].configure(bd=entries_bd, font=main_font)
                elif obj[0].winfo_class() in ['Label', 'Radiobutton']:
                    obj[0].configure(font=main_font)

                if len(obj) == 4:
                    obj[0].grid(row=row, column=obj[1], columnspan=obj[2], rowspan=obj[3])
                elif len(obj) == 3:
                    obj[0].grid(row=row, column=obj[1], columnspan=obj[2])
                else:
                    obj[0].grid(row=row, column=obj[1])


def finalise_window(window, center=True, topmost=False):

    window.update_idletasks()

    if center:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2
        window.geometry('+%d+%d' % (x, y))

    else:
        window.geometry('+%d+%d' % (0, 0))

    window.update_idletasks()

    window.lift()
    window.wm_attributes("-topmost", 1)
    if not topmost:
        window.after_idle(window.attributes, '-topmost', 0)

    window.deiconify()


def test_float_input(input_str, typing):

    if typing == '1':
        try:
            if float(input_str):
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_float_positive_input(input_str, typing):

    if typing == '1':
        try:
            if float(input_str) >= 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_int_input(input_str, typing):

    if typing == '1':
        try:
            if int(input_str):
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_int_positive_input(input_str, typing):

    if typing == '1':
        try:
            if int(input_str) >= 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_int_positive_non_zero_input(input_str, typing):

    if typing == '1':
        try:
            if int(input_str) > 0:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True


def test_phot_filter(input_str, typing):

    if typing == '1':
        try:
            if input_str in ['u', 'v', 'b', 'y', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']:
                return True
            else:
                return False
        except ValueError:
            return False

    else:
        return True
