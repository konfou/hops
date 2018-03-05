from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .hops_basics import *


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
    title_font = tuple(read_main_log('windows', 'title_font'))
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
                    if len(obj) == 5:
                        if obj[4] == 'title':
                            obj[0].configure(font=title_font)
                        else:
                            obj[0].configure(font=main_font)
                    else:
                        obj[0].configure(font=main_font)

                if len(obj) >= 4:
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
    # if not topmost:
    window.after_idle(window.attributes, '-topmost', 0)

    window.deiconify()

def fitting():

    fitting_directory = read_log('pipeline', 'fitting_directory')

    light_curve_file = read_log('fitting', 'light_curve_file')
    date = read_log('fitting', 'date')
    binning = read_log('fitting', 'binning')
    scatter = read_log('fitting', 'scatter')
    iterations = read_log('fitting', 'iterations')
    burn = read_log('fitting', 'burn')
    planet = read_log('fitting', 'planet')
    metallicity = read_log('fitting', 'metallicity')
    temperature = read_log('fitting', 'temperature')
    logg = read_log('fitting', 'logg')
    phot_filter = read_log('fitting', 'phot_filter')
    period = read_log('fitting', 'period')
    period_fit = read_log('fitting', 'period_fit')
    mid_time = read_log('fitting', 'mid_time')
    mid_time_fit = read_log('fitting', 'mid_time_fit')
    rp_over_rs = read_log('fitting', 'rp_over_rs')
    rp_over_rs_fit = read_log('fitting', 'rp_over_rs_fit')
    sma_over_rs = read_log('fitting', 'sma_over_rs')
    sma_over_rs_fit = read_log('fitting', 'sma_over_rs_fit')
    inclination = read_log('fitting', 'inclination')
    inclination_fit = read_log('fitting', 'inclination_fit')
    eccentricity = read_log('fitting', 'eccentricity')
    eccentricity_fit = read_log('fitting', 'eccentricity_fit')
    periastron = read_log('fitting', 'periastron')
    periastron_fit = read_log('fitting', 'periastron_fit')

    limb_darkening_coefficients = plc.clablimb('claret', logg, temperature, metallicity, phot_filter)

    light_curve = np.loadtxt(light_curve_file, unpack=True)

    if binning > 1:
        start = len(light_curve[0]) - (len(light_curve[0]) // binning) * binning
        light_curve_0 = np.mean(np.reshape(light_curve[0][start:],
                                           (light_curve[0].size // binning, binning)), 1)
        light_curve_1 = np.mean(np.reshape(light_curve[1][start:],
                                           (light_curve[1].size // binning, binning)), 1)
    else:
        light_curve_0 = light_curve[0]
        light_curve_1 = light_curve[1]

    light_curve_0 = light_curve_0[np.where(~np.isnan(light_curve_1))]
    light_curve_1 = light_curve_1[np.where(~np.isnan(light_curve_1))]

    moving_average = np.zeros_like(light_curve_0)
    for i in range(-5, 6):
        moving_average += np.roll(light_curve_1, i)
    for i in range(6):
        moving_average[i] = np.sum(light_curve_1[0:11])
        moving_average[-1 - i] = np.sum(light_curve_1[-11:])
    moving_average /= 11.0

    test = np.median(np.abs(light_curve_1 - moving_average))

    flag = np.where((np.abs(light_curve_1 - moving_average) < scatter * test))

    light_curve_0 = light_curve_0[flag]
    light_curve_1 = light_curve_1[flag]

    if not os.path.isdir(fitting_directory):
        os.mkdir(fitting_directory)
    else:
        fi = 2
        while os.path.isdir('{0}_{1}'.format(fitting_directory, str(fi))):
            fi += 1
        fitting_directory = '{0}_{1}'.format(fitting_directory, str(fi))
        os.mkdir(fitting_directory)

    if period_fit:
        period_fit = [period + period_fit[0], period + period_fit[1]]
    else:
        period_fit = False
    if mid_time_fit:
        mid_time_fit = [mid_time + mid_time_fit[0], mid_time + mid_time_fit[1]]
    else:
        mid_time_fit = False
    if rp_over_rs_fit:
        rp_over_rs_fit = [rp_over_rs * rp_over_rs_fit[0], rp_over_rs * rp_over_rs_fit[1]]
    else:
        rp_over_rs_fit = False
    if sma_over_rs_fit:
        sma_over_rs_fit = [sma_over_rs * sma_over_rs_fit[0], sma_over_rs * sma_over_rs_fit[1]]
    else:
        sma_over_rs_fit = False
    if inclination_fit:
        inclination_fit = [inclination + inclination_fit[0], inclination + inclination_fit[1]]
    else:
        inclination_fit = False
    if eccentricity_fit:
        eccentricity_fit = [eccentricity + eccentricity_fit[0], eccentricity + eccentricity_fit[1]]
    else:
        eccentricity_fit = False
    if periastron_fit:
        periastron_fit = [periastron + periastron_fit[0], periastron + periastron_fit[1]]
    else:
        periastron_fit = False

    mcmc_fit = plc.TransitAndPolyFitting([[light_curve_0, light_curve_1, np.ones_like(light_curve_1) *
                                         np.std(0.5 * (light_curve_1[:-1] - light_curve_1[1:]))]],
                                         method='claret',
                                         limb_darkening_coefficients=limb_darkening_coefficients,
                                         rp_over_rs=rp_over_rs,
                                         period=period,
                                         sma_over_rs=sma_over_rs,
                                         eccentricity=eccentricity,
                                         inclination=inclination,
                                         periastron=periastron,
                                         mid_time=mid_time,
                                         fit_rp_over_rs=rp_over_rs_fit,
                                         iterations=iterations,
                                         walkers=50,
                                         burn=burn,
                                         fit_first_order=True,
                                         fit_second_order=True,
                                         fit_period=period_fit,
                                         fit_sma_over_rs=sma_over_rs_fit,
                                         fit_eccentricity=eccentricity_fit,
                                         fit_inclination=inclination_fit,
                                         fit_periastron=periastron_fit,
                                         fit_mid_time=mid_time_fit,
                                         precision=3,
                                         exp_time=0,
                                         time_factor=1,
                                         counter=False,
                                         counter_window='FITTING'
                                         )

    mcmc_fit.run_mcmc()
    mcmc_fit.save_all(os.path.join(fitting_directory, 'results.pickle'))
    mcmc_fit.save_results(os.path.join(fitting_directory, 'results.txt'))
    mcmc_fit.plot_corner(os.path.join(fitting_directory, 'correlations.pdf'))
    mcmc_fit.plot_traces(os.path.join(fitting_directory, 'traces.pdf'))
    mcmc_fit.plot_models(os.path.join(fitting_directory, 'model.pdf'), planet, [date])
    figure = mcmc_fit.plot_detrended_models(os.path.join(fitting_directory, 'detrended_model.pdf'),
                                            planet, [date], return_plot=True)
    figure[0].savefig(os.path.join(fitting_directory, 'set_1_detrended_model_300dpi.jpg'),
                      dpi=300, transparent=True)
    figure[0].savefig(os.path.join(fitting_directory, 'set_1_detrended_model_600dpi.jpg'),
                      dpi=600, transparent=True)
    figure[0].savefig(os.path.join(fitting_directory, 'set_1_detrended_model_900dpi.jpg'),
                      dpi=900, transparent=True)
    mcmc_fit.save_models(os.path.join(fitting_directory, 'model.txt'))
    mcmc_fit.save_detrended_models(os.path.join(fitting_directory, 'detrended_model.txt'))

    shutil.copy('log.yaml', '{0}{1}log.yaml'.format(fitting_directory, os.sep))

    roott = Tk()

    exit_var_2 = BooleanVar(value=False)

    def break_and_exit():
        exit_var_2.set(True)

    initialise_window(roott, exit_command=break_and_exit)

    canvas = FigureCanvasTkAgg(figure[0], roott)
    canvas.get_tk_widget().pack()
    NavigationToolbar2TkAgg(canvas, roott)

    finalise_window(roott, topmost=True)

    while not exit_var_2.get():
        roott.update()

    roott.destroy()
