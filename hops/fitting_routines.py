from hops_basics import *

from transit_fitting import *
from clablimb import *
from transit import *


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

    light_curve = np.loadtxt(light_curve_file, unpack=True)

    if binning > 1:
        start = np.mod(len(light_curve[0]), binning)
        light_curve_0 = np.mean(np.reshape(light_curve[0][start:],
                                           (light_curve[0].size / binning, binning)), 1)
        light_curve_1 = np.mean(np.reshape(light_curve[1][start:],
                                           (light_curve[1].size / binning, binning)), 1)
    else:
        light_curve_0 = light_curve[0]
        light_curve_1 = light_curve[1]

    light_curve_0 = light_curve_0[np.where(~np.isnan(light_curve_1))]
    light_curve_1 = light_curve_1[np.where(~np.isnan(light_curve_1))]

    test = np.median(np.abs(light_curve_1[:-1] - light_curve_1[1:]))
    light_curve_0 = np.insert(np.append(light_curve_0, light_curve_0[-1]), 0, light_curve_0[0])
    light_curve_1 = np.insert(np.append(light_curve_1, light_curve_1[-1]), 0, light_curve_1[0])

    flag = np.where((np.abs(light_curve_1[1:-1] - light_curve_1[:-2]) < scatter * test) |
                    (np.abs(light_curve_1[1:-1] - light_curve_1[2:]) < scatter * test))

    root = Tk()

    def do_nothing():
        pass

    initialise_window(root, window_name='Fitting preview', exit_command=do_nothing)

    f = Figure(figsize=(7, 2 * 3))
    f.patch.set_facecolor('white')
    ax1 = f.add_subplot(211)
    ax2 = f.add_subplot(212)
    f.subplots_adjust(hspace=0, top=0.65)
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()

    limb_darkening_coefficients = clablimb(logg, temperature, metallicity, phot_filter)

    data_delta_t = light_curve_0[1:-1][flag] - light_curve_0[0]

    def mcmc_f(inputs, detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

        if inputs:
            detrend = detrend_zero * (1 + detrend_one * data_delta_t + detrend_two * data_delta_t * data_delta_t)
            transit_model = transit(limb_darkening_coefficients, model_rp_over_rs, period,
                                    sma_over_rs, eccentricity, inclination,
                                    periastron, mid_time + model_mid_time, time_array=light_curve_0[1:-1][flag])

            return detrend * transit_model

    def independent_f(detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):

        detrend = detrend_zero * (1 + detrend_one * data_delta_t + detrend_two * data_delta_t * data_delta_t)
        transit_model = transit(limb_darkening_coefficients, model_rp_over_rs, period,
                                sma_over_rs, eccentricity, inclination,
                                periastron, mid_time + model_mid_time, time_array=light_curve_0[1:-1][flag])

        return detrend, transit_model

    popt, pcov = curve_fit(mcmc_f, True, light_curve_1[1:-1][flag],
                           p0=[np.mean(light_curve_1[1:-1][flag]), 1, 1, rp_over_rs, 0])

    fit_detrend, fit_transit_model = independent_f(*popt)

    predicted_transit_model = transit(limb_darkening_coefficients, rp_over_rs, period,
                                      sma_over_rs, eccentricity, inclination,
                                      periastron, mid_time, time_array=light_curve_0[1:-1][flag])

    new_mid_time = (mid_time
                    + round((np.mean(light_curve_0) - mid_time) / period) * period
                    + popt[-1])

    phase = (light_curve_0 - new_mid_time) / period

    ax1.plot(phase[1:-1], light_curve_1[1:-1], 'ro', ms=3, mec='r')
    ax1.plot(phase[1:-1][flag], light_curve_1[1:-1][flag], 'ko', ms=3)
    ax1.plot(phase[1:-1][flag], fit_detrend * fit_transit_model, 'r-')
    ax1.set_yticks(ax1.get_yticks()[1:])
    ax1.tick_params(labelbottom='off')
    ax1.set_ylabel(r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=20)

    ax2.plot(phase[1:-1][flag], light_curve_1[1:-1][flag] / fit_detrend, 'ko', ms=3)
    ax2.plot(phase[1:-1][flag], fit_transit_model, 'r-')
    ax2.plot(phase[1:-1][flag], predicted_transit_model, 'c-')
    ax2.set_ylabel(r'$\mathrm{normalised} \ \mathrm{flux}$', fontsize=20)
    ax2.set_xlabel(r'$\mathrm{phase}$', fontsize=20)

    finalise_window(root, topmost=True)

    root.update()

    if askyesno('Fitting preview',
                'This is a quick fitting preview. Do you want to proceed to MCMC fitting?', parent=root):

        root.destroy()

        light_curve_0 = light_curve_0[1:-1][flag]
        light_curve_1 = light_curve_1[1:-1][flag]

        if not os.path.isdir(fitting_directory):
            os.mkdir(fitting_directory)
        else:
            fi = 2
            while os.path.isdir(fitting_directory + '_' + str(fi)):
                fi += 1
            fitting_directory = fitting_directory + '_' + str(fi)
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

        mcmc_fit = TransitAndPolyFitting([[light_curve_0, light_curve_1, np.ones_like(light_curve_1) *
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
                                         walkers=300,
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
                                         counter_window=True
                                         )

        mcmc_fit.run_mcmc()
        mcmc_fit.save_all(os.path.join(fitting_directory, 'results.pickle'))
        mcmc_fit.save_results(os.path.join(fitting_directory, 'results.txt'))
        mcmc_fit.plot_corner(os.path.join(fitting_directory, 'correlations.pdf'))
        mcmc_fit.plot_traces(os.path.join(fitting_directory, 'traces.pdf'))
        mcmc_fit.plot_models(os.path.join(fitting_directory, 'model.pdf'), planet, [date])
        mcmc_fit.plot_detrended_models(os.path.join(fitting_directory, 'detrended_model.pdf'), planet, [date])
        mcmc_fit.save_models(os.path.join(fitting_directory, 'model.txt'))
        mcmc_fit.save_detrended_models(os.path.join(fitting_directory, 'detrended_model.txt'))

        shutil.copy('log.yaml', '{0}{1}log.yaml'.format(fitting_directory, os.sep))

        # root.destroy()

        # root = Tk()
        #
        # exit_var_2 = BooleanVar(value=False)
        #
        # def break_and_exit():
        #     exit_var_2.set(True)
        #
        # initialise_window(root, window_name=fitting_window, exit_command=break_and_exit)
        #
        # f = Figure(figsize=(7, 2 * 2))
        # f.patch.set_facecolor('white')
        # ax1 = f.add_subplot(211)
        # ax2 = f.add_subplot(212)
        # f.subplots_adjust(hspace=0)
        # canvas = FigureCanvasTkAgg(f, master=root)
        # canvas.show()
        # canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        #
        # toolbar = NavigationToolbar2TkAgg(canvas, root)
        # toolbar.update()
        #
        # data_delta_t = light_curve_0 - light_curve_0[0]
        #
        # def independent_f(detrend_zero, detrend_one, detrend_two, model_rp_over_rs, model_mid_time):
        #
        #     detrend = detrend_zero * (1 + detrend_one * data_delta_t + detrend_two * data_delta_t * data_delta_t)
        #     transit_model = pylightcurve.transit(limb_darkening_coefficients, model_rp_over_rs, period,
        #                                          sma_over_rs, eccentricity, inclination,
        #                                          periastron, mid_time + model_mid_time,
        #                                          time_array=light_curve_0)
        #
        #     return detrend, transit_model
        #
        # fitting_results =
        #
        # fit_detrend, fit_transit_model = independent_f(*popt)
        #
        # new_mid_time = (mid_time
        #                 + round((np.mean(light_curve_0) - mid_time) / period) * period
        #                 + popt[-1])
        #
        # phase = (light_curve_0 - new_mid_time) / period
        #
        # ax1.plot(phase, light_curve_1, 'ko', ms=3)
        # ax1.plot(phase, fit_detrend * fit_transit_model, 'r-')
        # ax1.set_yticks(ax1.get_yticks()[1:])
        # ax1.tick_params(labelbottom='off')
        # ax1.set_ylabel(r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=20)
        #
        # ax2.plot(phase[1:-1][flag], light_curve_1[1:-1][flag] / fit_detrend, 'ko', ms=3)
        # ax2.plot(phase[1:-1][flag], fit_transit_model, 'r-')
        # ax2.plot(phase[1:-1][flag], predicted_transit_model, 'c-')
        # ax2.set_ylabel(r'$\mathrm{normalised} \ \mathrm{flux}$', fontsize=20)
        # ax2.set_xlabel(r'$\mathrm{phase}$', fontsize=20)
        #
        # finalise_window(root, topmost=True)
        #
        # while not exit_var_2.get():
        #     root.update()
        #
        # root.destroy()

    else:
        root.destroy()
