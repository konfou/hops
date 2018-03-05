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


def photometry():

    # get variables

    reduction_directory = read_log('pipeline', 'reduction_directory')
    light_curve_aperture_file = read_log('pipeline', 'light_curve_aperture_file')
    photometry_directory = read_log('pipeline', 'photometry_directory')
    photometry_file = read_log('pipeline', 'photometry_file')
    light_curve_gauss_file = read_log('pipeline', 'light_curve_gauss_file')
    results_figure = read_log('pipeline', 'results_figure')
    fov_figure = read_log('pipeline', 'fov_figure')
    mean_key = read_log('pipeline_keywords', 'mean_key')
    std_key = read_log('pipeline_keywords', 'std_key')
    align_x0_key = read_log('pipeline_keywords', 'align_x0_key')
    align_y0_key = read_log('pipeline_keywords', 'align_y0_key')
    align_u0_key = read_log('pipeline_keywords', 'align_u0_key')
    exposure_time_key = read_log('pipeline_keywords', 'exposure_time_key')
    observation_date_key = read_log('pipeline_keywords', 'observation_date_key')
    observation_time_key = read_log('pipeline_keywords', 'observation_time_key')
    star_std = read_log('alignment', 'star_std')
    search_window_std = read_log('alignment', 'search_window_std')
    target_ra_dec = read_log('photometry', 'target_ra_dec')
    sky_inner_aperture = read_log('photometry', 'sky_inner_aperture')
    sky_outer_aperture = read_log('photometry', 'sky_outer_aperture')
    max_comparisons = read_log('photometry', 'max_comparisons')
    targets_r_position = [read_log('photometry', 'target_r_position')]
    targets_u_position = [read_log('photometry', 'target_u_position')]
    targets_aperture = [read_log('photometry', 'target_aperture')]
    for comparison in range(max_comparisons):
        targets_r_position.append(read_log('photometry', 'comparison_{0}_r_position'.format(comparison)))
        targets_u_position.append(read_log('photometry', 'comparison_{0}_u_position'.format(comparison)))
        targets_aperture.append(read_log('photometry', 'comparison_{0}_aperture'.format(comparison)))

    science = glob.glob('{0}{1}*.f*t*'.format(reduction_directory, os.sep))
    science.sort()

    root = Tk()

    exit_var = BooleanVar(value=False)

    def break_and_exit():
        exit_var.set(True)

    initialise_window(root, exit_command=break_and_exit)

    label1 = Label(root, text='PHOTOMETRY')
    label2 = Label(root, text='FILE:')
    label3 = Label(root, text=' ')
    label4 = Label(root, text='COMPLETE:')
    label5 = Label(root, text=' ',)
    label6 = Label(root, text='TIME LEFT:')
    label7 = Label(root, text=' ')

    setup_window(root, [
        [[label1, 1, 2]],
        [[label2, 1], [label3, 2]],
        [[label4, 1], [label5, 2]],
        [[label6, 1], [label7, 2]],
    ])

    targets_files = []
    targets_hjd = []
    targets_x_position = []
    targets_y_position = []
    targets_x_std = []
    targets_y_std = []
    targets_gauss_flux = []
    targets_gauss_sky = []
    targets_aperture_flux = []
    targets_aperture_sky = []

    # for each science_file
    percent = 0
    lt0 = time.time()
    for counter, science_file in enumerate(science):

        fits = pf.open(science_file)

        if fits[1].header[align_x0_key]:

            targets_files.append(science_file)

            # calculate heliocentric julian date

            if observation_date_key == observation_time_key:
                local_time = ' '.join(fits[1].header[observation_date_key].split('T'))
                if counter == 1:
                    write_log('fitting', fits[1].header[observation_date_key].split('T')[0], 'date')
            else:
                local_time = ' '.join([fits[1].header[observation_date_key], fits[1].header[observation_time_key]])
                if counter == 1:
                    write_log('fitting', fits[1].header[observation_date_key], 'date')

            julian_date = (ephem.julian_date(float(ephem.Date(local_time))) +
                           fits[1].header[exposure_time_key] / (2.0 * 60.0 * 60.0 * 24.0))

            ra_target, dec_target = target_ra_dec.split()

            heliocentric_julian_date = jd_to_hjd(ra_target, dec_target, julian_date)

            targets_hjd.append(heliocentric_julian_date)

            # calculate gauss position, flux and sky and aperture flux and sky

            ref_x_position = fits[1].header[align_x0_key]
            ref_y_position = fits[1].header[align_y0_key]
            ref_u_position = fits[1].header[align_u0_key]

            for target in range(max_comparisons + 1):

                if targets_aperture[target] > 0:

                    norm, floor, x_mean, y_mean, x_std, y_std = \
                        fit_2d_gauss(fits[1].data,
                                     predicted_x_mean=(ref_x_position + targets_r_position[target] *
                                                       np.cos(ref_u_position + targets_u_position[target])),
                                     predicted_y_mean=(ref_y_position + targets_r_position[target] *
                                                       np.sin(ref_u_position + targets_u_position[target])),
                                     search_window=search_window_std * star_std)

                    targets_x_position.append(x_mean)
                    targets_y_position.append(y_mean)
                    targets_x_std.append(x_std)
                    targets_y_std.append(y_std)
                    targets_gauss_flux.append(2 * np.pi * norm * x_std * y_std)
                    targets_gauss_sky.append(floor)

                    flux_area = fits[1].data[int(y_mean) - targets_aperture[target]:
                                             int(y_mean) + targets_aperture[target] + 1,
                                             int(x_mean) - targets_aperture[target]:
                                             int(x_mean) + targets_aperture[target] + 1]
                    flux_pixels = (2 * targets_aperture[target] + 1) ** 2
                    flux = np.sum(flux_area)

                    sky_area_1 = int(sky_inner_aperture * targets_aperture[target])
                    sky_area_2 = int(sky_outer_aperture * targets_aperture[target])
                    fits[1].data[int(y_mean) - sky_area_1:int(y_mean) + sky_area_1 + 1,
                                 int(x_mean) - sky_area_1:int(x_mean) + sky_area_1 + 1] = 0
                    sky_area = fits[1].data[int(y_mean) - sky_area_2:int(y_mean) + sky_area_2 + 1,
                                            int(x_mean) - sky_area_2:int(x_mean) + sky_area_2 + 1]
                    sky_area = sky_area[np.where((sky_area > 0) &
                                                 (sky_area < fits[1].header[mean_key] + 3 * fits[1].header[std_key]))]
                    sky = np.sum(sky_area)
                    sky_pixels = sky_area.size

                    targets_aperture_flux.append(flux - flux_pixels * sky / sky_pixels)
                    targets_aperture_sky.append(sky / sky_pixels)

        # counter

        new_percent = round(100 * (counter + 1) / float(len(science)), 1)
        if new_percent != percent:
            lt1 = time.time()
            rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
            hours = rm_time / 3600.0
            minutes = (hours - int(hours)) * 60
            seconds = (minutes - int(minutes)) * 60
            label3.configure(text='     {0}     '.format(science_file.split(os.sep)[-1]))
            label5.configure(text='     {0}%    '.format(new_percent))
            label7.configure(text='     %dh %02dm %02ds     ' % (int(hours), int(minutes), int(seconds)))
            percent = new_percent

        if counter == 0:
            finalise_window(root, topmost=True)

        root.update()

        if exit_var.get():
            break

        if counter + 1 == len(science):
            write_log('pipeline', True, 'photometry_complete')

    root.destroy()

    if not exit_var.get():

        # save results, create photometry directory and move results there

        measurements_number = len(targets_files)
        targets_number = len(targets_x_position) // len(targets_files)
        comparisons_number = len(targets_x_position) // len(targets_files) - 1

        targets_hjd = np.array(targets_hjd)
        targets_x_position = np.swapaxes(np.reshape(targets_x_position, (measurements_number, targets_number)), 0, 1)
        targets_y_position = np.swapaxes(np.reshape(targets_y_position, (measurements_number, targets_number)), 0, 1)
        targets_x_std = np.swapaxes(np.reshape(targets_x_std, (measurements_number, targets_number)), 0, 1)
        targets_y_std = np.swapaxes(np.reshape(targets_y_std, (measurements_number, targets_number)), 0, 1)
        targets_gauss_flux = np.swapaxes(np.reshape(targets_gauss_flux, (measurements_number, targets_number)), 0, 1)
        targets_gauss_sky = np.swapaxes(np.reshape(targets_gauss_sky, (measurements_number, targets_number)), 0, 1)
        targets_aperture_flux = np.swapaxes(np.reshape(targets_aperture_flux,
                                                       (measurements_number, targets_number)), 0, 1)
        targets_aperture_sky = np.swapaxes(np.reshape(targets_aperture_sky,
                                                      (measurements_number, targets_number)), 0, 1)

        targets_results = [targets_hjd] + (list(targets_x_position) + list(targets_y_position) + list(targets_x_std) +
                                           list(targets_y_std) + list(targets_gauss_flux) + list(targets_gauss_sky) +
                                           list(targets_aperture_flux) + list(targets_aperture_sky))

        np.savetxt(photometry_file,
                   np.swapaxes(targets_results, 0, 1))

        np.savetxt(light_curve_gauss_file,
                   np.swapaxes([targets_hjd, targets_gauss_flux[0] / np.sum(targets_gauss_flux[1:], 0)], 0, 1))

        np.savetxt(light_curve_aperture_file,
                   np.swapaxes([targets_hjd, targets_aperture_flux[0] / np.sum(targets_aperture_flux[1:], 0)], 0, 1))

        if not os.path.isdir(photometry_directory):
            os.mkdir(photometry_directory)
        else:
            fi = 2
            while os.path.isdir('{0}_{1}'.format(photometry_directory, str(fi))):
                fi += 1
            photometry_directory = '{0}_{1}'.format(photometry_directory, str(fi))
            os.mkdir(photometry_directory)

        root = Tk()

        if comparisons_number > 1:
            f = Figure()
            f.set_figwidth(7)
            f.set_figheight(0.8 * root.winfo_screenheight() / f.get_dpi())
            ax = f.add_subplot(comparisons_number + 1, 1, 1)
        else:
            f = Figure()
            ax = f.add_subplot(1, 1, 1)

        exit_var_2 = BooleanVar(value=False)

        def break_and_exit():
            exit_var_2.set(True)

        initialise_window(root, exit_command=break_and_exit)

        f.patch.set_facecolor('white')
        canvas = FigureCanvasTkAgg(f, root)
        canvas.get_tk_widget().pack()
        NavigationToolbar2TkAgg(canvas, root)

        ax.plot(targets_hjd - targets_hjd[0], targets_aperture_flux[0] / np.sum(targets_aperture_flux[1:], 0)
                / np.median(targets_aperture_flux[0] / np.sum(targets_aperture_flux[1:], 0)), 'ko', ms=3)
        ax.plot(targets_hjd - targets_hjd[0], targets_gauss_flux[0] / np.sum(targets_gauss_flux[1:], 0)
                / np.median(targets_gauss_flux[0] / np.sum(targets_gauss_flux[1:], 0)), 'ro', ms=3, mec='r')
        ax.tick_params(labelbottom='off')
        ax.set_title(r'$\mathrm{Target}$')

        if comparisons_number > 1:
            for comp in range(comparisons_number):
                test_aperture_flux = list(targets_aperture_flux[1:])
                test_gauss_flux = list(targets_gauss_flux[1:])
                del test_aperture_flux[comp]
                del test_gauss_flux[comp]
                ax = f.add_subplot(comparisons_number + 1, 1, comp + 2)
                ax.plot(targets_hjd - targets_hjd[0], targets_aperture_flux[1:][comp] / np.sum(test_aperture_flux, 0)
                        / np.median(targets_aperture_flux[1:][comp] / np.sum(test_aperture_flux, 0)), 'ko', ms=3)
                ax.plot(targets_hjd - targets_hjd[0], targets_gauss_flux[1:][comp] / np.sum(test_gauss_flux, 0)
                        / np.median(targets_gauss_flux[1:][comp] / np.sum(test_gauss_flux, 0)),
                        'ro', ms=3, mec='r')
                ax.tick_params(labelbottom='off')
                ax.set_title(r'${0}{1}{2}$'.format('\mathrm{', 'Comparison \, {0}'.format(comp + 1), '}'))

        ax.tick_params(labelbottom='on')
        ax.set_xlabel(r'$\mathrm{\Delta t} \ \mathrm{[days]}$', fontsize=20)
        f.text(0.03, 0.5, r'$\mathrm{relative} \ \mathrm{flux}$', fontsize=20,
               ha='center', va='center', rotation='vertical')

        f.set_figheight(2 * (comparisons_number + 1))
        f.savefig(results_figure, dpi=200)

        shutil.move(fov_figure, '{0}/{1}'.format(photometry_directory, fov_figure))
        shutil.move(results_figure, '{0}/{1}'.format(photometry_directory, results_figure))
        shutil.move(photometry_file, '{0}/{1}'.format(photometry_directory, photometry_file))
        shutil.move(light_curve_gauss_file, '{0}/{1}'.format(photometry_directory, light_curve_gauss_file))
        shutil.move(light_curve_aperture_file, '{0}/{1}'.format(photometry_directory, light_curve_aperture_file))
        shutil.copy('log.yaml', '{0}/log.yaml'.format(photometry_directory))

        finalise_window(root, topmost=True)

        while not exit_var_2.get():
            root.update()

        root.destroy()
