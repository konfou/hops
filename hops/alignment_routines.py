from hops_basics import *


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
    window.after_idle(window.attributes, '-topmost', 0)

    window.deiconify()


def alignment():

    if read_log('pipeline', 'alignment_complete'):
        if not askyesno('Overwrite alignment', 'Alignment has been completed, do you want to run again?'):
            return 0

    write_log('pipeline', False, 'alignment_complete')

    # get variables

    reduction_directory = read_log('pipeline', 'reduction_directory')
    mean_key = read_log('pipeline_keywords', 'mean_key')
    std_key = read_log('pipeline_keywords', 'std_key')
    align_star_area_key = read_log('pipeline_keywords', 'align_star_area_key')
    align_x0_key = read_log('pipeline_keywords', 'align_x0_key')
    align_y0_key = read_log('pipeline_keywords', 'align_y0_key')
    align_u0_key = read_log('pipeline_keywords', 'align_u0_key')
    frame_low_std = read_log('windows', 'frame_low_std')
    frame_upper_std = read_log('windows', 'frame_upper_std')
    burn_limit = read_log('alignment', 'burn_limit')
    star_std = read_log('alignment', 'star_std')
    search_window_std = read_log('alignment', 'search_window_std')
    shift_tolerance = read_log('alignment', 'shift_tolerance_p')

    science = glob.glob('{0}{1}*.f*t*'.format(reduction_directory, os.sep))
    science.sort()

    fits = pf.open(science[0], memmap=False)

    shift_tolerance = int(max(len(fits[1].data), len(fits[1].data[0])) * (shift_tolerance / 100.0))
    y_length, x_length = fits[1].data.shape

    centroids = []
    std_limit = 5.0
    while len(centroids) == 0 and std_limit >= 1.0:
        centroids = tools.find_centroids(fits[1].data, mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                         std_limit=std_limit, burn_limit=7.0 * burn_limit / 8, star_std=2)
        std_limit -= 1.0
    calibration_stars = [[-pp[3], pp[1], pp[2], pp[0]] for pp in centroids]

    new_star_std = []
    for calibration_centroid in calibration_stars[:int(0.1 * len(calibration_stars))]:
        norm, floor, x_mean, y_mean, x_sigma, y_sigma = \
            tools.fit_2d_gauss_point(fits[1].data,
                                     predicted_x_mean=calibration_centroid[1],
                                     predicted_y_mean=calibration_centroid[2],
                                     search_window=2 * star_std)
        if not np.isnan(x_mean * y_mean):
            new_star_std.append(x_sigma)
    star_std = max(1, int(np.median(new_star_std)) - 1)
    write_log('alignment', star_std, 'star_std')

    centroids = []
    std_limit = 5.0
    while len(centroids) == 0 and std_limit >= 1.0:
        centroids = tools.find_centroids(fits[1].data, mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                         std_limit=std_limit, burn_limit=7.0 * burn_limit / 8, star_std=2 * star_std)
        std_limit -= 1.0
    calibration_stars = [[-pp[3], pp[1], pp[2], pp[0]] for pp in centroids]

    x_ref_position = np.nan
    y_ref_position = np.nan
    while np.isnan(x_ref_position * y_ref_position):
        norm, floor, x_mean, y_mean, x_sigma, y_sigma = tools.fit_2d_gauss_point(
            fits[1].data,
            predicted_x_mean=calibration_stars[0][1], predicted_y_mean=calibration_stars[0][2],
            search_window=2 * star_std)
        x_ref_position, y_ref_position = x_mean, y_mean
        del calibration_stars[0]

    centroids = tools.find_centroids(fits[1].data, mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                     std_limit=3.0, burn_limit=7.0 * burn_limit / 8, star_std=star_std)
    calibration_stars = [[-pp[3], pp[1], pp[2], pp[0]] for pp in centroids]
    calibration_stars.sort()

    # take the rest as calibration stars and calculate their polar coordinates relatively to the first
    calibration_stars_polar = []
    for calibration_star in calibration_stars[:100]:
        norm, floor, x_mean, y_mean, x_sigma, y_sigma = \
            tools.fit_2d_gauss_point(fits[1].data,
                                     predicted_x_mean=calibration_star[1],
                                     predicted_y_mean=calibration_star[2],
                                     search_window=2 * star_std)
        x_position, y_position = x_mean, y_mean
        r_position, u_position = tools.cartesian_to_polar(x_position, y_position, x_ref_position, y_ref_position)
        if not np.isnan(x_position * y_position):
            calibration_stars_polar.append([r_position, u_position])

    x0, y0, u0, comparisons = x_ref_position, y_ref_position, 0, calibration_stars_polar
    fits.close()

    # set the looking window and angular step

    root = Tk()

    exit_var = BooleanVar(value=False)

    def break_and_exit():
        exit_var.set(True)

    initialise_window(root, exit_command=break_and_exit)

    f = Figure()
    f.patch.set_facecolor('white')
    ax = f.add_subplot(111)
    ax.axis('off')
    canvas = FigureCanvasTkAgg(f, root)
    canvas.get_tk_widget().pack()

    fits = pf.open(science[0], memmap=False)
    ax.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
              vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
              vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])
    fits.close()
    circle = mpatches.Circle((x0, y0), 2 * search_window_std * star_std, ec='r', fill=False)
    ax.add_patch(circle)
    for ii in comparisons:
        circle = mpatches.Circle((x0 + ii[0] * np.cos(u0 + ii[1]), y0 + ii[0] * np.sin(u0 + ii[1])),
                                 2 * search_window_std * star_std, ec='w', fill=False)
        ax.add_patch(circle)

    frame1 = Frame(root)
    frame1.pack()

    label1 = Label(frame1, text='ALIGNMENT')
    label2 = Label(frame1, text='FILE:')
    label3 = Label(frame1, text=' ')
    label4 = Label(frame1, text='COMPLETE:')
    label5 = Label(frame1, text=' ')
    label6 = Label(frame1, text='TIME LEFT:')
    label7 = Label(frame1, text=' ')

    setup_window(frame1, [
        [[label1, 1, 2]],
        [[label2, 1], [label3, 2]],
        [[label4, 1], [label5, 2]],
        [[label6, 1], [label7, 2]],
    ])

    # for each science_file
    percent = 0
    skip_time = 0
    lt0 = time.time()
    for counter, science_file in enumerate(science):

        fits = pf.open(science_file, mode='update')

        # fast detection test
        centroids = tools.find_centroids(fits[1].data,
                                         x_low=int(x0 - shift_tolerance), x_upper=int(x0 + shift_tolerance + 1),
                                         y_low=int(y0 - shift_tolerance), y_upper=int(y0 + shift_tolerance + 1),
                                         x_centre=int(x0), y_centre=int(y0),
                                         mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                         std_limit=2.0, burn_limit=2 * burn_limit, star_std=star_std)

        if len(centroids) > 0:

            tests = []

            for ref_star in centroids:

                max_x = ref_star[1]
                max_y = ref_star[2]

                test = 0

                for comp in comparisons:

                    check_x = int(max_x + comp[0] * np.cos(u0 + comp[1]))
                    check_y = int(max_y + comp[0] * np.sin(u0 + comp[1]))
                    if 0 < check_x < x_length and 0 < check_y < y_length:
                        check_sum = np.sum(fits[1].data[check_y - star_std:check_y + star_std + 1,
                                           check_x - star_std:check_x + star_std + 1])
                        check_lim = (fits[1].header[mean_key] + fits[1].header[std_key]) * ((2 * star_std + 1) ** 2)
                        if check_sum > check_lim:
                            test += 1
                        else:
                            test -= 1

                    if test > 0.1 * len(comparisons):
                        break

                tests.append([test, max_x, max_y])

                if test > 0.1 * len(comparisons):
                    break

            tests.sort()
            test, max_x, max_y = tests[-1]
            if test < 0.1 * len(comparisons):
                stars_detected = False
                delta_skip_time = time.time()
                label3.configure(text='     ' + science_file.split(os.sep)[-1] + '     ')
                label5.configure(text='     -%    ')
                label7.configure(text='     -h -m -s     ')
                ax.cla()
                ax.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
                          vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
                          vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])
                ax.axis('off')
                canvas.show()
                root.update()
                skip_frame = askyesno('Alignment',
                                      'Stars not found close to their previous positions.\n'
                                      'Do you want to skip this frame?',
                                      parent=root)
            else:
                stars_detected = True
                delta_skip_time = time.time()
                skip_frame = False
                x0 = max_x
                y0 = max_y

        else:
            stars_detected = False
            delta_skip_time = time.time()
            label3.configure(text='     ' + science_file.split(os.sep)[-1] + '     ')
            label5.configure(text='     -%    ')
            label7.configure(text='     -h -m -s     ')
            ax.cla()
            ax.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
                      vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
                      vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])
            ax.axis('off')
            canvas.show()
            root.update()
            skip_frame = askyesno('Alignment',
                                  'Stars not found close to their previous positions.\n'
                                  'Do you want to skip this frame?',
                                  parent=root)

        if not stars_detected and not skip_frame:

            centroids = tools.find_centroids(fits[1].data, x_centre=int(x0), y_centre=int(y0),
                                             mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                             std_limit=2.0, burn_limit=2 * burn_limit, star_std=star_std)

            if len(centroids) > 0:

                tests = []

                for ref_star in centroids:

                    max_x = ref_star[1]
                    max_y = ref_star[2]

                    test = 0

                    for comp in comparisons:

                        check_x = int(max_x + comp[0] * np.cos(u0 + comp[1]))
                        check_y = int(max_y + comp[0] * np.sin(u0 + comp[1]))
                        if 0 < check_x < x_length and 0 < check_y < y_length:
                            check_sum = np.sum(fits[1].data[check_y - star_std:check_y + star_std + 1,
                                               check_x - star_std:check_x + star_std + 1])
                            check_lim = (fits[1].header[mean_key] + fits[1].header[std_key]) * ((2 * star_std + 1) ** 2)
                            if check_sum > check_lim:
                                test += 1
                            else:
                                test -= 1

                        if abs(test) > 0.1 * len(comparisons):
                            break

                    tests.append([test, max_x, max_y])

                    if test > 0.1 * len(comparisons):
                        break

                tests.sort()
                test, max_x, max_y = tests[-1]
                if test < 0.1 * len(comparisons):
                    stars_detected = False
                else:
                    stars_detected = True
                    x0 = max_x
                    y0 = max_y
                    delta_skip_time = time.time() - delta_skip_time
                    skip_time += delta_skip_time

            else:
                stars_detected = False

            # if the wider search works then update the values, otherwise continue to the rotation search
            if not stars_detected:

                ustep = np.arcsin(float(star_std) / comparisons[-1][0])
                centroids = tools.find_centroids(fits[1].data,
                                                 x_centre=x0, y_centre=y0,
                                                 mean=fits[1].header[mean_key], std=fits[1].header[std_key],
                                                 std_limit=2.0, burn_limit=2 * burn_limit, star_std=star_std,
                                                 flux_order=True)

                tests = []

                angles = np.array([np.pi, 0])
                for ff in np.arange(ustep, np.pi, ustep):
                    angles = np.append(angles, np.pi - ff)
                    angles = np.append(angles, np.pi + ff)

                for rotation in angles:

                    for ref_star in centroids:

                        max_x = ref_star[1]
                        max_y = ref_star[2]

                        test = 0

                        for comp in comparisons:

                            check_x = int(max_x + comp[0] * np.cos(rotation + comp[1]))
                            check_y = int(max_y + comp[0] * np.sin(rotation + comp[1]))
                            if 0 < check_x < x_length and 0 < check_y < y_length:
                                check_sum = np.sum(fits[1].data[check_y - star_std:check_y + star_std + 1,
                                                   check_x - star_std:check_x + star_std + 1])
                                check_lim = (fits[1].header[mean_key] +
                                             fits[1].header[std_key]) * ((2 * star_std + 1) ** 2)
                                if check_sum > check_lim:
                                    test += 1
                                else:
                                    test -= 1

                            if abs(test) > 0.1 * len(comparisons):
                                break

                        tests.append([test, max_x, max_y, rotation])

                        if test > 0.1 * len(comparisons):
                            break

                    if test > 0.1 * len(comparisons):
                            break

                tests.sort()
                test, max_x, max_y, rotation = tests[-1]
                if test < 0.1 * len(comparisons):
                    ax.cla()
                    ax.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
                              vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
                              vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])
                    ax.axis('off')

                    circle = mpatches.Circle((max_x, max_y), 2 * search_window_std * star_std, ec='r', fill=False)
                    ax.add_patch(circle)
                    for ii in comparisons:
                        circle = mpatches.Circle((max_x + ii[0] * np.cos(rotation + ii[1]),
                                                  max_y + ii[0] * np.sin(rotation + ii[1])),
                                                 2 * search_window_std * star_std, ec='w', fill=False)
                        ax.add_patch(circle)

                    ax.axis('off')
                    canvas.show()
                    root.update()

                    if not askyesno('Alignment',
                                    'Less than 10% of the stars detected.\nDo you accept this alignment?',
                                    parent=root):
                        stars_detected = False
                        delta_skip_time = time.time() - delta_skip_time
                        skip_time += delta_skip_time
                    else:
                        stars_detected = True
                        x0 = max_x
                        y0 = max_y
                        u0 = rotation
                        comparisons = comparisons[:test]
                        delta_skip_time = time.time() - delta_skip_time
                        skip_time += delta_skip_time
                else:
                    stars_detected = True
                    x0 = max_x
                    y0 = max_y
                    u0 = rotation
                    delta_skip_time = time.time() - delta_skip_time
                    skip_time += delta_skip_time

        if stars_detected:

            norm, floor, x_mean, y_mean, x_sigma, y_sigma = \
                tools.fit_2d_gauss(fits[1].data, predicted_x_mean=x0, predicted_y_mean=y0,
                                   search_window=2 * star_std)

            x0, y0 = x_mean, y_mean

            fits[1].header.set(align_x0_key, x0)
            fits[1].header.set(align_y0_key, y0)
            fits[1].header.set(align_u0_key, u0)
            fits[1].header.set(align_star_area_key, np.sqrt(x_sigma ** 2 + y_sigma ** 2))

            ax.cla()
            ax.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
                      vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
                      vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])
            ax.axis('off')

            circle = mpatches.Circle((x0, y0), 2 * search_window_std * star_std, ec='r', fill=False)
            ax.add_patch(circle)
            for ii in comparisons[:int(0.1 * len(comparisons))]:
                circle = mpatches.Circle((x0 + ii[0] * np.cos(u0 + ii[1]), y0 + ii[0] * np.sin(u0 + ii[1])),
                                         2 * search_window_std * star_std, ec='w', fill=False)
                ax.add_patch(circle)

            canvas.draw()

        else:

            fits[1].header.set(align_x0_key, False)
            fits[1].header.set(align_y0_key, False)
            fits[1].header.set(align_u0_key, False)
            fits[1].header.set(align_star_area_key, False)

        fits.flush()
        fits.close()

        # counter
        new_percent = round(100 * (counter + 1) / float(len(science)), 1)
        if new_percent != percent:
            lt1 = time.time()
            rm_time = (100 - new_percent) * (lt1 - lt0 - skip_time) / new_percent
            hours = rm_time / 3600.0
            minutes = (hours - int(hours)) * 60
            seconds = (minutes - int(minutes)) * 60
            label3.configure(text='     ' + science_file.split(os.sep)[-1] + '     ')
            label5.configure(text='     {0}%    '.format(new_percent))
            label7.configure(text='     %dh %02dm %02ds     ' % (int(hours), int(minutes), int(seconds)))
            percent = new_percent

        if counter == 0:
            finalise_window(root, topmost=True)

        root.update()

        if exit_var.get():
            break

        if counter + 1 == len(science):
            write_log('pipeline', True, 'alignment_complete')

    root.destroy()
