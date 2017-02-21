from hops_basics import *


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
    calibration_stars_number = read_log('alignment', 'calibration_stars_number')
    shift_tolerance = read_log('alignment', 'shift_tolerance')

    science = glob.glob('{0}{1}*.f*t*'.format(reduction_directory, os.sep))
    science.sort()

    fits = pf.open(science[0], memmap=False)

    centroids = []
    y_length, x_length = fits[0].data.shape

    frame_part = 10.5
    star_std_update = False
    while len(centroids) < 1 + calibration_stars_number and frame_part >= 2.5:

        std_limit = 9.0
        while len(centroids) < 1 + calibration_stars_number and std_limit >= 3:

            frame_limit = int(min(x_length / frame_part, y_length / frame_part))
            centroids = tools.find_centroids(fits[0].data,
                                             x_low=x_length / 2 - frame_limit, x_upper=x_length / 2 + frame_limit,
                                             y_low=y_length / 2 - frame_limit, y_upper=y_length / 2 + frame_limit,
                                             x_centre=x_length / 2, y_centre=y_length / 2,
                                             mean=fits[0].header[mean_key], std=fits[0].header[std_key],
                                             std_limit=std_limit, burn_limit=7.0 * burn_limit / 8, star_std=star_std)

            if len(centroids) > 0 and not star_std_update:

                calibration_centroids = [[-pp[3], pp[1], pp[2], pp[0]] for pp in centroids]

                new_star_std = []

                for calibration_centroid in calibration_centroids:
                    norm, floor, x_mean, y_mean, x_sigma, y_sigma = \
                        tools.fit_2d_gauss(fits[0].data,
                                           predicted_x_mean=calibration_centroid[1],
                                           predicted_y_mean=calibration_centroid[2],
                                           search_window=search_window_std * star_std)
                    new_star_std.append(x_sigma)

                star_std = max(1, int(np.mean(new_star_std)) - 1)

                centroids = tools.find_centroids(fits[0].data,
                                                 x_low=x_length / 2 - frame_limit, x_upper=x_length / 2 + frame_limit,
                                                 y_low=y_length / 2 - frame_limit, y_upper=y_length / 2 + frame_limit,
                                                 x_centre=x_length / 2, y_centre=y_length / 2,
                                                 mean=fits[0].header[mean_key], std=fits[0].header[std_key],
                                                 std_limit=std_limit, burn_limit=7.0 * burn_limit / 8,
                                                 star_std=star_std)

                write_log('alignment', star_std, 'star_std')
                star_std_update = True

            std_limit -= 2.0

        frame_part -= 2.0

    # reverse and short by brightness
    calibration_stars = [[-pp[3], pp[1], pp[2], pp[0]] for pp in centroids]
    calibration_stars.sort()
    calibration_stars = calibration_stars[:calibration_stars_number + 1]

    # take the first as reference star
    norm, floor, x_mean, y_mean, x_sigma, y_sigma = \
        tools.fit_2d_gauss(fits[0].data,
                           predicted_x_mean=calibration_stars[0][1],
                           predicted_y_mean=calibration_stars[0][2],
                           search_window=search_window_std * star_std)
    x_ref_position, y_ref_position = x_mean, y_mean
    del calibration_stars[0]

    # take the rest as calibration stars and calculate their polar coordinates relatively to the first
    calibration_stars_polar = []
    for calibration_star in calibration_stars:
        norm, floor, x_mean, y_mean, x_sigma, y_sigma = \
            tools.fit_2d_gauss(fits[0].data,
                               predicted_x_mean=calibration_star[1],
                               predicted_y_mean=calibration_star[2],
                               search_window=search_window_std * star_std)
        x_position, y_position = x_mean, y_mean
        r_position, u_position = tools.cartesian_to_polar(x_position, y_position, x_ref_position, y_ref_position)
        calibration_stars_polar.append([r_position, u_position])

    calibration_stars_polar.sort()

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
    ax.imshow(fits[0].data, origin='lower', cmap=cm.Greys_r,
              vmin=fits[0].header[mean_key] + frame_low_std * fits[0].header[std_key],
              vmax=fits[0].header[mean_key] + frame_upper_std * fits[0].header[std_key])
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

        ax.cla()
        ax.imshow(fits[0].data, origin='lower', cmap=cm.Greys_r,
                  vmin=fits[0].header[mean_key] + frame_low_std * fits[0].header[std_key],
                  vmax=fits[0].header[mean_key] + frame_upper_std * fits[0].header[std_key])
        ax.axis('off')

        # fast detection test
        centroids = tools.find_centroids(fits[0].data,
                                         x_low=int(x0 - shift_tolerance), x_upper=int(x0 + shift_tolerance + 1),
                                         y_low=int(y0 - shift_tolerance), y_upper=int(y0 + shift_tolerance + 1),
                                         x_centre=int(x0), y_centre=int(y0),
                                         mean=fits[0].header[mean_key], std=fits[0].header[std_key],
                                         std_limit=1.5, burn_limit=burn_limit, star_std=star_std)

        if len(centroids) > 0:
            stars_detected = True

            max_x = centroids[0][1]
            max_y = centroids[0][2]

            for comp in comparisons:
                check_x = int(max_x + comp[0] * np.cos(u0 + comp[1]))
                check_y = int(max_y + comp[0] * np.sin(u0 + comp[1]))
                check_sum = np.sum(fits[0].data[check_y - star_std:check_y + star_std + 1,
                                   check_x - star_std:check_x + star_std + 1])
                check_lim = (fits[0].header[mean_key] + fits[0].header[std_key]) * ((2 * star_std + 1) ** 2)
                if check_sum < check_lim:
                    stars_detected = False
                    break
        else:
            stars_detected = False

        # if the fast detection test works then update the values, otherwise continue to the wider search
        if stars_detected:
            x0 = max_x
            y0 = max_y

        else:

            delta_skip_time = time.time()

            label3.configure(text='     ' + science_file.split(os.sep)[-1] + '     ')
            label5.configure(text='     -%    ')
            label7.configure(text='     -h -m -s     ')

            canvas.show()
            root.update()

            if not askyesno('Alignment',
                            'Stars not found close to their previous positions.\nDo you want to skip this frame?',
                            parent=root):

                centroids = tools.find_centroids(fits[0].data,
                                                 x_centre=x0, y_centre=y0,
                                                 mean=fits[0].header[mean_key], std=fits[0].header[std_key],
                                                 std_limit=1.5, burn_limit=2 * burn_limit, star_std=star_std)

                for ref_star in centroids:

                    stars_detected = True

                    max_x = ref_star[1]
                    max_y = ref_star[2]

                    for comp in comparisons:
                        check_x = int(max_x + comp[0] * np.cos(u0 + comp[1]))
                        check_y = int(max_y + comp[0] * np.sin(u0 + comp[1]))
                        check_sum = np.sum(fits[0].data[check_y - star_std:check_y + star_std + 1,
                                           check_x - star_std:check_x + star_std + 1])
                        check_lim = (fits[0].header[mean_key] + fits[0].header[std_key]) * ((2 * star_std + 1) ** 2)
                        if check_sum < check_lim:
                            stars_detected = False
                            break

                    if stars_detected:
                        x0 = max_x
                        y0 = max_y
                        break

                # if the wider search works then update the values, otherwise continue to the rotation search
                if not stars_detected:

                    ustep = np.arcsin(float(star_std) / comparisons[-1][0])
                    centroids = tools.find_centroids(fits[0].data,
                                                     x_centre=x0, y_centre=y0,
                                                     mean=fits[0].header[mean_key], std=fits[0].header[std_key],
                                                     std_limit=1.5, burn_limit=2 * burn_limit, star_std=star_std,
                                                     flux_order=True)

                    for ref_star in centroids:

                        max_x = ref_star[1]
                        max_y = ref_star[2]

                        for rotation in np.arange(0, 2.0 * np.pi, ustep):

                            stars_detected = True

                            for comp in comparisons:
                                check_x = int(max_x + comp[0] * np.cos(rotation + comp[1]))
                                check_y = int(max_y + comp[0] * np.sin(rotation + comp[1]))
                                check_sum = np.sum(fits[0].data[check_y - star_std:check_y + star_std + 1,
                                                   check_x - star_std:check_x + star_std + 1])
                                check_lim = (fits[0].header[mean_key] +
                                             fits[0].header[std_key]) * ((2 * star_std + 1) ** 2)
                                if check_sum < check_lim:
                                    stars_detected = False
                                    break

                            if stars_detected:
                                x0 = max_x
                                y0 = max_y
                                u0 = rotation
                                break

                        if stars_detected:
                            break

            delta_skip_time = time.time() - delta_skip_time
            skip_time += delta_skip_time

        if stars_detected:

            norm, floor, x_mean, y_mean, x_sigma, y_sigma = \
                tools.fit_2d_gauss(fits[0].data,
                                   predicted_x_mean=x0, predicted_y_mean=y0,
                                   search_window=search_window_std * star_std)

            x0, y0 = x_mean, y_mean

            fits[0].header.set(align_x0_key, x0)
            fits[0].header.set(align_y0_key, y0)
            fits[0].header.set(align_u0_key, u0)
            fits[0].header.set(align_star_area_key, np.sqrt(x_sigma ** 2 + y_sigma ** 2))
            fits.flush()
            fits.close()

            circle = mpatches.Circle((x0, y0), 2 * search_window_std * star_std, ec='r', fill=False)
            ax.add_patch(circle)
            for ii in comparisons:
                circle = mpatches.Circle((x0 + ii[0] * np.cos(u0 + ii[1]), y0 + ii[0] * np.sin(u0 + ii[1])),
                                         2 * search_window_std * star_std, ec='w', fill=False)
                ax.add_patch(circle)

            canvas.draw()

        else:

            fits[0].header.set(align_x0_key, False)
            fits[0].header.set(align_y0_key, False)
            fits[0].header.set(align_u0_key, False)
            fits[0].header.set(align_star_area_key, False)
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
