from hops_basics import *


def reduction():

    if read_log('pipeline', 'reduction_complete'):
        if not askyesno('Overwrite files', 'Reduction has been completed, do you want to run again?'):
            return 0

    write_log('pipeline', False, 'reduction_complete')
    write_log('pipeline', False, 'alignment_complete')

    # get variables

    observation_files = read_log('pipeline', 'observation_files')
    reduction_directory = read_log('pipeline', 'reduction_directory')
    reduction_prefix = read_log('pipeline', 'reduction_prefix')
    exposure_time_key = read_log('pipeline_keywords', 'exposure_time_key')
    mean_key = read_log('pipeline_keywords', 'mean_key')
    std_key = read_log('pipeline_keywords', 'std_key')
    observation_date_key = read_log('pipeline_keywords', 'observation_date_key')
    observation_time_key = read_log('pipeline_keywords', 'observation_time_key')
    frame_low_std = read_log('windows', 'frame_low_std')
    frame_upper_std = read_log('windows', 'frame_upper_std')
    bias_files = read_log('reduction', 'bias_files')
    dark_files = read_log('reduction', 'dark_files')
    flat_files = read_log('reduction', 'flat_files')
    bin_to = read_log('reduction', 'bin_to')
    master_bias_method = read_log('reduction', 'master_bias_method')
    master_dark_method = read_log('reduction', 'master_dark_method')
    master_flat_method = read_log('reduction', 'master_flat_method')

    # check if reduction directory exists

    if os.path.isdir(reduction_directory):
        shutil.rmtree(reduction_directory)

    os.mkdir(reduction_directory)

    # create master bias

    bias_frames = []
    if len(bias_files) > 0:
        for bias_file in glob.glob('*{0}*.f*t*'.format(bias_files)):
            fits = pf.open(bias_file, memmap=False)
            bias_frames.append(fits[0].data)
            fits.close()

    if len(bias_frames) > 0:
        if master_bias_method == 'median':
            master_bias = np.median(bias_frames, 0)
        elif master_bias_method == 'mean':
            master_bias = np.mean(bias_frames, 0)
        else:
            master_bias = np.median(bias_frames, 0)
    else:
        master_bias = 0

    # create master dark

    dark_frames = []
    if len(str(dark_files)) > 0:
        for dark_file in glob.glob('*{0}*.f*t*'.format(dark_files)):
            fits = pf.open(dark_file, memmap=False)
            dark_frames.append((fits[0].data - master_bias) / fits[0].header[exposure_time_key])
            fits.close()

    if len(dark_frames) > 0:
        if master_dark_method == 'median':
            master_dark = np.median(dark_frames, 0)
        elif master_dark_method == 'mean':
            master_dark = np.mean(dark_frames, 0)
        else:
            master_dark = np.median(dark_frames, 0)
    else:
        master_dark = 0

    # create master flat

    flat_frames = []
    if len(str(flat_files)) > 0:
        for flat_file in glob.glob('*{0}*.f*t*'.format(flat_files)):
            fits = pf.open(flat_file, memmap=False)
            flat_frames.append(fits[0].data - master_bias - fits[0].header[exposure_time_key] * master_dark)
            fits.close()

    if len(flat_frames) > 0:
        if master_flat_method == 'median':
            master_flat = np.median(flat_frames, 0)
        elif master_flat_method == 'mean':
            master_flat = np.mean(flat_frames, 0)
        else:
            master_flat = np.mean(flat_frames, 0)
        master_flat = master_flat / np.median(master_flat)
    else:
        master_flat = 1

    # setup counter window

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

    frame1 = Frame(root)
    frame1.pack()

    label1 = Label(frame1, text='REDUCTION')
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

    # correct each observation_files file

    observation_files = glob.glob('*{0}*.f*t*'.format(observation_files))
    observation_files.sort()
    percent = 0
    lt0 = time.time()
    for counter, science_file in enumerate(observation_files):

        # correct it with master bias_files, master dark_files and master flat_files

        fits = pf.open(science_file, memmap=False)
        fits[0].data = (fits[0].data - master_bias - fits[0].header[exposure_time_key] * master_dark) / master_flat

        norm, floor, mean, std = tools.fit_distribution1d_gaussian(fits[0].data, binning=fits[0].data.size / bin_to)

        if np.isnan(norm):
            mean = np.mean(fits[0].data)
            std = np.std(fits[0].data)

        fits[0].header.set(mean_key, mean)
        fits[0].header.set(std_key, std)

        # write the new fits file

        if observation_date_key == observation_time_key:
                local_time = fits[0].header[observation_date_key]
                local_time = local_time.replace('-', '_').replace('T', '_').replace(':', '_') + '_'
        else:
                local_time = (fits[0].header[observation_date_key].replace('-', '_') + '_' +
                              fits[0].header[observation_time_key].replace(':', '_') + '_')

        hdu = pf.PrimaryHDU(header=fits[0].header, data=fits[0].data)
        hdu.writeto('{0}{1}{2}{3}{4}'.format(reduction_directory,
                                             os.sep, reduction_prefix, local_time, science_file.split(os.sep)[-1]))

        if counter == 0:
            ax.cla()
            ax.imshow(fits[0].data[::2, ::2], origin='lower', cmap=cm.Greys_r,
                      vmin=fits[0].header[mean_key] + frame_low_std * fits[0].header[std_key],
                      vmax=fits[0].header[mean_key] + frame_upper_std * fits[0].header[std_key])
            ax.axis('off')

            canvas.show()

        fits.close()

        # counter

        new_percent = round(100 * (counter + 1) / float(len(observation_files)), 1)
        if new_percent != percent:
            lt1 = time.time()
            rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
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

        if counter + 1 == len(observation_files):
            write_log('pipeline', True, 'reduction_complete')

    root.destroy()
