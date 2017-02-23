__version__ = 1.2

print 'Loading... Please wait for the main window to appear.'

from hops_basics import *
import reduction_routines as rdr
import alignment_routines as alr
import photometry_routines as phr
import fitting_routines as ftr

import exodata
import exodata.astroquantities as aq
from oec import *


def reduction_alignment_window():

    # #########
    # create and initialise the window
    # #########

    root = Tk()
    initialise_window(root)

    # get variables from log and set as tk variables those to be modified

    directory = StringVar(root, value=read_main_log('pipeline', 'directory'))
    directory_short = StringVar(root, value=read_main_log('pipeline', 'directory_short'))
    observation_files = StringVar(root, value=read_main_log('pipeline', 'observation_files'))
    bias_files = StringVar(root, value=read_main_log('reduction', 'bias_files'))
    dark_files = StringVar(root, value=read_main_log('reduction', 'dark_files'))
    flat_files = StringVar(root, value=read_main_log('reduction', 'flat_files'))
    exposure_time_key = StringVar(root, value=read_main_log('pipeline_keywords', 'exposure_time_key'))
    observation_date_key = StringVar(root, value=read_main_log('pipeline_keywords', 'observation_date_key'))
    observation_time_key = StringVar(root, value=read_main_log('pipeline_keywords', 'observation_time_key'))
    target_ra_dec = StringVar(root, value=read_main_log('photometry', 'target_ra_dec'))

    # set progress variables, useful for updating the window

    open_root2 = BooleanVar(root, value=False)
    open_root3 = BooleanVar(root, value=False)
    update_directory = BooleanVar(root, value=False)
    running = BooleanVar(root, value=False)

    # create widgets

    directory_label = Label(root, text='Directory')
    directory_entry = Button(root, textvariable=directory_short)

    observation_files_label = Label(root, text='Name identifier for observation files')
    observation_files_entry = Entry(root, textvariable=observation_files)
    observation_files_test = Label(root, text=' ')

    bias_files_label = Label(root, text='Name identifier for bias files')
    bias_files_entry = Entry(root, textvariable=bias_files)
    bias_files_test = Label(root, text=' ')

    dark_files_label = Label(root, text='Name identifier for dark files')
    dark_files_entry = Entry(root, textvariable=dark_files)
    dark_files_test = Label(root, text=' ')

    flat_files_label = Label(root, text='Name identifier for flat files')
    flat_files_entry = Entry(root, textvariable=flat_files)
    flat_files_test = Label(root, text=' ')

    show_files_button = Button(root, text='Show files')

    exposure_time_key_label = Label(root, text='Exposure time header keyword')
    exposure_time_key_entry = Entry(root, textvariable=exposure_time_key)
    exposure_time_key_test = Label(root, text=' ')

    observation_date_key_label = Label(root, text='Observation date header keyword')
    observation_date_key_entry = Entry(root, textvariable=observation_date_key)
    observation_date_key_test = Label(root, text=' ')

    observation_time_key_label = Label(root, text='Observation time header keyword')
    observation_time_key_entry = Entry(root, textvariable=observation_time_key)
    observation_time_key_test = Label(root, text=' ')

    target_ra_dec_label = Label(root, text='Target RA DEC\n(hh:mm:ss +/-dd:mm:ss)')
    target_ra_dec_entry = Entry(root, textvariable=target_ra_dec)
    target_ra_dec_test = Label(root, text=' ')

    show_header_button = Button(root, text='Show header')

    reduction_alignment_button = Button(root, text='REDUCTION & ALIGNMENT')

    # define the function that updates the window

    def update_window(entry):

        if not entry:
            return 0

        if running.get():

            directory_entry['state'] = DISABLED
            observation_files_entry['state'] = DISABLED
            bias_files_entry['state'] = DISABLED
            dark_files_entry['state'] = DISABLED
            flat_files_entry['state'] = DISABLED
            show_files_button['state'] = DISABLED
            exposure_time_key_entry['state'] = DISABLED
            observation_date_key_entry['state'] = DISABLED
            observation_time_key_entry['state'] = DISABLED
            target_ra_dec_entry['state'] = DISABLED
            show_header_button['state'] = DISABLED
            reduction_alignment_button['state'] = DISABLED

        elif not os.path.isdir(directory.get()):

            directory_entry['state'] = NORMAL
            observation_files_entry['state'] = DISABLED
            bias_files_entry['state'] = DISABLED
            dark_files_entry['state'] = DISABLED
            flat_files_entry['state'] = DISABLED
            show_files_button['state'] = DISABLED
            exposure_time_key_entry['state'] = DISABLED
            observation_date_key_entry['state'] = DISABLED
            observation_time_key_entry['state'] = DISABLED
            target_ra_dec_entry['state'] = DISABLED
            show_header_button['state'] = DISABLED
            reduction_alignment_button['state'] = DISABLED

        else:

            directory_entry['state'] = NORMAL

            if update_directory.get():

                os.chdir(directory.get())
                copy_main_log()

                observation_files.set(read_log('pipeline', 'observation_files'))
                bias_files.set(read_log('reduction', 'bias_files'))
                dark_files.set(read_log('reduction', 'dark_files'))
                flat_files.set(read_log('reduction', 'flat_files'))
                exposure_time_key.set(read_log('pipeline_keywords', 'exposure_time_key'))
                observation_date_key.set(read_log('pipeline_keywords', 'observation_date_key'))
                observation_time_key.set(read_log('pipeline_keywords', 'observation_time_key'))
                target_ra_dec.set(read_log('photometry', 'target_ra_dec'))

            observation_files_entry['state'] = NORMAL
            bias_files_entry['state'] = NORMAL
            dark_files_entry['state'] = NORMAL
            flat_files_entry['state'] = NORMAL

            if open_root2.get():
                show_files_button['state'] = DISABLED
            else:
                show_files_button['state'] = NORMAL

            check_science = test_file_number(observation_files_entry.get())
            observation_files_test.configure(text=check_science[1])

            check_bias = test_file_number(bias_files_entry.get())
            bias_files_test.configure(text=check_bias[1])

            check_dark = test_file_number(dark_files_entry.get())
            dark_files_test.configure(text=check_dark[1])

            check_flat = test_file_number(flat_files_entry.get())
            flat_files_test.configure(text=check_flat[1])

            if not check_science[0]:

                exposure_time_key_entry['state'] = DISABLED
                observation_date_key_entry['state'] = DISABLED
                observation_time_key_entry['state'] = DISABLED
                target_ra_dec_entry['state'] = DISABLED
                show_header_button['state'] = DISABLED
                reduction_alignment_button['state'] = DISABLED

            else:

                target_ra_dec_entry['state'] = NORMAL
                exposure_time_key_entry['state'] = NORMAL
                observation_date_key_entry['state'] = NORMAL
                observation_time_key_entry['state'] = NORMAL

                if open_root3.get():
                    show_header_button['state'] = DISABLED
                else:
                    show_header_button['state'] = NORMAL

                check_ra_dec = test_coordinates(target_ra_dec_entry.get())
                target_ra_dec_test.configure(text=check_ra_dec[1])

                check_exposure_time = test_fits_keyword(observation_files_entry.get(),
                                                        exposure_time_key_entry.get())
                exposure_time_key_test.configure(text=check_exposure_time[1])

                check_observation_date = test_fits_keyword(observation_files_entry.get(),
                                                           observation_date_key_entry.get())
                observation_date_key_test.configure(text=check_observation_date[1])

                check_observation_time = test_fits_keyword(observation_files_entry.get(),
                                                           observation_time_key_entry.get())
                observation_time_key_test.configure(text=check_observation_time[1])

                if (check_ra_dec[0] and check_exposure_time[0] and
                        check_observation_date[0] and check_observation_time[0]
                        and not open_root2.get() and not open_root3.get()):

                    reduction_alignment_button['state'] = NORMAL

                else:

                    reduction_alignment_button['state'] = DISABLED

        root.update()

    update_directory = BooleanVar(root, value=True)
    update_window('a')
    update_directory = BooleanVar(root, value=False)

    # define actions for the different buttons, including calls to the function that updates the window

    def choose_directory():

        new_directory = tkFileDialog.askdirectory()

        if len(new_directory) > 0:
            directory.set(new_directory)
            directory_short.set('/'.join(new_directory.split('/')[-2:]))
            write_main_log('pipeline', directory.get(), 'directory')
            write_main_log('pipeline', directory_short.get(), 'directory_short')

            update_directory.set(True)
            update_window('a')
            update_directory.set(False)

    def show_content():

        open_root2.set(True)
        update_window('a')

        root2 = Tk()
        root2.geometry('{0}x{1}'.format(root2.winfo_screenwidth() / 3, int(0.8 * root2.winfo_screenheight())))

        def root2_close():
            root2.destroy()
            open_root2.set(False)
            update_window('a')

        initialise_window(root2, window_name='Files list', exit_command=root2_close)

        scrollbar = Scrollbar(root2)
        scrollbar.pack(side=RIGHT, fill=Y)

        files_list = Listbox(root2, yscrollcommand=scrollbar.set, font='Courier')
        files_list.insert(END, '  List of files in your directory:')
        files_list.insert(END, '  ')

        for ii in glob.glob('*.f*t*'):
            files_list.insert(END, '  ' + str(ii).split(os.sep)[-1])

        files_list.pack(side=LEFT, fill=BOTH, expand=True)
        scrollbar.config(command=files_list.yview)

        finalise_window(root2, center=False, topmost=True)
        root2.mainloop()

    def show_header():

        open_root3.set(True)
        update_window('a')

        root3 = Tk()
        root3.geometry('{0}x{1}'.format(root3.winfo_screenwidth() / 3, int(0.8 * root3.winfo_screenheight())))

        def root3_close():
            root3.destroy()
            open_root3.set(False)
            update_window('a')

        initialise_window(root3, window_name='Fits file header', exit_command=root3_close)

        scrollbar = Scrollbar(root3)
        scrollbar.pack(side=RIGHT, fill=Y)

        science_header = pf.open(glob.glob('*{0}*.f*t*'.format(observation_files_entry.get()))[0])[0].header

        header_list = Listbox(root3, yscrollcommand=scrollbar.set, font='Courier')
        header_list.insert(END, '  Keywords:      Values:')
        header_list.insert(END, '  ')

        for ii in science_header:
            if ii != '':
                header_list.insert(END, '  ' + str(ii[:10]) +
                                   ' ' * (15 - len(str(ii[:10]))) + str(science_header[ii]))

        header_list.pack(side=LEFT, fill=BOTH, expand=True)
        scrollbar.config(command=header_list.yview)

        finalise_window(root3, center=False, topmost=True)
        root3.mainloop()

    def reduction_alignment():

        running.set(True)
        update_window('a')

        write_log('pipeline', directory.get(), 'directory')
        write_log('pipeline', directory_short.get(), 'directory_short')
        write_log('pipeline', observation_files.get(), 'observation_files')
        write_log('reduction', bias_files.get(), 'bias_files')
        write_log('reduction', dark_files.get(), 'dark_files')
        write_log('reduction', flat_files.get(), 'flat_files')
        write_log('photometry', target_ra_dec.get(), 'target_ra_dec')
        write_log('pipeline_keywords', exposure_time_key.get(), 'exposure_time_key')
        write_log('pipeline_keywords', observation_date_key.get(), 'observation_date_key')
        write_log('pipeline_keywords', observation_time_key.get(), 'observation_time_key')

        rdr.reduction()

        if read_log('pipeline', 'reduction_complete'):
            alr.alignment()

        if read_log('pipeline', 'reduction_complete') and read_log('pipeline', 'alignment_complete'):
            root.destroy()
        else:
            running.set(False)
            update_window('a')

    # connect actions to widgets

    directory_entry['command'] = choose_directory
    observation_files_entry.bind(sequence='<KeyRelease>', func=update_window)
    bias_files_entry.bind(sequence='<KeyRelease>', func=update_window)
    dark_files_entry.bind(sequence='<KeyRelease>', func=update_window)
    flat_files_entry.bind(sequence='<KeyRelease>', func=update_window)
    show_files_button['command'] = show_content
    exposure_time_key_entry.bind(sequence='<KeyRelease>', func=update_window)
    observation_date_key_entry.bind(sequence='<KeyRelease>', func=update_window)
    observation_time_key_entry.bind(sequence='<KeyRelease>', func=update_window)
    target_ra_dec_entry.bind(sequence='<KeyRelease>', func=update_window)
    show_header_button['command'] = show_header
    reduction_alignment_button['command'] = reduction_alignment

    # setup window

    photo = PhotoImage(file=holomon_logo)
    logo_label = Label(root, image=photo)
    window_label = Label(root, text='Reduction & Alignment')
    created_by_label = Label(root, text=read_main_log('windows', 'created_by').replace(',', '\n'))

    setup_window(root, [
        [],
        [[logo_label, 0, 1, 6], [window_label, 1, 3, 1, 'title']],
        [],
        [[directory_label, 1], [directory_entry, 2]],
        [],
        [[observation_files_label, 1], [observation_files_entry, 2], [observation_files_test, 3]],
        [[bias_files_label, 1], [bias_files_entry, 2], [bias_files_test, 3]],
        [[created_by_label, 0, 1, 3], [dark_files_label, 1], [dark_files_entry, 2], [dark_files_test, 3]],
        [[flat_files_label, 1], [flat_files_entry, 2], [flat_files_test, 3]],
        [],
        [[show_files_button, 2]],
        [],
        [[target_ra_dec_label, 1], [target_ra_dec_entry, 2], [target_ra_dec_test, 3]],
        [[exposure_time_key_label, 1], [exposure_time_key_entry, 2], [exposure_time_key_test, 3]],
        [[observation_date_key_label, 1], [observation_date_key_entry, 2], [observation_date_key_test, 3]],
        [[observation_time_key_label, 1], [observation_time_key_entry, 2], [observation_time_key_test, 3]],
        [],
        [[show_header_button, 2]],
        [],
        [[reduction_alignment_button, 1, 3]],
        [],
    ])

    # finalise and show window

    finalise_window(root)
    root.mainloop()


def photometry_window():

    # #########
    # create and initialise window
    # #########

    root = Tk()
    root4 = Tk()

    def phot_exit():
        root.destroy()
        root4.destroy()
        os._exit(0)

    def root4_close():

        open_root4.set(False)
        update_window('a')

        root4.withdraw()

    initialise_window(root, exit_command=phot_exit)
    initialise_window(root4, exit_command=root4_close)

    # get variables from log and set as tk variables those to be modified

    observation_files = read_log('pipeline', 'observation_files')
    reduction_directory = read_log('pipeline', 'reduction_directory')
    fov_figure = read_log('pipeline', 'fov_figure')
    mean_key = read_log('pipeline_keywords', 'mean_key')
    std_key = read_log('pipeline_keywords', 'std_key')
    align_x0_key = read_log('pipeline_keywords', 'align_x0_key')
    align_y0_key = read_log('pipeline_keywords', 'align_y0_key')
    frame_low_std = read_log('windows', 'frame_low_std')
    frame_upper_std = read_log('windows', 'frame_upper_std')
    burn_limit = read_log('alignment', 'burn_limit')
    star_std = read_log('alignment', 'star_std')
    search_window_std = read_log('alignment', 'search_window_std')
    max_comparisons = read_log('photometry', 'max_comparisons')
    max_targets = max_comparisons + 1

    targets_x_position = [DoubleVar(root, value=read_log('photometry', 'target_x_position'))]
    for comparison in range(max_comparisons):
        targets_x_position.append(
            DoubleVar(root, value=read_log('photometry', 'comparison_{0}_x_position'.format(comparison + 1))))
        if not targets_x_position[-1].get():
            targets_x_position[-1].set(0)

    targets_y_position = [DoubleVar(root, value=read_log('photometry', 'target_y_position'))]
    for comparison in range(max_comparisons):
        targets_y_position.append(
            DoubleVar(root, value=read_log('photometry', 'comparison_{0}_y_position'.format(comparison + 1))))
        if not targets_y_position[-1].get():
            targets_y_position[-1].set(0)

    targets_aperture = [IntVar(root, value=read_log('photometry', 'target_aperture'))]
    for comparison in range(max_comparisons):
        targets_aperture.append(
            IntVar(root, value=read_log('photometry', 'comparison_{0}_aperture'.format(comparison + 1))))
        if not targets_aperture[-1].get():
            targets_aperture[-1].set(0)

    # set progress variables, useful for updating the tk windows

    open_root4 = BooleanVar(root, value=False)
    finalised_root4 = BooleanVar(root, value=False)
    targets_indication = IntVar(root, value=0)
    click_test_x = DoubleVar(root, value=-100)
    click_test_y = DoubleVar(root, value=-100)
    click_off_axis = BooleanVar(root, value=False)
    running = BooleanVar(root, value=False)

    # create the plot in the additional window

    fits = pf.open(glob.glob('{0}{1}*{2}*.f*t*'.format(reduction_directory, os.sep, observation_files))[0],
                   memmap=False)
    f = Figure()
    f.patch.set_facecolor('white')
    ax = f.add_subplot(111)
    canvas = FigureCanvasTkAgg(f, root4)
    canvas.get_tk_widget().pack()
    NavigationToolbar2TkAgg(canvas, root4)

    ax.imshow(fits[0].data, origin='lower', cmap=cm.Greys_r,
              vmin=fits[0].header[mean_key] + frame_low_std * fits[0].header[std_key],
              vmax=fits[0].header[mean_key] + frame_upper_std * fits[0].header[std_key])

    targets_box = [mpatches.Rectangle((targets_x_position[0].get() - targets_aperture[0].get(),
                                       targets_y_position[0].get() - targets_aperture[0].get()),
                                      2 * targets_aperture[0].get() + 1, 2 * targets_aperture[0].get() + 1,
                                      ec='r', fill=False)]
    for comparison in range(max_comparisons):
        targets_box.append(mpatches.Rectangle((targets_x_position[comparison + 1].get() -
                                               targets_aperture[comparison + 1].get(),
                                               targets_y_position[comparison + 1].get() -
                                               targets_aperture[comparison + 1].get()),
                                              2 * targets_aperture[comparison + 1].get() + 1,
                                              2 * targets_aperture[comparison + 1].get() + 1,
                                              ec='#07fefc', fill=False))

    for target in range(max_targets):
        ax.add_patch(targets_box[target])

    targets_text = [ax.text(targets_x_position[0].get(),
                            targets_y_position[0].get() - targets_aperture[0].get() - 1, 'T',
                            color='r', fontsize=20, va='top')]

    for comparison in range(max_comparisons):
        targets_text.append(ax.text(targets_x_position[comparison + 1].get()
                                    + targets_aperture[comparison + 1].get() + 1,
                                    targets_y_position[comparison + 1].get()
                                    - targets_aperture[comparison + 1].get() - 1,
                                    'C{0}'.format(comparison + 1), color='#07fefc', fontsize=20, va='top'))

    # create widgets

    position_label = Label(root, text='     Position     ')

    box_semi_length_label = Label(root, text='Box semi-length')

    targets_indication_entry = [Radiobutton(root, text='      Target           ', variable=targets_indication, value=0)]
    for comparison in range(max_comparisons):
        targets_indication_entry.append(Radiobutton(root, text='Comparison {0}     '.format(comparison + 1),
                                                    variable=targets_indication, value=comparison + 1))

    targets_x_position_label = [Label(root, textvar=targets_x_position[0])]
    for comparison in range(max_comparisons):
        targets_x_position_label.append(Label(root, textvar=targets_x_position[comparison + 1]))

    targets_y_position_label = [Label(root, textvar=targets_y_position[0])]
    for comparison in range(max_comparisons):
        targets_y_position_label.append(Label(root, textvar=targets_y_position[comparison + 1]))

    targets_aperture_entry = [Entry(root, textvar=targets_aperture[0], validate='key')]
    for comparison in range(max_comparisons):
        targets_aperture_entry.append(Entry(root, textvar=targets_aperture[comparison + 1], validate='key'))

    for target in range(max_targets):
        targets_aperture_entry[target]['validatecommand'] = \
            (targets_aperture_entry[target].register(test_int_positive_non_zero_input), '%P', '%d')

    show_fov_button = Button(root, text='Show FOV')

    photometry_button = Button(root, text='PHOTOMETRY')

    proceed_to_fitting_button = Button(root, text='PROCEED TO FITTING')

    # define the function that updates the window

    def update_window(event):

        if running.get():

            for i_target in range(max_targets):
                targets_indication_entry[i_target]['state'] = DISABLED
                targets_aperture_entry[i_target]['state'] = DISABLED

            show_fov_button['state'] = DISABLED
            photometry_button['state'] = DISABLED
            proceed_to_fitting_button['state'] = DISABLED

        else:

            if open_root4.get():
                show_fov_button['state'] = DISABLED
            else:
                show_fov_button['state'] = NORMAL

            for i_target in range(max_targets):
                targets_indication_entry[i_target]['state'] = NORMAL

            try:
                if event.inaxes is not None:

                    click_off_axis.set(False)

                    if (event.xdata, event.ydata) == (click_test_x.get(), click_test_y.get()):

                        frame_limit = 2 * search_window_std * star_std
                        centroids = tools.find_centroids(fits[0].data,
                                                         x_low=int(event.xdata - frame_limit),
                                                         x_upper=int(event.xdata + frame_limit),
                                                         y_low=int(event.ydata - frame_limit),
                                                         y_upper=int(event.ydata + frame_limit),
                                                         x_centre=int(event.xdata), y_centre=int(event.ydata),
                                                         mean=fits[0].header[mean_key], std=fits[0].header[std_key],
                                                         std_limit=3.0, burn_limit=burn_limit, star_std=star_std)

                        if centroids.size > 0:

                            norm, floor, x_mean, y_mean, x_sigma, y_sigma = \
                                tools.fit_2d_gauss(fits[0].data,
                                                   predicted_x_mean=centroids[0][1], predicted_y_mean=centroids[0][2],
                                                   search_window=search_window_std * star_std)

                            if np.sqrt((x_mean - event.xdata) ** 2 + (y_mean - event.ydata) ** 2) < 3 * x_sigma:

                                targets_x_position[targets_indication.get()].set(round(x_mean, 1))
                                targets_y_position[targets_indication.get()].set(round(y_mean, 1))
                                targets_aperture[targets_indication.get()].set(abs(int(4 * x_sigma)))

                        click_test_x.set(-100)
                        click_test_y.set(-100)

                    else:
                        click_test_x.set(event.xdata)
                        click_test_y.set(event.ydata)

                else:

                    click_test_x.set(-100)
                    click_test_y.set(-100)

                    if click_off_axis.get():

                        targets_x_position[targets_indication.get()].set(0)
                        targets_y_position[targets_indication.get()].set(0)
                        targets_aperture[targets_indication.get()].set(0)

                        click_off_axis.set(False)

                    else:
                        click_off_axis.set(True)

            except AttributeError:
                pass

            try:

                for i_target in range(max_targets):

                    if 0 in [targets_x_position[i_target].get(), targets_y_position[i_target].get()]:

                        targets_box[i_target].set_xy((-10000, -10000))

                        targets_text[i_target].set_x(-10000)
                        targets_text[i_target].set_y(-10000)

                        targets_aperture_entry[i_target]['state'] = DISABLED

                    else:

                        targets_box[i_target].set_xy((targets_x_position[i_target].get() -
                                                      targets_aperture[i_target].get(),
                                                      targets_y_position[i_target].get() -
                                                      targets_aperture[i_target].get()))

                        targets_box[i_target].set_width(2 * targets_aperture[i_target].get() + 1)
                        targets_box[i_target].set_height(2 * targets_aperture[i_target].get() + 1)

                        targets_text[i_target].set_x(targets_x_position[i_target].get() +
                                                     targets_aperture[i_target].get() + 1)
                        targets_text[i_target].set_y(targets_y_position[i_target].get() -
                                                     targets_aperture[i_target].get() - 1)

                        targets_aperture_entry[i_target]['state'] = NORMAL

            except ValueError:
                photometry_button['state'] = DISABLED

            if 0 in [targets_x_position[0].get(), targets_y_position[0].get(),
                     targets_x_position[1].get(), targets_y_position[1].get()]:
                photometry_button['state'] = DISABLED

            elif open_root4.get():
                photometry_button['state'] = DISABLED

            else:
                photometry_button['state'] = NORMAL

            if read_log('pipeline', 'photometry_complete') and not open_root4.get():
                proceed_to_fitting_button['state'] = NORMAL

            else:
                proceed_to_fitting_button['state'] = DISABLED

        canvas.draw()
        root.update()
        root4.update()

    update_window('a')

    # define actions for the different buttons, including calls to the function that updates the window

    def photometry():

        running.set(True)
        update_window('a')

        write_log('photometry', targets_x_position[0].get(), 'target_x_position')
        write_log('photometry', targets_y_position[0].get(), 'target_y_position')
        write_log('photometry', targets_aperture[0].get(), 'target_aperture')
        target_polar = tools.cartesian_to_polar(targets_x_position[0].get(), targets_y_position[0].get(),
                                                fits[0].header[align_x0_key], fits[0].header[align_y0_key])
        write_log('photometry', float(target_polar[0]), 'target_r_position')
        write_log('photometry', float(target_polar[1]), 'target_u_position')

        for i_comparison in range(max_comparisons):
            write_log('photometry', targets_x_position[i_comparison + 1].get(),
                      'comparison_{0}_x_position'.format(i_comparison + 1))
            write_log('photometry', targets_y_position[i_comparison + 1].get(),
                      'comparison_{0}_y_position'.format(i_comparison + 1))
            write_log('photometry', targets_aperture[i_comparison + 1].get(),
                      'comparison_{0}_aperture'.format(i_comparison + 1))

            if 0 not in [targets_x_position[i_comparison + 1].get(), targets_y_position[i_comparison + 1].get()]:

                target_polar = tools.cartesian_to_polar(targets_x_position[i_comparison + 1].get(),
                                                        targets_y_position[i_comparison + 1].get(),
                                                        fits[0].header[align_x0_key], fits[0].header[align_y0_key])

            else:

                target_polar = [0, 0]

            write_log('photometry', float(target_polar[0]), 'comparison_{0}_r_position'.format(i_comparison + 1))
            write_log('photometry', float(target_polar[1]), 'comparison_{0}_u_position'.format(i_comparison + 1))

        f.savefig(fov_figure, dpi=200)
        phr.photometry()

        running.set(False)
        update_window('a')

    def show_fov():

        if not finalised_root4.get():
            finalised_root4.set(True)
            open_root4.set(True)
            update_window('a')
            finalise_window(root4, center=False, topmost=True)

        else:
            open_root4.set(True)
            update_window('a')
            root4.deiconify()

    def proceed_to_fitting():
        root.destroy()
        root4.destroy()

    # connect actions to widgets

    f.canvas.callbacks.connect('button_press_event', update_window)
    show_fov_button['command'] = show_fov
    photometry_button['command'] = photometry
    proceed_to_fitting_button['command'] = proceed_to_fitting

    for target in range(max_targets):
        targets_aperture_entry[target].bind(sequence='<KeyRelease>', func=update_window)

    # setup window

    photo = PhotoImage(file=holomon_logo)
    logo_label = Label(root, image=photo)
    window_label = Label(root, text='Photometry')
    created_by_label = Label(root, text=read_main_log('windows', 'created_by').replace(',', '\n'))

    setup_list = [
        [],
        [[logo_label, 0, 1, 6], [window_label, 1, 4, 1, 'title']],
        [],
        [[position_label, 2, 2], [box_semi_length_label, 4]],
        []
    ]

    for target in range(max_targets):

        if target == 2:
            setup_list.append([[created_by_label, 0, 1, 3],
                               [targets_indication_entry[target], 1], [targets_x_position_label[target], 2],
                               [targets_y_position_label[target], 3], [targets_aperture_entry[target], 4]])
        else:
            setup_list.append([[targets_indication_entry[target], 1], [targets_x_position_label[target], 2],
                               [targets_y_position_label[target], 3], [targets_aperture_entry[target], 4]])

    setup_list.append([])
    setup_list.append([[show_fov_button, 4]])
    setup_list.append([])
    setup_list.append([[photometry_button, 1, 4]])
    setup_list.append([])
    setup_list.append([[proceed_to_fitting_button, 1, 4]])
    setup_list.append([])

    setup_window(root, setup_list)

    # finalise and show window

    finalise_window(root)

    root.mainloop()
    root4.mainloop()


def fitting_window():

    # #########
    # create and initialise fitting window
    # #########

    root = Tk()
    initialise_window(root)

    # get variables from log and set as tk variables those to be modified

    catalogue = exodata.OECDatabase(oec(), stream=True)

    light_curve_file = StringVar(value=read_log('fitting', 'light_curve_file'))
    light_curve_file_short = StringVar(value=os.path.split(light_curve_file.get())[1])
    light_curve_dir_short = StringVar(value=os.path.split(os.path.split(light_curve_file.get())[0])[1])
    planet_search = StringVar(value=read_log('fitting', 'planet_search'))
    planet = StringVar(value=read_log('fitting', 'planet'))
    binning = IntVar(value=read_log('fitting', 'binning'))
    scatter = DoubleVar(value=read_log('fitting', 'scatter'))
    iterations = IntVar(value=read_log('fitting', 'iterations'))
    burn = IntVar(value=read_log('fitting', 'burn'))
    metallicity = DoubleVar(value=read_log('fitting', 'metallicity'))
    temperature = DoubleVar(value=read_log('fitting', 'temperature'))
    logg = DoubleVar(value=read_log('fitting', 'logg'))
    phot_filter = StringVar(value=read_log('fitting', 'phot_filter'))
    period = DoubleVar(value=read_log('fitting', 'period'))
    mid_time = DoubleVar(value=read_log('fitting', 'mid_time'))
    rp_over_rs = DoubleVar(value=read_log('fitting', 'rp_over_rs'))
    sma_over_rs = DoubleVar(value=read_log('fitting', 'sma_over_rs'))
    inclination = DoubleVar(value=read_log('fitting', 'inclination'))
    eccentricity = DoubleVar(value=read_log('fitting', 'eccentricity'))
    periastron = DoubleVar(value=read_log('fitting', 'periastron'))

    # set progress variables, useful for updating the window

    update_planet = BooleanVar(root, value=False)
    running = BooleanVar(root, value=False)
    go_back = BooleanVar(root, value=False)

    # create widgets

    light_curve_file_label = Label(root, text='Light-curve file')
    light_curve_dir_label = Label(root, textvar=light_curve_dir_short)
    light_curve_file_entry = Button(root, textvar=light_curve_file_short)

    binning_label = Label(root, text='Binning')
    binning_entry = Entry(root, textvariable=binning, validate='key')
    binning_entry['validatecommand'] = (binning_entry.register(test_int_positive_non_zero_input), '%P', '%d')

    scatter_label = Label(root, text='Scatter limit')
    scatter_entry = Entry(root, textvariable=scatter, validate='key')
    scatter_entry['validatecommand'] = (scatter_entry.register(test_float_positive_input), '%P', '%d')

    iterations_label = Label(root, text='Iterations')
    iterations_entry = Entry(root, textvariable=iterations, validate='key')
    iterations_entry['validatecommand'] = (iterations_entry.register(test_int_positive_non_zero_input), '%P', '%d')

    burn_label = Label(root, text='Burned iterations')
    burn_entry = Entry(root, textvariable=burn, validate='key')
    burn_entry['validatecommand'] = (burn_entry.register(test_int_positive_non_zero_input), '%P', '%d')

    metallicity_label = Label(root, text='Stellar metallicity [Fe/H, dex]')
    metallicity_entry = Entry(root, textvariable=metallicity, validate='key')
    metallicity_entry['validatecommand'] = (metallicity_entry.register(test_float_input), '%P', '%d')

    temperature_label = Label(root, text='Stellar temperature [K]')
    temperature_entry = Entry(root, textvariable=temperature, validate='key')
    temperature_entry['validatecommand'] = (temperature_entry.register(test_float_positive_input), '%P', '%d')

    logg_label = Label(root, text='Stellar log(g) [cm/s^2]')
    logg_entry = Entry(root, textvariable=logg, validate='key')
    logg_entry['validatecommand'] = (logg_entry.register(test_float_positive_input), '%P', '%d')

    phot_filter_label = Label(root, text='Filter')
    phot_filter_entry = Entry(root, textvariable=phot_filter, validate='key')
    phot_filter_entry['validatecommand'] = (phot_filter_entry.register(test_phot_filter), '%P', '%d')

    period_label = Label(root, text='Period [days]')
    period_entry = Entry(root, textvariable=period, validate='key')
    period_entry['validatecommand'] = (period_entry.register(test_float_positive_input), '%P', '%d')

    mid_time_label = Label(root, text='Mid-time [days, HJD]')
    mid_time_entry = Entry(root, textvariable=mid_time, validate='key')
    mid_time_entry['validatecommand'] = (mid_time_entry.register(test_float_positive_input), '%P', '%d')

    rp_over_rs_label = Label(root, text='Rp/Rs')
    rp_over_rs_entry = Entry(root, textvariable=rp_over_rs, validate='key')
    rp_over_rs_entry['validatecommand'] = (rp_over_rs_entry.register(test_float_positive_input), '%P', '%d')

    sma_over_rs_label = Label(root, text='a/Rs')
    sma_over_rs_entry = Entry(root, textvariable=sma_over_rs, validate='key')
    sma_over_rs_entry['validatecommand'] = (sma_over_rs_entry.register(test_float_positive_input), '%P', '%d')

    inclination_label = Label(root, text='Inclination [deg]')
    inclination_entry = Entry(root, textvariable=inclination, validate='key')
    inclination_entry['validatecommand'] = (inclination_entry.register(test_float_positive_input), '%P', '%d')

    eccentricity_label = Label(root, text='Eccentricity')
    eccentricity_entry = Entry(root, textvariable=eccentricity, validate='key')
    eccentricity_entry['validatecommand'] = (eccentricity_entry.register(test_float_positive_input), '%P', '%d')

    periastron_label = Label(root, text='Periastron [deg]')
    periastron_entry = Entry(root, textvariable=periastron, validate='key')
    periastron_entry['validatecommand'] = (periastron_entry.register(test_float_positive_input), '%P', '%d')

    planet_label = Label(root, text='Planet')
    combostyle = ttk.Style()
    combostyle.theme_create('combostyle', parent='alt',
                            settings={'TCombobox': {'configure':
                                                    {'selectbackground': 'white',
                                                     'fieldbackground': 'white',
                                                     'background': 'white'}}})
    combostyle.theme_use('combostyle')
    planet_entry = ttk.Combobox(root, textvariable=planet, state='readonly')
    planet_search_entry = Entry(root, textvariable=planet_search)

    return_to_photometry_button = Button(root, text='RETURN TO PHOTOMETRY')

    fitting_button = Button(root, text='FITTING')

    exit_hops_button = Button(root, text='EXIT')

    # define the function that updates the window

    def update_window(entry):

        if not entry:
            return 0

        if planet.get() == 'Choose Planet':

            catalogue_planets = []

            ra_dec_target = read_log('photometry', 'target_ra_dec')
            ra_dec_target = ephem.Equatorial(ra_dec_target.split()[0], ra_dec_target.split()[1])
            ra_target, dec_target = ra_dec_target.ra * 180 / np.pi, ra_dec_target.dec * 180 / np.pi

            for catalogue_planet in catalogue.planets:
                if not np.isnan(catalogue_planet.system.dec):
                    catalogue_planets.append([np.sqrt((catalogue_planet.system.dec.deg - dec_target) ** 2
                                                      + (catalogue_planet.system.ra.deg - ra_target) ** 2),
                                              catalogue_planet.name])
            catalogue_planets.sort()

            planet.set(catalogue_planets[0][1])
            planet_search.set(catalogue_planets[0][1])

            catalogue_planet = catalogue.searchPlanet(planet.get())

            if not np.isnan(catalogue_planet.star.Z):
                metallicity.set(catalogue_planet.star.Z)
            else:
                metallicity.set(0.0)
            if not np.isnan(catalogue_planet.star.T):
                temperature.set(float(catalogue_planet.star.T))
            else:
                temperature.set(0.0)
            if not np.isnan(catalogue_planet.star.calcLogg()):
                logg.set(round(float(catalogue_planet.star.calcLogg()), 3))
            else:
                logg.set(0.0)
            if not np.isnan(catalogue_planet.P):
                period.set(float(catalogue_planet.P))
            else:
                period.set(0.0)
            if not np.isnan(catalogue_planet.transittime):
                mid_time.set(float(catalogue_planet.transittime))
            else:
                mid_time.set(0.0)
            if not np.isnan(catalogue_planet.R) and not np.isnan(catalogue_planet.star.R):
                rp_over_rs.set(
                    round(float(catalogue_planet.R.rescale(aq.m) / catalogue_planet.star.R.rescale(aq.m)), 5))
            else:
                rp_over_rs.set(0.0)
            if not np.isnan(catalogue_planet.calcSMA()) and not np.isnan(catalogue_planet.star.R):
                sma_over_rs.set(round(float(catalogue_planet.calcSMA().rescale(aq.m) /
                                            catalogue_planet.star.R.rescale(aq.m)), 5))
            else:
                sma_over_rs.set(0.0)
            if not np.isnan(catalogue_planet.i):
                inclination.set(float(catalogue_planet.i))
            else:
                inclination.set(0.0)
            if not np.isnan(catalogue_planet.e):
                eccentricity.set(float(catalogue_planet.e))
            else:
                eccentricity.set(0.0)
            if not np.isnan(catalogue_planet.periastron):
                periastron.set(float(catalogue_planet.periastron))
            else:
                periastron.set(0.0)

        if running.get():

            light_curve_file_entry['state'] = DISABLED
            binning_entry['state'] = DISABLED
            scatter_entry['state'] = DISABLED
            iterations_entry['state'] = DISABLED
            burn_entry['state'] = DISABLED
            metallicity_entry['state'] = DISABLED
            temperature_entry['state'] = DISABLED
            logg_entry['state'] = DISABLED
            phot_filter_entry['state'] = DISABLED
            period_entry['state'] = DISABLED
            mid_time_entry['state'] = DISABLED
            rp_over_rs_entry['state'] = DISABLED
            sma_over_rs_entry['state'] = DISABLED
            inclination_entry['state'] = DISABLED
            eccentricity_entry['state'] = DISABLED
            periastron_entry['state'] = DISABLED
            planet_entry['state'] = DISABLED
            planet_search_entry['state'] = DISABLED
            return_to_photometry_button['state'] = DISABLED
            fitting_button['state'] = DISABLED
            exit_hops_button['state'] = DISABLED

        elif not os.path.isfile(light_curve_file.get()):

            light_curve_file_entry['state'] = NORMAL
            binning_entry['state'] = DISABLED
            scatter_entry['state'] = DISABLED
            iterations_entry['state'] = DISABLED
            burn_entry['state'] = DISABLED
            metallicity_entry['state'] = DISABLED
            temperature_entry['state'] = DISABLED
            logg_entry['state'] = DISABLED
            phot_filter_entry['state'] = DISABLED
            period_entry['state'] = DISABLED
            mid_time_entry['state'] = DISABLED
            rp_over_rs_entry['state'] = DISABLED
            sma_over_rs_entry['state'] = DISABLED
            inclination_entry['state'] = DISABLED
            eccentricity_entry['state'] = DISABLED
            periastron_entry['state'] = DISABLED
            planet_entry['state'] = DISABLED
            planet_search_entry['state'] = DISABLED
            return_to_photometry_button['state'] = DISABLED
            fitting_button['state'] = DISABLED
            exit_hops_button['state'] = NORMAL

        else:

            light_curve_file_entry['state'] = NORMAL
            binning_entry['state'] = NORMAL
            scatter_entry['state'] = NORMAL
            iterations_entry['state'] = NORMAL
            burn_entry['state'] = NORMAL
            metallicity_entry['state'] = NORMAL
            temperature_entry['state'] = NORMAL
            logg_entry['state'] = NORMAL
            phot_filter_entry['state'] = NORMAL
            period_entry['state'] = NORMAL
            mid_time_entry['state'] = NORMAL
            rp_over_rs_entry['state'] = NORMAL
            sma_over_rs_entry['state'] = NORMAL
            inclination_entry['state'] = NORMAL
            eccentricity_entry['state'] = NORMAL
            periastron_entry['state'] = NORMAL
            planet_entry['state'] = NORMAL
            planet_search_entry['state'] = NORMAL

            if isinstance(catalogue.searchPlanet(planet_search.get()), list):
                planet_entry['values'] = tuple([ppp.name for ppp in catalogue.searchPlanet(planet_search.get())])
            elif catalogue.searchPlanet(planet_search.get()):
                planet_entry['values'] = tuple([catalogue.searchPlanet(planet_search.get()).name])
            else:
                planet_entry['values'] = tuple([])

            if update_planet.get():

                obj_planet = catalogue.searchPlanet(planet.get())

                if not np.isnan(obj_planet.star.Z):
                    metallicity.set(obj_planet.star.Z)
                else:
                    metallicity.set(0.0)
                if not np.isnan(obj_planet.star.T):
                    temperature.set(float(obj_planet.star.T))
                else:
                    temperature.set(0.0)
                if not np.isnan(obj_planet.star.calcLogg()):
                    logg.set(round(float(obj_planet.star.calcLogg()), 3))
                else:
                    logg.set(0.0)
                if not np.isnan(obj_planet.P):
                    period.set(float(obj_planet.P))
                else:
                    period.set(0.0)
                if not np.isnan(obj_planet.transittime):
                    mid_time.set(float(obj_planet.transittime))
                else:
                    mid_time.set(0.0)
                if not np.isnan(obj_planet.R) and not np.isnan(obj_planet.star.R):
                    rp_over_rs.set(round(float(obj_planet.R.rescale(aq.m) / obj_planet.star.R.rescale(aq.m)), 5))
                else:
                    rp_over_rs.set(0.0)
                if not np.isnan(obj_planet.calcSMA()) and not np.isnan(obj_planet.star.R):
                    sma_over_rs.set(round(float(obj_planet.calcSMA().rescale(aq.m) /
                                                obj_planet.star.R.rescale(aq.m)), 5))
                else:
                    sma_over_rs.set(0.0)
                if not np.isnan(obj_planet.i):
                    inclination.set(float(obj_planet.i))
                else:
                    inclination.set(0.0)
                if not np.isnan(obj_planet.e):
                    eccentricity.set(float(obj_planet.e))
                else:
                    eccentricity.set(0.0)
                if not np.isnan(obj_planet.periastron):
                    periastron.set(float(obj_planet.periastron))
                else:
                    periastron.set(0.0)

            enable_buttons = True

            if not os.path.isfile(light_curve_file.get()):
                enable_buttons = False

            for input_entry in [binning_entry, scatter_entry, iterations_entry, burn_entry, phot_filter_entry,
                                metallicity_entry, temperature_entry, logg_entry, period_entry, mid_time_entry,
                                rp_over_rs_entry, sma_over_rs_entry, inclination_entry, eccentricity_entry,
                                periastron_entry]:

                if len(str(input_entry.get())) == 0:
                    enable_buttons = False

            if enable_buttons:
                fitting_button['state'] = NORMAL

            else:
                fitting_button['state'] = DISABLED

            return_to_photometry_button['state'] = NORMAL
            exit_hops_button['state'] = NORMAL

        planet_entry.selection_clear()

    update_window('a')

    # define actions for the different buttons, including calls to the function that updates the window

    def choose_light_curve_file():

        new_light_curve_file = tkFileDialog.askopenfilename()

        if len(new_light_curve_file) > 0:
            light_curve_file.set(new_light_curve_file)
            light_curve_dir_short.set(os.path.split(os.path.split(light_curve_file.get())[0])[1])
            light_curve_file_short.set(os.path.split(light_curve_file.get())[1])
            update_window('a')

    def choose_planet(entry):

        if not entry:
            return 0

        update_planet.set(True)
        update_window('a')
        update_planet.set(False)

    def return_to_photometry():

        write_log('fitting', light_curve_file.get(), 'light_curve_file')
        write_log('fitting', light_curve_file_short.get(), 'light_curve_file_short')
        write_log('fitting', planet_search.get(), 'planet_search')
        write_log('fitting', planet.get(), 'planet')
        write_log('fitting', binning.get(), 'binning')
        write_log('fitting', scatter.get(), 'scatter')
        write_log('fitting', iterations.get(), 'iterations')
        write_log('fitting', burn.get(), 'burn')
        write_log('fitting', metallicity.get(), 'metallicity')
        write_log('fitting', temperature.get(), 'temperature')
        write_log('fitting', logg.get(), 'logg')
        write_log('fitting', phot_filter.get(), 'phot_filter')
        write_log('fitting', period.get(), 'period')
        write_log('fitting', mid_time.get(), 'mid_time')
        write_log('fitting', rp_over_rs.get(), 'rp_over_rs')
        write_log('fitting', sma_over_rs.get(), 'sma_over_rs')
        write_log('fitting', inclination.get(), 'inclination')
        write_log('fitting', eccentricity.get(), 'eccentricity')
        write_log('fitting', periastron.get(), 'periastron')

        go_back.set(True)
        root.destroy()

    def fitting():

        running.set(True)
        update_window('a')

        write_log('fitting', light_curve_file.get(), 'light_curve_file')
        write_log('fitting', light_curve_file_short.get(), 'light_curve_file_short')
        write_log('fitting', planet_search.get(), 'planet_search')
        write_log('fitting', planet.get(), 'planet')
        write_log('fitting', binning.get(), 'binning')
        write_log('fitting', scatter.get(), 'scatter')
        write_log('fitting', iterations.get(), 'iterations')
        write_log('fitting', burn.get(), 'burn')
        write_log('fitting', metallicity.get(), 'metallicity')
        write_log('fitting', temperature.get(), 'temperature')
        write_log('fitting', logg.get(), 'logg')
        write_log('fitting', phot_filter.get(), 'phot_filter')
        write_log('fitting', period.get(), 'period')
        write_log('fitting', mid_time.get(), 'mid_time')
        write_log('fitting', rp_over_rs.get(), 'rp_over_rs')
        write_log('fitting', sma_over_rs.get(), 'sma_over_rs')
        write_log('fitting', inclination.get(), 'inclination')
        write_log('fitting', eccentricity.get(), 'eccentricity')
        write_log('fitting', periastron.get(), 'periastron')

        ftr.fitting()

        running.set(False)
        update_window('a')

    def exit_hops():
        root.destroy()

    light_curve_file_entry['command'] = choose_light_curve_file
    planet_entry.bind('<<ComboboxSelected>>', choose_planet)
    planet_search_entry.bind(sequence='<KeyRelease>', func=update_window)
    binning_entry.bind(sequence='<KeyRelease>', func=update_window)
    scatter_entry.bind(sequence='<KeyRelease>', func=update_window)
    iterations_entry.bind(sequence='<KeyRelease>', func=update_window)
    burn_entry.bind(sequence='<KeyRelease>', func=update_window)
    phot_filter_entry.bind(sequence='<KeyRelease>', func=update_window)
    metallicity_entry.bind(sequence='<KeyRelease>', func=update_window)
    temperature_entry.bind(sequence='<KeyRelease>', func=update_window)
    logg_entry.bind(sequence='<KeyRelease>', func=update_window)
    period_entry.bind(sequence='<KeyRelease>', func=update_window)
    mid_time_entry.bind(sequence='<KeyRelease>', func=update_window)
    rp_over_rs_entry.bind(sequence='<KeyRelease>', func=update_window)
    sma_over_rs_entry.bind(sequence='<KeyRelease>', func=update_window)
    inclination_entry.bind(sequence='<KeyRelease>', func=update_window)
    eccentricity_entry.bind(sequence='<KeyRelease>', func=update_window)
    periastron_entry.bind(sequence='<KeyRelease>', func=update_window)
    return_to_photometry_button['command'] = return_to_photometry
    fitting_button['command'] = fitting
    exit_hops_button['command'] = exit_hops

    # setup window

    photo = PhotoImage(file=holomon_logo)
    logo_label = Label(root, image=photo)
    window_label = Label(root, text='Fitting')
    created_by_label = Label(root, text=read_main_log('windows', 'created_by').replace(',', '\n'))

    setup_window(root, [
        [],
        [[logo_label, 0, 1, 6], [window_label, 1, 4, 1, 'title']],
        [],
        [[light_curve_file_label, 1, 1, 2], [light_curve_dir_label, 2],
         [planet_label, 3, 1, 2], [planet_search_entry, 4]],
        [[light_curve_file_entry, 2], [planet_entry, 4]],
        [],
        [],
        [[created_by_label, 0, 1, 3],
         [binning_label, 1], [binning_entry, 2], [iterations_label, 3], [iterations_entry, 4]],
        [[scatter_label, 1], [scatter_entry, 2], [burn_label, 3], [burn_entry, 4]],
        [],
        [[phot_filter_label, 1], [phot_filter_entry, 2], [period_label, 3], [period_entry, 4]],
        [[metallicity_label, 1], [metallicity_entry, 2], [mid_time_label, 3], [mid_time_entry, 4]],
        [[temperature_label, 1], [temperature_entry, 2], [rp_over_rs_label, 3], [rp_over_rs_entry, 4]],
        [[logg_label, 1], [logg_entry, 2], [sma_over_rs_label, 3], [sma_over_rs_entry, 4]],
        [[inclination_label, 3], [inclination_entry, 4]],
        [[eccentricity_label, 3], [eccentricity_entry, 4]],
        [[periastron_label, 3], [periastron_entry, 4]],
        [],
        [[fitting_button, 1, 4]],
        [],
        [[return_to_photometry_button, 1, 4]],
        [],
        [[exit_hops_button, 1, 4]],
        [],
    ])

    # finalise and show  window

    finalise_window(root)
    root.mainloop()

    return go_back.get()


if __name__ == '__main__':
    reduction_alignment_window()
    run = 1
    while run == 1:
        photometry_window()
        run = fitting_window()