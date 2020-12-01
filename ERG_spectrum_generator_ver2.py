""" This code aims at creating a simple spectrum plotting GUI """
# RISH, Kyoto University
# Ian 2020.11

# ----------------------------------------------
import os
import numpy as np
import matplotlib.pyplot as plt
import spacepy.pycdf as cdf
import tkinter as tk
from tkinter import Tk
from tkinter import ttk
from tkinter import Menu
from tkinter.filedialog import askopenfilename
from tkinter import messagebox as msg
import datetime
# -----------------------------------------------

# -----------------------------------------------
FILE_DATA = {}  # cdf file dictionary
FILE_NAME = ''
POSITION = []   # positon of time interval
PARAMETER = 1022    # for each time break, the lasting period of waveform data
TIME_ERROR = np.timedelta64(9, 'h')  # time deviation between UTC and JTC
METHOD_TMP = 0   # method switch
YLOG_tmp = 0    # ylog switch
MGNT = np.zeros(0)  # MGF wave magnitude
MGF_EPOCH = np.zeros(0)  # MGF epoch
Q = 1.60217662e-19  # charge of electron
M = 9.10938356e-31  # mass of electron
# -----------------------------------------------

# main window
win = Tk()  # create the main window
win.title('ERG Spectrogram Generator')  # name

# tab window
tab_control = ttk.Notebook(win)  # tab controller <-- 2 tabs are needed
tab1 = ttk.Frame(tab_control)   # tab announcement
tab_control.add(tab1, text='Spectrum Plot')   # tab initialization <-- for single wave field
tab2 = ttk.Frame(tab_control)   # tab announcement
tab_control.add(tab2, text='Others')   # tab initialization <-- Bx, By, Bz, Ex, Ey
tab_control.pack(expand=1, fill="both")    # pack to make visible

# tab1 frame
data_select = ttk.LabelFrame(tab1, text=' Data Select ')    # data select frame
data_select.grid(column=0, row=0, padx=8, pady=8, rowspan=2, sticky='WENS')
FFT_setting = ttk.LabelFrame(tab1, text=' FFT Setting ')    # FFT setting frame
FFT_setting.grid(column=1, row=0, padx=8, pady=8, sticky='WENS')
plot_setting = ttk.LabelFrame(tab1, text=' Plot Setting ')    # plot setting frame
plot_setting.grid(column=1, row=1, padx=8, pady=8, sticky='WENS')


# Select WFC file
def _select_wfc():
    global FILE_NAME, FILE_DATA, MGNT, MGF_EPOCH
    file_select = Tk()
    file_select.withdraw()  # hide window
    filename = askopenfilename(initialdir=os.getcwd(),
                               title=' Select WFC file', filetypes=[('CDF Files', '*.cdf')])    # open file
    f_opening_label['text'] = 'opening :    '
    f_label['text'] = os.path.split(filename)[-1]
    FILE_NAME = os.path.split(filename)[-1][:-4]
    file_data = cdf.CDF(filename)   # use pycdf to read cdf file
    FILE_DATA = file_data
    waveform_combo['value'], time_combo['value'], time_offset_combo['value'] = \
        list(file_data), list(file_data), list(file_data)  # select data
    # renew widgets
    split_button.configure(state='!disabled')  # split button
    split_segments['value'] = 0  # split segment
    MGNT, MGF_EPOCH = np.zeros(0), np.zeros(0)  # MGF file <- a new mgf file needs to be chosen


# Select MGF file <- optimization shall be made
def _select_mgf():
    global MGNT, MGF_EPOCH
    mgf_file = Tk()
    mgf_file.withdraw()  # hide window
    mgf_path = askopenfilename(initialdir=os.getcwd(),
                               title=' Select MGF file', filetypes=[('CDF Files', '*.cdf')])    # open file
    mgf_name = os.path.split(mgf_path)[-1][:-4]
    mgf_data = cdf.CDF(mgf_path)
    MGNT = mgf_data['magt_8sec'][:]  # choose magnitude file
    MGF_EPOCH =mgf_data['epoch_8sec'][:]  # choose epoch file


# data split
def _split():
    split_button.configure(state='disabled')  # disable split button
    global POSITION
    time_tmp = time_combo.get()  # epoch name in WFC data
    time = FILE_DATA[time_tmp][:]
    epoch_stamp = [time[i].timestamp() for i in range(time.shape[0])]   # convert epoch to timestamp
    epoch_interval = np.where(np.diff(epoch_stamp) > 10)    # find interval of time
    position = np.insert(epoch_interval[0], 0, -1)   # insert initial time
    position = np.insert(position, position.shape[0], time.shape[0])    # insert end time
    number_segment = position.shape[0]
    split_segments['value'] = list(range(1, number_segment))  # assign segment number to split combobox
    POSITION = position


# plot spectrogram
def _spectrum():
    wave_tmp, split_tmp, time_offset_tmp = waveform_combo.get(), split_segments.get(), time_offset_combo.get()
    split_tmp = int(split_tmp)
    time_offset = FILE_DATA[time_offset_tmp][:]  # read time offsets
    start = (POSITION[split_tmp - 1] + 1) * time_offset.shape[0] + PARAMETER  # calculate start point of one segment
    end = (POSITION[split_tmp]) * time_offset.shape[0]  # end point of one segment
    # wave
    wave = FILE_DATA[wave_tmp][:]  # obtain waveform data
    shp = wave.shape
    wave_array = wave.reshape(shp[0] * shp[1], )  # reshape field data to one dimensional array
    wave_segment = wave_array[start:end]  # chosen wave segment
    # time
    time_tmp = time_combo.get()
    time = FILE_DATA[time_tmp][:]   # read epoch data
    position = POSITION
    time_segment_start = position[split_tmp - 1] + 1  # start point of epoch segment
    time_segment_end = position[split_tmp]  # end point of epoch segment
    loc1, loc2 = 0, 0  # find position in MGF epoch <- for reading mgf magnitude data
    for loc1 in range(MGF_EPOCH.shape[0]):
        if time[time_segment_start] <= MGF_EPOCH[loc1]: break
    for loc2 in range(MGF_EPOCH.shape[0]):
        if time[time_segment_end-1] <= MGF_EPOCH[loc2]: break
    mgf_time = MGF_EPOCH[loc1-1:loc2+1]  # mgf time segment
    mgf_mgnt = MGNT[loc1-1:loc2+1]  # mgf magnitude segment
    # make sure mgf file is not missed
    if mgf_mgnt.shape[0] < 1:
        msg.showwarning('Warning', 'Please select MGF data!')
        return
    # calculate fce
    fce = Q / (M * 2 * np.pi) * mgf_mgnt * 1e-9
    fce_half = fce / 2
    # time segment
    time_segment_datetime = time[time_segment_start:time_segment_end]
    time_segment_stamp = [time_segment_datetime[i].timestamp() + time_offset[j] / 1000  # insert time offsets into epoch
                          for i in range(time_segment_datetime.shape[0]) for j in range(time_offset.shape[0])]
    time_segment = time_segment_stamp[PARAMETER:]  # corresponding time segment <- a parameter value is important here
    # FFT
    fs, nfft, noverlap = fs_entry.get(), nfft_entry.get(), noverlap_entry.get()  # obtain FFT setting
    fs, nfft, noverlap = int(fs), int(nfft), int(noverlap)
    # figure
    ylim_d, ylim_u = ylim_entry_1.get(), ylim_entry_2.get()  # obtain plot setting
    clim_d, clim_u = clim_entry_1.get(), clim_entry_2.get()
    fig_w, fig_h = figsize_entry_1.get(), figsize_entry_2.get()
    ylim_u, ylim_d, clim_u, clim_d, fig_w, fig_h = int(ylim_u), int(ylim_d), \
                                                   int(clim_u), int(clim_d), int(fig_w), int(fig_h)
    # plot
    method = method_button['text']  # method switch
    ylog = ylog_button['text']  # ylog switch
    plt.figure(figsize=(fig_w, fig_h))  # new figure
    cmap = cmap_combo.get()
    cmap_index = {'autumn': plt.autumn,
                  'bone': plt.bone,
                  'copper': plt.copper,
                  'cool': plt.cool,
                  'flag': plt.flag,
                  'gray': plt.gray,
                  'hot': plt.hot,
                  'hsv': plt.hsv,
                  'inferno': plt.inferno,
                  'jet': plt.jet,
                  'magma': plt.magma,
                  'nipy_spectral': plt.nipy_spectral,
                  'pink': plt.pink,
                  'plasma': plt.plasma,
                  'prism': plt.prism,
                  'spring': plt.spring,
                  'summer': plt.summer,
                  'virdis': plt.viridis,
                  'winter': plt.winter}  # colormap type selection dictionary
    cmap_index[cmap]()  # choose one colormap type
    spec_data, spec_freq, spec_time, spec_img \
        = plt.specgram(wave_segment, NFFT=nfft, Fs=fs, noverlap=noverlap, scale='dB')  # specgram method
    ax_x, ax_y = plt.gca().get_position().x0, plt.gca().get_position().y0  # left_bottom corner position of current axes
    ratio = spec_time[-1] / (mgf_mgnt.shape[0] - 1)  # for drawing fce line in specgram method
    mgf_time_tick = np.arange(0.1, spec_time[-1], ratio)
    mgf_time_tick = np.insert(mgf_time_tick, mgf_time_tick.shape[0], spec_time[-1])
    if method == 'specgram':  # specgram method
        time_tick = []
        for time in range(spec_time.shape[0]):
            # find real time corresponding to FFT position
            time_trans = datetime.datetime.fromtimestamp(time_segment[int((time + 0.5) * (nfft - noverlap))])
            time_tick.append(time_trans.strftime('%H:%M:%S')+'\n' +
                             time_trans.strftime('%f')+'\n'+time_trans.strftime('%Y.%m.%d'))    # create time ticks
        plt.plot(mgf_time_tick, fce, color="r", linestyle="--", label="electron cyclotron freq.")
        plt.plot(mgf_time_tick, fce_half, color="r", linestyle="--", label="half electron cyclotron freq.")
        plt.xticks(list(spec_time), time_tick)    # set time ticks
        plt.locator_params(axis='x', nbins=10)    # arrange time ticks
        plt.figtext(ax_x - 0.05, ax_y - 0.05, 'Time:\nms:\nDate:')  # label
        if ylog == 'log': plt.yscale('log')  # ylog form
    else:  # pcolormesh method
        plt.clf()  # clear current specgram figure
        time_tick = []
        for time in range(spec_time.shape[0]):  # find real time corresponding to FFT position
            time_tick.append(time_segment[int((time + 0.5) * (nfft - noverlap))])
        time_tick = np.array(time_tick)
        time_tick = time_tick * 1e9     # increase precision
        time_tick_final = time_tick.astype('datetime64[ns]') + TIME_ERROR  # convert float to numpy.datetime64 type
        spec_data_db = 10 * np.log10(spec_data)  # convert data scale to dB form, consistent with specgram method
        plt.pcolormesh(time_tick_final, spec_freq, spec_data_db)  # pcolormesh
        plt.plot(mgf_time, fce, color="r", linestyle="--", label="electron cyclotron freq.")
        plt.plot(mgf_time, fce_half, color="r", linestyle="--", label="half electron cyclotron freq.")
        plt.figtext(ax_x - 0.05, ax_y - 0.02, 'Time:')  # label
        if ylog == 'log': plt.yscale('log')
    cbar = plt.colorbar()
    cbar.set_label('dB', rotation=270)
    plt.clim(clim_d, clim_u)  # colorbar range
    plt.ylim(ylim_d, ylim_u)  # frequency range
    plt.title(FILE_NAME.upper() + '  ' + wave_tmp[0:2] + '  Spectrum')
    plt.ylabel('Frequency (Hz)')
    plt.xlabel('Time (UTC)')
    plt.show()


# plural spectrum
def _plural():
    split_tmp, time_offset_tmp, time_tmp = split_segments.get(), time_offset_combo.get(), time_combo.get()
    split_tmp = int(split_tmp)
    time_offset, time = FILE_DATA[time_offset_tmp][:], FILE_DATA[time_tmp][:]
    start = (POSITION[split_tmp - 1] + 1) * time_offset.shape[0] + PARAMETER
    end = (POSITION[split_tmp]) * time_offset.shape[0]
    position = POSITION
    time_segment_start = position[split_tmp - 1] + 1
    time_segment_end = position[split_tmp]
    time_segment_datetime = time[time_segment_start:time_segment_end]
    time_segment_stamp = [time_segment_datetime[i].timestamp() + time_offset[j] / 1000
                          for i in range(time_segment_datetime.shape[0]) for j in range(time_offset.shape[0])]
    time_segment = time_segment_stamp[PARAMETER:]  # corresponding time segment
    data_list = list(FILE_DATA)
    wave_data_indices = \
        [i for i, s in enumerate(data_list) if '_waveform' in s]  # find all waveforms in current WFC data
    wave_name = [data_list[j] for j in wave_data_indices]  # list containing waveform names
    loc1, loc2 = 0, 0
    for loc1 in range(MGF_EPOCH.shape[0]):
        if time[time_segment_start] <= MGF_EPOCH[loc1]: break
    for loc2 in range(MGF_EPOCH.shape[0]):
        if time[time_segment_end-1] <= MGF_EPOCH[loc2]: break
    mgf_time = MGF_EPOCH[loc1 - 1:loc2 + 1]
    mgf_mgnt = MGNT[loc1 - 1:loc2 + 1]
    fce = Q / (M * 2 * np.pi) * mgf_mgnt * 1e-9
    fce_half = fce / 2
    # FFT
    fs, nfft, noverlap = fs_entry.get(), nfft_entry.get(), noverlap_entry.get()
    fs, nfft, noverlap = int(fs), int(nfft), int(noverlap)
    # figure
    ylim_d, ylim_u = ylim_entry_1.get(), ylim_entry_2.get()
    clim_d, clim_u = clim_entry_1.get(), clim_entry_2.get()
    fig_w, fig_h = figsize_entry_1.get(), figsize_entry_2.get()
    ylim_u, ylim_d, clim_u, clim_d, fig_w, fig_h = int(ylim_u), int(ylim_d), \
                                                   int(clim_u), int(clim_d), int(fig_w), int(fig_h)
    # plot
    method = method_button['text']
    ylog = ylog_button['text']
    plt.figure(figsize=(fig_w, fig_h))
    cmap = cmap_combo.get()
    cmap_index = {'autumn': plt.autumn,
                  'bone': plt.bone,
                  'copper': plt.copper,
                  'cool': plt.cool,
                  'flag': plt.flag,
                  'gray': plt.gray,
                  'hot': plt.hot,
                  'hsv': plt.hsv,
                  'inferno': plt.inferno,
                  'jet': plt.jet,
                  'magma': plt.magma,
                  'nipy_spectral': plt.nipy_spectral,
                  'pink': plt.pink,
                  'plasma': plt.plasma,
                  'prism': plt.prism,
                  'spring': plt.spring,
                  'summer': plt.summer,
                  'virdis': plt.viridis,
                  'winter': plt.winter}
    cmap_index[cmap]()
    num = len(wave_name)  # number of waveform types
    spec_time_method1 = np.zeros(0)
    time_tick_method1 = []
    for index in range(num):    # make subplots
        time_tick = []
        wave_component = FILE_DATA[wave_name[index]][:]
        shp = wave_component.shape
        wave_component_array = wave_component.reshape(shp[0] * shp[1], )
        wave_component_segment = wave_component_array[start:end]
        plt.subplot(num, 1, index+1)
        spec_data_p, spec_freq_p, spec_time_p, spec_img_p = \
            plt.specgram(wave_component_segment, NFFT=nfft, Fs=fs, noverlap=noverlap, scale='dB')
        if method == 'specgram':
            ratio = spec_time_p[-1] / (mgf_mgnt.shape[0] - 1)
            mgf_time_tick = np.arange(0.1, spec_time_p[-1], ratio)
            mgf_time_tick = np.insert(mgf_time_tick, mgf_time_tick.shape[0], spec_time_p[-1])
            plt.plot(mgf_time_tick, fce, color="r", linestyle="--", label="electron cyclotron freq.")
            plt.plot(mgf_time_tick, fce_half, color="r", linestyle="--", label="half electron cyclotron freq.")
            if index == num-1:
                for time in range(spec_time_p.shape[0]):
                    time_trans = datetime.datetime.fromtimestamp(
                        time_segment[int((time + 0.5) * (nfft - noverlap))])
                    time_tick.append(time_trans.strftime('%H:%M:%S') + '\n' +
                                     time_trans.strftime('%f') + '\n' + time_trans.strftime('%Y.%m.%d'))
                spec_time_method1 = spec_time_p
                time_tick_method1 = time_tick
            if ylog == 'log': plt.yscale('log')
        else:
            for time in range(spec_time_p.shape[0]):
                time_tick.append(time_segment[int((time + 0.5) * (nfft - noverlap))])
            plt.cla()
            time_tick = np.array(time_tick)
            time_tick = time_tick * 1e9
            time_tick_final = time_tick.astype('datetime64[ns]') + TIME_ERROR
            spec_data_db = 10 * np.log10(spec_data_p)
            plt.pcolormesh(time_tick_final, spec_freq_p, spec_data_db)
            plt.plot(mgf_time, fce, color="r", linestyle="--", label="electron cyclotron freq.")
            plt.plot(mgf_time, fce_half, color="r", linestyle="--", label="half electron cyclotron freq.")
            if ylog == 'log': plt.yscale('log')
        cbar = plt.colorbar()
        cbar.set_label('dB', rotation=270)
        plt.clim(clim_d, clim_u)
        plt.title(FILE_NAME.upper() + '  ' + wave_name[index][0:2].upper() + '  Spectrum')
        plt.ylim(ylim_d, ylim_u)
        plt.ylabel('Frequency (Hz)')
        if index < num-1:
            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax_x, ax_y = plt.gca().get_position().x0, plt.gca().get_position().y0
    if method == 'specgram':
        plt.xticks(list(spec_time_method1), time_tick_method1)
        plt.locator_params(axis='x', nbins=10)
        plt.figtext(ax_x - 0.05, ax_y - 0.05, 'Time:\nms:\nDate:')
    else:
        plt.figtext(ax_x - 0.05, ax_y - 0.02, 'Time:')
    plt.xlabel('Time (UTC)')
    plt.show()


# Exit GUI cleanly
def _quit():
    win.quit()
    win.destroy()
    exit()


# About massage box
def _info():
    msg.showinfo('v1.0.0', 'Author: Ian, RISH, Kyoto University\n 2020.11')


# Method switch
def _method():
    global METHOD_TMP
    METHOD_TMP += 1
    if METHOD_TMP % 2:
        method_button.configure(text='specgram')
    else:
        method_button.configure(text='pcolormesh')


# Ylog switch
def _method1():
    global YLOG_tmp
    YLOG_tmp += 1
    if YLOG_tmp % 2:
        ylog_button.configure(text='linear')
    else:
        ylog_button.configure(text='log')


# default parameters setting
def _default():
    time_combo.set('epoch')
    time_offset_combo.set('time_offsets')
    fs_entry.delete(0, tk.END)
    fs_entry.insert(0, '65536')
    nfft_entry.delete(0, tk.END)
    nfft_entry.insert(0, '8192')
    noverlap_entry.delete(0, tk.END)
    noverlap_entry.insert(0, '4096')
    ylim_entry_1.delete(0, tk.END)
    ylim_entry_1.insert(0, '1')
    ylim_entry_2.delete(0, tk.END)
    ylim_entry_2.insert(0, '3000')
    clim_entry_1.delete(0, tk.END)
    clim_entry_1.insert(0, '-10')
    clim_entry_2.delete(0, tk.END)
    clim_entry_2.insert(0, '10')
    figsize_entry_1.delete(0, tk.END)
    figsize_entry_1.insert(0, '20')
    figsize_entry_2.delete(0, tk.END)
    figsize_entry_2.insert(0, '7')
    cmap_combo.set('jet')


# Creating a menu bar
menu_bar = Menu(win)
win.config(menu=menu_bar)

# Adding menu items
file_menu = Menu(menu_bar, tearoff=0)   # file menu
file_menu.add_command(label='New...', command=_select_wfc)  # open file
file_menu.add_command(label='Mgf Data', command=_select_mgf)
file_menu.add_command(label='Default', command=_default)
file_menu.add_command(label="Exit", command=_quit)  # quit
menu_bar.add_cascade(label="File", menu=file_menu)  # connect file menu to menu bar
help_menu = Menu(menu_bar, tearoff=0)   # help menu
help_menu.add_command(label='About', command=_info)  # information
menu_bar.add_cascade(label="Help", menu=help_menu)

# data select frame
# file name
f_opening_label = ttk.Label(data_select, text='Opening ...', width=30)
f_opening_label.grid(column=0, row=0, sticky='W')
f_label = ttk.Label(data_select, text='', width=30)
f_label.grid(column=0, row=1, sticky='W')

# select waveform data
waveform_label = ttk.Label(data_select, text='Waveform data:').grid(column=0, row=2, columnspan=2, sticky='W')
content1 = tk.StringVar()
waveform_combo = ttk.Combobox(data_select, width=25, textvariable=content1, state='readonly')
waveform_combo.grid(column=0, row=3, columnspan=2, sticky='WENS')

# select time data
epoch_label = ttk.Label(data_select, text='Epoch data:').grid(column=0, row=4, sticky='W')
content2 = tk.StringVar()
time_combo = ttk.Combobox(data_select, width=25, textvariable=content2, state='readonly')
time_combo.grid(column=0, row=5, columnspan=2, sticky='WENS')

# Select time offset
time_offset_label = ttk.Label(data_select, text='Time offset:').grid(column=0, row=6, sticky='W')
content3 = tk.StringVar()
time_offset_combo = ttk.Combobox(data_select, width=25, textvariable=content3, state='readonly')
time_offset_combo.grid(column=0, row=7, columnspan=2, sticky='WENS')

# Segment splitter
split_button = ttk.Button(data_select, text='Split Data:', command=_split)
split_button.grid(column=0, row=8, pady=10, sticky='WENS')
segment = tk.IntVar()
split_segments = ttk.Combobox(data_select, width=10, textvariable=segment, state='readonly')
split_segments.grid(column=1, row=8)

# Method select
method_button = ttk.Button(data_select, text='Specgram', command=_method)
method_button.grid(column=1, row=0, rowspan=2, sticky='WENS')


# FFT setting frame
# fs
fs_label = ttk.Label(FFT_setting, text='Fs: ').grid(column=0, row=0, sticky='WE')
fs = tk.DoubleVar()
fs_entry = ttk.Entry(FFT_setting, textvariable=fs)
fs_entry.grid(column=1, row=0, sticky='WENS')

# nfft
nfft_label = ttk.Label(FFT_setting, text='NFFT: ').grid(column=0, row=1, sticky='WE')
nfft = tk.DoubleVar()
nfft_entry = ttk.Entry(FFT_setting, textvariable=nfft)
nfft_entry.grid(column=1, row=1, sticky='WENS')

# overlap
noverlap_label = ttk.Label(FFT_setting, text='noverlap: ').grid(column=0, row=2, sticky='WE')
noverlap = tk.DoubleVar()
noverlap_entry = ttk.Entry(FFT_setting,textvariable=noverlap)
noverlap_entry.grid(column=1, row=2, sticky='WENS')

# window form
window_label = ttk.Label(FFT_setting, text='window: ').grid(column=0, row=3, sticky='WE')
window_form = tk.StringVar()
window_form.set('Hanning')
window_combo = ttk.Combobox(FFT_setting, textvariable=window_form, state='readonly')
window_combo.grid(column=1, row=3, sticky='WENS')


# plot setting frame
# ylim
ylim_label = ttk.Label(plot_setting, text='Y range: ').grid(column=0, row=0, sticky='W')
ylog_button = ttk.Button(plot_setting, text='linear', command=_method1)
ylog_button.grid(column=1, row=0)
ylim_1 = tk.IntVar()
ylim_entry_1 = ttk.Entry(plot_setting, width=6, textvariable=ylim_1)
ylim_entry_1.grid(column=2, row=0)
ylim_mid = ttk.Label(plot_setting, text=' - ').grid(column=3, row=0)
ylim_2 = tk.IntVar()
ylim_entry_2 = ttk.Entry(plot_setting, width=6, textvariable=ylim_2)
ylim_entry_2.grid(column=4, row=0)

# clim
clim_label = ttk.Label(plot_setting, text='Colorbar range: ').grid(column=0, row=1, sticky='W', columnspan=2)
clim_1 = tk.IntVar()
clim_entry_1 = ttk.Entry(plot_setting, width=6, textvariable=clim_1)
clim_entry_1.grid(column=2, row=1)
clim_mid = ttk.Label(plot_setting, text=' - ').grid(column=3, row=1)
clim_2 = tk.IntVar()
clim_entry_2 = ttk.Entry(plot_setting, width=6, textvariable=clim_2)
clim_entry_2.grid(column=4, row=1)

# figsize
figsize_label = ttk.Label(plot_setting, text='Figure size: ').grid(column=0, row=2, sticky='W', columnspan=2)
figsize_1 = tk.IntVar()
figsize_entry_1 = ttk.Entry(plot_setting, width=6, textvariable=figsize_1)
figsize_entry_1.grid(column=2, row=2)
figsize_mid = ttk.Label(plot_setting, text=' x ').grid(column=3, row=2)
figsize_2 = tk.IntVar()
figsize_entry_2 = ttk.Entry(plot_setting, width=6, textvariable=figsize_2)
figsize_entry_2.grid(column=4, row=2)

# colormap type
cmap_label = ttk.Label(plot_setting, text='Colormap type: ').grid(column=0, row=3, sticky='W', columnspan=2)
cmap_type = tk.StringVar()
cmap_combo = ttk.Combobox(plot_setting, textvariable=cmap_type, state='readonly', width=13)
cmap_combo['value'] = ('autumn', 'bone', 'copper', 'cool', 'flag', 'gray',
                       'hot', 'hsv', 'inferno', 'jet', 'magma', 'nipy_spectral',
                       'pink', 'plasma', 'prism', 'spring', 'summer', 'virdis', 'winter')
cmap_combo.grid(column=2, row=3, columnspan=3)
cmap_combo.current(0)


# plot frame
plot_button = ttk.Button(tab1, text='Plural Plot', command=_plural)
plot_button.grid(column=2, row=0, ipadx=10, ipady=10, sticky='WENS')
plot_button = ttk.Button(tab1, text='Single Plot', command=_spectrum)
plot_button.grid(column=2, row=1, ipadx=10, ipady=10, sticky='WENS')

# erg icon
win.iconbitmap('erg.ico')

win.mainloop()  # starts the window
