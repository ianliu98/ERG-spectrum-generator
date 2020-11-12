""" This code aims at creating a simple spectrum plotting GUI """

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
from matplotlib.colors import LogNorm
# -----------------------------------------------

# -----------------------------------------------
FILE_DATA = {}  # cdf file dictionary
FILE_NAME = ''
TIME_STAMP = []  # list to store float timestamp from datetime.datetime type
POSITION = []   # positon of time interval
PARAMETER = 1022    # for each time break, the lasting period of waveform data
TIME_ERROR = np.timedelta64(9, 'h')  # time deviation between UTC and JTC
METHOD_TMP = 0   # method selecting use
Q = 1.60217662e-19  # charge of electron
M = 9.10938356e-31  # mass of electron
# -----------------------------------------------

win = Tk()  # create the main window
win.title('ERG Spectrogram Generator')  # name

tab_control = ttk.Notebook(win)  # tab controller <-- 2 tabs are needed

tab1 = ttk.Frame(tab_control)   # tab announcement
tab_control.add(tab1, text='Spectrum Plot')   # tab initialization <-- for single wave field
tab2 = ttk.Frame(tab_control)   # tab announcement
tab_control.add(tab2, text='Others')   # tab initialization <-- Bx, By, Bz, Ex, Ey

tab_control.pack(expand=1, fill="both")    # pack to make visible

data_select = ttk.LabelFrame(tab1, text=' Data Select ')    # data select frame
data_select.grid(column=0, row=0, padx=8, pady=8, rowspan=2, sticky='WENS')
FFT_setting = ttk.LabelFrame(tab1, text=' FFT Setting ')    # FFT setting frame
FFT_setting.grid(column=1, row=0, padx=8, pady=8, sticky='WENS')
plot_setting = ttk.LabelFrame(tab1, text=' Plot Setting ')    # plot setting frame
plot_setting.grid(column=1, row=1, padx=8, pady=8, sticky='WENS')


# Selecting file
def _select():
    global FILE_NAME
    file_select = Tk()
    file_select.withdraw()  # hide window
    filename = askopenfilename(initialdir=os.getcwd(),
                               title=' Select file', filetypes=[('CDF Files', '*.cdf')])    # open file
    f_label['text'] = 'opening :    ' + os.path.split(filename)[-1]
    FILE_NAME = os.path.split(filename)[-1][:-4]
    file_data = cdf.CDF(filename)
    global FILE_DATA
    FILE_DATA = file_data
    waveform_combo['value'], time_combo['value'], time_offset_combo['value'] = \
        list(file_data), list(file_data), list(file_data)  # select data
    split_button.configure(state='!disabled')  # make split button able after open file every time


# data split
def _split():
    split_button.configure(state='disabled')  # make split button disabled after split
    global TIME_STAMP, POSITION
    time_tmp, time_offset_tmp = time_combo.get(), time_offset_combo.get()  # epoch and time_offsets name in WFC data
    time = FILE_DATA[time_tmp][:]
    time_offset = FILE_DATA[time_offset_tmp][:]
    time_stamp = (time[i].timestamp() + time_offset[j] / 1000   # group epoch and time_offsets
                  for i in range(time.shape[0]) for j in range(time_offset.shape[0]))  # convert datetime.datetime type to timestamp
    time_stamp_list = list(time_stamp)
    del time, time_offset
    position = np.where(np.diff(time_stamp_list) > 10)  # find time intervals
    position = position[0]
    position = np.insert(position, 0, 0)
    position = np.insert(position, position.shape[0], time_stamp_list[-1])
    number_segment = position.shape[0]
    split_segments['value'] = list(range(1, number_segment)) # assign segment number to split combobox
    TIME_STAMP, POSITION = time_stamp_list, position


# plot spectrogram
def _spectrum():
    global TIME_STAMP, POSITION, PARAMETER
    # wave
    wave_tmp, split_tmp = waveform_combo.get(), split_segments.get()  
    wave = FILE_DATA[wave_tmp][:]  # obtain waveform data
    shp = wave.shape
    wave_array = wave.reshape(shp[0] * shp[1], )  # reshape field data to one dimensional array
    split_tmp = int(split_tmp)
    wave_segment = wave_array[POSITION[split_tmp-1]+PARAMETER:POSITION[split_tmp]]  # chosed wave segment
    time_segment = TIME_STAMP[POSITION[split_tmp-1]+PARAMETER:POSITION[split_tmp]]  # corresponding time segment
    fce = Q/(M * 2 * np.pi) * wave_segment * 1e-12
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
    method = method_button['text']  # obtain value of method button
    plt.figure(figsize=(fig_w, fig_h))
    cmap = cmap_combo.get()
    cmap_index = {'autumn': plt.autumn,  # colormap type selection dictionary
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
    cmap_index[cmap]()  # choose one cmap type 
    spec_data, spec_freq, spec_time, spec_img = plt.specgram(wave_segment, NFFT=nfft, Fs=fs, noverlap=noverlap)  # specgram
    ax_x, ax_y = plt.gca().get_position().x0, plt.gca().get_position().y0  # left_bottom cornor position of current axes
    fce_fft = []
    if method == 'specgram':  # specgram method
        time_tick = []
        for time in range(spec_time.shape[0]):
            time_trans = datetime.datetime.fromtimestamp(time_segment[int((time + 0.5) * (nfft - noverlap))])   # find real time corresponding to FFT position
            fce_fft.append(fce[int((time + 0.5) * (nfft - noverlap))])
            time_tick.append(time_trans.strftime('%H:%M:%S')+'\n' +
                             time_trans.strftime('%f')+'\n'+time_trans.strftime('%Y.%m.%d'))    # create time ticks
        fce_fft = np.array(fce_fft)
        # print(fce_fft)
        plt.plot(spec_time, fce_fft, color="k", linestyle="--", label="electron cyclotron freq.")
        plt.plot(spec_time, fce_fft/2, color="k", linestyle="--", label="half electron cyclotron freq.")
        plt.xticks(list(spec_time), time_tick)    # set time ticks
        plt.locator_params(axis='x', nbins=10)
        plt.figtext(ax_x - 0.05, ax_y - 0.05, 'Time:\nms:\nDate:')
        plt.colorbar()
        plt.clim(clim_d, clim_u)
    else:  # pcolormesh method
        plt.clf()
        time_tick = []
        fce_fft = []
        for time in range(spec_time.shape[0]):
            time_tick.append(time_segment[int((time + 0.5) * (nfft - noverlap))])   # find real time corresponding to FFT position
            fce_fft.append(fce[int((time + 0.5) * (nfft - noverlap))])
        fce_fft = np.array(fce_fft)
        time_tick = np.array(time_tick)
        time_tick = time_tick * 1e9     # increase precision
        time_tick_final = time_tick.astype('datetime64[ns]') + TIME_ERROR  # convert float to numpy.datetime64 type
        plt.pcolormesh(time_tick_final, spec_freq, spec_data, norm=LogNorm())  # pcolormesh
        # plt.plot(spec_time, fce_fft, color="k", linestyle="--", label="electron cyclotron freq.")
        # plt.plot(spec_time, fce_fft/2, color="k", linestyle="--", label="half electron cyclotron freq.")
        plt.figtext(ax_x - 0.05, ax_y - 0.02, 'Time:')
        plt.colorbar()
        plt.clim(10 ** clim_d, 10 ** clim_u)
    plt.ylim(ylim_d, ylim_u)
    plt.title(FILE_NAME.upper() + '  ' + wave_tmp[0:2] + '  Spectrum')
    plt.ylabel('Frequency (Hz)')
    plt.xlabel('Time (UTC)')
    plt.show()


# plural spectrum
def _plural():
    global FILE_DATA, POSITION, PARAMETER, TIME_STAMP
    split_tmp = split_segments.get()
    split_tmp = int(split_tmp)
    start, end = POSITION[split_tmp-1] + PARAMETER, POSITION[split_tmp]
    data_list = list(FILE_DATA)
    wave_data_indices = [i for i, s in enumerate(data_list) if '_waveform' in s]  # find all waveforms in current WFC data
    wave_name = [data_list[j] for j in wave_data_indices]
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
    num = len(wave_name)
    spec_time_method1 = np.zeros(0)
    time_tick_method1 = []
    for index in range(num):    # make subplots
        time_tick = []
        wave_component = FILE_DATA[wave_name[index]][:]
        shp = wave_component.shape
        wave_component_array = wave_component.reshape(shp[0] * shp[1], )
        wave_component_segment = wave_component_array[start:end]
        # fce = Q / (M * 2 * np.pi) * wave_component_segment * 1e-12
        # fce_half = fce / 2
        time_component_segment = TIME_STAMP[start:end]
        plt.subplot(num, 1, index+1)
        spec_data_p, spec_freq_p, spec_time_p, spec_img_p = \
            plt.specgram(wave_component_segment, NFFT=nfft, Fs=fs, noverlap=noverlap)
        if method == 'specgram':
            if index == num-1:
                for time in range(spec_time_p.shape[0]):
                    time_trans = datetime.datetime.fromtimestamp(
                        time_component_segment[int((time + 0.5) * (nfft - noverlap))])
                    time_tick.append(time_trans.strftime('%H:%M:%S') + '\n' +
                                     time_trans.strftime('%f') + '\n' + time_trans.strftime('%Y.%m.%d'))
                spec_time_method1 = spec_time_p
                time_tick_method1 = time_tick
            plt.colorbar()
            plt.clim(clim_d, clim_u)
        else:
            for time in range(spec_time_p.shape[0]):
                time_tick.append(time_component_segment[int((time + 0.5) * (nfft - noverlap))])
            plt.cla()
            time_tick = np.array(time_tick)
            time_tick = time_tick * 1e9
            time_tick_final = time_tick.astype('datetime64[ns]') + TIME_ERROR
            plt.pcolormesh(time_tick_final, spec_freq_p, spec_data_p, norm=LogNorm())
            plt.colorbar()
            plt.clim(10 ** clim_d, 10 ** clim_u)
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
    msg.showinfo('v1.0.0', 'Author: Ian, RISH, Kyoto University\n 2020.11.09')


# Method select
def _method():
    global METHOD_TMP
    METHOD_TMP += 1
    if METHOD_TMP % 2:
        method_button.configure(text='specgram')
    else:
        method_button.configure(text='pcolormesh')


# Creating a menu bar
menu_bar = Menu(win)
win.config(menu=menu_bar)

# Adding menu items
file_menu = Menu(menu_bar, tearoff=0)   # file menu
file_menu.add_command(label='New...', command=_select)  # open file
file_menu.add_command(label="Exit", command=_quit)  # quit
menu_bar.add_cascade(label="File", menu=file_menu)  # connect file menu to menu bar
help_menu = Menu(menu_bar, tearoff=0)   # help menu
help_menu.add_command(label='About', command=_info)  # information
menu_bar.add_cascade(label="Help", menu=help_menu)


# file name
f_label = ttk.Label(data_select, text='Opening ...', width=30)
f_label.grid(column=0, row=0, sticky='W')

# select waveform data
waveform_label = ttk.Label(data_select, text='Select a waveform data:').grid(column=0, row=1, columnspan=2, sticky='W')
content1 = tk.StringVar()
waveform_combo = ttk.Combobox(data_select, width=25, textvariable=content1, state='readonly')
waveform_combo.grid(column=0, row=2, columnspan=2, sticky='WENS')

# select time data
epoch_label = ttk.Label(data_select, text='Select time data:').grid(column=0, row=3, sticky='W')
content2 = tk.StringVar()
time_combo = ttk.Combobox(data_select, width=25, textvariable=content2, state='readonly')
time_combo.grid(column=0, row=4, columnspan=2, sticky='WENS')

# Select time offset
time_offset_label = ttk.Label(data_select, text='Select time offset:').grid(column=0, row=5, sticky='W')
content3 = tk.StringVar()
time_offset_combo = ttk.Combobox(data_select, width=25, textvariable=content3, state='readonly')
time_offset_combo.grid(column=0, row=6, columnspan=2, sticky='WENS')

# Segment splitter
split_button = ttk.Button(data_select, text='Split Data:', command=_split)
split_button.grid(column=0, row=7, pady=10, sticky='WENS')
segment = tk.IntVar()
split_segments = ttk.Combobox(data_select, width=10, textvariable=segment, state='readonly')
split_segments.grid(column=1, row=7)

# Method select
method_button = ttk.Button(data_select, text='Specgram', command=_method)
method_button.grid(column=1, row=0, sticky='WENS')


# FFT setting
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


# plot setting
# ylim
ylim_label = ttk.Label(plot_setting, text='Y range: ').grid(column=0, row=0, sticky='W')
ylim_1 = tk.IntVar()
ylim_entry_1 = ttk.Entry(plot_setting, width=6, textvariable=ylim_1)
ylim_entry_1.grid(column=1, row=0)
ylim_mid = ttk.Label(plot_setting, text=' - ').grid(column=2, row=0)
ylim_2 = tk.IntVar()
ylim_entry_2 = ttk.Entry(plot_setting, width=6, textvariable=ylim_2)
ylim_entry_2.grid(column=3, row=0)
# clim
clim_label = ttk.Label(plot_setting, text='Colorbar range(10x): ').grid(column=0, row=1, sticky='W')
clim_1 = tk.IntVar()
clim_entry_1 = ttk.Entry(plot_setting, width=6, textvariable=clim_1)
clim_entry_1.grid(column=1, row=1)
clim_mid = ttk.Label(plot_setting, text=' - ').grid(column=2, row=1)
clim_2 = tk.IntVar()
clim_entry_2 = ttk.Entry(plot_setting, width=6, textvariable=clim_2)
clim_entry_2.grid(column=3, row=1)
# figsize
figsize_label = ttk.Label(plot_setting, text='Figure size: ').grid(column=0, row=2, sticky='W')
figsize_1 = tk.IntVar()
figsize_entry_1 = ttk.Entry(plot_setting, width=6, textvariable=figsize_1)
figsize_entry_1.grid(column=1, row=2)
figsize_mid = ttk.Label(plot_setting, text=' x ').grid(column=2, row=2)
figsize_2 = tk.IntVar()
figsize_entry_2 = ttk.Entry(plot_setting, width=6, textvariable=figsize_2)
figsize_entry_2.grid(column=3, row=2)
# colormap type
cmap_label = ttk.Label(plot_setting, text='Colormap type: ').grid(column=0, row=3, sticky='W')
cmap_type = tk.StringVar()
cmap_combo = ttk.Combobox(plot_setting, textvariable=cmap_type, state='readonly', width=13)
cmap_combo['value'] = ('autumn', 'bone', 'copper', 'cool', 'flag', 'gray',
                       'hot', 'hsv', 'inferno', 'jet', 'magma', 'nipy_spectral',
                       'pink', 'plasma', 'prism', 'spring', 'summer', 'virdis', 'winter')
cmap_combo.grid(column=1, row=3, columnspan=3)
cmap_combo.current(0)


# plot
plot_button = ttk.Button(tab1, text='Plural Plot', command=_plural)
plot_button.grid(column=2, row=0, ipadx=10, ipady=10, sticky='WENS')
plot_button = ttk.Button(tab1, text='Single Plot', command=_spectrum)
plot_button.grid(column=2, row=1, ipadx=10, ipady=10, sticky='WENS')

win.iconbitmap('erg.ico')

win.mainloop()  # starts the window

