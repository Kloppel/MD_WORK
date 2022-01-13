# coding=utf-8

import matplotlib.pyplot as plt
import numpy as np

def parse_rmsd_file(filepath):

    rmsds_dict = {}
    traj_list = []
    f = open(filepath,'r')

    for line in f:
        entries = line.split()
        if entries:
            if entries[0] == 'frame':
                for entry in entries:
                    if entry != 'frame':
                        rmsds_dict[entry] = []
                        traj_list.append(entry)
            else:
                for i,traj in enumerate(traj_list):
                    if entries[i+1] != 'NA':
                        rmsds_dict[traj].append(float(entries[i+1]))
                    # else:
                    #     # rmsds_dict[traj].append(None)
                    #     continue
    return rmsds_dict


if __name__ == '__main__':

    # flu_oxi = '/user/jdragelj/projects/cco_fluorescein/rmsd_data/fluO_RMSD.dat'
    # rmsd_flu_oxi = parse_rmsd_file(flu_oxi)
    # flu_red = '/user/jdragelj/projects/cco_fluorescein/rmsd_data/fluR_RMSD.dat'
    # rmsd_flu_red = parse_rmsd_file(flu_red)
    #
    #
    # # wt = '/user/jdragelj/wt_RMSD.dat'
    # # rmsd_wt = parse_rmsd_file(wt)
    # # print rmsd.keys()
    #
    # outfolder = '/user/jdragelj/Desktop/'
    #
    # time_at_frames = np.arange(0, 645, 4)
    # frame_scale = []
    # for time in time_at_frames:
    #     time = time / 20.0
    #     frame_scale.append(time)
    #
    # x_lim = [0, 33]
    # x_ticks = np.arange(0, 34, 1)
    # labelsx = np.arange(0, 34, 3)
    # string_labels_x = []
    # for i in np.arange(0, 34, 1):
    #     if i not in labelsx:
    #         string_labels_x.append('')
    #     else:
    #         string_labels_x.append(i)
    #
    # y1_lim = [0, 1.6]
    #
    # plt.rcParams.update({'font.size': 20})
    # fig, ax1 = plt.subplots()
    #
    #
    # y1 = rmsd_flu_oxi['flu110'][::4]
    # y2 = rmsd_flu_oxi['flu000'][:645][::4]
    # y3 = rmsd_flu_oxi['flu11u'][::4]
    # y4 = rmsd_flu_oxi['flu100'][::4]
    #
    # ax1.plot(frame_scale, y1, color='black', lw=1.5)
    # ax1.plot(frame_scale, y2, color='green', lw=1.5)
    # ax1.plot(frame_scale, y3, color='blue', lw=1.5)
    # ax1.plot(frame_scale, y4, color='magenta', lw=1.5)
    #
    # ax1.set_ylim(y1_lim[0], y1_lim[1])
    # ax1.set_xticklabels(string_labels_x)
    # ax1.tick_params(direction='out', top='off', right='off')
    # plt.xlim(x_lim[0], x_lim[1])
    # plt.xticks(x_ticks)
    #
    # plt.savefig(outfolder + 'RMSD_oxi_dominant.png', dpi=1200)
    #
    # fig, ax2 = plt.subplots()
    #
    # y5 = rmsd_flu_red['flu11d'][::4]
    # y6 = rmsd_flu_red['flu11u'][::4]
    # y7 = rmsd_flu_red['flu000'][:645][::4]
    #
    # ax2.plot(frame_scale, y5, color='blue', lw=1.5)
    # ax2.plot(frame_scale, y6, color='red', lw=1.5)
    # ax2.plot(frame_scale, y7, color='green', lw=1.5)
    #
    # ax2.set_ylim(y1_lim[0], y1_lim[1])
    # ax2.set_xticklabels(string_labels_x)
    # ax2.tick_params(direction='out', top='off', right='off')
    # plt.xlim(x_lim[0], x_lim[1])
    # plt.xticks(x_ticks)
    #
    # plt.savefig(outfolder + 'RMSD_red_dominant.png', dpi=1200)



    wt = '/user/jdragelj/projects/cco_fluorescein/rmsd_data/wt_RMSD.dat'
    rmsd_wt = parse_rmsd_file(wt)
    outfolder = '/user/jdragelj/Desktop/'


    y1_lim = [0, 1.4]
    y1_ticks = np.arange(0, 1.5, 0.2)
    time_at_frames = np.arange(0, 405, 4)
    frame_scale = []
    for time in time_at_frames:
        time = time / 20.0
        frame_scale.append(time)

    x_lim = [0, 21]
    x_ticks = np.arange(0, 22, 1)
    labelsx = np.arange(0, 22, 3)
    string_labels_x = []
    for i in np.arange(0, 22, 1):
        if i not in labelsx:
            string_labels_x.append('')
        else:
            string_labels_x.append(i)

    fig = plt.figure()
    ax3 = fig.add_subplot(111)

    yx1 = rmsd_wt['wt11O'][::4]
    yx2 = rmsd_wt['wt11R'][::4]
    yx3 = rmsd_wt['wt01R'][::4]
    yx4 = rmsd_wt['wt10O'][::4]

    # ax3.plot(frame_scale, yx1, color='black', lw=1.5)
    # ax3.plot(frame_scale, yx2, color='red', lw=1.5)
    # ax3.plot(frame_scale, yx3, color='cyan', lw=1.5)
    # ax3.plot(frame_scale, yx4, color='olive', lw=1.5)



    # ax3.plot(frame_scale, yx1, color='black', lw=1.0, label=label)
    ax3.plot(frame_scale, yx1, color='black', lw=1.0)
    ax3.plot(frame_scale, yx2, color='red', lw=1.0)
    ax3.plot(frame_scale, yx3, color='cyan', lw=1.0)
    ax3.plot(frame_scale, yx4, color='olive', lw=1.0)

    # ax3.plot(frame_scale, yx1, color='black')
    # ax3.plot(frame_scale, yx2, color='red')
    # ax3.plot(frame_scale, yx3, color='cyan')
    # ax3.plot(frame_scale, yx4, color='olive')

    ax3.set_ylim(y1_lim[0], y1_lim[1])
    ax3.set_xticklabels(string_labels_x)
    ax3.tick_params(direction='out', top='off', right='off')
    plt.xlim(x_lim[0], x_lim[1])
    ax3.set_yticks(y1_ticks)
    plt.xticks(x_ticks)
    # plt.rcParams["figure.figsize"] = (1, 3)
    fig.set_size_inches(4.0,2.5)
    plt.rcParams.update({'font.size': 12})
    plt.rc('font', family='Times New Roman')
    # plt.xlabel(r'time[ns]')
    # plt.xlabel(r'time[ns], $^\lambda$')
    # plt.xlabel(r'O/His73$^{+}$/His526$^{\epsilon}$/Flu$^-$')
    # plt.xlabel(r'O /His73$\mathregular{^+}$/His526$\mathregular{^\epsilon}$/Flu$\mathregular{^-}$')
    # plt.xlabel(r'O His73$\mathregular{^+}$His526$\mathregular{^\epsilon}$Flu$\mathregular{^-}$')
    # plt.ylabel(r'RMSD[$\AA$]')
    plt.tight_layout()
    # plt.legend(loc='lower right')

    plt.savefig(outfolder + 'RMSD_wt_dominant.png', dpi=300)







    # flu = '/user/jdragelj/Desktop/flu_pnas_rmsd/all_cut.dat'
    # rmsd_xtra = parse_rmsd_file(flu)
    #
    # # flu_oxi = '/user/jdragelj/fluO_RMSD.dat'
    # # rmsd_flu_oxi = parse_rmsd_file(flu_oxi)
    # #
    # # wt = '/user/jdragelj/wt_RMSD.dat'
    # # rmsd_wt = parse_rmsd_file(wt)
    # # print rmsd.keys()
    #
    # outfolder = '/user/jdragelj/Desktop/'
    #
    #
    #
    # y1_lim = [0, 1.6]
    # y1_ticks = np.arange(0, 1.7, 0.2)
    # time_at_frames = np.arange(0, 945, 2)
    # frame_scale = []
    # for time in time_at_frames:
    #     time = time / 20.0
    #     frame_scale.append(time)
    #
    # x_lim = [0, 48]
    # x_ticks = np.arange(0, 49, 1)
    # labelsx = np.arange(0, 49, 4)
    # string_labels_x = []
    # for i in np.arange(0, 49, 1):
    #     if i not in labelsx:
    #         string_labels_x.append('')
    #     else:
    #         string_labels_x.append(i)
    #
    # fig = plt.figure()
    # ax3 = fig.add_subplot(111)
    #
    # yx1 = rmsd['mol0']
    # yx2 = rmsd['mol2']
    # yx3 = rmsd['mol3']
    #
    # # ax3.plot(frame_scale, yx1, color='black', lw=1.0, label=label)
    # ax3.plot(frame_scale, yx1, color='black', lw=0.5)
    # ax3.plot(frame_scale, yx2, color='red', lw=0.5)
    # ax3.plot(frame_scale, yx3, color='cyan', lw=0.5)
    #
    # # ax3.plot(frame_scale, yx1, color='black')
    # # ax3.plot(frame_scale, yx2, color='red')
    # # ax3.plot(frame_scale, yx3, color='cyan')
    # # ax3.plot(frame_scale, yx4, color='olive')
    #
    # ax3.set_ylim(y1_lim[0], y1_lim[1])
    # ax3.set_xticklabels(string_labels_x)
    # ax3.tick_params(direction='out', top='off', right='off')
    # plt.xlim(x_lim[0], x_lim[1])
    # ax3.set_yticks(y1_ticks)
    # plt.xticks(x_ticks)
    # # plt.rcParams["figure.figsize"] = (1, 3)
    # fig.set_size_inches(4.0,2.5)
    # plt.rcParams.update({'font.size': 12})
    # plt.rc('font', family='Times New Roman')
    # # plt.xlabel(r'time[ns]')
    # # plt.xlabel(r'time[ns], $^\lambda$')
    # # plt.xlabel(r'O/His73$^{+}$/His526$^{\epsilon}$/Flu$^-$')
    # # plt.xlabel(r'O /His73$\mathregular{^+}$/His526$\mathregular{^\epsilon}$/Flu$\mathregular{^-}$')
    # # plt.xlabel(r'O His73$\mathregular{^+}$His526$\mathregular{^\epsilon}$Flu$\mathregular{^-}$')
    # # plt.ylabel(r'RMSD[$\AA$]')
    # plt.tight_layout()
    # # plt.legend(loc='lower right')
    #
    # plt.savefig(outfolder + 'RMSD_flu_pnas.png', dpi=1200)







    flu_oxi = '/user/jdragelj/projects/cco_fluorescein/rmsd_data/fluO_RMSD.dat'
    rmsd_flu_oxi = parse_rmsd_file(flu_oxi)
    flu_red = '/user/jdragelj/projects/cco_fluorescein/rmsd_data/fluR_RMSD.dat'
    rmsd_flu_red = parse_rmsd_file(flu_red)

    flu_ani = '/user/jdragelj/Desktop/flu_pnas_rmsd/all_cut.dat'
    rmsd_xtra_ani = parse_rmsd_file(flu_ani)

    outfolder = '/user/jdragelj/Desktop/'

    y1_lim = [0, 1.6]
    y1_ticks = np.arange(0, 1.7, 0.2)
    time_at_frames = np.arange(0, 945, 4)
    frame_scale = []
    for time in time_at_frames:
        time = time / 20.0
        frame_scale.append(time)


    time_at_frames1 = np.arange(0, 645, 4)
    frame_scale1 = []
    for time1 in time_at_frames1:
        time1 = time1 / 20.0
        frame_scale1.append(time1)


    x_lim = [0, 48]
    x_ticks = np.arange(0, 49, 1)
    labelsx = np.arange(0, 49, 4)
    string_labels_x = []
    for i in np.arange(0, 49, 1):
        if i not in labelsx:
            string_labels_x.append('')
        else:
            string_labels_x.append(i)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)


    # y1 = rmsd_flu_oxi['flu110'][::4]
    yx2 = rmsd_xtra_ani['mol2'][::2]
    y2 = rmsd_flu_oxi['flu000'][:645][::4]
    y3 = rmsd_flu_oxi['flu11u'][::4]
    y4 = rmsd_flu_oxi['flu100'][::4]

    # ax1.plot(frame_scale, y1, color='black', lw=1.5)
    ax1.plot(frame_scale, yx2, color='black', lw=1.0)
    ax1.plot(frame_scale1, y2, color='green', lw=1.0)
    ax1.plot(frame_scale1, y3, color='blue', lw=1.0)
    ax1.plot(frame_scale1, y4, color='magenta', lw=1.0)

    ax1.set_ylim(y1_lim[0], y1_lim[1])
    ax1.set_xticklabels(string_labels_x)
    ax1.tick_params(direction='out', top='off', right='off')
    plt.xlim(x_lim[0], x_lim[1])
    plt.xticks(x_ticks)

    fig.set_size_inches(4.0,2.5)
    plt.rcParams.update({'font.size': 12})
    plt.rc('font', family='Times New Roman')
    plt.tight_layout()
    plt.savefig(outfolder + 'RMSD_oxi_dominant.png', dpi=1200)


    fig, ax2 = plt.subplots()
    # ax2= fig.add_subplot(111)

    # yx3 = rmsd_xtra_ani['mol3'][::2]
    yx1 = rmsd_xtra_ani['mol0'][::2]
    y5 = rmsd_flu_red['flu11d'][::4]
    y6 = rmsd_flu_red['flu11u'][::4]
    y7 = rmsd_flu_red['flu000'][:645][::4]

    ax2.plot(frame_scale1, y5, color='red', lw=1.0)
    ax2.plot(frame_scale1, y6, color='blue', lw=1.0)
    ax2.plot(frame_scale1, y7, color='green', lw=1.0)
    ax2.plot(frame_scale, yx1, color='magenta', lw=1.0)
    # ax2.plot(frame_scale, yx3, color='blue', lw=1.0)

    ax2.set_ylim(y1_lim[0], y1_lim[1])
    ax2.set_xticklabels(string_labels_x)
    ax2.tick_params(direction='out', top='off', right='off')
    plt.xlim(x_lim[0], x_lim[1])
    plt.xticks(x_ticks)

    fig.set_size_inches(4.0,2.5)
    plt.rcParams.update({'font.size': 12})
    plt.rc('font', family='Times New Roman')
    plt.tight_layout()

    plt.savefig(outfolder + 'RMSD_red_dominant.png', dpi=1200)
