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

    # rmsd = '/user/jdragelj/Pm_state.dat'
    # rmsd_data = parse_rmsd_file(rmsd)
    # outfolder = '/user/jdragelj/Desktop/'
    # y1_lim = [0, 1.5]
    # y1_ticks = np.arange(0, 1.6, 0.1)
    # time_at_frames = np.arange(0, 501, 4)
    # frame_scale = []
    # for time in time_at_frames:
    #     time = time / 20.0
    #     frame_scale.append(time)
    # x_lim = [0, 26]
    # x_ticks = np.arange(0, 26, 1)
    # labelsx = np.arange(0, 26, 5)
    # string_labels_x = []
    # for i in np.arange(0, 26, 1):
    #     if i not in labelsx:
    #         string_labels_x.append('')
    #     else:
    #         string_labels_x.append(i)
    # fig = plt.figure()
    # ax3 = fig.add_subplot(111)
    # yx1 = rmsd_data['mol26'][::4]
    # yx2 = rmsd_data['mol27'][::4]
    # yx3 = rmsd_data['mol28'][::4]
    # yx4 = rmsd_data['mol29'][::4]
    # yx5 = rmsd_data['mol30'][::4]
    # yx6 = rmsd_data['mol31'][::4]
    # yx7 = rmsd_data['mol32'][::4]
    # ax3.plot(frame_scale, yx2, color='green', lw=0.75)
    # ax3.plot(frame_scale, yx3, color='yellow', lw=0.75)
    # ax3.plot(frame_scale, yx5, color='magenta', lw=0.75)
    # ax3.plot(frame_scale, yx6, color='orange', lw=0.75)
    # ax3.plot(frame_scale, yx7, color='black', lw=0.75)
    # ax3.plot(frame_scale, yx1, color='blue', lw=0.75)
    # ax3.plot(frame_scale, yx4, color='red', lw=0.75)
    # ax3.set_ylim(y1_lim[0], y1_lim[1])
    # ax3.set_xticklabels(string_labels_x)
    # ax3.tick_params(direction='out', top='off', right='off')
    # plt.xlim(x_lim[0], x_lim[1])
    # ax3.set_yticks(y1_ticks)
    # plt.xticks(x_ticks)
    # fig.set_size_inches(4.0,2.5)
    # plt.rcParams.update({'font.size': 10})
    # plt.rc('font', family='Times New Roman')
    # plt.tight_layout()
    # plt.savefig(outfolder + 'RMSD_pm.png', dpi=300)

    outfolder = '/home/users/j/jdragelj/projects/cco_pls/rmsd_PRDa/'

    # rmsd = '/user/jdragelj/PF_state.dat'
    rmsd = outfolder + 'PF_state_fix.dat'
    rmsd_data = parse_rmsd_file(rmsd)

    y1_lim = [0, 1.5]
    y1_ticks = np.arange(0, 1.6, 0.1)
    time_at_frames = np.arange(0, 501, 4)
    frame_scale = []
    for time in time_at_frames:
        time = time / 20.0
        frame_scale.append(time)
    x_lim = [0, 26]
    x_ticks = np.arange(0, 26, 1)
    labelsx = np.arange(0, 26, 5)
    string_labels_x = []
    for i in np.arange(0, 26, 1):
        if i not in labelsx:
            string_labels_x.append('')
        else:
            string_labels_x.append(i)
    fig = plt.figure()
    ax3 = fig.add_subplot(111)
    #####
    yxprdah1 = rmsd_data['prdah1'][::4]
    yxprdah2 = rmsd_data['prdah2'][::4]
    yxprda = rmsd_data['prda-'][::4]
    #####
    yx1 = rmsd_data['mol43'][::4]
    yx2 = rmsd_data['mol44'][::4]
    yx3 = rmsd_data['mol45'][::4]
    yx4 = rmsd_data['mol46'][::4]
    yx5 = rmsd_data['mol47'][::4]
    yx6 = rmsd_data['mol48'][::4]
    yx7 = rmsd_data['mol51'][::4]
    ####-
    ax3.plot(frame_scale, yxprdah1, color='grey', lw=0.75)
    ax3.plot(frame_scale, yxprdah2, color='lime', lw=0.75)
    ax3.plot(frame_scale, yxprda, color='brown', lw=0.75)
    #####
    ax3.plot(frame_scale, yx2, color='green', lw=0.75)
    ax3.plot(frame_scale, yx3, color='yellow', lw=0.75)
    ax3.plot(frame_scale, yx5, color='magenta', lw=0.75)
    ax3.plot(frame_scale, yx6, color='orange', lw=0.75)
    ax3.plot(frame_scale, yx7, color='black', lw=0.75)
    ax3.plot(frame_scale, yx1, color='blue', lw=0.75)
    ax3.plot(frame_scale, yx4, color='red', lw=0.75)
    ax3.set_ylim(y1_lim[0], y1_lim[1])
    ax3.set_xticklabels(string_labels_x)
    ax3.tick_params(direction='out', top='off', right='off')
    plt.xlim(x_lim[0], x_lim[1])
    ax3.set_yticks(y1_ticks)
    plt.xticks(x_ticks)
    fig.set_size_inches(4.0,2.5)
    plt.rcParams.update({'font.size': 10})
    plt.rc('font', family='Times New Roman')
    plt.tight_layout()
    plt.savefig(outfolder + 'RMSD_pf.png', dpi=300)


    # rmsd = '/user/jdragelj/Pr_state.dat'
    # rmsd_data_o = parse_rmsd_file(rmsd)
    # rmsd_data = {}
    # for key, values in rmsd_data_o.iteritems():
    #     new_values = []
    #     for i, value in enumerate(values):
    #         if i < 501:
    #             new_values.append(value)
    #     rmsd_data[key]=new_values
    # outfolder = '/user/jdragelj/Desktop/'
    # y1_lim = [0, 1.5]
    # y1_ticks = np.arange(0, 1.6, 0.1)
    # time_at_frames = np.arange(0, 501, 4)
    # frame_scale = []
    # for time in time_at_frames:
    #     time = time / 20.0
    #     frame_scale.append(time)
    # x_lim = [0, 26]
    # x_ticks = np.arange(0, 26, 1)
    # labelsx = np.arange(0, 26, 5)
    # string_labels_x = []
    # for i in np.arange(0, 26, 1):
    #     if i not in labelsx:
    #         string_labels_x.append('')
    #     else:
    #         string_labels_x.append(i)
    # fig = plt.figure()
    # ax3 = fig.add_subplot(111)
    # yx1 = rmsd_data['mol33'][::4]
    # yx2 = rmsd_data['mol34'][::4]
    # yx3 = rmsd_data['mol35'][::4]
    # yx4 = rmsd_data['mol36'][::4]
    # yx5 = rmsd_data['mol37'][::4]
    # yx6 = rmsd_data['mol38'][::4]
    # yx7 = rmsd_data['mol39'][::4]
    # yx8 = rmsd_data['mol40'][::4]
    # yx9 = rmsd_data['mol41'][::4]
    # yx10 = rmsd_data['mol42'][::4]
    # rmsd = '/user/jdragelj/Pr_extra_rmsd.dat'
    # rmsd_data_o = parse_rmsd_file(rmsd)
    # rmsd_data = {}
    # for key, values in rmsd_data_o.iteritems():
    #     new_values = []
    #     for i, value in enumerate(values):
    #         if i < 501:
    #             new_values.append(value)
    #     rmsd_data[key]=new_values
    # yx11 = rmsd_data['mol6'][::4]
    # ax3.plot(frame_scale, yx2, color='green', lw=0.75)
    # ax3.plot(frame_scale, yx3, color='yellow', lw=0.75)
    # ax3.plot(frame_scale, yx5, color='magenta', lw=0.75)
    # ax3.plot(frame_scale, yx6, color='black', lw=0.75)
    # # ax3.plot(frame_scale, yx7, color='aqua', lw=0.75)
    # time_at_frames1 = np.arange(0, 301, 4)
    # frame_scale1 = []
    # for time in time_at_frames1:
    #     time = time / 20.0
    #     frame_scale1.append(time)
    # ax3.plot(frame_scale1, yx8, color='brown', lw=0.75)
    # ax3.plot(frame_scale1, yx9, color='grey', lw=0.75)
    # ax3.plot(frame_scale1, yx10, color='lime', lw=0.75)
    # ax3.plot(frame_scale, yx11, color='orange', lw=0.75)
    # ax3.plot(frame_scale, yx1, color='blue', lw=0.75)
    # ax3.plot(frame_scale, yx4, color='red', lw=0.75)
    # ax3.set_ylim(y1_lim[0], y1_lim[1])
    # ax3.set_xticklabels(string_labels_x)
    # ax3.tick_params(direction='out', top='off', right='off')
    # plt.xlim(x_lim[0], x_lim[1])
    # ax3.set_yticks(y1_ticks)
    # plt.xticks(x_ticks)
    # fig.set_size_inches(4.0,2.5)
    # plt.rcParams.update({'font.size': 10})
    # plt.rc('font', family='Times New Roman')
    # plt.tight_layout()
    # plt.savefig(outfolder + 'RMSD_pr.png', dpi=300)
    #
    #
    # rmsd = '/user/jdragelj/F_state.dat'
    # rmsd_data = parse_rmsd_file(rmsd)
    # outfolder = '/user/jdragelj/Desktop/'
    # y1_lim = [0, 1.5]
    # y1_ticks = np.arange(0, 1.6, 0.1)
    # time_at_frames = np.arange(0, 501, 4)
    # frame_scale = []
    # for time in time_at_frames:
    #     time = time / 20.0
    #     frame_scale.append(time)
    # x_lim = [0, 26]
    # x_ticks = np.arange(0, 26, 1)
    # labelsx = np.arange(0, 26, 5)
    # string_labels_x = []
    # for i in np.arange(0, 26, 1):
    #     if i not in labelsx:
    #         string_labels_x.append('')
    #     else:
    #         string_labels_x.append(i)
    # fig = plt.figure()
    # ax3 = fig.add_subplot(111)
    # yx1 = rmsd_data['mol17'][::4]
    # yx2 = rmsd_data['mol18'][::4]
    # yx3 = rmsd_data['mol19'][::4]
    # yx4 = rmsd_data['mol20'][::4]
    # yx5 = rmsd_data['mol21'][::4]
    # yx6 = rmsd_data['mol22'][::4]
    # yx7 = rmsd_data['mol25'][::4]
    # ax3.plot(frame_scale, yx2, color='green', lw=0.75)
    # ax3.plot(frame_scale, yx3, color='yellow', lw=0.75)
    # ax3.plot(frame_scale, yx5, color='magenta', lw=0.75)
    # ax3.plot(frame_scale, yx6, color='orange', lw=0.75)
    # ax3.plot(frame_scale, yx7, color='black', lw=0.75)
    # ax3.plot(frame_scale, yx1, color='blue', lw=0.75)
    # ax3.plot(frame_scale, yx4, color='red', lw=0.75)
    # ax3.set_ylim(y1_lim[0], y1_lim[1])
    # ax3.set_xticklabels(string_labels_x)
    # ax3.tick_params(direction='out', top='off', right='off')
    # plt.xlim(x_lim[0], x_lim[1])
    # ax3.set_yticks(y1_ticks)
    # plt.xticks(x_ticks)
    # fig.set_size_inches(4.0,2.5)
    # plt.rcParams.update({'font.size': 10})
    # plt.rc('font', family='Times New Roman')
    # plt.tight_layout()
    # plt.savefig(outfolder + 'RMSD_f.png', dpi=300)








