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


    folder = '/work/jdragelj/projects/cco_pls/MD/'
    rmsd = folder + 'Pm.dat'
    rmsd_data = parse_rmsd_file(rmsd)
    outfolder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/new_models/PRD_H_PRA_a3/PLOTS/'
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
    # fig = plt.figure()
    # ax3 = fig.add_subplot(111)
    # yx1 = rmsd_data['mol0'][::4]
    # yx2 = rmsd_data['mol1'][::4]
    # yx3 = rmsd_data['mol2'][::4]
    # ax3.plot(frame_scale, yx1, color='black', lw=0.75)
    # ax3.plot(frame_scale, yx2, color='red', lw=0.75)
    # ax3.plot(frame_scale, yx3, color='blue', lw=0.75)
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
    #
    #
    # rmsd = folder + 'PF.dat'
    # rmsd_data = parse_rmsd_file(rmsd)
    # # y1_lim = [0, 1.6]
    # # y1_ticks = np.arange(0, 1.7, 0.1)
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
    # yx1 = rmsd_data['mol21'][::4]
    # yx2 = rmsd_data['mol22'][::4]
    # yx3 = rmsd_data['mol23'][::4]
    # ax3.plot(frame_scale, yx1, color='black', lw=0.75)
    # ax3.plot(frame_scale, yx2, color='red', lw=0.75)
    # ax3.plot(frame_scale, yx3, color='blue', lw=0.75)
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
    # plt.savefig(outfolder + 'RMSD_pf.png', dpi=300)
    #
    #
    # rmsd = folder + 'Pr.dat'
    # rmsd_data = parse_rmsd_file(rmsd)
    # # y1_lim = [0, 1.6]
    # # y1_ticks = np.arange(0, 1.7, 0.1)
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
    # yx1 = rmsd_data['mol25'][::4]
    # yx2 = rmsd_data['mol27'][::4]
    # yx3 = rmsd_data['mol28'][::4]
    # ax3.plot(frame_scale, yx1, color='black', lw=0.75)
    # ax3.plot(frame_scale, yx2, color='red', lw=0.75)
    # ax3.plot(frame_scale, yx3, color='blue', lw=0.75)
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


    rmsd = folder + 'F.dat'
    rmsd_data = parse_rmsd_file(rmsd)
    # y1_lim = [0, 1.6]
    # y1_ticks = np.arange(0, 1.7, 0.1)
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
    yx1 = rmsd_data['mol4'][::4]
    yx2 = rmsd_data['mol5'][::4]
    yx3 = rmsd_data['mol6'][::4]
    ax3.plot(frame_scale, yx1, color='black', lw=0.75)
    ax3.plot(frame_scale, yx2, color='red', lw=0.75)
    ax3.plot(frame_scale, yx3, color='blue', lw=0.75)
    ax3.set_ylim(y1_lim[0], y1_lim[1])
    ax3.set_xticklabels(string_labels_x)
    ax3.tick_params(direction='out', top='off', right='off')
    plt.xlim(x_lim[0], x_lim[1])
    ax3.set_yticks(y1_ticks)
    plt.xticks(x_ticks)
    fig.set_size_inches(4.0,2.5)

    # from matplotlib import rc
    # rc('text', usetex=True)
    # fig.text(2.0, 0.1, 'Time / ns', fontsize=14, ha='center')
    # fig.text(0.1, 0.125, 'RMSD / \AA', fontsize=14, rotation='horizontal', va='center')
    # fig.suptitle('F', x=0.46, y=0.98, fontsize=16)
    # leg = fig.legend(bbox_to_anchor=(0.8, 0.5), loc='center left')
    # plt.subplots_adjust(right=0.8)

    font = {'family': 'calibri',
            'size': 10}

    plt.rc('font', **font)

    # plt.rcParams.update({'font.size': 10})
    # plt.rc('font', family='Calibri')
    plt.tight_layout()

    plt.savefig('/home/users/j/jdragelj/Desktop/RMSD_test.png', dpi=300)
    # plt.savefig(outfolder + 'RMSD_f.png', dpi=300)

    plt.show()









