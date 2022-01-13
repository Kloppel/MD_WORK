# coding=utf-8

import matplotlib
matplotlib.use('Agg')

import kbp2
import re
import numpy as np
import pickle
import os
# from kbp2.workspace_jd import plot_tools


# coding=utf-8

# This is storage for general utility external_scripts -> for example transforming data into file types etc

# coding=utf-8

import os
import re
import matplotlib.pyplot as plt
import numpy as np
import kbp2


def get_pka_frame(folder, resi, frame, frame_folder_name):
    '''
    :param folder: folder where the titrations have been done; 'done/' has been added as it is commonly used with pka_manager
    :param resi: residue that has been titrated
    :param frame: frame number
    :param frame_folder_name: name of the folder of the titrated frame
    :return: pka value of the residue in the chosen frame
    '''

    result_file = folder + 'done/%s%i/results.dat' % (frame_folder_name,frame)
    for line in open(result_file, 'r'):
        residue_pka = line.strip()
        entries = re.split(r'[:]', residue_pka)
        residue = entries[0]
        pka = entries[1]
        if residue == resi:
            return pka

def get_pkas_trajectory_residue(result_folder, residue, frames, frame_folder_name):
    '''
    :param result_folder: folder where the titrations have been done
    :param residue: residue that has been titrated
    :param frames: a list of frame numbers for ploting
    :return:
    '''
    results = []
    for frame in frames:
        resi_pka = get_pka_frame(result_folder, residue, frame, frame_folder_name)
        results.append(resi_pka)
    return results

def get_charmm_interation_energy_trajectory(result_folder, frames):
    '''
    :param result_folder: folder where the calculations has been performed
    :param frames: list of frame numbers
    :return:
    '''
    results = []
    for frame in frames:
            energies = kbp2.charmm.get_charmm_energy(result_folder + '/%i/' % frame, inter=True)
            results.append(energies['ELEC'][0])
    return results

def plot_pkas_frames(result_folder, out_folder, out_filename, frames, frame_scale, frame_folder_name, \
                         x_ticks=None, x_lim=None, y_lim=None, y_ticks=None, residue_selection = None):


    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    if out_folder[-1] != '/':
        out_folder = out_folder+'/'

    ### get pkas for residues in question ###
    results_dict = {}
    for resi in residue_selection:
        results_dict[resi] = get_pkas_trajectory_residue(result_folder, resi, frames, frame_folder_name)

    ###############################
    #### PLOT selected residues ###
    ###############################
    if len(residue_selection) > 8:
        raise AssertionError('Too many residues to plot together!')

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    colors = ['blue', 'black', 'red', 'magenta', 'cyan', 'yellow', 'grey', 'orange']
    # linestyles = ['-', '-', '-']
    # marker = ['x', 'o', '.']

    for k, resi in enumerate(residue_selection):
        color = colors[k]
        # color = 'black'
        y = np.array(results_dict[resi], dtype = 'float_')
        y_mean = sum(y)/len(y)
        ax1.plot(frame_scale, y, label='%s: %.2f, max %.2f' % (resi, y_mean, max(y)), lw=0.75, color=color)
        # ax1.plot(frame_scale, y, lw=0.75, color=color)

    # scale, limits, x, y, range
    if y_lim is not None:
        plt.ylim(y_lim[0],y_lim[1])
    if y_ticks is not None:
        plt.yticks(y_ticks)
    if x_lim is not None:
        plt.xlim(x_lim[0],x_lim[1])
    if x_ticks is not None:
        plt.xticks(x_ticks)

    # for frame, value in zip(frame_scale, y):
        # print frame, value

    ax1.tick_params(direction = 'out', top='off', right='off', length=3.0)
    ax1.spines['top'].set_linewidth(0.5)
    ax1.spines['left'].set_linewidth(0.5)
    ax1.spines['right'].set_linewidth(0.5)
    ax1.spines['bottom'].set_linewidth(0.5)
    # ax1.legend(loc='upper right')
    # fig.set_size_inches(3.2, 2.7)
    # plt.plot([2.6, 2.15])
    plt.rcParams.update({'font.size': 12})
    plt.savefig(out_folder + '%s_pka.png' % (out_filename+'_final'))





if __name__ == '__main__':

    # for state in ['Pm', 'PF', 'Pr', 'F']:
    # for state in ['PF', 'Pr', 'F']:
    # for state in ['PF']:
    #     # folder_titrations = '/mnt/fu-scratch//jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_xtra/prd/10ang_cavities_09_voda/' % state
    #     folder_titrations = '/mnt/fu-scratch//jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_xtra/prd/10ang_cavities_08_voda/' % state
    #     out_folder_general = '/user/jdragelj/projects/cco_pls/PRD/titr_plots/discussion/'
    #     trajectories = ['md_%s_prd-' % state, 'md_%s_prdh1' % state, 'md_%s_prdh2' % state]
    #     for trajectory in trajectories:
    #         if state == 'Pr':
    #             frames = np.arange(0,601,4)
    #             time_at_frames = np.arange(0,601,4)
    #             x_lim = [0, 30]
    #             x_ticks = np.arange(0.0, 31.0, 1.0)
    #         else:
    #             frames = np.arange(0,501,4)
    #             time_at_frames = np.arange(0,501,4)
    #             x_lim = [0, 25]
    #             x_ticks = np.arange(0.0, 26.0, 1.0)
    #         residue_selection = ['PRD-3_EHEM']
    #         for trajectory in trajectories:
    #             # out_folder = out_folder_general + '%s/' % state
    #             out_folder = out_folder_general + '%s_08cav/' % state
    #             result_folder = folder_titrations + trajectory + '/'
    #             out_file = trajectory
    #             frame_scale = []
    #             for time in time_at_frames:
    #                 time = time/20.0
    #                 frame_scale.append(time)
    #             frame_folder_name = 'frame'
    #             y_lim = [-10,20]
    #             y_ticks = np.arange(-10, 21, 1.0)
    #             plot_tools.plot_pkas_frames(result_folder, out_folder, out_file, frames, frame_scale,
    #                                         frame_folder_name,x_ticks, x_lim, y_lim, y_ticks, residue_selection)

    # # for state in ['Pm', 'PF', 'Pr', 'F']:
    # for state in ['PF']:
    #     # folder_titrations = '/mnt/fu-scratch//jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_pra_a3/10ang_cavities_09_voda/' % state
    #     folder_titrations = '/mnt/fu-scratch//jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_pra_a3/10ang_cavities_08_voda/' % state
    #     out_folder_general = '/user/jdragelj/projects/cco_pls/PRAa3/titr_plots/discussion/'
    #     trajectories = ['md_pra-', 'md_prah1', 'md_prah2']
    #     for trajectory in trajectories:
    #         if state == 'Pr':
    #             frames = np.arange(0, 601, 4)
    #             time_at_frames = np.arange(0, 601, 4)
    #             x_lim = [0, 30]
    #             x_ticks = np.arange(0.0, 31.0, 1.0)
    #         else:
    #             frames = np.arange(0, 501, 4)
    #             time_at_frames = np.arange(0, 501, 4)
    #             x_lim = [0, 25]
    #             x_ticks = np.arange(0.0, 26.0, 1.0)
    #         residue_selection = ['PRA-1_EHEM']
    #         for trajectory in trajectories:
    #             out_folder = out_folder_general + '%s_08cav/' % state
    #             # out_folder = out_folder_general + '%s/' % state
    #             result_folder = folder_titrations + trajectory + '/'
    #             out_file = trajectory
    #             frame_scale = []
    #             for time in time_at_frames:
    #                 time = time/20.0
    #                 frame_scale.append(time)
    #             frame_folder_name = 'frame'
    #             y_lim = [-10,20]
    #             y_ticks = np.arange(-10, 21, 1.0)
    #             plot_tools.plot_pkas_frames(result_folder, out_folder, out_file, frames, frame_scale,
    #                                         frame_folder_name,x_ticks, x_lim, y_lim, y_ticks, residue_selection)

    # folder_titrations = '/mnt/fu-scratch//jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/10ang_cavities_09_voda/'
    # out_folder_general = '/user/jdragelj/projects/cco_pls/PRAa/titr_plots/discussion/'
    # trajectories = ['md_pra-', 'md_prah1', 'md_prah2']
    # for trajectory in trajectories:
    #     frames = np.arange(0,301,4)
    #     time_at_frames = np.arange(0,301,4)
    #     x_lim = [0, 15]
    #     x_ticks = np.arange(0.0, 16.0, 1.0)
    #     residue_selection = ['PRA-1_EHEM']
    #     for trajectory in trajectories:
    #         out_folder = out_folder_general
    #         result_folder = folder_titrations + trajectory + '/'
    #         out_file = trajectory
    #         frame_scale = []
    #         for time in time_at_frames:
    #             time = time/20.0
    #             frame_scale.append(time)
    #         frame_folder_name = 'frame'
    #         y_lim = [-10,20]
    #         y_ticks = np.arange(-10, 21, 1.0)
    #         plot_tools.plot_pkas_frames(result_folder, out_folder, out_file, frames, frame_scale,
    #                                     frame_folder_name,x_ticks, x_lim, y_lim, y_ticks, residue_selection)
    #
    #



    # for state in ['PF', 'F']:
    #     folder_titrations = '/mnt/fu-scratch//jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_pra_a3/10ang_cavities_09_voda/' % state
    #     out_folder_general = '/user/jdragelj/projects/cco_pls/PRAa3/titr_plots/discussion/'
    #     frames = np.arange(0, 301, 4)
    #     time_at_frames = np.arange(0, 301, 4)
    #     x_lim = [0, 15]
    #     x_ticks = np.arange(0.0, 16.0, 1.0)
    #     residue_selection = ['PRA-1_EHEM', 'ASP-407_ACHA']
    #     out_folder = out_folder_general + '%s/' % state
    #     result_folder = folder_titrations + 'md_%s_prah_asp_dep/' % state
    #     out_file = 'md_%s_prah_asp_dep' % state
    #     frame_scale = []
    #     for time in time_at_frames:
    #         time = time/20.0
    #         frame_scale.append(time)
    #     frame_folder_name = 'frame'
    #     y_lim = [-10,20]
    #     y_ticks = np.arange(-10, 21, 1.0)
    #     plot_tools.plot_pkas_frames(result_folder, out_folder, out_file, frames, frame_scale,
    #                                 frame_folder_name,x_ticks, x_lim, y_lim, y_ticks, residue_selection)

    # for state in ['PF', 'F']:
    #     folder_titrations = '/mnt/fu-scratch//jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_pra_a3/10ang_cavities_09_voda/' % state
    #     out_folder_general = '/user/jdragelj/projects/cco_pls/PRAa3/titr_plots/discussion/'
    #     frames = np.arange(200, 301, 1)
    #     time_at_frames = np.arange(200, 301, 1)
    #     x_lim = [10, 15]
    #     x_ticks = np.arange(10.0, 16.0, 1.0)
    #     residue_selection = ['PRA-1_EHEM', 'ASP-407_ACHA']
    #     out_folder = out_folder_general + '%s/' % state
    #     result_folder = folder_titrations + 'md_%s_prah_asp_dep/' % state
    #     out_file = 'md_%s_prah_asp_dep_detail' % state
    #     frame_scale = []
    #     for time in time_at_frames:
    #         time = time/20.0
    #         frame_scale.append(time)
    #     frame_folder_name = 'frame'
    #     y_lim = [-10,21]
    #     y_ticks = np.arange(-10, 22, 1.0)
    #     plot_tools.plot_pkas_frames(result_folder, out_folder, out_file, frames, frame_scale,
    #                                 frame_folder_name,x_ticks, x_lim, y_lim, y_ticks, residue_selection)
    #
    # for state in ['PF', 'F']:
    #     folder_titrations = '/mnt/fu-scratch//jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_md_crystal/10ang_cavities_09_voda/' % state
    #     out_folder_general = '/user/jdragelj/projects/cco_pls/crystal/titr_plots/discussion/'
    #     frames = np.arange(0, 501, 4)
    #     time_at_frames = np.arange(0, 501, 4)
    #     x_lim = [0, 25]
    #     x_ticks = np.arange(0.0, 26.0, 5.0)
    #     residue_selection = ['PRA-1_EHEM', 'ASP-407_ACHA']
    #     out_folder = out_folder_general + '%s/' % state
    #     result_folder = folder_titrations
    #     out_file = 'md_crystal_%s_pra_asp' % state
    #     frame_scale = []
    #     for time in time_at_frames:
    #         time = time/20.0
    #         frame_scale.append(time)
    #     frame_folder_name = 'frame'
    #     y_lim = [-10,21]
    #     y_ticks = np.arange(-10, 22, 1.0)
    #     plot_tools.plot_pkas_frames(result_folder, out_folder, out_file, frames, frame_scale,
    #                                 frame_folder_name,x_ticks, x_lim, y_lim, y_ticks, residue_selection)

    ######################
    ### plot pkas mix #### publication plots
    ######################
    import matplotlib.pyplot as plt

    outfolder = '/home/users/j/jdragelj/Desktop/'

    title1 = 'pra_hemea'
    frames1 = np.arange(0, 301, 2)
    time_at_frames1 = np.arange(0, 301, 2)
    frame_scale1 = []
    for time in time_at_frames1:
        time = time / 20.0
        frame_scale1.append(time)

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/10ang_cavities_09_voda/md_pra-/'
    pkas_pra_dep_a = get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_GHEM', frames1, 'frame')
    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/10ang_cavities_09_voda/md_prah1/'
    pkas_pra_h1_a = get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_GHEM', frames1, 'frame')
    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/10ang_cavities_09_voda/md_prah2/'
    pkas_pra_h2_a = get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_GHEM', frames1, 'frame')


    x_lim = [0, 15]
    x_ticks = np.arange(0, 16, 1)
    labelsx = np.arange(0, 16, 3)
    string_labels_x = []
    for i in np.arange(0, 16, 1):
        if i not in labelsx:
            string_labels_x.append('')
        else:
            string_labels_x.append(i)

    y1 = pkas_pra_dep_a
    y2 = pkas_pra_h1_a
    y3 = pkas_pra_h2_a

    y1_lim = [-10, 5]
    y1_ticks = np.arange(-10, 6, 1)
    labelsy = np.arange(-10, 6, 2)
    string_labels_y = []
    for i in np.arange(-10, 6, 1):
        if i not in labelsy:
            string_labels_y.append('')
        else:
            string_labels_y.append(i)

    plt.rcParams.update({'font.size': 20})
    fig, ax1 = plt.subplots()

    ax1.plot(frame_scale1, y1, color='purple', lw=1.5)
    ax1.plot(frame_scale1, y2, color='red', lw=1.5)
    ax1.plot(frame_scale1, y3, color='green', lw=1.5)

    ax1.set_ylim(y1_lim[0], y1_lim[1])
    ax1.set_yticks(y1_ticks)
    ax1.set_xticklabels(string_labels_x)
    ax1.set_yticklabels(string_labels_y)
    ax1.tick_params(direction='out', top='off', right='off')
    plt.xlim(x_lim[0], x_lim[1])
    plt.xticks(x_ticks)

    plt.savefig(outfolder + 'PRA_a_pkas.png')

    #
    ###### PRD

    outfolder = '/home/users/j/jdragelj/Desktop/'

    title1 = 'prd_hemea3'
    frames1 = np.arange(0, 501, 2)
    time_at_frames1 = np.arange(0, 501, 2)
    frame_scale1 = []
    for time in time_at_frames1:
        time = time / 20.0
        frame_scale1.append(time)

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_cavities_09_voda/md_Pr_prd-_2/'
    pkas_prd_dep_Pr = get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames1, 'frame')

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_xtra/prd/10ang_cavities_09_voda/md_Pm_prd-/'
    pkas_prd_dep_Pm = get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames1, 'frame')

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_xtra/prd/10ang_cavities_09_voda/md_PF_prd-/'
    pkas_prd_dep_PF = get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames1, 'frame')

    pka_results_folder = '/mnt/fu-scratch//jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_xtra/prd/10ang_cavities_09_voda/md_F_prd-/'
    pkas_prd_dep_F = get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames1, 'frame')

    x_lim = [0, 25]
    x_ticks = np.arange(0, 26, 1)
    labelsx = np.arange(0, 26, 5)
    string_labels_x = []
    for i in np.arange(0, 26, 1):
        if i not in labelsx:
            string_labels_x.append('')
        else:
            string_labels_x.append(i)

    y1 = pkas_prd_dep_Pr
    y2 = pkas_prd_dep_Pm
    y3 = pkas_prd_dep_PF
    y4 = pkas_prd_dep_F

    y1_lim = [-6, 6]
    y1_ticks = np.arange(-6, 7, 1)
    labelsy = np.arange(-6, 7, 2)
    string_labels_y = []
    for i in np.arange(-6, 7, 1):
        if i not in labelsy:
            string_labels_y.append('')
        else:
            string_labels_y.append(i)

    plt.rcParams.update({'font.size': 20})
    fig, ax1 = plt.subplots()

    ax1.plot(frame_scale1, y2, color='orange', lw=1.5)
    ax1.plot(frame_scale1, y3, color='purple', lw=1.5)
    ax1.plot(frame_scale1, y4, color='black', lw=1.5)
    ax1.plot(frame_scale1, y1, color='green', lw=1.5)


    ax1.set_ylim(y1_lim[0], y1_lim[1])
    ax1.set_yticks(y1_ticks)
    ax1.set_xticklabels(string_labels_x)
    ax1.set_yticklabels(string_labels_y)
    ax1.tick_params(direction='out', top='off', right='off')
    plt.xlim(x_lim[0], x_lim[1])
    plt.xticks(x_ticks)

    plt.savefig(outfolder + 'PRD_a3_dominant_pkas.png')
    exit()


    outfolder = '/home/users/j/jdragelj/Desktop/'

    title1 = 'prd_hemea3'
    frames1 = np.arange(0, 501, 4)
    time_at_frames1 = np.arange(0, 501, 4)
    frame_scale1 = []
    for time in time_at_frames1:
        time = time / 20.0
        frame_scale1.append(time)

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_cavities_09_voda/md_Pr_prdh1/'
    pkas_prd_h1_Pr = get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames1, 'frame')
    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_cavities_09_voda/md_Pr_prdh2/'
    pkas_prd_h2_Pr = get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames1, 'frame')

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_xtra/prd/10ang_cavities_09_voda/md_PF_prdh2/'
    pkas_prd_h2_PF = get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames1, 'frame')

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_xtra/prd/10ang_cavities_09_voda/md_F_prdh1/'
    pkas_prd_h1_F = get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames1, 'frame')

    x_lim = [0, 25]
    x_ticks = np.arange(0, 26, 1)
    labelsx = np.arange(0, 26, 5)
    string_labels_x = []
    for i in np.arange(0, 26, 1):
        if i not in labelsx:
            string_labels_x.append('')
        else:
            string_labels_x.append(i)

    y1 = pkas_prd_h1_Pr
    y2 = pkas_prd_h2_Pr
    y3 = pkas_prd_h2_PF
    y4 = pkas_prd_h1_F

    y1_lim = [0, 20]
    y1_ticks = np.arange(0, 21, 1)
    labelsy = np.arange(0, 21, 2)
    string_labels_y = []
    for i in np.arange(0, 21, 1):
        if i not in labelsy:
            string_labels_y.append('')
        else:
            string_labels_y.append(i)

    plt.rcParams.update({'font.size': 20})
    fig, ax1 = plt.subplots()

    ax1.plot(frame_scale1, y2, color='skyblue', lw=1.5)
    ax1.plot(frame_scale1, y3, color='blue', lw=1.5)
    ax1.plot(frame_scale1, y4, color='cyan', lw=1.5)
    ax1.plot(frame_scale1, y1, color='red', lw=1.5)


    ax1.set_ylim(y1_lim[0], y1_lim[1])
    ax1.set_yticks(y1_ticks)
    ax1.set_xticklabels(string_labels_x)
    ax1.set_yticklabels(string_labels_y)
    ax1.tick_params(direction='out', top='off', right='off')
    plt.xlim(x_lim[0], x_lim[1])
    plt.xticks(x_ticks)

    plt.savefig(outfolder + 'PRD_a3_other_pkas.png')
    plt.show()

    exit()

    ###### PRA

    # outfolder = '/user/jdragelj/projects/cco_pls/PRAa3/titr_plots/pub/'
    outfolder = '/home/users/j/jdragelj/Desktop/'

    title1 = 'pra_hemea3'
    frames1 = np.arange(0, 501, 4)
    time_at_frames1 = np.arange(0, 501, 4)
    frame_scale1 = []
    for time in time_at_frames1:
        time = time / 20.0
        frame_scale1.append(time)

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra_a3/10ang_cavities_09_voda/md_pra-/'
    # pkas_pra_dep_Pr = plot_tools.get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames1, 'frame')
    pkas_pra_dep_Pr = get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames1, 'frame')

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_pra_a3/10ang_cavities_09_voda/md_pra-/'
    # pkas_pra_dep_Pm = plot_tools.get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames1, 'frame')
    pkas_pra_dep_Pm = get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames1, 'frame')

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_pra_a3/10ang_cavities_09_voda/md_pra-/'
    # pkas_pra_dep_PF = plot_tools.get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames1, 'frame')
    pkas_pra_dep_PF = get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames1, 'frame')

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_pra_a3/10ang_cavities_09_voda/md_pra-/'
    # pkas_pra_dep_F = plot_tools.get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames1, 'frame')
    pkas_pra_dep_F = get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames1, 'frame')

    x_lim = [0, 25]
    x_ticks = np.arange(0, 26, 1)
    labelsx = np.arange(0, 26, 5)
    string_labels_x = []
    for i in np.arange(0, 26, 1):
        if i not in labelsx:
            string_labels_x.append('')
        else:
            string_labels_x.append(i)

    y1 = pkas_pra_dep_Pr
    y2 = pkas_pra_dep_Pm
    y3 = pkas_pra_dep_PF
    y4 = pkas_pra_dep_F

    y1_lim = [-6, 6]
    y1_ticks = np.arange(-6, 7, 1)
    labelsy = np.arange(-6, 7, 2)
    string_labels_y = []
    for i in np.arange(-6, 7, 1):
        if i not in labelsy:
            string_labels_y.append('')
        else:
            string_labels_y.append(i)

    plt.rcParams.update({'font.size': 20})
    fig, ax1 = plt.subplots()

    ax1.plot(frame_scale1, y2, color='green', lw=1.5)
    ax1.plot(frame_scale1, y3, color='purple', lw=1.5)
    ax1.plot(frame_scale1, y4, color='red', lw=1.5)
    ax1.plot(frame_scale1, y1, color='black', lw=1.5)


    ax1.set_ylim(y1_lim[0], y1_lim[1])
    ax1.set_yticks(y1_ticks)
    ax1.set_xticklabels(string_labels_x)
    ax1.set_yticklabels(string_labels_y)
    ax1.tick_params(direction='out', top='off', right='off')
    plt.xlim(x_lim[0], x_lim[1])
    plt.xticks(x_ticks)

    plt.savefig(outfolder + 'PRA_a3_dominant_pkas.png')
    plt.show()

    # exit()

    outfolder = '/home/users/j/jdragelj/Desktop/'

    title1 = 'pra_hemea3'
    frames1 = np.arange(0, 501, 4)
    time_at_frames1 = np.arange(0, 501, 4)
    frame_scale1 = []
    for time in time_at_frames1:
        time = time / 20.0
        frame_scale1.append(time)

    pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra_a3/10ang_cavities_09_voda/md_prah2/'
    pkas_pra_h2_Pr = get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames1, 'frame')


    x_lim = [0, 25]
    x_ticks = np.arange(0, 26, 1)
    labelsx = np.arange(0, 26, 5)
    string_labels_x = []
    for i in np.arange(0, 26, 1):
        if i not in labelsx:
            string_labels_x.append('')
        else:
            string_labels_x.append(i)

    y1 = pkas_pra_h2_Pr

    y1_lim = [0, 20]
    y1_ticks = np.arange(0, 21, 1)
    labelsy = np.arange(0, 21, 2)
    string_labels_y = []
    for i in np.arange(0, 21, 1):
        if i not in labelsy:
            string_labels_y.append('')
        else:
            string_labels_y.append(i)

    plt.rcParams.update({'font.size': 20})
    fig, ax1 = plt.subplots()

    ax1.plot(frame_scale1, y1, color='orange', lw=1.5)

    ax1.set_ylim(y1_lim[0], y1_lim[1])
    ax1.set_yticks(y1_ticks)
    ax1.set_xticklabels(string_labels_x)
    ax1.set_yticklabels(string_labels_y)
    ax1.tick_params(direction='out', top='off', right='off')
    plt.xlim(x_lim[0], x_lim[1])
    plt.xticks(x_ticks)

    plt.savefig(outfolder + 'PRA_a3_other_pkas.png')
    plt.show()