# coding=utf-8

import os
import re

import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.rc({'font': "Times New Roman"})

import numpy as np
import kbp2
from workspace_jd import tools, salt_bridge_distances
# from kbp2.workspace_jd.cco_fluorescin import interaction_ener


#todo: Redo functions, make them general

def get_pka_frame(folder, resi, frame, frame_folder_name):
    '''
    :param folder:
    :param resi:
    :param frame:
    :return:
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
    :param result_folder:
    :param residue:
    :param frames:
    :return:
    '''
    results = []
    for frame in frames:
        resi_pka = get_pka_frame(result_folder, residue, frame, frame_folder_name)
        results.append(resi_pka)
    return results

def plot_pkas_frames(result_folder, out_folder, out_file, title, frames, frame_scale,
                     frame_folder_name, residue_selection = None, residue_labels=None):
    """0
    :param out_folder:
    :param out_file:
    :param title:
    :param frames:
    :param frame_scale:
    :param residue_selection:
    :return:
    """

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    if out_folder[-1] != '/':
        out_folder = out_folder+'/'

    results_dict = {}
    for resi in residue_selection:
        results_dict[resi] = get_pkas_trajectory_residue(result_folder, resi, frames, frame_folder_name)

    ####################
    #### PLOT single ###
    ####################
    x = frame_scale
    print max(x), 'ns'
    for j, residue in enumerate(residue_selection):

        residue_folder = out_folder + 'resi_' + residue + '/'
        if not os.path.exists(residue_folder):
            os.mkdir(residue_folder)
        if residue_folder[-1] != '/':
            residue_folder = residue_folder+'/'

        # y - pKas
        y = np.array(results_dict[residue], dtype = 'float_')
        # average value
        y_mean = [np.mean(y) for i in x]
        plt.figure()
        # plt.ylim(0,15)
        plt.plot(x,y, label=residue_labels[j], marker='.', color='black')
        plt.plot(x,y_mean, label='Avg. %.2f' % y_mean[0], linestyle=':', color='black')
        plt.title(title)
        plt.xlabel('t(ns)')
        plt.ylabel('pKa')
        plt.legend(loc='lower right')
        plt.savefig(residue_folder + '%s_pka.png' % (out_file+'_'+residue))
        plt.close()



    #################
    #### PLOT all ###
    #################
    if len(residue_selection) > 1:
        if len(residue_selection) <= 3:
            x = frame_scale
            plt.figure()
            linestyles = ['-','-','-']
            marker = ['x','o','.']
            for k, resi in enumerate(residue_selection):
                y = np.array(results_dict[resi], dtype = 'float_')
                y_mean = [np.mean(y) for i in x]
                style = linestyles[k]
                plt.plot(x,y, label=residue_labels[k], marker=marker[k], linestyle=style, color='k')
                plt.plot(x,y_mean,  label='Avg. %.2f' % y_mean[0], linestyle='--', color='k')


            plt.xlabel('t(ns)')
            plt.ylabel('pKa')

            # plt.ylim(0,15)

            plt.legend(loc='upper left',prop={'size':10})
            plt.title(title)
            plt.rcParams.update({'font.size': 17})
            plt.savefig(out_folder + '%s_pka.png' % (out_file+'_final'))

        else:
            from termcolor import colored
            colors_dict = {'b': 'blue', \
                'g': 'green', \
                'r': 'red', \
                'c': 'cyan', \
                'm': 'magenta', \
                'y': 'yellow', \
                'k': 'grey', \
                'grey': 'white'}

            x = frame_scale
            plt.figure()
            for k, resi in enumerate(residue_selection):
                y = np.array(results_dict[resi], dtype = 'float_')
                y_mean = [np.mean(y) for i in x]


                colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'grey']
                color=colors[k]
                plt.plot(x,y, color=color)
                plt.plot(x,y_mean,  label= resi + ': %.2f' % y_mean[0], color =color, linestyle='-')

                print colored(resi, colors_dict[color]), colored(round(y_mean[0], 2), colors_dict[color])

            # plt.ylim(-10,20)
            yticks = np.arange(-10.0, 21.0, 1.0)
            plt.yticks(yticks)

            plt.xlabel('t(ns)')
            plt.ylabel('pKa')

            plt.legend(loc='upper left', prop={'size':7})
            plt.title(title)
            plt.rcParams.update({'font.size': 17})
            plt.savefig(out_folder + '%s_pka.png' % (out_file+'_final'))



def plot_pkas_frames_pub(result_folder, out_folder, out_file, title, frames, frame_scale, frame_folder_name, \
                         x_ticks=None, x_lim=None, y_lim=None, y_ticks=None, residue_selection = None, \
                         residue_labels=None):
    """
    :param out_folder:
    :param out_file:
    :param title:
    :param frames:
    :param frame_scale:
    :param residue_selection:
    :return:
    """

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    if out_folder[-1] != '/':
        out_folder = out_folder+'/'

    results_dict = {}
    for resi in residue_selection:
        results_dict[resi] = get_pkas_trajectory_residue(result_folder, resi, frames, frame_folder_name)

    ####################
    #### PLOT single ###
    ####################
    x = frame_scale
    for j, residue in enumerate(residue_selection):
        residue_folder = out_folder + 'resi_' + residue + '/'
        if not os.path.exists(residue_folder):
            os.mkdir(residue_folder)
        if residue_folder[-1] != '/':
            residue_folder = residue_folder+'/'
        # y - pKas
        y = np.array(results_dict[residue], dtype = 'float_')
        # average value
        y_mean = [np.mean(y) for i in x]
        plt.figure()
        # scale, limits, x, y, range
        if y_lim is not None:
            plt.ylim(y_lim[0],y_lim[1])
        if y_ticks is not None:
            plt.yticks(y_ticks)
        if x_lim is not None:
            plt.xlim(x_lim[0],x_lim[1])
        if x_ticks is not None:
            plt.xticks(x_ticks)
        plt.plot(x,y, label=residue_labels[j], marker='.', color='black')
        plt.plot(x,y_mean, label='Avg. %.2f' % y_mean[0], linestyle=':', color='black')
        plt.title(title)
        plt.savefig(residue_folder + '%s_pka.png' % (out_file+'_'+residue))
        plt.close()

    #################
    #### PLOT all ###
    #################
    if len(residue_selection) > 1:
        if len(residue_selection) > 3:
            raise AssertionError('Too many residues to plot together!')
        x = frame_scale
        plt.figure()
        # linestyles = ['--','-','_']
        linestyles = ['-','-','-']
        marker = ['x','.',' ']
        # colors = ['k','k','grey']
        colors = ['k','grey','k']
        for k, resi in enumerate(residue_selection):
            y = np.array(results_dict[resi], dtype = 'float_')
            y_mean = [np.mean(y) for i in x]
            style = linestyles[k]
            plt.plot(x,y, label=residue_labels[k], marker=marker[k], linestyle=linestyles[k], color=colors[k])
            # plt.plot(x,y, label=residue_labels[k], marker=marker[k], linestyle=linestyles[k], color='k')
            # plt.plot(x,y, label=residue_labels[k], marker='.', linestyle=style, color='k')
            # plt.plot(x,y_mean,  label='Avg. %.2f' % y_mean[0], linestyle=':', color='k')
            # plt.plot(x,y_mean,  label='Avg. %.2f' % y_mean[0], linestyle='-', color='k')
            print resi, np.mean(y), marker[k]

        # scale, limits, x, y, range
        if y_lim is not None:
            plt.ylim(y_lim[0],y_lim[1])
        if y_ticks is not None:
            plt.yticks(y_ticks)
        if x_lim is not None:
            plt.xlim(x_lim[0],x_lim[1])
        if x_ticks is not None:
            plt.xticks(x_ticks)

        plt.rcParams.update({'font.size': 20})
        plt.savefig(out_folder + '%s_pka.png' % (out_file+'_final'))

def make_csv_many_trajectories(property, source_folder, out_folder, trajectories, frames, frame_folder_name, residue_selection):
    '''
    :param property:
    :param source_folder:
    :param out_folder:
    :param trajectories:
    :param frames:
    :param residue_selection:
    :return:
    '''

    results_dict = {}
    for residue in residue_selection:
        results_dict[residue]={}
        for trajectory in trajectories:
            results_dict[residue][trajectory] = []
            result_folder = source_folder + trajectory + '/'
            if property == 'pka':
                results_dict[residue][trajectory] = get_pkas_trajectory_residue(result_folder, residue, frames, frame_folder_name)
            # elif property == 'sasa':
            #     results_dict[residue][trajectory] = get_sasa_trajectory_residue(result_folder, residue, frames)
            elif property == 'charm_inter_elec':
                error =  'Not done yet!'
                raise AssertionError(error)

    for residue in residue_selection:
        csv_filename = out_folder + residue + '_%s.csv' % property
        comma_file = open(csv_filename, 'w')
        comma_file.write(tools.string_creation([' '] + list(frames)))
        comma_file.write('\n')
        for trajectory in trajectories:
            comma_file.write(tools.string_creation([trajectory] + results_dict[residue][trajectory]))
            comma_file.write('\n')

def frames_plot_two_properties_pub(md_name, prop1_name, prop2_name, prop1_data, prop2_data, time_scale, outfolder,
                                   x_ticks, x_lim, y1_lim, y1_ticks, y2_lim, y2_ticks, string_labels_x,
                                   string_labels_y1, string_labels_y2):

    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    if outfolder[-1] != '/':
        outfolder = outfolder+'/'

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # ax2 = fig.add_subplot(111)
    ax2 = ax1.twinx()


    # ax1.plot(time_scale, prop1_data, color='black', lw=0.75)
    ax1.plot(time_scale, prop1_data, color='red', lw=0.75)
    # ax2.plot(time_scale, prop2_data, linestyle=':', color='black', marker='.', lw=1.0, markersize=1)
    # ax2.plot(time_scale, prop2_data, color='green', lw=0.75)
    ax2.plot(time_scale, prop2_data, color='blue', lw=0.75)

    # ax1.set_xlabel('time')

    # ax1.set_ylabel(prop1_name, color='g')
    # ax2.set_ylabel(prop2_name, color='b')
    # ax1.legend(loc='upper left')
    # ax2.legend(loc='upper right')

    ax1.set_ylim(y1_lim[0], y1_lim[1])
    ax1.set_yticks(y1_ticks)
    ax2.set_ylim(y2_lim[0], y2_lim[1])
    ax2.set_yticks(y2_ticks)

    ax1.set_xticklabels(string_labels_x)
    ax2.set_xticklabels(string_labels_x)

    ax1.set_yticklabels(string_labels_y1)
    ax2.set_yticklabels(string_labels_y2)

    # ax1.tick_params(direction='in')
    # ax1.tick_params(direction='out')
    # ax2.tick_params(direction='in')
    # ax2.tick_params(direction='out')

    ax1.tick_params(direction='out',top='off', length=3.0)
    ax2.tick_params(direction = 'out', top='off', bottom='off', left='off', length=3.0)

    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_linewidth(0.5)
    ax2.spines['left'].set_linewidth(0.5)
    ax2.spines['right'].set_linewidth(0.5)
    ax2.spines['bottom'].set_linewidth(0.5)
    fig.set_size_inches(3.2, 2.7)
    plt.rcParams.update({'font.size': 10})
    plt.tight_layout()

    plt.xlim(x_lim[0], x_lim[1])
    plt.xticks(x_ticks)

    plt.savefig(outfolder + '%s_%s_%s.png' % (md_name, prop1_name, prop2_name), dpi=300)
    # plt.savefig(outfolder + '%s_%s_%s.png' % (md_name, prop1_name, prop2_name))
    # plt.show()


def frames_plot_three_properties_pub(md_name, prop1_name, prop2_name, prop1_data, prop2_data, prop3_data, time_scale, outfolder,
                                   x_ticks, x_lim, y1_lim, y1_ticks, y2_lim, y2_ticks, string_labels_x,
                                   string_labels_y1, string_labels_y2):

    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    if outfolder[-1] != '/':
        outfolder = outfolder+'/'

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # ax2 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    # ax1.plot(time_scale, prop1_data, color='black', lw=0.75)
    # ax2.plot(time_scale, prop2_data, linestyle=':', color='black', marker='.', lw=1.0, markersize=1)
    # ax2.plot(time_scale, prop3_data, linestyle=':', color='grey', marker='.', lw=1.0, markersize=1)

    ax1.plot(time_scale, prop1_data, color='black', lw=0.75)
    ax2.plot(time_scale, prop2_data, linestyle='-', color='magenta', marker='x', lw=0.5, markersize=2.0)
    ax2.plot(time_scale, prop3_data, linestyle='-', color='green', marker='.', lw=0.5, markersize=3.0)

    for i,j,k in zip(time_scale, prop2_data, prop3_data):
        if i*20 in [80,100,120,140,160,180,200,300]:
            print i, j, k
        if i*20 in [300, 330, 440]:
            print i, j, k



    # ax1.set_xlabel('time')

    # ax1.set_ylabel(prop1_name, color='g')
    # ax2.set_ylabel(prop2_name, color='b')
    # ax1.legend(loc='upper left')
    # ax2.legend(loc='upper right')

    ax1.set_ylim(y1_lim[0], y1_lim[1])
    ax1.set_yticks(y1_ticks)
    ax2.set_ylim(y2_lim[0], y2_lim[1])
    ax2.set_yticks(y2_ticks)

    ax1.set_xticklabels(string_labels_x)
    ax2.set_xticklabels(string_labels_x)

    ax1.set_yticklabels(string_labels_y1)
    ax2.set_yticklabels(string_labels_y2)

    # ax1.tick_params(direction='in')
    # ax1.tick_params(direction='out')
    # ax2.tick_params(direction='in')
    # ax2.tick_params(direction='out')

    ax1.tick_params(direction='out',top='off', length=3.0)
    ax2.tick_params(direction = 'out', top='off', bottom='off', left='off', length=3.0)

    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(True)
    # ax1.spines['left'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    # ax2.spines['top'].set_linewidth(0.5)
    ax1.spines['left'].set_linewidth(0.5)
    # ax2.spines['left'].set_linewidth(0.5)
    ax2.spines['right'].set_linewidth(0.5)
    ax2.spines['bottom'].set_linewidth(0.5)

    ax2.spines['right'].set_color('grey')
    # ax2.yaxis.label.set_color('black')
    ax2.tick_params(axis='y', colors='black')



    fig.set_size_inches(3.2, 2.7)
    plt.rcParams.update({'font.size': 10})
    plt.tight_layout()

    plt.xlim(x_lim[0], x_lim[1])
    plt.xticks(x_ticks)

    plt.savefig(outfolder + '%s_%s_%s.png' % (md_name, prop1_name, prop2_name), dpi=300)
    # plt.savefig(outfolder + '%s_%s_%s.png' % (md_name, prop1_name, prop2_name))
    # plt.show()


if __name__ == '__main__':

    # for state in ['Pm', 'PF', 'Pr', 'F']:
    #     for prot in ['prd-']:
    #
    #         # out_folder_general = '/user/jdragelj/projects/cco_pls/PRD/salt-bridges_double/%s/' % state
    #         out_folder_general = '/home/users/j/jdragelj/projects/cco_pls/sb_plots/new/'
    #
    #         frames_folder_general = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/new_models/PRD_H_PRA_a3/MD/%s_prd-/' % state
    #
    #         # frames_folder_general = '/mnt/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_%s_xtra6/' % (state, state, prot)
    #         # frames_folder_general = '/mnt/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_PF_crystal/'
    #         # y1_lim = [-10, 20]
    #         # y1_ticks = np.arange(-10, 21, 1)
    #         y1_lim = [-10, 8]
    #         y1_ticks = np.arange(-10, 9, 1)
    #         y2_lim = [4, 8]
    #         y2_ticks = np.arange(4, 8.1, 0.1)
    #
    #         # if state != 'Pr':
    #         if state != 'X':
    #             x_lim = [0, 25]
    #             x_ticks = np.arange(0, 26, 1)
    #         else:
    #             x_lim = [0, 30]
    #             x_ticks = np.arange(0, 31, 1)
    #
    #         frames = np.arange(0, 501, 10)
    #         time_at_frames1 = np.arange(0, 501, 10)
    #
    #         frames_folder1 = frames_folder_general + 'frames_voda/'
    #         # frames_folder1 = frames_folder_general + 'frames_voda_1/'
    #         frame_scale1 = []
    #         for time in time_at_frames1:
    #             time = time / 20.0
    #             frame_scale1.append(time)
    #
    #         # if state != 'Pr':
    #         if state != 'X':
    #             labelsx = np.arange(0, 26, 5)
    #             string_labels_x = []
    #             for i in np.arange(0, 26, 1):
    #                 if i not in labelsx:
    #                     string_labels_x.append('')
    #                 else:
    #                     string_labels_x.append(i)
    #         else:
    #             # labelsx = np.arange(0, 31, 2)
    #             # string_labels_x = []
    #             # for i in np.arange(0, 31, 2):
    #             #     if i not in labelsx:
    #             #         string_labels_x.append('')
    #             #     else:
    #             #         string_labels_x.append(i)
    #
    #             labelsx = np.arange(0, 31, 5)
    #             string_labels_x = []
    #             for i in np.arange(0, 31, 1):
    #                 if i not in labelsx:
    #                     string_labels_x.append('')
    #                 else:
    #                     string_labels_x.append(i)
    #
    #         labelsy1 = np.arange(-10, 9, 2)
    #         string_labels_y1 = []
    #         for i in np.arange(-10, 9, 1):
    #             if i not in labelsy1:
    #                 string_labels_y1.append('')
    #             else:
    #                 string_labels_y1.append(i)
    #
    #         # labelsy1 = np.arange(-10, 21, 2)
    #         # string_labels_y1 = []
    #         # for i in np.arange(-10, 21, 1):
    #         #     if i not in labelsy1:
    #         #         string_labels_y1.append('')
    #         #     else:
    #         #         string_labels_y1.append(i)
    #
    #         labelsy2 = np.arange(4.0, 9.0, 1.0)
    #         string_labels_y2 = []
    #         # for i in np.arange(4.0, 8.1, 0.1):
    #         for i in y2_ticks:
    #             # if i not in labelsy2:
    #             #     string_labels_y2.append('')
    #             # else:
    #             #     print i
    #             #     string_labels_y2.append(int(i))
    #
    #             found = False
    #             for k in labelsy2:
    #                 # print abs(round(i - k,1))
    #                 if abs(round(i - k,1)) == 0.0:
    #                     string_labels_y2.append(round(i,0))
    #                     found = True
    #             if not found:
    #                     string_labels_y2.append('')
    #
    #
    #         # pka_results_folder = '/mnt/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_xtra/prd/10ang_cavities_09_voda/md_%s_%s/' % (state, state, prot)
    #         pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/new_models/PRD_H_PRA_a3/titrate/md_%s_prd-/' % state
    #         # pka_results_folder = '/mnt/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_md_crystal/10ang_cavities_09_voda/'
    #
    #         pkas = get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames, 'frame')
    #         # pkas = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEMEA3', frames, 'frame')
    #
    #         name = 'md_%s_%s' % (state, prot)
    #
    #         # sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder1, frames, 'PRD-3_EHEM',
    #         #                                                          'ARG-481_ACHA', name, out_folder_general,
    #         #                                                          double='PRD-3_EHEM')
    #         # # sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder1, frames, 'PRD-3_EHEMEA3', 'ARG-481_ACHAIN', name, out_folder_general, double='PRD-3_EHEMEA3')
    #         #
    #         # sb_dist_min = []
    #         # sb_dist_other = []
    #         # for dist_tuple in sb_dist:
    #         #     sb_dist_min.append(dist_tuple[0])
    #         #     sb_dist_other.append(dist_tuple[1])
    #         #
    #         # csv_file = '/home/users/j/jdragelj/Desktop/sb_plots/data.dat'
    #         # f = open(csv_file, 'w')
    #         # # f.write('time,d1,d2 \n')
    #         # for time, pka, d1, d2 in zip(frame_scale1, pkas, sb_dist_min, sb_dist_other):
    #         #     # print time, pka, d1, d2
    #         #     f.write('%s,%s,%s \n' % (str(time), str(d1), str(d2)))
    #
    #         # frames_plot_three_properties_pub('md_%s_%s' % (state, prot), 'pka', 'sb_dist', pkas,
    #         #                                                        sb_dist_min, sb_dist_other, frame_scale1,
    #         #                                                        out_folder_general, x_ticks, x_lim, y1_lim, y1_ticks,
    #         #                                                        y2_lim,
    #         #                                                        y2_ticks, string_labels_x, string_labels_y1,
    #         #                                                        string_labels_y2)
    #
    #         sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder1, frames, 'PRD-3_EHEM',
    #                                                                  'ARG-481_ACHA', name, out_folder_general)
    #
    #         print min(sb_dist)
    #
    #         frames_plot_two_properties_pub(md_name='md_%s_%s' % (state, prot), prop1_name='pka', prop2_name='sb_dist', prop1_data=pkas,
    #                                          prop2_data=sb_dist, time_scale=frame_scale1,
    #                                          outfolder=out_folder_general, x_ticks=x_ticks, x_lim=x_lim, y1_lim=y1_lim, y1_ticks=y1_ticks,
    #                                          y2_lim=y2_lim,
    #                                          y2_ticks=y2_ticks, string_labels_x=string_labels_x, string_labels_y1=string_labels_y1,
    #                                          string_labels_y2=string_labels_y2)
    #

    # for state in ['F', 'PF', 'Pm', 'Pr']:
    #     for prot in ['prah1', 'prah2', 'pra-']:
    #todo: keep this, this is asp data
    for state in ['F']:
        for prot in ['prah2']:

            out_folder_general = '/home/users/j/jdragelj/Desktop/%s/' % state
            if state != 'Pr':
                frames_folder_general = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/%s_state/md_%s/' % (state,prot)
            else:
                frames_folder_general = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/md_%s/' % (prot)

            # y1_lim = [-10, 20]
            # y1_ticks = np.arange(-10, 21, 1)
            y1_lim = [4, 20]
            y1_ticks = np.arange(4, 21, 1)
            y2_lim = [2, 8]
            y2_ticks = np.arange(2, 9, 1)

            if state != 'Pr':
                x_lim = [0, 25]
                x_ticks = np.arange(0, 26, 1)
            else:
                x_lim = [0, 31]
                x_ticks = np.arange(0, 31, 2)

            if state != 'Pr':
                frames = np.arange(0, 501, 4)
                time_at_frames1 = np.arange(0, 501, 4)
            else:
                frames = np.arange(0, 601, 4)
                time_at_frames1 = np.arange(0, 601, 4)

            frames_folder1 = frames_folder_general + 'frames_voda/'
            frame_scale1 = []
            for time in time_at_frames1:
                time = time / 20.0
                frame_scale1.append(time)

            if state != 'Pr':
                labelsx = np.arange(0, 26, 5)
                string_labels_x = []
                for i in np.arange(0, 26, 1):
                    if i not in labelsx:
                        string_labels_x.append('')
                    else:
                        string_labels_x.append(i)
            else:
                labelsx = np.arange(0, 31, 2)
                string_labels_x = []
                for i in np.arange(0, 31, 2):
                    if i not in labelsx:
                        string_labels_x.append('')
                    else:
                        string_labels_x.append(i)


            # labelsy1 = np.arange(-10, 21, 2)
            # string_labels_y1 = []
            # for i in np.arange(-10, 21, 1):
            #     if i not in labelsy1:
            #         string_labels_y1.append('')
            #     else:
            #         string_labels_y1.append(i)

            labelsy1 = np.arange(4, 21, 2)
            string_labels_y1 = []
            for i in np.arange(4, 21, 1):
                if i not in labelsy1:
                    string_labels_y1.append('')
                else:
                    string_labels_y1.append(i)

            labelsy2 = np.arange(2, 9, 1)
            string_labels_y2 = []
            for i in np.arange(2, 9, 1):
                if i not in labelsy2:
                    string_labels_y2.append('')
                else:
                    string_labels_y2.append(i)

            pka_results_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_pra_a3/10ang_cavities_09_voda/md_%s/' % (state, prot)

            pkas = get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames, 'frame')

            name = 'md_%s_%s' % (state, prot)

            sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder1, frames, 'PRA-1_EHEM', 'ASP-407_ACHA', name, out_folder_general)

            # sb_dist = sb_dist[2:]

            frames_plot_two_properties_pub('md_%s_%s' % (state, prot), 'pka', 'sb_dist', pkas, sb_dist, frame_scale1,
                                                                 out_folder_general, x_ticks, x_lim, y1_lim, y1_ticks, y2_lim,
                                                                 y2_ticks, string_labels_x, string_labels_y1, string_labels_y2)