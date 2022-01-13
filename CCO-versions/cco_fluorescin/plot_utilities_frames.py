# coding=utf-8

import os
import re

import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.rc({'font': "Times New Roman"})

import numpy as np
import kbp2
from workspace_jd import tools
from kbp2.workspace_jd.cco_fluorescin import interaction_ener


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


def get_sasa_trajectory_residue(result_folder, residue, frames):
    '''
    :param result_folder:
    :param residue:
    :param frames:
    :return:
    '''
    results_sasa = []
    for frame in frames:
        sasa = kbp2.charmm.get_charmm_sasa(result_folder + '/%i/' % frame)
        results_sasa.append(sasa)
    return results_sasa

def get_charmm_interation_energy_trajectory(result_folder, frames):
    '''
    :param result_folder:
    :param frames:
    :return:
    '''
    results = []
    for frame in frames:
            energies = kbp2.charmm.get_charmm_energy(result_folder + '/%i/' % frame, inter=True)
            results.append(energies['ELEC'][0])
    return results


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
            elif property == 'sasa':
                results_dict[residue][trajectory] = get_sasa_trajectory_residue(result_folder, residue, frames)
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

def plot_binding_energy(out_folder, out_file, title, result_folder, frame_list, time_scale, ordered_components, y_lim_set=None):
    '''

    :param out_folder:
    :param out_file:
    :param title:
    :param result_folder:
    :param frame_list:
    :param time_scale:
    :param ordered_components:
    :return:
    '''

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    if out_folder[-1] != '/':
        out_folder = out_folder+'/'


    if result_folder[-1] != '/':
        result_folder = result_folder+'/'

    results_ene = []
    for frame in frame_list:
        bind_ene = interaction_ener.get_binding_energy(result_folder+'%i/'%frame, ordered_components)
        results_ene.append(bind_ene)

    x=time_scale
    y = np.array(results_ene, dtype = 'float_')
    y_mean = [np.mean(y) for i in x]
    plt.figure()
    plt.plot(x,y, marker='o', color='black')
    plt.plot(x,y_mean,  label='Avg. %.2f' % y_mean[0], linestyle=':', color='black')
    plt.title(title)
    if y_lim_set is not None:
        plt.ylim(y_lim_set[0], y_lim_set[1])
    plt.xlabel('t(ns)')
    plt.ylabel('Binding energy')
    plt.legend(loc='upper right')
    plt.savefig(out_folder + '%s_elec_charmm.png' % out_file)
    # plt.show

    return y_mean[0]


def plot_total_energy(out_folder, out_file, title, result_folder, frame_list, time_scale, ordered_components, y_lim_set=None, correction=None):
    '''

    :param out_folder:
    :param out_file:
    :param title:
    :param result_folder:
    :param frame_list:
    :param time_scale:
    :param ordered_components:
    :return:
    '''

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    if out_folder[-1] != '/':
        out_folder = out_folder+'/'


    if result_folder[-1] != '/':
        result_folder = result_folder+'/'

    results_ene = []
    for frame in frame_list:
        bind_ene = interaction_ener.get_binding_energy(result_folder+'%i/'%frame, ordered_components, KJ=False)
        if correction is not None:
            bind_ene = bind_ene - correction
        results_ene.append(bind_ene)

    x=time_scale
    y = np.array(results_ene, dtype = 'float_')
    y_mean = [np.mean(y) for i in x]

    plt.figure()
    plt.plot(x,y, marker='o', color='black')
    plt.plot(x,y_mean,  label='Avg. %.2f' % y_mean[0], linestyle=':', color='black')
    plt.title(title)
    if y_lim_set is not None:
        plt.ylim(y_lim_set[0], y_lim_set[1])
    plt.xlabel('t(ns)')
    plt.ylabel('Total energy (kcal/mol)')
    plt.legend(loc='upper right')
    plt.savefig(out_folder + '%s_elec_charmm.png' % out_file)

    return y_mean[0]


def plot_energy_frames(out_folder, out_file, title, result_folder, frame_list, time_scale, method=None):
    '''
    :param out_folder:
    :param out_file:
    :param title:
    :param result_folder:
    :param frame_list:
    :param time_scale:
    :param method:
    :return:
    '''

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    if out_folder[-1] != '/':
        out_folder = out_folder+'/'

    results_ene = []
    for frame in frame_list:
        if method == 'apbs':
            result_file = result_folder + '/%i/interaction_energy.dat' % frame
            for line in open(result_file, 'r'):
                if 'KJ/mol' in line:
                    line = line.split()
                    energy = line[0]
                    results_ene.append(energy)
        elif method == 'charm_inter_elec':
            # energies = kbp2.charmm.get_charmm_energy(result_folder + '/%i/' % frame, inter=True)
            energies = kbp2.charmm.get_charmm_energy(result_folder + '/%i/' % frame, inter=True, KJ=False)
            results_ene.append(energies['ELEC'][0])


    x=time_scale
    y = np.array(results_ene, dtype = 'float_')
    y_mean = [np.mean(y) for i in x]
    plt.figure()
    plt.plot(x,y, marker='o', color='black')
    plt.plot(x,y_mean,  label='Avg. %.2f' % y_mean[0], linestyle=':', color='black')
    plt.title(title)
    plt.xlabel('t(ns)')
    plt.ylabel('Elec. Interaction Energy (kcal/mol)')
    # plt.ylabel('Elec. Interaction Energy (KJ/mol)')
    plt.legend(loc='upper right')
    plt.ylim(-20,5,1)
    plt.savefig(out_folder + '%s_elec_charmm.png' % out_file)

def plot_sasa_frames(out_folder, out_file, title, result_folder, frame_list, time_scale, x_ticks=None, x_lim=None,
                     y_lim=None, y_ticks=None):

    '''
    :param out_folder:
    :param out_file:
    :param title:
    :param result_folder:
    :param frame_list:
    :param time_scale:
    :return:
    '''

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    if out_folder[-1] != '/':
        out_folder = out_folder+'/'

    plt.rcParams.update({'font.size': 15})

    results_ene = []
    for frame in frame_list:
        sasa = kbp2.charmm.get_charmm_sasa(result_folder + '/%i/' % frame)
        results_ene.append(sasa)

    x=time_scale
    y = np.array(results_ene, dtype = 'float_')
    y_mean = [np.mean(y) for i in x]

    plt.figure()
    # plt.plot(x,y, marker='o', color='black')
    plt.plot(x,y, color='black', lw=1.5)
    plt.xlim()
    # plt.plot(x,y_mean,  label='Avg. %.2f' % y_mean[0], linestyle=':', color='black')
    plt.title(title)
    # plt.xlabel('t(ns)')
    # plt.ylabel('SASA')
    # plt.legend(loc='upper right')

    plt.xlim(x_lim[0], x_lim[1])
    plt.ylim(y_lim[0], y_lim[1])
    plt.xticks(x_ticks)
    plt.yticks(y_ticks)

    plt.savefig(out_folder + '%s_sasa.png' % out_file)

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


    ax1.plot(time_scale, prop1_data, color='black', lw=0.75)
    # ax2.plot(time_scale, prop2_data, linestyle=':', color='black', marker='.', lw=1.0, markersize=1)
    ax2.plot(time_scale, prop2_data, color='green', lw=0.75)

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
    # plt.show()HoxG_Ntag_NiSI_3RGW_small_rigid_longer_equi


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
HoxG_Ntag_NiSI_3RGW_small_rigid_longer_equi
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