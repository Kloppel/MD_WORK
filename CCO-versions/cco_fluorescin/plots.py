# coding=utf-8

import kbp2
from kbp2.workspace_jd.cco_fluorescin import plot_utilities_frames
import numpy as np
import matplotlib.pyplot as plt
import pickle


if __name__ == '__main__':


    #########################
    ## plot charmm energy ###
    #########################
    # global_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/charmm_binding_energy/4_80/'
    # global_output = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/charmm_binding_energy/4_80/'
    #
    # # selected_mds_after_interaction_energy = {'hsd73_hsd526_flu-':[180],\
    # #                                      'hsd73_hsd526_red_fluh':[270],\
    # #                                      'hsd73_hse526_flu-':[180],\
    # #                                      'hsd73_hsp526_flu-':[270, 90],\
    # #                                      'hsd73_hsp526_red_flu2':[270],\
    # #                                      'hsd73_hsp526_red_flu-':[0]}
    #
    # # selected_mds_after_interaction_energy = {'hsd73_hsp526_red_flu2':[270]}
    # # selected_mds_after_interaction_energy = {'hsd73_hse526_flu-':[180]}
    # # selected_mds_after_interaction_energy = {'hsd73_hsd526_red_flu-':[90]}
    # selected_mds_after_interaction_energy = {'hsd73_hse526_flu2':[0,90,180,270]}
    # mds_selection_reorganized = {}
    # mds_selection_reorganized[0] = []
    # mds_selection_reorganized[90] = []
    # mds_selection_reorganized[180] = []
    # mds_selection_reorganized[270] = []
    # for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
    #     for position in position_list:
    #         if pattern not in mds_selection_reorganized[position]:
    #             mds_selection_reorganized[position].append(pattern)
    #
    # for position, pattern_list in mds_selection_reorganized.iteritems():
    #     for pattern in pattern_list:
    # global_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/interaction_energy_charmm/charmm/'
    # global_output = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/interaction_energy_charmm/charmm/results/'
    # positions = [0,90,180,270]
    # for position in positions:
    #     result_folder = global_workfolder + str(position) + '/'
    #     subfolders = os.listdir(result_folder)
    #     for pattern in subfolders:
    #         if '526' in pattern:
    #             out_folder = global_output + str(position) + '/'
    #             # out_folder = global_output + str(position) + '/' + pattern + '/'
    #             result_folder = global_workfolder + str(position) + '/' + pattern
    #             # frames = np.arange(44,205,2)
    #             frames = np.arange(24,105,2)
    #             new_pos, new_pat = pattern_conversion_flu(position, pattern)
    #             title = new_pat
    #             out_file = title + new_pos + ' * ' +' (%s)' % (pattern)
    #             # time_at_frames = np.arange(0,161,2)
    #             time_at_frames = np.arange(0,81,2)
    #             frame_scale = []
    #             for time in time_at_frames:
    #                 time = time/20.0
    #                 frame_scale.append(time)
    #             # residue_selection =  ['HSP-526_ACHA', 'FLU-1_FLUR']
    #             # simple_resi_list = convert_residues(residue_selection)
    #             frames_plot.plot_energy_frames(out_folder, out_file, title, result_folder, frames, frame_scale, method='charm_inter_elec')
    #     #         ordered_components = ['total', 'flu', 'protein']
    #     #         frames_plot.plot_binding_energy(out_folder, out_file, title, result_folder, frames, frame_scale, ordered_components)
# ######################################################################################################################################################################################################################################
# ######################################################################################################################################################################################################################################
#

    # #############################################
    # ### TOTAL ENERGY CHARMM WITH GBVM AND PBC ###
    # #############################################
    #
    # import os
    # # global_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/energy_calculations/charmm_total_energy_notermini/'
    # global_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/energy_calculations/new_charmm_total_notermini/'
    # global_output = '/user/jdragelj/projects/cco_fluorescein/total_energy/'
    # positions = [90]
    # # frames = np.arange(0,445,2)
    # frames = np.arange(44,105,2)
    #
    # dict_data = {}
    # for posi in positions:
    #     result_folder = global_workfolder + str(posi) + '/'
    #     subfolders = os.listdir(result_folder)
    #     dict_data[posi] = {}
    #
    #     if posi == 90:
    #         subfolders = ['hsd73_hse526_red_flu-','hsd73_hse526_flu-']
    #     # if posi == 270:
    #     #     subfolders = ['hsd73_hsd526_flu-']
    #     # if posi == 90:
    #     #     subfolders = ['hsd73_hse526_flu-']
    #
    #     for subf in subfolders:
    #         if '526' in subf:
    #             out_folder = global_output + str(posi) + '/'
    #             result_folder = global_workfolder + str(posi) + '/' + subf + '/'
    #             title = subf
    #             out_file = title + str(posi) + ' * ' +' (%s)' % (subf)
    #             # time_at_frames = np.arange(0,445,2)
    #             time_at_frames = np.arange(44,105,2)
    #             frame_scale = []
    #             for time in time_at_frames:
    #                 time = time/20.0
    #                 frame_scale.append(time)
    #             # ordered_components = ['total', 'flu', 'protein']
    #             ordered_components = ['total']
    #             y_lim_set = [-350,350]
    #             # y_lim_set = None
    #             correction = -12200
    #             # correction = None
    #             mean_ene = plot_utilities_frames.plot_total_energy(out_folder, out_file, title, result_folder, frames, frame_scale, ordered_components, y_lim_set=y_lim_set, correction=correction)
    #             # mean_ene = frames_plot.plot_binding_energy(out_folder, out_file, title, result_folder, frames, frame_scale, ordered_components)
    #             dict_data[posi][subf] = []
    #             dict_data[posi][subf].append(float(mean_ene))

    ### comma file writeout
    # if len(positions) == 4:
    #     csv_filename_o  = global_output + 'O_energy.csv'
    #     comma_fileo = open(csv_filename_o, 'w')
    #     comma_fileo.write('\n')
    #     csv_filename_r  = global_output + 'R_energy.csv'
    #     comma_filer = open(csv_filename_r, 'w')
    #     comma_filer.write('\n')
    #     for subf in subfolders:
    #         if '526' in subf:
    #             where = check_protons(subf)
    #             name = str(where) + '_' + subf[6:]
    #             if 'red' in subf:
    #                 comma_filer.write(tools.string_creation([name] + dict_data[0][subf] + dict_data[90][subf] + dict_data[180][subf]+ dict_data[270][subf]))
    #                 comma_filer.write('\n')
    #             else:
    #                 comma_fileo.write(tools.string_creation([name] + dict_data[0][subf] + dict_data[90][subf] + dict_data[180][subf]+ dict_data[270][subf]))
    #                 comma_fileo.write('\n')


    ### todo:recheck what this is again!
    # # file = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/180/hsd73_hse526_flu-/md_ex_10/output/frames_water/out.out'
    # # file = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/270/hsd73_hsp526_red_flu2/md_ex_15/output/frames_water/out.out'
    # # file = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsd73_hsd526_red_flu-/md/output/frames_water/out.out'
    # # file = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/0/hsd73_hse526_flu2/md/output/frames_water/out.out'
    # f = open(file,'r')
    # energies = []
    # energies2 = []
    # i = 0
    # j = 0
    # for line in f:
    #     if i > 405:
    #         continue
    #     if j > 405:
    #         continue
    #     else:
    #         if 'INTE EXTERN>' in line:
    #             data = line.split()
    #             if i==j:
    #                 if i >= 44:
    #                     energies2.append(float(data[2]))
    #                 j += 2
    #             i += 1
    #
    #
    # frames = np.arange(44,405,2)
    # frames_time = np.arange(0,361,2)
    # time_scale = []
    # for time in frames_time:
    #     time = time/20.0
    #     time_scale.append(time)
    # plt.figure()
    # plt.plot(time_scale,energies2, marker='.', color='black')
    # x = time_scale
    # y = energies2
    # y_mean = [np.mean(y) for i in x]
    # plt.plot(x,y_mean,  label='Avg. %.2f' % y_mean[0], linestyle=':', color='k')
    # plt.legend(loc='upper left',prop={'size':10})
    # plt.ylim(-45,-20)
    # plt.xlabel('t(ns)')
    # plt.ylabel('Interaction energy FLU / protein-membrane-water')
    # plt.savefig('/user/jdragelj/Desktop/180.png')
    # plt.show()



######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
###########################################################################################################################################################
############################################ PKA ##########################################################################################################
###########################################################################################################################################################


 #    #################
 #    ## plots pkas ###
 #    #################
 #    folder_titrations = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/'
 #
 #    out_folder_general = '/user/jdragelj/projects/cco_fluorescein/titration_plots/discussion/with_flu/'
 #    # out_folder_general = '/user/jdragelj/projects/cco_fluorescein/titration_plots/publication/with_flu/'
 #
 #    print('Newst_trajcs_July2017')
 #    selected_mds_after_interaction_energy = {'hsp73_hsp526_red_flu-_2': [90], \
 #                                             'hsp73_hsp526_red_fluh': [90], \
 #                                             'hsp73_hsp526_red_flu2': [90], \
 # \
 #                                             'hsp73_hsp526_fluh': [90], \
 #                                             'hsp73_hsp526_flu2': [90], \
 #                                             'hsp73_hsp526_flu-_2': [90], \
 # \
 #                                             'hsp73_hse526_flu-': [90], \
 #                                             'hsp73_hse526_fluh': [90], \
 #                                             'hsp73_hse526_flu2': [90], \
 # \
 #                                             'hsp73_hse526_red_flu-': [90], \
 #                                             'hsp73_hse526_red_fluh': [90], \
 #                                             'hsp73_hse526_red_flu2': [90], \
 # \
 #                                             # 'hsd73_hse526_red_flu-':[90], \
 #                                             # 'hsd73_hse526_flu-':[90],\
 #
 #                                             }
 #
 #    selected_mds_after_interaction_energy = {'hsd73_hse526_red_flu-':[90], \
 #                                             'hsd73_hse526_flu-':[90],\
 #                                             }
 #
 #
 #
 #    mds_selection_reorganized = {}
 #    mds_selection_reorganized[0] = []
 #    mds_selection_reorganized[90] = []
 #    mds_selection_reorganized[180] = []
 #    mds_selection_reorganized[270] = []
 #    for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
 #        for position in position_list:
 #            if pattern not in mds_selection_reorganized[position]:
 #                mds_selection_reorganized[position].append(pattern)
 #
 #    for position, pattern_list in mds_selection_reorganized.iteritems():
 #        for pattern in pattern_list:
 #            out_folder = out_folder_general + str(position) + '/'
 #            result_folder = folder_titrations + str(position) + '/' + pattern + '/'
 #
 #            frames = np.arange(0,645,2)
 #            time_at_frames = np.arange(0,645,2)
 #
 #
 #            title = pattern
 #            out_file = pattern
 #
 #
 #            frame_scale = []
 #            for time in time_at_frames:
 #                time = time/20.0
 #                frame_scale.append(time)
 #
 #            residue_selection =  ['HSP-526_ACHA', 'FLU-1_FLUR', 'HSP-73_BCHA']
 #            simple_resi_list = residue_selection
 #
 #            print "----New----"
 #            print position, pattern
 #            frame_folder_name = 'frame'
 #
 #            # y_lim = [0,17]
 #            # y_ticks = np.arange(0.0, 18.0, 1.0)
 #            # x_lim = [0,20]
 #            # x_ticks = np.arange(0.0, 21.0, 5.0)
 #            # plot_utilities_frames.plot_pkas_frames_pub(result_folder, out_folder, out_file, title, frames, frame_scale,
 #            #                                            frame_folder_name,x_ticks, x_lim, y_lim, y_ticks, residue_selection,
 #            #                                            simple_resi_list)
 #            plot_utilities_frames.plot_pkas_frames(result_folder, out_folder, out_file, title, frames, frame_scale,
 #                                                   frame_folder_name, residue_selection, simple_resi_list)
 #
 #            # pattern_list = []
 #            # pattern_list.append(pattern)
 #            # source_folder = folder_titrations + str(position) + '/'
 #            # if not os.path.exists(out_folder + pattern + '/'):
 #            #     os.mkdir(out_folder + pattern + '/')
 #            # plot_utilities_frames.make_csv_many_trajectories('pka', source_folder, out_folder + pattern + '/', pattern_list, frames, frame_folder_name, residue_selection)


    # #########################
    # ### plots pkas mutant ###
    # #########################
    #
    # folder_titrations = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/mutant_no_flu/'
    # out_folder_general = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/titration_plots/publication/mutant/'
    # # out_folder_general = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/titration_plots/histidines_noflu/'
    # patterns = os.listdir(folder_titrations)
    #
    # frames = np.arange(2,103,1)
    # time_at_frames = np.arange(0,101,1)
    #
    # for pattern in patterns:
    #     out_folder = out_folder_general + '/'
    #     result_folder = folder_titrations + pattern + '/'
    #     title = pattern
    #     out_file = title
    #
    #     frame_scale = []
    #     for time in time_at_frames:
    #         time = time/10.0
    #         frame_scale.append(time)
    #
    #     residue_selection =  ['HSP-526_ACHA', 'HSP-73_BCHA']
    #     simple_resi_list = convert_residues(residue_selection)
    #     print pattern
    #     frame_folder_name = 'frame'
    #
    #     y_lim = [-5,15]
    #     y_ticks = np.arange(-5.0, 15.5, 1.0)
    #
    #     x_lim = [0,10]
    #     x_ticks = np.arange(0.0, 11.0, 2.0)
    #
    #
    #
    #     plot_utilities_frames.plot_pkas_frames_pub(result_folder, out_folder, out_file, title, frames, frame_scale,
    #                                        frame_folder_name,x_ticks, x_lim, y_lim, y_ticks, residue_selection,
    #                                        simple_resi_list)

        # plot_utilities_frames.plot_pkas_frames_pub(result_folder, out_folder, out_file, title, frames, frame_scale, frame_folder_name, residue_selection, simple_resi_list)

        # plot_utilities_frames.plot_pkas_frames(result_folder, out_folder, out_file, title, frames, frame_scale, frame_folder_name, residue_selection, simple_resi_list)
        # pattern_list = []
        # pattern_list.append(pattern)
        # source_folder = folder_titrations + '/'
        # if not os.path.exists(out_folder + pattern + '/'):
        #     os.mkdir(out_folder + pattern + '/')
        # plot_utilities_frames.make_csv_many_trajectories('pka', source_folder, out_folder + pattern + '/', pattern_list, frames, frame_folder_name, residue_selection)


    # ############################
    # ### plots pkas wild type ###
    # ############################
    # import os
    #
    # folder_titrations = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/'
    # out_folder_general = '/user/jdragelj/projects/cco_fluorescein/titration_plots/discussion/wild_type/'
    #
    # patterns = os.listdir(folder_titrations)
    #
    # frames = np.arange(0,405,2)
    # time_at_frames = np.arange(0,405,2)
    #
    # for pattern in patterns:
    #     out_folder = out_folder_general + '/'
    #     result_folder = folder_titrations + pattern + '/'
    #     title = pattern
    #     out_file = title
    #
    #     frame_scale = []
    #     for time in time_at_frames:
    #         time = time/20.0
    #         frame_scale.append(time)
    #     residue_selection =  ['HSP-526_ACHA', 'HSP-73_BCHA']
    #     simple_resi_list = residue_selection
    #     print pattern
    #     frame_folder_name = 'frame'
    #
    #     y_lim = [0,14]
    #     y_ticks = np.arange(0.0, 14.5, 1.0)
    #
    #     x_lim = [0,20]
    #     x_ticks = np.arange(0.0, 21.0, 2.0)
    #
    #
    #     plot_utilities_frames.plot_pkas_frames_pub(result_folder, out_folder, out_file, title, frames, frame_scale,
    #                                                frame_folder_name,x_ticks, x_lim, y_lim, y_ticks, residue_selection,
    #                                                simple_resi_list)
    #
    #
    #     plot_utilities_frames.plot_pkas_frames(result_folder, out_folder, out_file, title, frames, frame_scale, frame_folder_name, residue_selection, simple_resi_list)
    #
    #     pattern_list = []
    #     pattern_list.append(pattern)
    #     source_folder = folder_titrations + '/'
    #     if not os.path.exists(out_folder + pattern + '/'):
    #         os.mkdir(out_folder + pattern + '/')
    #     plot_utilities_frames.make_csv_many_trajectories('pka', source_folder, out_folder + pattern + '/', pattern_list, frames, frame_folder_name, residue_selection)

######################################################################################################################################################################################################################################
    ######################
    ### plot pkas mix #### publication plots
    ######################

    outfolder = '/user/jdragelj/projects/cco_fluorescein/titration_plots/publication/with_flu/90/'

    title1 = 'hsp73_hsp526_flu-_2'
    frames1 = np.arange(4, 645, 4)
    time_at_frames1 = np.arange(0, 641, 4)
    frame_scale1 = []
    for time in time_at_frames1:
        time = time / 20.0
        frame_scale1.append(time)

    pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/' + title1 + '/'
    pkas_flu1 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'FLU-1_FLUR', frames1, 'frame')
    pkas_731 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames1, 'frame')
    pkas_5261 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames1, 'frame')

    title5 = 'hsd73_hse526_flu-'
    frames5 = np.arange(4, 645, 4)
    time_at_frames5 = np.arange(0, 641, 4)
    frame_scale5 = []
    for time in time_at_frames5:
        time = time / 20.0
        frame_scale5.append(time)

    pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/' + title5 + '/'
    pkas_flu5 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'FLU-1_FLUR', frames5, 'frame')
    pkas_735 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames5, 'frame')
    pkas_5265 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames5, 'frame')

    title6 = 'hsp73_hsp526_fluh'
    frames6 = np.arange(4, 645,4)
    time_at_frames6 = np.arange(0, 641,4)
    frame_scale6 = []
    for time in time_at_frames6:
        time = time / 20.0
        frame_scale6.append(time)

    pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/' + title6 + '/'
    pkas_flu6 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'FLU-1_FLUR', frames6, 'frame')
    pkas_736 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames6, 'frame')
    pkas_5266 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames6, 'frame')

    title7 = 'hsp73_hse526_flu-'
    frames7 = np.arange(4, 645,4)
    time_at_frames7 = np.arange(0, 641, 4)
    frame_scale7 = []
    for time in time_at_frames7:
        time = time / 20.0
        frame_scale7.append(time)

    pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/' + title7 + '/'
    pkas_flu7 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'FLU-1_FLUR', frames7, 'frame')
    pkas_737 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames7, 'frame')
    pkas_5267 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames7, 'frame')

    title2 = 'hsp73_hsp526_red_fluh'
    frames2 = np.arange(4, 645,4)
    time_at_frames2 = np.arange(0, 641, 4)
    frame_scale2 = []
    for time in time_at_frames2:
        time = time / 20.0
        frame_scale2.append(time)

    pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/' + title2 + '/'
    pkas_flu2 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'FLU-1_FLUR', frames2, 'frame')
    pkas_732 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames2, 'frame')
    pkas_5262 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames2, 'frame')

    title3 = 'hsp73_hsp526_red_flu2'
    frames3 = np.arange(4, 645, 4)
    time_at_frames3 = np.arange(0, 641, 4)
    frame_scale3 = []
    for time in time_at_frames3:
        time = time / 20.0
        frame_scale3.append(time)

    pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/' + title3 + '/'
    pkas_flu3 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'FLU-1_FLUR', frames3, 'frame')
    pkas_733 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames3, 'frame')
    pkas_5263 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames3, 'frame')

    title4 = 'hsd73_hse526_red_flu-'
    frames4 = np.arange(4, 645, 4)
    time_at_frames4 = np.arange(0, 641, 4)
    frame_scale4 = []
    for time in time_at_frames4:
        time = time / 20.0
        frame_scale4.append(time)

    pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/' + title4 + '/'
    pkas_flu4 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'FLU-1_FLUR', frames4, 'frame')
    pkas_734 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames4, 'frame')
    pkas_5264 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames4, 'frame')

    for residue in ['HSP-526_ACHA', 'HSP-73_BCHA', 'FLU-1_FLUR']:

        x_lim = [0, 32]
        x_ticks = np.arange(0, 33, 1)
        labelsx = np.arange(0, 33, 4)
        string_labels_x = []
        for i in np.arange(0, 33, 1):
            if i not in labelsx:
                string_labels_x.append('')
            else:
                string_labels_x.append(i)

        if residue == 'HSP-526_ACHA':
            y1 = pkas_5261
            y2 = pkas_5262
            y3 = pkas_5263
            y4 = pkas_5264
            y5 = pkas_5265
            y6 = pkas_5266
            y7 = pkas_5267

            y1_lim = [-1, 11]
            y1_ticks = np.arange(-1, 12, 1)
            labelsy = np.arange(-1, 12, 2)
            string_labels_y = []
            for i in np.arange(-1, 12, 1):
                if i not in labelsy:
                    string_labels_y.append('')
                else:
                    string_labels_y.append(i)

        if residue == 'HSP-73_BCHA':
            y1 = pkas_731
            y2 = pkas_732
            y3 = pkas_733
            y4 = pkas_734
            y5 = pkas_735
            y6 = pkas_736
            y7 = pkas_737

            y1_lim = [-2, 16]
            y1_ticks = np.arange(-2, 17, 1)
            labelsy = np.arange(-2, 17, 3)
            string_labels_y = []
            for i in np.arange(-2, 17, 1):
                if i not in labelsy:
                    string_labels_y.append('')
                else:
                    string_labels_y.append(i)

        if residue == 'FLU-1_FLUR':
            y1 = pkas_flu1
            y2 = pkas_flu2
            y3 = pkas_flu3
            y4 = pkas_flu4
            y5 = pkas_flu5
            y6 = pkas_flu6
            y7 = pkas_flu7

            y1_lim = [3, 9]
            y1_ticks = np.arange(3, 10, 1)
            labelsy = np.arange(3, 10, 1)
            string_labels_y = []
            for i in np.arange(3, 10, 1):
                if i not in labelsy:
                    string_labels_y.append('')
                else:
                    string_labels_y.append(i)


        plt.rcParams.update({'font.size': 20})
        fig, ax1 = plt.subplots()

        ax1.plot(frame_scale1, y1, color='black', lw=2.0)
        ax1.plot(frame_scale5, y5, color='green', lw=1.5)
        ax1.plot(frame_scale6, y6, color='blue', lw=1.5)
        ax1.plot(frame_scale7, y7, color='magenta', lw=1.5)

        # ax1.plot(frame_scale1, y1, color='black', lw=1.5, linestyle='-')
        # ax1.plot(frame_scale5, y5, color='gray', lw=1.5, linestyle='-')
        # ax1.plot(frame_scale6, y6, color='black', lw=1.5, linestyle='--')
        # ax1.plot(frame_scale7, y7, color='black', lw=1.5, linestyle='-.')

        ax1.set_ylim(y1_lim[0], y1_lim[1])
        ax1.set_yticks(y1_ticks)
        ax1.set_xticklabels(string_labels_x)
        ax1.set_yticklabels(string_labels_y)
        ax1.tick_params(direction='out', top='off', right='off')
        plt.xlim(x_lim[0], x_lim[1])
        plt.xticks(x_ticks)

        plt.savefig(outfolder + '%s_oxi_dominant.png' % residue)
        # plt.show()

        fig, ax2 = plt.subplots()

        ax2.plot(frame_scale2, y2, color='blue', lw=1.5)
        ax2.plot(frame_scale3, y3, color='red', lw=2.0)
        ax2.plot(frame_scale4, y4, color='green', lw=1.5)

        # ax2.plot(frame_scale2, y2, color='black', lw=1.5, linestyle='--')
        # ax2.plot(frame_scale3, y3, color='black', lw=1.5, linestyle=':')
        # ax2.plot(frame_scale4, y4, color='black', lw=1.5, linestyle='-')

        ax2.set_ylim(y1_lim[0], y1_lim[1])
        ax2.set_yticks(y1_ticks)
        ax2.set_xticklabels(string_labels_x)
        ax2.set_yticklabels(string_labels_y)
        ax2.tick_params(direction='out', top='off', right='off')
        plt.xlim(x_lim[0], x_lim[1])
        plt.xticks(x_ticks)

        plt.savefig(outfolder + '%s_red_dominant.png' % residue)
        # plt.show()

        if residue == 'FLU-1_FLUR':
            y1_lim = [3, 7]
            y1_ticks = np.arange(3, 8, 1)
            labelsy = np.arange(3, 8, 1)
            string_labels_y = []
            for i in np.arange(3, 8, 1):
                if i not in labelsy:
                    string_labels_y.append('')
                else:
                    string_labels_y.append(i)

        if residue == 'HSP-73_BCHA':
            y1_lim = [2, 16]
            y1_ticks = np.arange(2, 17, 1)
            labelsy = np.arange(2, 17, 2)
            string_labels_y = []
            for i in np.arange(2, 17, 1):
                if i not in labelsy:
                    string_labels_y.append('')
                else:
                    string_labels_y.append(i)

        if residue == 'HSP-526_ACHA':
            y1_lim = [4, 10]
            y1_ticks = np.arange(4, 11, 1)
            labelsy = np.arange(4, 11, 1)
            string_labels_y = []
            for i in np.arange(4, 11, 1):
                if i not in labelsy:
                    string_labels_y.append('')
                else:
                    string_labels_y.append(i)

        fig = plt.figure()
        ax3 = fig.add_subplot(111)

        ax3.plot(frame_scale3, y3, color='black', lw=0.75, linestyle=(0, (3, 2)), marker='.', markersize=1.5)
        ax3.plot(frame_scale1, y1, color='black', lw=0.75)


        ax3.spines['top'].set_linewidth(0.5)
        ax3.spines['left'].set_linewidth(0.5)
        ax3.spines['right'].set_linewidth(0.5)
        ax3.spines['bottom'].set_linewidth(0.5)
        fig.set_size_inches(3.2, 2.7)
        plt.rcParams.update({'font.size': 10})
        # plt.tight_layout()
        ax3.set_ylim(y1_lim[0], y1_lim[1])
        ax3.set_yticks(y1_ticks)
        ax3.set_xticklabels(string_labels_x)
        ax3.set_yticklabels(string_labels_y)
        ax3.tick_params(direction='out', top='off', right='off', length=3.0)
        plt.xlim(x_lim[0], x_lim[1])
        plt.xticks(x_ticks)
        plt.plot([2.6, 2.15])

        plt.savefig(outfolder + '%s_oxi_red_dominant.png' % residue, dpi=300)
        # plt.show()


######################################################################################################################################################################################################################################    #pub#
    # outfolder = '/user/jdragelj/projects/cco_fluorescein/titration_plots/publication/wild_type/'
    #
    # title1 = 'hsp73_hsp526_red'
    # frames1 = np.arange(4, 405, 4)
    # time_at_frames1 = np.arange(0, 401, 4)
    # frame_scale1 = []
    # for time in time_at_frames1:
    #     time = time / 20.0
    #     frame_scale1.append(time)
    #
    # pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/' + title1 + '/'
    # pkas_731 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames1,
    #                                                              'frame')
    # pkas_5261 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames1,
    #                                                               'frame')
    #
    # title2 = 'hsp73_hsp526'
    # frames2 = np.arange(4, 405, 4)
    # time_at_frames2 = np.arange(0, 401, 4)
    # frame_scale2 = []
    # for time in time_at_frames2:
    #     time = time / 20.0
    #     frame_scale2.append(time)
    #
    # pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/' + title2 + '/'
    # pkas_732 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames2,
    #                                                              'frame')
    # pkas_5262 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames2,
    #                                                               'frame')
    #
    # title3 = 'hsd73_hsp526_red'
    # frames3 = np.arange(4, 405, 4)
    # time_at_frames3 = np.arange(0, 401, 4)
    # frame_scale3 = []
    # for time in time_at_frames3:
    #     time = time / 20.0
    #     frame_scale3.append(time)
    #
    # pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/' + title3 + '/'
    # pkas_733 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames3,
    #                                                              'frame')
    # pkas_5263 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames3,
    #                                                               'frame')
    #
    # title4 = 'hsp73_hse526'
    # frames4 = np.arange(4, 405, 4)
    # time_at_frames4= np.arange(0, 401, 4)
    # frame_scale4 = []
    # for time in time_at_frames4:
    #     time = time / 20.0
    #     frame_scale4.append(time)
    #
    # pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/' + title4 + '/'
    # pkas_734 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames4,
    #                                                              'frame')
    # pkas_5264 = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-526_ACHA', frames4,
    #                                                               'frame')
    #
    #
    #
    # for residue in ['HSP-526_ACHA', 'HSP-73_BCHA']:
    #
    #     x_lim = [0, 20]
    #     x_ticks = np.arange(0, 21, 1)
    #     labelsx = np.arange(0, 21, 2)
    #     string_labels_x = []
    #     for i in np.arange(0, 21, 1):
    #         if i not in labelsx:
    #             string_labels_x.append('')
    #         else:
    #             string_labels_x.append(i)
    #
    #     if residue == 'HSP-526_ACHA':
    #         y1 = pkas_5261
    #         y2 = pkas_5262
    #         y3 = pkas_5263
    #         y4 = pkas_5264
    #
    #         y1_lim = [3, 9]
    #         y1_ticks = np.arange(3, 10, 1)
    #         labelsy = np.arange(3, 10, 1)
    #         string_labels_y = []
    #         for i in np.arange(3, 10, 1):
    #             if i not in labelsy:
    #                 string_labels_y.append('')
    #             else:
    #                 string_labels_y.append(i)
    #
    #
    #     if residue == 'HSP-73_BCHA':
    #         y1 = pkas_731
    #         y2 = pkas_732
    #         y3 = pkas_733
    #         y4 = pkas_734
    #
    #         y1_lim = [-2, 14]
    #         y1_ticks = np.arange(-2, 15, 1)
    #         labelsy = np.arange(-2, 15, 2)
    #         string_labels_y = []
    #         for i in np.arange(-2, 15, 1):
    #             if i not in labelsy:
    #                 string_labels_y.append('')
    #             else:
    #                 string_labels_y.append(i)
    #
    #     plt.rcParams.update({'font.size': 20})
    #     fig, ax1 = plt.subplots()
    #
    #     ax1.plot(frame_scale1, y1, color='black', lw=1.5)
    #     ax1.plot(frame_scale2, y2, color='red', lw=1.5)
    #     ax1.plot(frame_scale3, y3, color='cyan', lw=1.5)
    #     ax1.plot(frame_scale4, y4, color='olive', lw=1.5)
    #
    #     ax1.set_ylim(y1_lim[0], y1_lim[1])
    #     ax1.set_yticks(y1_ticks)
    #     ax1.set_xticklabels(string_labels_x)
    #     ax1.set_yticklabels(string_labels_y)
    #     ax1.tick_params(direction='out', top='off', right='off')
    #     plt.xlim(x_lim[0], x_lim[1])
    #     plt.xticks(x_ticks)
    #
    #     plt.savefig(outfolder + '%s_dominant.png' % residue)
    #     # plt.show()
    #
    #     fig, ax3 = plt.subplots()
    #
    #     ax3.plot(frame_scale3, y1, color='black', lw=2.0, linestyle='--', marker='.')
    #     ax3.plot(frame_scale1, y2, color='black', lw=1.5)
    #
    #     ax3.set_ylim(y1_lim[0], y1_lim[1])
    #     ax3.set_yticks(y1_ticks)
    #     ax3.set_xticklabels(string_labels_x)
    #     ax3.set_yticklabels(string_labels_y)
    #     ax3.tick_params(direction='out', top='off', right='off')
    #     plt.xlim(x_lim[0], x_lim[1])
    #     plt.xticks(x_ticks)
    #
    #     plt.savefig(outfolder + '%s_oxi_red_dominant.png' % residue)
    #     # plt.show()
    #
    #
    #
    #

# ######################################################################################################################################################################################################################################
# ######################################################################################################################################################################################################################################

    # ###########################
    # ### protonation curve #####
    # ###########################
    #
    # sourcefolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/'
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/titration_plots/protonation_curves/'
    # # workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/titration_curves_5ns_all_frames/each_frame/'
    # frames = np.arange(44,105,2)
    # frame_folder_name = 'frames_1'
    # positions = [0,90,180,270]
    #
    # # selected_mds_after_interaction_energy = {'hsd73_hsd526_flu-':[180],\
    # #                                      'hsd73_hsd526_red_fluh':[270],\
    # #                                      'hsd73_hsp526_flu-':[270],\
    # #                                      'hsd73_hsp526_red_flu2':[270],\
    # #                                      'hsd73_hse526_flu-':[180],\
    # #                                      'hsd73_hsp526_red_flu-':[0]}
    #
    # # selected_mds_after_interaction_energy = {'hsd73_hse526_red_flu-':[90],\
    # #                                      'hsd73_hsp526_red_flu-':[180],\
    # #                                      'hsd73_hse526_flu-':[90],\
    # #                                      'hsd73_hse526_fluh':[180]}
    #
    # selected_mds_after_interaction_energy = {'hsd73_hse526_red_flu-':[90],\
    #                                      'hsd73_hse526_flu-':[90]}
    #
    # # selected_mds_after_interaction_energy = {'hsd73_hsd526_flu-':[180],\
    # #                                      'hsd73_hsd526_red_fluh':[270],\
    # #                                      'hsd73_hsp526_flu-':[270],\
    # #                                      'hsd73_hsp526_red_flu2':[270],\
    # #                                      'hsd73_hse526_flu-':[180],\
    # #                                      'hsd73_hsp526_red_flu-':[0],\
    # #                                      'hsd73_hsd526_flu2':[90],\
    # #                                      'hsd73_hsd526_fluh':[180], \
    # #                                      'hsd73_hsp526_red_fluh':[180],\
    # #                                      'hsd73_hsp526_red_fluh':[0]}
    #
    # mds_selection_reorganized = {}
    # mds_selection_reorganized[0] = []
    # mds_selection_reorganized[90] = []
    # mds_selection_reorganized[180] = []
    # mds_selection_reorganized[270] = []
    # for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
    #     for position in position_list:
    #         if pattern not in mds_selection_reorganized[position]:
    #             mds_selection_reorganized[position].append(pattern)
    # # residue_list = ['FLU-1_FLUR', 'HSP-526_ACHA']
    # residue_list = ['FLU-1_FLUR']
    # for position, pattern_list in mds_selection_reorganized.iteritems():
    #     for pattern in pattern_list:
    #         print pattern
    #         sub_workfolder = workfolder +  str(position) + '/'
    #         if not os.path.exists(sub_workfolder):
    #             os.mkdir(sub_workfolder)
    #         sub_sub_workfolder = sub_workfolder + '/' + pattern + '/'
    #         sub_folder = sourcefolder + str(position) + '/' + pattern + '/'
    #         if not os.path.exists(sub_sub_workfolder):
    #             os.mkdir(sub_sub_workfolder)
    #         combined_results1 = kbp2.kbp_results.FrameKbpResults()
    #         # frame_range = np.arange(44,105,2)
    #         frame_range = np.arange(44,445,2)
    #         for frame in frame_range:
    #             file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
    #             frame_result = pickle.load( open( file, "rb" ) )
    #             combined_results1.add_task(frame_result.kbp_results[0], frame_nr=frame)
    #         combined_results1.average_over_frames(frame_range=frame_range)
    #         combined_results = kbp2.kbp_results.FrameKbpResults()
    #         combined_results.add_task(combined_results1.avg_results)
    #         combined_results.combine_frames_karlsberg(cpus=3)
    #         for residue_descr in residue_list:
    #             residue_nr = combined_results.descr.get_residue_index(residue_descr)
    #             deprot_curve = combined_results.combined_results.deprot_curves[residue_nr]
    #             ph_values = combined_results.descr.ph_values
    #             plt.figure()
    #             final_pkas = combined_results.combined_results.pkas
    #             pka = final_pkas[residue_descr]
    #             dep_refined = deprot_curve[26:39]
    #             ph_refined = ph_values[26:39]
    #             plt.plot(ph_refined, dep_refined, linestyle = '_', color = 'k')
    #             resname, resid, segname = re.split('[-_]', residue_descr)
    #             if resname == 'HSP':
    #                 resname = 'HIS'
    #             plt.title("Deprotonation curve for " + resname + ' pKa %.2f' % pka)
    #             plt.ylim([-0.01, 1.01])
    #             plt.savefig('%s%s_titr_cur.png' % (sub_sub_workfolder, residue_descr))
    #             plt.xticks(np.arange(3.0, 9.0, 0.5))
    #             plt.close()
    #         plt.figure()
    #         pkas_string = ''
    #         if len(residue_list) > 3:
    #             raise AssertionError("Too many residues for now")
    #         for residue_descr in residue_list:
    #             residue_nr = combined_results.descr.get_residue_index(residue_descr)
    #             deprot_curve = combined_results.combined_results.deprot_curves[residue_nr]
    #             ph_values = combined_results.descr.ph_values
    #             final_pkas = combined_results.combined_results.pkas
    #             resname, resid, segname = re.split('[-_]', residue_descr)
    #             if resname == 'HSP':
    #                 resname = 'HIS'
    #             pka = final_pkas[residue_descr]
    #             pkas_string = ' ' + resname + ' %.2f' % pka
    #             if residue_descr == 'FLU-1_FLUR':
    #                 color = 'r'
    #             elif residue_descr == 'HSP-526_ACHA':
    #                 color = 'b'
    #             dep_refined = deprot_curve[26:39]
    #             ph_refined = ph_values[26:39]
    #             plt.plot(ph_refined, dep_refined, '%sx-' % color, label = pkas_string)
    #             # plt.plot(ph_values, deprot_curve, '%sx-' % color, label = pkas_string)
    #             plt.xticks(np.arange(3.0, 9.0, 0.5))
    #             # plt.xticks(np.arange(min(ph_values), max(ph_values)+1, 1.0))
    #             plt.ylim([-0.01, 1.01])
    #             # plt.show()
    #         plt.xticks(np.arange(3.0, 9.0, 0.5))
    #         new_pos, new_pattern = pattern_conversion_flu(position, pattern)
    #         plt.title(new_pattern)
    #         # plt.title(pattern)
    #         plt.legend(loc = 'upper left')
    #         plt.savefig('%s%s.png' % (sub_sub_workfolder, pattern))


######################################################################################################################################################################################################################################
    ############################
    #### protonation curve many trajectories final ##### publication plots
    ############################

    sourcefolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/'
    workfolder = '/user/jdragelj/projects/cco_fluorescein/titration_plots/publication/'

    frames = np.arange(344,645,2)
    frame_range = [344,645]

    # frame_folder_name = 'frames_1'
    positions = [90]

 #    selected_mds_after_interaction_energy = {'hsp73_hsp526_red_fluh': [90], \
 #                                             'hsp73_hsp526_red_flu2': [90], \
 # \
 #                                             'hsp73_hsp526_fluh': [90], \
 #                                             'hsp73_hsp526_flu-_2': [90], \
 #                                             'hsp73_hse526_flu-': [90], \
 # \
 #                                             'hsd73_hse526_red_flu-': [90], \
 #                                             'hsd73_hse526_flu-': [90], \
 #                                             }


    selected_mds_after_interaction_energy = {'hsp73_hsp526_red_flu2': [90], \
 \

                                             'hsp73_hsp526_flu-_2': [90], \
                                             }

    mds_selection_reorganized = {}
    mds_selection_reorganized[0] = []
    mds_selection_reorganized[90] = []
    mds_selection_reorganized[180] = []
    mds_selection_reorganized[270] = []
    for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
        for position in position_list:
            if pattern not in mds_selection_reorganized[position]:
                mds_selection_reorganized[position].append(pattern)

    residue_list = ['FLU-1_FLUR', 'HSP-73_BCHA', 'HSE-526_ACHA']
    fig, ax3 = plt.subplots()

    types = ['kx-', 'kx--']
    ti = 0
    pka_diff = []
    list_of_comb_results = []
    #oxi
    combined_results_average = kbp2.kbp_results.FrameKbpResults()
    for position, pattern_list in mds_selection_reorganized.iteritems():
        for pattern in pattern_list:
            if 'red' not in pattern:
                sub_workfolder = workfolder +  str(position) + '/'
                sub_sub_workfolder = sub_workfolder + '/' + pattern + '/'
                sub_folder = sourcefolder + str(position) + '/' + pattern + '/'

                combined_results = kbp2.kbp_results.FrameKbpResults()

                for i, frame in enumerate(frames):
                    file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
                    frame_result = pickle.load( open( file, "rb" ) )
                    combined_results.add_task(frame_result.kbp_results[0], frame_nr=frame)
                combined_results.average_over_frames(frame_range)
                combined_results_average.add_task(combined_results.avg_results)
    combined_results_average.combine_frames_karlsberg(cpus=3)
    list_of_comb_results.append(combined_results_average)

    #red
    combined_results_average = kbp2.kbp_results.FrameKbpResults()
    for position, pattern_list in mds_selection_reorganized.iteritems():
        for pattern in pattern_list:
            if 'red' in pattern:
                sub_workfolder = workfolder +  str(position) + '/'
                sub_sub_workfolder = sub_workfolder + '/' + pattern + '/'
                sub_folder = sourcefolder + str(position) + '/' + pattern + '/'
                combined_results = kbp2.kbp_results.FrameKbpResults()
                for i, frame in enumerate(frames):
                    file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
                    frame_result = pickle.load( open( file, "rb" ) )
                    combined_results.add_task(frame_result.kbp_results[0], frame_nr=frame)
                combined_results.average_over_frames(frame_range)
                combined_results_average.add_task(combined_results.avg_results)
    combined_results_average.combine_frames_karlsberg(cpus=3)
    list_of_comb_results.append(combined_results_average)


    print len(list_of_comb_results)
    for i, combined_results_average in enumerate(list_of_comb_results):
        residue_descr = residue_list[0]
        residue_nr = combined_results_average.descr.get_residue_index(residue_descr)
        deprot_curve = combined_results_average.combined_results.deprot_curves[residue_nr]
        ph_values = combined_results_average.descr.ph_values
        final_pkas = combined_results_average.combined_results.pkas
        pka = final_pkas[residue_descr]
        pka = round(pka,2)
        pka_diff.append(pka)

        dep_refined = deprot_curve[26:39]
        ph_refined = ph_values[26:39]

        if i == 0:
            ax3.plot(ph_refined, dep_refined, linestyle='--', color = 'black', lw=1.5)
            print ph_refined
            f = open('/user/jdragelj/projects/cco_fluorescein/titration_plots/publication/titration_data_oxi2.dat', 'w')
        else:
            ax3.plot(ph_refined, dep_refined, color = 'black', lw=1.5)
            f = open('/user/jdragelj/projects/cco_fluorescein/titration_plots/publication/titration_data_red2.dat', 'w')

        for ph, occ in zip(ph_refined, dep_refined):
            to_write = '%.1f;%.8f\n' %(ph, occ)
            f.write(to_write)
        f.close()


        ti +=1

    labelsy = np.arange(0.0, 1.01, 0.2)
    string_labels_y = []
    for i in np.arange(0.0, 1.01, 0.1):
        if i not in labelsy:
            string_labels_y.append('')
        else:
            string_labels_y.append(i)

    plt.rcParams.update({'font.size': 10})
    pk_dif = abs(pka_diff[0]-pka_diff[1])
    print pk_dif
    plt.ylim([-0.00, 1.00])
    time = (len(frames)-1)/10

    plt.xticks(np.arange(3.0, 9.5, 1))
    plt.yticks(np.arange(0.0, 1.1, 0.1))
    plt.xlim(3,9)
    ax3.set_yticklabels(string_labels_y)



    ax3.spines['top'].set_linewidth(0.5)
    ax3.spines['left'].set_linewidth(0.5)
    ax3.spines['right'].set_linewidth(0.5)
    ax3.spines['bottom'].set_linewidth(0.5)
    fig.set_size_inches(3.2, 2.7)
    plt.rcParams.update({'font.size': 10})
    plt.tight_layout()

    ax3.tick_params(direction='out', top='off', right='off')

    plt.savefig('%spka_shift2_%ins.png' % (workfolder,time), dpi=300)
    # plt.show()


# ######################################################################################################################################################################################################################################
# ######################################################################################################################################################################################################################################
#     #########################
#     ## plots salt-bridges ###
#     #########################
#
#     frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/'
#     # out_folder_general = '/user/jdragelj/projects/cco_fluorescein/salt_bridges/fluO5/'
#     # out_folder_general = '/user/jdragelj/projects/cco_fluorescein/salt_bridges/fluO6/'
#     out_folder_general = '/user/jdragelj/projects/cco_fluorescein/salt_bridges/90/'
#
#
#     print('his526 - flu O5')
#     selected_mds_after_interaction_energy = {
#                                              'hsp73_hsp526_red_fluh': [90], \
#                                              'hsp73_hsp526_red_flu2': [90], \
#                                              'hsp73_hsp526_fluh': [90], \
#                                              'hsp73_hsp526_flu-_2': [90], \
#                                              # 'hsp73_hse526_flu-': [90], \
#                                              }
#
#     residue_2 = 'FLU-1_FLUR'
#     residue_1 = 'HSP-526_ACHA'
#
#     mds_selection_reorganized = {}
#     mds_selection_reorganized[0] = []
#     mds_selection_reorganized[90] = []
#     mds_selection_reorganized[180] = []
#     mds_selection_reorganized[270] = []
#     for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
#         for position in position_list:
#             if pattern not in mds_selection_reorganized[position]:
#                 mds_selection_reorganized[position].append(pattern)
#
#     for position, pattern_list in mds_selection_reorganized.iteritems():
#         for pattern in pattern_list:
#
#             if 'flu2' in pattern:
#                 residue_2 = 'FLH-1_FLUR'
#             if 'flu-' in pattern:
#                 residue_2 = 'FLX-1_FLUR'
#             if 'fluh' in pattern:
#                 residue_2 = 'FLU-1_FLUR'
#
#             out_folder = out_folder_general + str(position) + '/'
#             frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_10/frames_1/'
#
#             frames = np.arange(0, 645, 2)
#             time_at_frames = np.arange(0, 645, 2)
#             title = str(position) + pattern
#             out_file = pattern
#
#             frame_scale = []
#             for time in time_at_frames:
#                 time = time / 20.0
#                 frame_scale.append(time)
#
#             print "--------"
#             print position, pattern
#             frame_folder_name = 'frame'
#
#             y_lim = [2, 13]
#             y_ticks = np.arange(2.0, 13.0, 1.0)
#
#             x_lim = [0, 32]
#             x_ticks = np.arange(0.0, 33.0, 2.0)
#
#             salt_bridge_distances.plot_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file,
#                                                                     title, frames, frame_scale, frames_folder,
#                                                                     x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim,
#                                                                     y_ticks=y_ticks)
#
#     print('73 - 78')
#     selected_mds_after_interaction_energy = {'hsp73_hsp526_red_fluh': [90], \
#                                              'hsp73_hsp526_red_flu2': [90], \
#                                              'hsp73_hse526_red_fluh': [90], \
#  \
#                                              'hsp73_hsp526_fluh': [90], \
#                                              'hsp73_hsp526_flu-_2': [90], \
#                                              'hsp73_hse526_flu-': [90], \
#                                              }
#
#     residue_1 = 'HSP-73_BCHA'
#     residue_2 = 'GLU-78_BCHA'
#
#     mds_selection_reorganized = {}
#     mds_selection_reorganized[0] = []
#     mds_selection_reorganized[90] = []
#     mds_selection_reorganized[180] = []
#     mds_selection_reorganized[270] = []
#     for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
#         for position in position_list:
#             if pattern not in mds_selection_reorganized[position]:
#                 mds_selection_reorganized[position].append(pattern)
#
#     for position, pattern_list in mds_selection_reorganized.iteritems():
#         for pattern in pattern_list:
#             out_folder = out_folder_general + str(position) + '/'
#             frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_10/frames_1/'
#
#             frames = np.arange(0, 645, 2)
#             time_at_frames = np.arange(0, 645, 2)
#             title = str(position) + pattern
#             out_file = pattern
#
#             frame_scale = []
#             for time in time_at_frames:
#                 time = time / 20.0
#                 frame_scale.append(time)
#
#             print "--------"
#             print position, pattern
#             frame_folder_name = 'frame'
#
#             y_lim = [2, 13]
#             y_ticks = np.arange(2.0, 13.0, 1.0)
#
#             x_lim = [0, 32]
#             x_ticks = np.arange(0.0, 33.0, 2.0)
#
#             salt_bridge_distances.plot_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file, title, frames, frame_scale, frames_folder,
#                          x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim, y_ticks=y_ticks)
#
#
#     print('526 - 525')
#     selected_mds_after_interaction_energy = {'hsp73_hsp526_red_fluh': [90], \
#                                              'hsp73_hsp526_red_flu2': [90], \
#                                              'hsp73_hse526_red_fluh': [90], \
#  \
#                                              'hsp73_hsp526_fluh': [90], \
#                                              'hsp73_hsp526_flu-_2': [90], \
#                                              'hsp73_hse526_flu-': [90], \
#                                              }
#
#     # selected_mds_after_interaction_energy = {'hsd73_hse526_flu-': [90], \
#     #                                          'hsd73_hse526_red_flu-': [90]}
#     residue_1 = 'HSP-526_ACHA'
#     residue_2 = 'GLU-525_ACHA'
#
#     mds_selection_reorganized = {}
#     mds_selection_reorganized[0] = []
#     mds_selection_reorganized[90] = []
#     mds_selection_reorganized[180] = []
#     mds_selection_reorganized[270] = []
#     for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
#         for position in position_list:
#             if pattern not in mds_selection_reorganized[position]:
#                 mds_selection_reorganized[position].append(pattern)
#
#     for position, pattern_list in mds_selection_reorganized.iteritems():
#         for pattern in pattern_list:
#
#             if 'hsp526' in pattern:
#                 residue_1 = 'HSP-526_ACHA'
#             if 'hse526' in pattern:
#                 residue_1 = 'HSE-526_ACHA'
#
#             out_folder = out_folder_general + str(position) + '/'
#             frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_10/frames_1/'
#             if pattern == 'hsd73_hse526_flu-':
#                 frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_17_ex_12/frames_1/'
#             if pattern == 'hsd73_hse526_red_flu-':
#                 frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_17_ex_10/frames_1/'
#
#             frames = np.arange(0, 645, 2)
#             time_at_frames = np.arange(0, 645, 2)
#             title = str(position) + pattern
#             out_file = pattern
#
#             frame_scale = []
#             for time in time_at_frames:
#                 time = time / 20.0
#                 frame_scale.append(time)
#
#             print "--------"
#             print position, pattern
#             frame_folder_name = 'frame'
#
#             y_lim = [2, 13]
#             y_ticks = np.arange(2.0, 13.0, 1.0)
#
#             x_lim = [0, 32]
#             x_ticks = np.arange(0.0, 33.0, 2.0)
#
#             salt_bridge_distances.plot_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file,
#                                                                     title, frames, frame_scale, frames_folder,
#                                                                     x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim,
#                                                                     y_ticks=y_ticks)
#
#     print('526 - flu')
#     selected_mds_after_interaction_energy = {'hsp73_hsp526_red_fluh': [90], \
#                                              'hsp73_hsp526_red_flu2': [90], \
#                                              'hsp73_hse526_red_fluh': [90], \
#  \
#                                              'hsp73_hsp526_fluh': [90], \
#                                              'hsp73_hsp526_flu-_2': [90], \
#                                              'hsp73_hse526_flu-': [90], \
#                                              }
#     residue_1 = 'HSP-526_ACHA'
#     residue_2 = 'FLU-1_FLUR'
#
#     mds_selection_reorganized = {}
#     mds_selection_reorganized[0] = []
#     mds_selection_reorganized[90] = []
#     mds_selection_reorganized[180] = []
#     mds_selection_reorganized[270] = []
#     for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
#         for position in position_list:
#             if pattern not in mds_selection_reorganized[position]:
#                 mds_selection_reorganized[position].append(pattern)
#
#     for position, pattern_list in mds_selection_reorganized.iteritems():
#         for pattern in pattern_list:
#
#             if 'hsp526' in pattern:
#                 residue_1 = 'HSP-526_ACHA'
#             if 'hse526' in pattern:
#                 residue_1 = 'HSE-526_ACHA'
#
#             if 'flu2' in pattern:
#                 residue_2 = 'FLH-1_FLUR'
#             if 'flu-' in pattern:
#                 residue_2 = 'FLX-1_FLUR'
#             if 'fluh' in pattern:
#                 residue_2 = 'FLU-1_FLUR'
#
#             out_folder = out_folder_general + str(position) + '/'
#             frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_10/frames_1/'
#
#             frames = np.arange(0, 645, 2)
#             time_at_frames = np.arange(0, 645, 2)
#             title = str(position) + pattern
#             out_file = pattern
#
#             frame_scale = []
#             for time in time_at_frames:
#                 time = time / 20.0
#                 frame_scale.append(time)
#
#             print "--------"
#             print position, pattern
#             frame_folder_name = 'frame'
#
#             y_lim = [2, 13]
#             y_ticks = np.arange(2.0, 13.0, 1.0)
#
#             x_lim = [0, 32]
#             x_ticks = np.arange(0.0, 33.0, 2.0)
#
#             salt_bridge_distances.plot_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file,
#                                                                     title, frames, frame_scale, frames_folder,
#                                                                     x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim,
#                                                                     y_ticks=y_ticks)
#
#
#     print('arg62 - flu')
#     selected_mds_after_interaction_energy = {'hsp73_hsp526_red_fluh': [90], \
#                                              'hsp73_hsp526_red_flu2': [90], \
#                                              'hsp73_hse526_red_fluh': [90], \
# \
#                                              'hsp73_hsp526_fluh': [90], \
#                                              'hsp73_hsp526_flu-_2': [90], \
#                                              'hsp73_hse526_flu-': [90], \
#                                              }
#
#     selected_mds_after_interaction_energy = {'hsd73_hse526_flu-': [90], \
#                                              'hsd73_hse526_red_flu-': [90]}
#
#     residue_2 = 'FLU-1_FLUR'
#     residue_1 = 'ARG-62_BCHA'
#
#     mds_selection_reorganized = {}
#     mds_selection_reorganized[0] = []
#     mds_selection_reorganized[90] = []
#     mds_selection_reorganized[180] = []
#     mds_selection_reorganized[270] = []
#     for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
#         for position in position_list:
#             if pattern not in mds_selection_reorganized[position]:
#                 mds_selection_reorganized[position].append(pattern)
#
#     for position, pattern_list in mds_selection_reorganized.iteritems():
#         for pattern in pattern_list:
#
#             if 'flu2' in pattern:
#                 residue_2 = 'FLH-1_FLUR'
#             if 'flu-' in pattern:
#                 residue_2 = 'FLX-1_FLUR'
#             if 'fluh' in pattern:
#                 residue_2 = 'FLU-1_FLUR'
#
#             out_folder = out_folder_general + str(position) + '/'
#             frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_10/frames_1/'
#             if pattern == 'hsd73_hse526_flu-':
#                 frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_17_ex_12/frames_1/'
#             if pattern == 'hsd73_hse526_red_flu-':
#                 frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_17_ex_10/frames_1/'
#
#             frames = np.arange(0, 645, 2)
#             time_at_frames = np.arange(0, 645, 2)
#             title = str(position) + pattern
#             out_file = pattern
#
#             frame_scale = []
#             for time in time_at_frames:
#                 time = time / 20.0
#                 frame_scale.append(time)
#
#             print "--------"
#             print position, pattern
#             frame_folder_name = 'frame'
#
#             y_lim = [2, 13]
#             y_ticks = np.arange(2.0, 13.0, 1.0)
#
#             x_lim = [0, 32]
#             x_ticks = np.arange(0.0, 33.0, 2.0)
#
#             salt_bridge_distances.plot_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file,
#                                                                     title, frames, frame_scale, frames_folder,
#                                                                     x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim,
#                                                                     y_ticks=y_ticks)
#
#     print('526 - 528')
#     selected_mds_after_interaction_energy = {'hsp73_hsp526_red_fluh': [90], \
#                                              'hsp73_hsp526_red_flu2': [90], \
#                                              'hsp73_hse526_red_fluh': [90], \
# \
#                                              'hsp73_hsp526_fluh': [90], \
#                                              'hsp73_hsp526_flu-_2': [90], \
#                                              'hsp73_hse526_flu-': [90], \
#                                              }
#
#     selected_mds_after_interaction_energy = {'hsd73_hse526_flu-': [90], \
#                                              'hsd73_hse526_red_flu-': [90]}
#
#
#
#     residue_1 = 'HSP-526_ACHA'
#     residue_2 = 'ASP-528_ACHA'
#
#     mds_selection_reorganized = {}
#     mds_selection_reorganized[0] = []
#     mds_selection_reorganized[90] = []
#     mds_selection_reorganized[180] = []
#     mds_selection_reorganized[270] = []
#     for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
#         for position in position_list:
#             if pattern not in mds_selection_reorganized[position]:
#                 mds_selection_reorganized[position].append(pattern)
#
#     for position, pattern_list in mds_selection_reorganized.iteritems():
#         for pattern in pattern_list:
#
#             if 'hsp526' in pattern:
#                 residue_1 = 'HSP-526_ACHA'
#             if 'hse526' in pattern:
#                 residue_1 = 'HSE-526_ACHA'
#
#             out_folder = out_folder_general + str(position) + '/'
#             frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_10/frames_1/'
#             if pattern == 'hsd73_hse526_flu-':
#                 frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_17_ex_12/frames_1/'
#             if pattern == 'hsd73_hse526_red_flu-':
#                 frames_folder = frames_folder_general + str(position) + '/' + pattern + '/md_ex_17_ex_10/frames_1/'
#
#             frames = np.arange(0, 645, 2)
#             time_at_frames = np.arange(0, 645, 2)
#             title = str(position) + pattern
#             out_file = pattern
#
#             frame_scale = []
#             for time in time_at_frames:
#                 time = time / 20.0
#                 frame_scale.append(time)
#
#             print "--------"
#             print position, pattern
#             frame_folder_name = 'frame'
#
#             y_lim = [2, 13]
#             y_ticks = np.arange(2.0, 13.0, 1.0)
#
#             x_lim = [0, 32]
#             x_ticks = np.arange(0.0, 33.0, 2.0)
#
#             salt_bridge_distances.plot_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file,
#                                                                     title, frames, frame_scale, frames_folder,
#                                                                     x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim,
#                                                                     y_ticks=y_ticks)
#
#     ########################
#     # plots salt-bridges ###
#     ########################
#
#     frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/wild_type_mds/'
#     out_folder_general = '/user/jdragelj/projects/cco_fluorescein/salt_bridges/wt/'
#
#
#     # print('525 - 526')
#
#     residue_1 = 'HSP-526_ACHA'
#     residue_2 = 'GLU-525_ACHA'
#
#     pattern_list = ['hsp73_hsp526', 'hsp73_hse526', 'hsp73_hsp526_red', 'hsd73_hsp526_red']
#
#     for pattern in pattern_list:
#         out_folder = out_folder_general + '/'
#         frames_folder = frames_folder_general + '/' + pattern + '/md_ex_10/frames_1/'
#
#         if 'hsp526' in pattern:
#             residue_1 = 'HSP-526_ACHA'
#         if 'hse526' in pattern:
#             residue_1 = 'HSE-526_ACHA'
#
#         frames = np.arange(0, 405, 2)
#         time_at_frames = np.arange(0, 405, 2)
#         title = pattern
#         out_file = pattern
#
#
#         frame_scale = []
#         for time in time_at_frames:
#             time = time / 20.0
#             frame_scale.append(time)
#
#         print "--------"
#         print pattern
#         frame_folder_name = 'frame'
#
#         y_lim = [2, 13]
#         y_ticks = np.arange(2.0, 13.0, 1.0)
#
#         x_lim = [0, 20]
#         x_ticks = np.arange(0.0, 21.0, 2.0)
#
#
#         salt_bridge_distances.plot_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file, title, frames, frame_scale, frames_folder,
#                      x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim, y_ticks=y_ticks)
#     print('526 - 528')
#
#     residue_1 = 'HSP-526_ACHA'
#     residue_2 = 'ASP-528_ACHA'
#
#     pattern_list = ['hsp73_hsp526', 'hsp73_hse526', 'hsp73_hsp526_red', 'hsd73_hsp526_red']
#
#     for pattern in pattern_list:
#         out_folder = out_folder_general + '/'
#         frames_folder = frames_folder_general + '/' + pattern + '/md_ex_10/frames_1/'
#
#         if 'hsp526' in pattern:
#             residue_1 = 'HSP-526_ACHA'
#         if 'hse526' in pattern:
#             residue_1 = 'HSE-526_ACHA'
#
#         frames = np.arange(0, 405, 2)
#         time_at_frames = np.arange(0, 405, 2)
#         title = pattern
#         out_file = pattern
#
#
#         frame_scale = []
#         for time in time_at_frames:
#             time = time / 20.0
#             frame_scale.append(time)
#
#         print "--------"
#         print pattern
#         frame_folder_name = 'frame'
#
#         y_lim = [2, 13]
#         y_ticks = np.arange(2.0, 13.0, 1.0)
#
#         x_lim = [0, 20]
#         x_ticks = np.arange(0.0, 21.0, 2.0)
#
#
#         salt_bridge_distances.plot_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file, title, frames, frame_scale, frames_folder,
#                      x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim, y_ticks=y_ticks)
#
#     print('73 - 78')
#
#     residue_1 = 'HSP-73_BCHA'
#     residue_2 = 'GLU-78_BCHA'
#
#     pattern_list = ['hsp73_hsp526', 'hsp73_hse526', 'hsp73_hsp526_red', 'hsd73_hsp526_red']
#
#     for pattern in pattern_list:
#         out_folder = out_folder_general + '/'
#         frames_folder = frames_folder_general + '/' + pattern + '/md_ex_10/frames_1/'
#
#         if 'hsp73' in pattern:
#             residue_1 = 'HSP-73_BCHA'
#         if 'hsd73' in pattern:
#             residue_1 = 'HSD-73_BCHA'
#
#         frames = np.arange(0, 405, 2)
#         time_at_frames = np.arange(0, 405, 2)
#         title = pattern
#         out_file = pattern
#
#
#         frame_scale = []
#         for time in time_at_frames:
#             time = time / 20.0
#             frame_scale.append(time)
#
#         print "--------"
#         print pattern
#         frame_folder_name = 'frame'
#
#         y_lim = [2, 13]
#         y_ticks = np.arange(2.0, 13.0, 1.0)
#
#         x_lim = [0, 20]
#         x_ticks = np.arange(0.0, 21.0, 2.0)
#
#
#         salt_bridge_distances.plot_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file, title, frames, frame_scale, frames_folder,
#                      x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim, y_ticks=y_ticks)

######################################################################################################################################################################################################################################
    # #########################
    # ## plots salt-bridges ### publication
    # #########################
    #
    # frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/wild_type_mds/'
    # out_folder_general = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/salt_bridges/wt/'
    #
    #
    # residue_1 = 'HSP-73_BCHA'
    # residue_2 = 'GLU-78_BCHA'
    #
    # out_folder = out_folder_general + '/'
    # out_file = 'wt-CcO'
    #
    # title1 = 'hsp73_hsp526'
    # frames1 = np.arange(4, 305, 2)
    # time_at_frames1 = np.arange(0, 301, 2)
    # frames_folder1 = frames_folder_general + '/' + title1 + '/md_ex_10/frames_1/'
    # frame_scale1 = []
    # for time in time_at_frames1:
    #     time = time / 20.0
    #     frame_scale1.append(time)
    #
    # title2 = 'hsp73_hsp526_red'
    # frames2 = np.arange(4, 205, 2)
    # time_at_frames2 = np.arange(0, 201, 2)
    # frames_folder2 = frames_folder_general + '/' + title2 + '/md/frames_1/'
    # frame_scale2 = []
    # for time in time_at_frames2:
    #     time = time / 20.0
    #     frame_scale2.append(time)
    #
    # y_lim = [2, 9]
    # y_ticks = np.arange(2.0, 10.0, 1.0)
    #
    # x_lim = [0, 15]
    # x_ticks = np.arange(0.0, 16.0, 2.0)
    #
    #
    # salt_bridge_distances.plot_many_trajs_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file, title1, title2, frames1,
    #                                              frames2, frame_scale1, frame_scale2, frames_folder1, frames_folder2,
    #                                              x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim, y_ticks=y_ticks)
    #
    #
    # #########################
    # ## plots salt-bridges ### publication
    # #########################
    #
    # frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/wild_type_mds/'
    # out_folder_general = '/user/jdragelj/projects/cco_fluorescein/salt_bridges/wt/'
    #
    # residue_1 = 'HSP-526_ACHA'
    # residue_2 = 'GLU-525_ACHA'
    #
    # out_folder = out_folder_general + '/'
    # out_file = 'wt-CcO'
    #
    # title1 = 'hsp73_hsp526'
    # frames1 = np.arange(4, 405, 2)
    # time_at_frames1 = np.arange(0, 401, 2)
    # frames_folder1 = frames_folder_general + title1 + '/md_ex_10/frames_1/'
    # frame_scale1 = []
    # for time in time_at_frames1:
    #     time = time / 20.0
    #     frame_scale1.append(time)
    #
    # title2 = 'hsp73_hsp526_red'
    # frames2 = np.arange(4, 405, 2)
    # time_at_frames2 = np.arange(0, 401, 2)
    # frames_folder2 = frames_folder_general + title2 + '/md_ex_10/frames_1/'
    # frame_scale2 = []
    # for time in time_at_frames2:
    #     time = time / 20.0
    #     frame_scale2.append(time)
    #
    # y_lim = [2, 12]
    # y_ticks = np.arange(2, 13, 1)
    #
    # x_lim = [0, 20]
    # x_ticks = np.arange(0, 21, 1)
    #
    # salt_bridge_distances.plot_many_trajs_salt_bridge_distances_frames(residue_1, residue_2, out_folder, out_file, title1, title2, frames1,
    #                                              frames2, frame_scale1, frame_scale2, frames_folder1, frames_folder2,
    #                                              x_ticks=x_ticks, x_lim=x_lim, y_lim=y_lim, y_ticks=y_ticks)
    #

# ######################################################################################################################################################################################################################################
# ######################################################################################################################################################################################################################################
    ### plots pka/salt-bridge distance his73 ### publication
    ############################################
    #
    # print 'plots pka/salt-bridge distance his73'
    #
    # out_folder_general = '/user/jdragelj/projects/cco_fluorescein/salt_bridges/90/'
    # frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/'
    #
    # y1_lim = [-2, 16]
    # y1_ticks = np.arange(-2, 17, 1)
    #
    # y2_lim = [2, 10]
    # y2_ticks = np.arange(2, 11, 1)
    #
    # x_lim = [0, 32]
    # x_ticks = np.arange(0, 33, 2)
    #
    #
    # title1 = 'hsp73_hsp526_flu-_2'
    # frames1 = np.arange(4, 645, 2)
    # time_at_frames1 = np.arange(0, 641, 2)
    # frames_folder1 = frames_folder_general + title1 + '/md_ex_10/frames_1/'
    # frame_scale1 = []
    # for time in time_at_frames1:
    #     time = time / 20.0
    #     frame_scale1.append(time)
    #
    # labelsx = np.arange(0, 33, 4)
    # string_labels_x = []
    # for i in np.arange(0, 33, 2):
    #     if i not in labelsx:
    #         string_labels_x.append('')
    #     else:
    #         string_labels_x.append(i)
    #
    # labelsy1 = np.arange(-2, 17, 2)
    # string_labels_y1 = []
    # for i in np.arange(-2, 17, 1):
    #     if i not in labelsy1:
    #         string_labels_y1.append('')
    #     else:
    #         string_labels_y1.append(i)
    #
    # labelsy2 = np.arange(2, 11, 1)
    # string_labels_y2 = []
    # for i in np.arange(2, 11, 1):
    #     if i not in labelsy2:
    #         string_labels_y2.append('')
    #     else:
    #         string_labels_y2.append(i)
    #
    # pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/' + title1 + '/'
    # pkas = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames1, 'frame')
    # sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder1, frames1, 'HSP-73_BCHA', 'GLU-78_BCHA', title1, out_folder_general)
    # sb_dist = sb_dist[2:]
    # plot_utilities_frames.frames_plot_two_properties_pub(title1, 'pka', 'sb_dist', pkas, sb_dist, frame_scale1,
    #                                                      out_folder_general, x_ticks, x_lim, y1_lim, y1_ticks, y2_lim,
    #                                                      y2_ticks, string_labels_x, string_labels_y1, string_labels_y2)
    #
    # title2 = 'hsp73_hsp526_red_flu2'
    # frames2 = np.arange(4, 645, 2)
    # time_at_frames2 = np.arange(0, 641, 2)
    # frames_folder2 = frames_folder_general + title2 + '/md_ex_10/frames_1/'
    # frame_scale2 = []
    # for time in time_at_frames2:
    #     time = time / 20.0
    #     frame_scale2.append(time)
    #
    # labelsx = np.arange(0, 33, 4)
    # string_labels_x = []
    # for i in np.arange(0, 33, 2):
    #     if i not in labelsx:
    #         string_labels_x.append('')
    #     else:
    #         string_labels_x.append(i)
    #
    # labelsy1 = np.arange(-2, 17, 2)
    # string_labels_y1 = []
    # for i in np.arange(-2, 17, 1):
    #     if i not in labelsy1:
    #         string_labels_y1.append('')
    #     else:
    #         string_labels_y1.append(i)
    #
    # labelsy2 = np.arange(2, 11, 1)
    # string_labels_y2 = []
    # for i in np.arange(2, 11, 1):
    #     if i not in labelsy2:
    #         string_labels_y2.append('')
    #     else:
    #         string_labels_y2.append(i)
    #
    # pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/' + title2 + '/'
    # pkas = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames2, 'frame')
    # sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder2, frames2, 'HSP-73_BCHA', 'GLU-78_BCHA', title2, out_folder_general)
    # sb_dist = sb_dist[2:]
    # plot_utilities_frames.frames_plot_two_properties_pub(title2, 'pka', 'sb_dist', pkas, sb_dist, frame_scale2,
    #                                                      out_folder_general, x_ticks, x_lim, y1_lim, y1_ticks, y2_lim,
    #                                                      y2_ticks, string_labels_x, string_labels_y1, string_labels_y2)

# # ######################################################################################################################################################################################################################################
#    #salt-bridges / pka wild type publication plots
#
#     out_folder_general = '/user/jdragelj/projects/cco_fluorescein/salt_bridges/wt/'
#     frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/wild_type_mds/'
#
#     y1_lim = [0, 14]
#     y1_ticks = np.arange(0, 15, 1)
#
#     y2_lim = [2, 10]
#     y2_ticks = np.arange(2, 11, 1)
#
#     x_lim = [0, 20]
#     x_ticks = np.arange(0, 21, 1)
#
#
#     title1 = 'hsp73_hsp526'
#     frames1 = np.arange(4, 405, 2)
#     time_at_frames1 = np.arange(0, 401, 2)
#     frames_folder1 = frames_folder_general + '/' + title1 + '/md_ex_10/frames_1/'
#     frame_scale1 = []
#     for time in time_at_frames1:
#         time = time / 20.0
#         frame_scale1.append(time)
#
#     labelsx = np.arange(0, 21, 4)
#     string_labels_x = []
#     for i in np.arange(0, 21, 1):
#         if i not in labelsx:
#             string_labels_x.append('')
#         else:
#             string_labels_x.append(i)
#
#     labelsy1 = np.arange(0, 15, 2)
#     string_labels_y1 = []
#     for i in np.arange(0, 15, 1):
#         if i not in labelsy1:
#             string_labels_y1.append('')
#         else:
#             string_labels_y1.append(i)
#
#     labelsy2 = np.arange(2, 11, 2)
#     string_labels_y2 = []
#     for i in np.arange(2, 11, 1):
#         if i not in labelsy2:
#             string_labels_y2.append('')
#         else:
#             string_labels_y2.append(i)
#
#     pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/' + title1 + '/'
#     pkas = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames1, 'frame')
#     sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder1, frames1, 'HSP-73_BCHA', 'GLU-78_BCHA', title1, out_folder_general)
#     sb_dist = sb_dist[2:]
#     plot_utilities_frames.frames_plot_two_properties_pub(title1, 'pka', 'sb_dist', pkas, sb_dist, frame_scale1,
#                                                          out_folder_general, x_ticks, x_lim, y1_lim, y1_ticks, y2_lim,
#                                                          y2_ticks, string_labels_x, string_labels_y1, string_labels_y2)
#
#
#     title2 = 'hsp73_hsp526_red'
#     frames2 = np.arange(4, 405, 2)
#     time_at_frames2 = np.arange(0, 401, 2)
#     frames_folder2 = frames_folder_general + '/' + title2 + '/md_ex_10/frames_1/'
#     frame_scale2 = []
#     for time in time_at_frames2:
#         time = time / 20.0
#         frame_scale2.append(time)
#
#     pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/' + title2 + '/'
#     pkas = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'HSP-73_BCHA', frames2, 'frame')
#     sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder2, frames2, 'HSP-73_BCHA', 'GLU-78_BCHA', title2, out_folder_general)
#     sb_dist = sb_dist[2:]
#     plot_utilities_frames.frames_plot_two_properties_pub(title2, 'pka', 'sb_dist', pkas, sb_dist, frame_scale2,
#                                                          out_folder_general, x_ticks, x_lim, y1_lim, y1_ticks, y2_lim,
#                                                          y2_ticks, string_labels_x, string_labels_y1, string_labels_y2)
# #
# # #####################################################################################################################################################################################################################################









