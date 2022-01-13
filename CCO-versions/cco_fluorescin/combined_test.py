# coding=utf-8

import os
import cPickle as pickle
import numpy as np
import kbp2
from kbp2.workspace_jd.cco_fluorescin import plots


results_list_oxi = {}
results_list_red = {}
workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/'

print('Newst_trajcs_July2017')
selected_mds_after_interaction_energy = {'hsp73_hsp526_red_flu-_2': [90], \
                                         'hsp73_hsp526_red_fluh': [90], \
                                         'hsp73_hsp526_red_flu2': [90], \

                                         'hsp73_hsp526_fluh': [90], \
                                         'hsp73_hsp526_flu2': [90], \
                                         'hsp73_hsp526_flu-_2': [90], \

                                         'hsp73_hse526_flu-': [90], \
                                         'hsp73_hse526_fluh': [90], \
                                         'hsp73_hse526_flu2': [90], \

                                         'hsp73_hse526_red_flu-': [90],\
                                         'hsp73_hse526_red_fluh': [90],\
                                         'hsp73_hse526_red_flu2': [90],\

                                         # 'hsd73_hse526_red_flu-':[90], \
                                         # 'hsd73_hse526_flu-':[90],\

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


# xs = [4, 44, 84, 124, 164, 204, 244, 284, 324, 364, 404, 444]
xs = np.arange(244, 445, 10)
log = open('/user/jdragelj/Desktop/fludiff_last10.dat', 'w')
for x in xs:
    result_folder = '/user/jdragelj/projects/cco_fluorescein/kbp2_evaluations/'
    for position, pattern_list in mds_selection_reorganized.iteritems():
        for pattern in pattern_list:
            sub_folder = workfolder + str(position) + '/' + pattern + '/'
            folder_result = kbp2.kbp_results.FrameKbpResults()

            frames = np.arange(x,445,2)
            frame_range = [x, 445]


            for frame in frames:
                file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
                frame_result = pickle.load( open( file, "rb" ) )
                if frame_result.kbp_results[0].conf_energy != 0.0:
                    frame_result.kbp_results[0].conf_energy = 0.0
                folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
            folder_result.average_over_frames(frame_range)
            if 'red' in pattern:
                results_list_red[(pattern, str(position))] = folder_result
            else:
                results_list_oxi[(pattern, str(position))] = folder_result

    if not os.path.exists(result_folder):
        os.mkdir(result_folder)

    shift = None

    combined_results_oxi = kbp2.kbp_results.FrameKbpResults()
    order_oxi = []
    for name_tuple, folder_result in results_list_oxi.iteritems():
        pattern, pos_str = name_tuple
        new_position, new_pattern = plots.pattern_conversion_flu(pos_str, pattern)
        traj_curr = pattern+'_'+pos_str
        order_oxi.append(traj_curr)
        combined_results_oxi.add_task(folder_result.avg_results)
    combined_results_oxi.combine_frames_karlsberg(cpus=3)

    nr_of_pacs = len(combined_results_oxi.combined_results.conformer_occs[0])

    #############
    folder = result_folder
    final_pkas = combined_results_oxi.combined_results.pkas
    for residue_descr, pka in final_pkas.iteritems():
        if residue_descr == 'FLU-1_FLUR':
            pka_oxi_flu = pka


    combined_results_red = kbp2.kbp_results.FrameKbpResults()
    order_red = []
    for name_tuple, folder_result in results_list_red.iteritems():
        pattern, pos_str = name_tuple
        new_position, new_pattern = plots.pattern_conversion_flu(pos_str, pattern)
        traj_curr = pattern+'_'+pos_str
        order_red.append(traj_curr)
        combined_results_red.add_task(folder_result.avg_results)
    combined_results_red.combine_frames_karlsberg(cpus=3)

    nr_of_pacs = len(combined_results_red.combined_results.conformer_occs[0])

    #############
    folder = result_folder
    final_pkas = combined_results_red.combined_results.pkas
    for residue_descr, pka in final_pkas.iteritems():
        if residue_descr == 'FLU-1_FLUR':
            pka_red_flu = pka

    shift = pka_red_flu - pka_oxi_flu
    print shift
    log.write(str(shift)+'\n')
log.close()


















    #####################################################

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/titrations/0_his73b_lys354_flu/hsd73b_lsn354_fluh/done/'
    #
    # results_to_combine = kbp2.kbp_results.FrameKbpResults()
    #
    # frame_range = [20, 22, 24]
    #
    # for frame in frame_range:
    #     file = workfolder + 'frame' + str(frame) + '/' + 'results.pkl'
    #     frame_result = pickle.load( open( file, "rb" ) )
    #     print frame_result.kbp_results[0]
    #     results_to_combine.add_task(frame_result.kbp_results[0], frame_nr=frame)
    #
    #
    # results_to_combine.average_over_frames(frame_range=[20,22,24])
    #
    # print results_to_combine.avg_results.conf_energy
    # print results_to_combine.avg_results.g

    # print results_to_combine.avg_results.
    # # results_to_combine.combine_frames_karlsberg()
    #
    # #printout
    # final_pkas = results_to_combine.combined_results.pkas
    # filename = workfolder + 'results_new.dat'
    # results_file = open(filename, 'w')
    #
    # for residue_descr, pka in final_pkas.iteritems():
    #     s = " %13s: %0.2f " % (residue_descr, pka)
    #     results_file.write(s)
    #     results_file.write("\n")
    #
    #
    # results_to_combine.average_over_frames(frame_range=[0,1])
    #
    # print results_to_combine.avg_results.conf_energy





