# coding=utf-8

import os
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import kbp2



###########################################################################################################################
results_list_oxi = {}
results_list_red = {}
workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/'
# workfolder = '/media/jdragelj/88bc6e86-bfdd-41a3-b3f0-9921ae91f5b7/jovan/backup/scratch/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/'
# workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/'
print workfolder
#
# print('Newst_trajcs_July2017_30ns')
# selected_mds_after_interaction_energy = {'hsp73_hsp526_red_flu-_2': [90], \
#                                          'hsp73_hsp526_red_fluh': [90], \
#                                          'hsp73_hsp526_red_flu2': [90], \
#
#                                          'hsp73_hsp526_fluh': [90], \
#                                          'hsp73_hsp526_flu2': [90], \
#                                          'hsp73_hsp526_flu-_2': [90], \
#
#                                          'hsp73_hse526_flu-': [90], \
#                                          'hsp73_hse526_fluh': [90], \
#                                          'hsp73_hse526_flu2': [90], \
#
#                                          'hsp73_hse526_red_flu-': [90],\
#                                          'hsp73_hse526_red_fluh': [90],\
#                                          'hsp73_hse526_red_flu2': [90],\
#
#                                          'hsd73_hse526_red_flu-': [90], \
#                                          'hsd73_hse526_flu-': [90], \
#                                           }
# #
selected_mds_after_interaction_energy = {'hsp73_hsp526_red_fluh': [90], \
                                         'hsp73_hsp526_red_flu2': [90], \

                                         'hsp73_hsp526_fluh': [90], \
                                         'hsp73_hsp526_flu-_2': [90], \
                                         'hsp73_hse526_flu-': [90], \

                                         'hsd73_hse526_red_flu-': [90], \
                                         'hsd73_hse526_flu-': [90], \
                                          }

# selected_mds_after_interaction_energy = {'hsp73_hsp526_red_flu2': [90], \
#                                          'hsp73_hsp526_flu-_2': [90]}


# print 'preliminary work!'
# selected_mds_after_interaction_energy = {'hsd73_hse526_flu-':[90],\
#                                          'hsd73_hse526_fluh':[180],\
#                                          'hsd73_hsp526_flu2':[270],\
#                                          'hsd73_hsp526_red_flu-':[180],\
#                                          'hsd73_hsp526_red_fluh':[180],\
#                                          'hsd73_hse526_red_flu-':[90]}

mds_selection_reorganized = {}
mds_selection_reorganized[0] = []
mds_selection_reorganized[90] = []
mds_selection_reorganized[180] = []
mds_selection_reorganized[270] = []
for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
    for position in position_list:
        if pattern not in mds_selection_reorganized[position]:
            mds_selection_reorganized[position].append(pattern)


# frames = np.arange(44,245,2)
# frame_range = [44, 245]
# frames = np.arange(144,345,2)
# frame_range = [144, 345]
# frames = np.arange(244,445,2)
# frame_range = [244, 445]
# frames = np.arange(244,345,2)
# frame_range = [244, 345]
# frames = np.arange(344,545,2)
# frame_range = [344, 545]
# frames = np.arange(344,445,2)
# frame_range = [344, 445]
# frames = np.arange(444,645,2)
# frame_range = [444, 645]

frames = np.arange(344,645,2)
frame_range = [344, 645]
print 'time interval:',(frame_range[0]-4)/20,'-',(frame_range[1]-4)/20,'ns'

# result_folder = '/user/jdragelj/projects/cco_fluorescein/kbp2_evaluations/'
result_folder = '/user/jdragelj/Desktop/'
for position, pattern_list in mds_selection_reorganized.iteritems():
    for pattern in pattern_list:
        sub_folder = workfolder + str(position) + '/' + pattern + '/'
        folder_result = kbp2.kbp_results.FrameKbpResults()



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


# ###########################################################################################################################

# #
# selected_mds_after_interaction_energy2 = {'hsd73_hse526_flu-':[90],\
#                                          'hsd73_hse526_fluh':[180],\
#                                          'hsd73_hsp526_flu2':[270],\
#                                          'hsd73_hsp526_red_flu-':[180],\
#                                          'hsd73_hsp526_red_fluh':[180],\
#                                          'hsd73_hse526_red_flu-':[90]}
# #

#
#
# mds_selection_reorganized2 = {}
# mds_selection_reorganized2[0] = []
# mds_selection_reorganized2[90] = []
# mds_selection_reorganized2[180] = []
# mds_selection_reorganized2[270] = []
# for pattern2, position_list2 in selected_mds_after_interaction_energy2.iteritems():
#     for position2 in position_list2:
#         if pattern2 not in mds_selection_reorganized2[position2]:
#             mds_selection_reorganized2[position2].append(pattern2)
#
#
# workfolder = '/media/jdragelj/88bc6e86-bfdd-41a3-b3f0-9921ae91f5b7/jovan/backup/scratch/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/'
# # result_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/kbp_evaluations/'
# result_folder = '/scratch/scratch/jdragelj/Desktop/'
#
# print '3ns trajectroires are on 4TB drive'
#
# frames = np.arange(244, 545, 2)
# frame_range = [244, 545]
#
# for position, pattern_list in mds_selection_reorganized2.iteritems():
#     for pattern in pattern_list:
#         sub_folder = workfolder + str(position) + '/' + pattern + '/'
#         folder_result = kbp2.kbp_results.FrameKbpResults()
#         print position, pattern
#         for frame in frames:
#
#             if pattern == 'hsd73_hse526_flu-':
#                 frames = np.arange(344, 645, 2)
#                 frame_range = [344, 645]
#
#             file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
#             frame_result = pickle.load( open( file, "rb" ) )
#             folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
#         folder_result.average_over_frames(frame_range)
#         if 'red' in pattern:
#             # traj_name = str(position) + '_' + pattern
#             results_list_red[(pattern, str(position))] = folder_result
#             # results_list_red[(pattern, str(position))] = folder_result
#         else:
#             # traj_name = str(position) + '_' + pattern
#             results_list_oxi[(pattern, str(position))] = folder_result
#             # results_list_oxi[(pattern, str(position))] = folder_result
#
# ##########################################################################################################################


if not os.path.exists(result_folder):
    os.mkdir(result_folder)

shift = None

print 'oxi'
combined_results_oxi = kbp2.kbp_results.FrameKbpResults()
order_oxi = []
for name_tuple, folder_result in results_list_oxi.iteritems():
    pattern, pos_str = name_tuple
    traj_curr = pattern+'_'+pos_str
    order_oxi.append(traj_curr)
    combined_results_oxi.add_task(folder_result.avg_results)
combined_results_oxi.combine_frames_karlsberg(cpus=3)

nr_of_pacs = len(combined_results_oxi.combined_results.conformer_occs[0])
plt.figure()
x =  combined_results_oxi.combined_results.descr.ph_values
types = ['-.', '--', '-']
ti = 0
for i in range(nr_of_pacs):
    y = combined_results_oxi.combined_results.conformer_occs[:, i]
    if nr_of_pacs <= 3:
        plt.plot(x, y, linestyle=types[ti],color='k')
        ti +=1
    else:
        plt.plot(x, y)
    for x1, y1 in zip(x, y):
        # if x1 == 7.0:
        if x1 == 8.0:
            if round(y1,2) >= 0.01:
                print order_oxi[i], x1, y1

plt.ylim([-0.01, 1.01])
plt.title('Trajectories in state O')
plt.legend(order_oxi, prop={'size':10})
plt.ylabel('occupancy', fontsize=20)
plt.xlabel('pH', fontsize=20)
# plt.plot([7.0,7.0,7.0,7.0,7.0,7.0], [0.0,0.2,0.4,0.6,0.8,1.0], linestyle = '-', color='k')
plt.plot([8.0,8.0,8.0,8.0,8.0,8.0], [0.0,0.2,0.4,0.6,0.8,1.0], linestyle = '-', color='k')
# plt.savefig(result_folder+'oxi.jpeg')
plt.show()




#############
folder = result_folder
final_pkas = combined_results_oxi.combined_results.pkas
# filename = folder + 'results_oxi.dat'
# results_file = open(filename, 'w')
for residue_descr, pka in final_pkas.iteritems():
    s = " %13s: %0.2f " % (residue_descr, pka)
    print s
    # results_file.write(s)
    # results_file.write("\n")

    if residue_descr == 'FLU-1_FLUR':
        pka_oxi_flu = pka

print '------------------------------------------------------------------------------------------'

print 'red'
combined_results_red = kbp2.kbp_results.FrameKbpResults()
order_red = []
for name_tuple, folder_result in results_list_red.iteritems():
    pattern, pos_str = name_tuple
    traj_curr = pattern+'_'+pos_str
    order_red.append(traj_curr)
    combined_results_red.add_task(folder_result.avg_results)
combined_results_red.combine_frames_karlsberg(cpus=3)

nr_of_pacs = len(combined_results_red.combined_results.conformer_occs[0])

plt.figure()
x =  combined_results_red.combined_results.descr.ph_values
types = ['-.', '--', '-']
ti = 0
for i in range(nr_of_pacs):
    y = combined_results_red.combined_results.conformer_occs[:, i]
    if nr_of_pacs <= 3:
        plt.plot(x, y, linestyle=types[ti],color='k')
        ti +=1
    else:
        plt.plot(x, y)
    for x1, y1 in zip(x, y):
        # if x1 == 7.0:
        if x1 == 8.0:
            if round(y1, 2) >= 0.01:
                print order_red[i], x1, y1

plt.ylim([-0.01, 1.01])
plt.title('Trajectories in state R')
plt.ylabel('occupancy', fontsize=20)
plt.xlabel('pH', fontsize=20)
plt.legend(order_red, prop={'size':10})
# plt.plot([7.0,7.0,7.0,7.0,7.0,7.0], [0.0,0.2,0.4,0.6,0.8,1.0], linestyle = '-', color='k')
plt.plot([8.0,8.0,8.0,8.0,8.0,8.0], [0.0,0.2,0.4,0.6,0.8,1.0], linestyle = '-', color='k')
# plt.savefig(result_folder+'red.jpeg')
plt.show()



#############
folder = result_folder
final_pkas = combined_results_red.combined_results.pkas
# filename = folder + 'results_red.dat'
# results_file = open(filename, 'w')
for residue_descr, pka in final_pkas.iteritems():
    s = " %13s: %0.2f " % (residue_descr, pka)
    print s
    # results_file.write(s)
    # results_file.write("\n")

    if residue_descr == 'FLU-1_FLUR':
        pka_red_flu = pka

shift = round(pka_red_flu, 2) - round(pka_oxi_flu,2)
print 'pKa shift: ', shift