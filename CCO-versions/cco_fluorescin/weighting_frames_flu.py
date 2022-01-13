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
# workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/'
# workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/Cua_hemea/'




selected_mds_after_interaction_energy = {'hsp73_hsp526_red_flu2': [90], \
                                         'hsp73_hsp526_flu-_2': [90]}


mds_selection_reorganized = {}
mds_selection_reorganized[0] = []
mds_selection_reorganized[90] = []
mds_selection_reorganized[180] = []
mds_selection_reorganized[270] = []
for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
    for position in position_list:
        if pattern not in mds_selection_reorganized[position]:
            mds_selection_reorganized[position].append(pattern)


result_folder = '/user/jdragelj/Desktop/'
for position, pattern_list in mds_selection_reorganized.iteritems():
    for pattern in pattern_list:
        sub_folder = workfolder + str(position) + '/' + pattern + '/'

        frames = np.arange(344,645,2)

        for frame in frames:
            folder_result = kbp2.kbp_results.FrameKbpResults()
            file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
            frame_result = pickle.load( open( file, "rb" ) )
            if frame_result.kbp_results[0].conf_energy != 0.0:
                frame_result.kbp_results[0].conf_energy = 0.0
            folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
            frame_range = [frame, frame]
            folder_result.average_over_frames(frame_range)

            if 'red' in pattern:
                results_list_red[(frame, pattern)] = folder_result
            else:
                results_list_oxi[(frame, pattern)] = folder_result


if not os.path.exists(result_folder):
    os.mkdir(result_folder)

shift = None

print 'oxi'
combined_results_oxi = kbp2.kbp_results.FrameKbpResults()
order_oxi = []
for name_tuple, folder_result in results_list_oxi.iteritems():
    frame, pattern= name_tuple
    traj_curr = str(frame) + '_' + pattern
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
        if x1 == 7.0:
            if y1 > 0.045:
                print order_oxi[i], x1, y1

plt.ylim([-0.01, 1.01])
plt.title('Trajectories in state O')
plt.legend(order_oxi, prop={'size':10})
plt.ylabel('occupancy', fontsize=20)
plt.xlabel('pH', fontsize=20)
plt.plot([7.0,7.0,7.0,7.0,7.0,7.0], [0.0,0.2,0.4,0.6,0.8,1.0], linestyle = '-', color='k')
plt.show()




#############
folder = result_folder
final_pkas = combined_results_oxi.combined_results.pkas
for residue_descr, pka in final_pkas.iteritems():
    s = " %13s: %0.2f " % (residue_descr, pka)
    print s
    if residue_descr == 'FLU-1_FLUR':
        pka_oxi_flu = pka

print '------------------------------------------------------------------------------------------'

print 'red'
combined_results_red = kbp2.kbp_results.FrameKbpResults()
order_red = []
for name_tuple, folder_result in results_list_red.iteritems():
    frame, pattern= name_tuple
    traj_curr = str(frame) + '_' + pattern
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
        if x1 == 7.0:
            if y1 > 0.1:
                print order_red[i], x1, y1

plt.ylim([-0.01, 1.01])
plt.title('Trajectories in state R')
plt.ylabel('occupancy', fontsize=20)
plt.xlabel('pH', fontsize=20)
plt.legend(order_red, prop={'size':10})
plt.plot([7.0,7.0,7.0,7.0,7.0,7.0], [0.0,0.2,0.4,0.6,0.8,1.0], linestyle = '-', color='k')
plt.show()



#############
folder = result_folder
final_pkas = combined_results_red.combined_results.pkas
for residue_descr, pka in final_pkas.iteritems():
    s = " %13s: %0.2f " % (residue_descr, pka)
    print s
    if residue_descr == 'FLU-1_FLUR':
        pka_red_flu = pka

# shift = pka_red_flu - pka_oxi_flu
# print 'pKa shift: ', shift