# coding=utf-8

import cPickle as pickle
import os

import numpy as np

import kbp2
from kbp2.workspace_jd.cco_fluorescin import plots


#####################
### COMBINED MDS ####
#####################

his_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/'
# his_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/wild_type/'
patterns = os.listdir(his_folder)

frames = np.arange(104,405,2)
frame_range = [104, 405]

results_list_red = {}
results_list_oxi = {}

for pattern in patterns:

    # if pattern not in ['hsp73_hsp526', 'hsp73_hse526', 'hsp73_hsp526_red', 'hsd73_hsp526_red']:
    #     continue

    sub_folder = his_folder + pattern + '/'
    folder_result = kbp2.kbp_results.FrameKbpResults()
    for frame in frames:
        file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
        frame_result = pickle.load( open( file, "rb" ) )
        folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
    folder_result.average_over_frames(frame_range)
    if 'red' in pattern:
        results_list_red[pattern] = folder_result
    else:
        results_list_oxi[pattern] = folder_result

###########################################################################################################################

from matplotlib import pyplot as plt
print 'oxi'
combined_results_oxi = kbp2.kbp_results.FrameKbpResults()
order_oxi = []
for name_tuple, folder_result in results_list_oxi.iteritems():
    pattern = name_tuple
    # new_pattern = plots.pattern_conversion_his(pattern)
    traj_curr = pattern
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
        if (i % 2 == 0):
            linest = '--'
        else:
            linest = '-'
        plt.plot(x, y, linestyle=linest)

plt.ylim([-0.01, 1.01])
plt.title('Trajectories in state O')
plt.legend(order_oxi, prop={'size':10})
plt.ylabel('occupancy', fontsize=20)
plt.xlabel('pH', fontsize=20)
plt.plot([7.0,7.0,7.0,7.0,7.0,7.0], [0.0,0.2,0.4,0.6,0.8,1.0], linestyle = '-', color='k')
# plt.show()



#############
folder = '/user/jdragelj/Desktop/'
final_pkas = combined_results_oxi.combined_results.pkas
# filename = folder + 'results_oxi.dat'
# results_file = open(filename, 'w')
for residue_descr, pka in final_pkas.iteritems():
    s = " %13s: %0.2f " % (residue_descr, pka)
    print s
    # results_file.write(s)
    # results_file.write("\n")

plt.show()


print '------------------------------------------------------------------------------------------'

print 'red'
combined_results_red = kbp2.kbp_results.FrameKbpResults()
order_red = []
for name_tuple, folder_result in results_list_red.iteritems():
    pattern = name_tuple
    # new_pattern = plots.pattern_conversion_his(pattern)
    traj_curr = pattern
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
        if (i % 2 == 0):
            linest = '--'
        else:
            linest = '-'
        plt.plot(x, y, linestyle=linest)

plt.ylim([-0.01, 1.01])
plt.title('Trajectories in state R')
plt.ylabel('occupancy', fontsize=20)
plt.xlabel('pH', fontsize=20)
plt.legend(order_red, prop={'size':10})
plt.plot([7.0,7.0,7.0,7.0,7.0,7.0], [0.0,0.2,0.4,0.6,0.8,1.0], linestyle = '-', color='k')
# plt.show()


#############
folder = '/user/jdragelj/Desktop/'
final_pkas = combined_results_red.combined_results.pkas
# filename = folder + 'results_red.dat'
# results_file = open(filename, 'w')
for residue_descr, pka in final_pkas.iteritems():
    s = " %13s: %0.2f " % (residue_descr, pka)
    print s
    # results_file.write(s)
    # results_file.write("\n")

plt.show()

