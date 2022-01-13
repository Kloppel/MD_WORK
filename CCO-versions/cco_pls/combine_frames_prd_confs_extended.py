# coding=utf-8

import cPickle as pickle
import numpy as np
import kbp2
from matplotlib import pyplot as plt

# REMARK: conditions for prior pKa computations
#  kbp2_settings.titr_residue_charges_0 = True
#  titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2_np.yaml'

# print 'as a security check tims script again'

##############
### INPUT ####
##############


mds = ['md_Pr_prd-']
titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_cavities_09_voda/'

residue_list = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM', 'ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', \
            'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', \
            'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']

# #########################################
# ### GETTING DATA IN THE RIGHT FORMAT ####
# #########################################
#
results = {}
for md in mds:
    frames = np.arange(20, 121, 2)
    frame_range = [20, 121]
    sub_folder = titration_folder + md + '/'
    # kbp result for the entire MD (empty)
    folder_result = kbp2.kbp_results.FrameKbpResults()
    for frame in frames:
        file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
        # taking kbp result for each frame
        frame_result = pickle.load( open( file, "rb" ) )
        # adding each kbp frame result to folder result
        folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
    # averaging energy terms
    folder_result.average_over_frames(frame_range)
    # stroing the result for the entire MD
    results['open-1-6'] = folder_result
    # results[md] = folder_result

for md in mds:
    frames = np.arange(14, 221, 2)
    frame_range = [140, 221]
    sub_folder = titration_folder + md + '/'
    # kbp result for the entire MD (empty)
    folder_result = kbp2.kbp_results.FrameKbpResults()
    for frame in frames:
        file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
        # taking kbp result for each frame
        frame_result = pickle.load(open(file, "rb"))
        # adding each kbp frame result to folder result
        folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
    # averaging energy terms
    folder_result.average_over_frames(frame_range)
    # stroing the result for the entire MD
    results['open-7-11'] = folder_result
    # results[md] = folder_result

for md in mds:
    frames = np.arange(240, 441, 2)
    frame_range = [240, 441]
    sub_folder = titration_folder + md + '/'
    # kbp result for the entire MD (empty)
    folder_result = kbp2.kbp_results.FrameKbpResults()
    for frame in frames:
        file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
        # taking kbp result for each frame
        frame_result = pickle.load( open( file, "rb" ) )
        # adding each kbp frame result to folder result
        folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
    # averaging energy terms
    folder_result.average_over_frames(frame_range)
    # stroing the result for the entire MD
    results['closed-12-22'] = folder_result

mds = ['md_Pr_prd-_2']
titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_cavities_09_voda/'
for md in mds:
    frames = np.arange(60, 221, 2)
    frame_range = [60, 221]
    sub_folder = titration_folder + md + '/'
    # kbp result for the entire MD (empty)
    folder_result = kbp2.kbp_results.FrameKbpResults()
    for frame in frames:
        file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
        # taking kbp result for each frame
        frame_result = pickle.load( open( file, "rb" ) )
        # adding each kbp frame result to folder result
        folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
    # averaging energy terms
    folder_result.average_over_frames(frame_range)
    # stroing the result for the entire MD
    results['new_open-3-11'] = folder_result

for md in mds:
    frames = np.arange(240, 341, 2)
    frame_range = [240, 341]
    sub_folder = titration_folder + md + '/'
    # kbp result for the entire MD (empty)
    folder_result = kbp2.kbp_results.FrameKbpResults()
    for frame in frames:
        file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
        # taking kbp result for each frame
        frame_result = pickle.load( open( file, "rb" ) )
        # adding each kbp frame result to folder result
        folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
    # averaging energy terms
    folder_result.average_over_frames(frame_range)
    # stroing the result for the entire MD
    results['new_open-12-17'] = folder_result

for md in mds:
    frames = np.arange(400, 501, 2)
    frame_range = [400, 501]
    sub_folder = titration_folder + md + '/'
    # kbp result for the entire MD (empty)
    folder_result = kbp2.kbp_results.FrameKbpResults()
    for frame in frames:
        file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
        # taking kbp result for each frame
        frame_result = pickle.load( open( file, "rb" ) )
        # adding each kbp frame result to folder result
        folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
    # averaging energy terms
    folder_result.average_over_frames(frame_range)
    # stroing the result for the entire MD
    results['new_open-20-25'] = folder_result




###########################
### COMPUTING THE PKAS ####
###########################

combined_results = kbp2.kbp_results.FrameKbpResults()
order = []
for md, folder_result in results.iteritems():
    order.append(md)
    combined_results.add_task(folder_result.avg_results)
combined_results.combine_frames_karlsberg(cpus=3)

############################
### RESULTS PRINT/PLOT  ####
############################

nr_of_pacs = len(combined_results.combined_results.conformer_occs[0])
plt.figure()
x =  combined_results.combined_results.descr.ph_values
for i in range(nr_of_pacs):
    y = combined_results.combined_results.conformer_occs[:, i]
    plt.plot(x, y)

    for x1, y1 in zip(x, y):
        if x1 == 7.0:
            print order[i], x1, y1

plt.ylim([-0.01, 1.01])
plt.legend(order, prop={'size':10})
plt.ylabel('MD occupancies', fontsize=20)
plt.xlabel('pH', fontsize=20)
plt.plot([7.0,7.0,7.0,7.0,7.0,7.0], [0.0,0.2,0.4,0.6,0.8,1.0])


final_pkas = combined_results.combined_results.pkas
for residue_descr, pka in final_pkas.iteritems():
    s = " %13s: %0.2f " % (residue_descr, pka)
    print s
plt.show()

