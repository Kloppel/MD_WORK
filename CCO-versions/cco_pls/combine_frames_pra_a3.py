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

# mds = ['md_pra-','md_prah1','md_prah2']
mds = ['md_pra-','md_prah1']
# frames = np.arange(20,501,2)
# frame_range = [20,501]
frames = np.arange(50,501,50)
frame_range = [50,501]

# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_pra_a3/10ang_cavities_09_voda/'
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra_a3/10ang_cavities_09_voda/'

# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_pra_a3/10ang_cavities_09_voda/'
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_pra_a3/10ang_cavities_08_voda/'
titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_pra_a3/10ang_no_cavities_voda/'
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_pra_a3/10ang_no_cavities_voda_tapbs01/'

# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_pra_a3/10ang_cavities_09_voda/'
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_pra_a3/10ang_cavities_08_voda/'
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_pra_a3/10ang_no_cavities_voda/'
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_pra_a3/10ang_no_cavities_voda_tapbs01/'

residue_list = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM', 'ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', \
            'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', \
            'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']


#########################################
### GETTING DATA IN THE RIGHT FORMAT ####
#########################################

results = {}
for md in mds:

    # frames = np.arange(20, 501, 2)
    # frame_range = [20, 501]
    frames = np.arange(50, 501, 50)
    frame_range = [50, 501]
    if md == 'md_prah2':
        # frames = np.arange(20, 141, 2)
        # frame_range = [20, 141]
        frames = np.arange(50, 141, 50)
        frame_range = [50, 141]

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
    results[md] = folder_result

# # titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/md_crystal/cavities_voda/'
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/md_crystal/10ang_voda/'
# sub_folder = titration_folder
# folder_result = kbp2.kbp_results.FrameKbpResults()
# for frame in frames:
#     file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
#     frame_result = pickle.load(open(file, "rb"))
#     folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
# folder_result.average_over_frames(frame_range)
# results['crystal'] = folder_result

# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_pra_a3/10ang_cavities_09_voda/md_PF_prah_asph/'
# sub_folder = titration_folder
# folder_result = kbp2.kbp_results.FrameKbpResults()
# frames = np.arange(20, 200, 2)
# frame_range = [20, 200]
# for frame in frames:
#     file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
#     frame_result = pickle.load(open(file, "rb"))
#     folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
# folder_result.average_over_frames(frame_range)
# results['prah_asph'] = folder_result

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

###############################
### PLOT TITRATION CURVES  ####
###############################
# residue_list = ['PRD-3_EHEM']
# # residue_list = ['PRA-1_GHEM']
# for residue in residue_list:
#     fig, ax3 = plt.subplots()
#     residue_descr = residue
#
#     residue_nr = combined_results.descr.get_residue_index(residue_descr)
#     deprot_curve = combined_results.combined_results.deprot_curves[residue_nr]
#     ph_values = combined_results.descr.ph_values
#     plt.title("Deprotonation curve for %s" %residue)
#     ax3.plot(ph_values, deprot_curve)
#     plt.show()





