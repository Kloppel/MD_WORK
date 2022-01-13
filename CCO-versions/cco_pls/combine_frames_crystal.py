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

# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_md_crystal/10ang_cavities_09_voda/'
# frames = np.arange(20,501,4)
# frame_range = [20,501]
# #
titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_md_crystal/10ang_cavities_09_voda/'
frames = np.arange(20,501,4)
frame_range = [20,501]
frames = np.arange(20,301,4)
frame_range = [20,301]
titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_pra_a3/10ang_cavities_09_voda/md_PF_prah_asp_dep/'
frames = np.arange(0,301,4)
frame_range = [0,301]
frames = np.arange(20,301,4)
frame_range = [20,301]
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_pra_a3/10ang_cavities_09_voda/md_PF_prah_asp_dep_t/'
# frames = np.arange(0,301,4)
# frame_range = [0,301]

# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_md_crystal/10ang_cavities_09_voda/'
# frames = np.arange(20,601,4)
# frame_range = [20,601]
# #
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_md_crystal/10ang_cavities_09_voda/'
# frames = np.arange(20,501,4)
# frame_range = [20,501]
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_pra_a3/10ang_cavities_09_voda/md_F_prah_asp_dep/'
# frames = np.arange(0,301,4)
# frame_range = [0,301]
# titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_pra_a3/10ang_cavities_09_voda/md_F_prah_asp_dep_t/'
# frames = np.arange(0,301,4)
# frame_range = [0,301]


# residue_list = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM', 'ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', \
#             'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', \
#             'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']

residue_list = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM']


#########################################
### GETTING DATA IN THE RIGHT FORMAT ####
#########################################

results = {}
sub_folder = titration_folder
folder_result = kbp2.kbp_results.FrameKbpResults()
for frame in frames:
    file = sub_folder + 'done/' + 'frame' + str(frame) + '/' + 'results.pkl'
    frame_result = pickle.load(open(file, "rb"))
    folder_result.add_task(frame_result.kbp_results[0], frame_nr=frame)
folder_result.average_over_frames(frame_range)
results['crystal'] = folder_result







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
plt.show()

final_pkas = combined_results.combined_results.pkas
# for residue_descr, pka in final_pkas.iteritems():
#     s = " %13s: %0.2f " % (residue_descr, pka)
#     print s

# print 'PRD-3_EHEM', final_pkas['PRD-3_EHEM']
# print 'PRA-1_EHEM', final_pkas['PRA-1_EHEM']
# print 'PRD-3_GHEM', final_pkas['PRD-3_GHEM']
# print 'PRA-1_GHEM', final_pkas['PRA-1_GHEM']

# print 'PRA-1_EHEM', final_pkas['PRA-1_EHEM']
print 'DPP-407_ACHA', final_pkas['DPP-407_ACHA']



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





