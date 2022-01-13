# coding=utf-8

import cPickle as pickle
import numpy as np
import kbp2
from matplotlib import pyplot as plt


residue_list = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM', 'ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', \
            'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', \
            'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']

############
### PRAa ###
############

combined_results = kbp2.kbp_results.FrameKbpResults()
file = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/Pr/praa/pra-/results.pkl'
frame_result = pickle.load( open( file, "rb" ) )
combined_results.add_task(frame_result.kbp_results[0])
file = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/Pr/praa/prah1/results.pkl'
frame_result = pickle.load( open( file, "rb" ) )
combined_results.add_task(frame_result.kbp_results[0])
file = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/Pr/praa/prah2/results.pkl'
frame_result = pickle.load( open( file, "rb" ) )
combined_results.add_task(frame_result.kbp_results[0])
combined_results.combine_frames_karlsberg(cpus=3)

final_pkas = combined_results.combined_results.pkas
for residue_descr, pka in final_pkas.iteritems():
    if residue_descr == 'PRA-1_GHEM':
        s = " %13s: %0.2f " % (residue_descr, pka)
        print s


#############
### PRDa3 ###
#############

combined_results = kbp2.kbp_results.FrameKbpResults()
file = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/Pr/prda3/prd-/results.pkl'
frame_result = pickle.load( open( file, "rb" ) )
combined_results.add_task(frame_result.kbp_results[0])
file = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/Pr/prda3/prdh1/results.pkl'
frame_result = pickle.load( open( file, "rb" ) )
combined_results.add_task(frame_result.kbp_results[0])
file = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/Pr/prda3/prdh2/results.pkl'
frame_result = pickle.load( open( file, "rb" ) )
combined_results.add_task(frame_result.kbp_results[0])
combined_results.combine_frames_karlsberg(cpus=3)

final_pkas = combined_results.combined_results.pkas
for residue_descr, pka in final_pkas.iteritems():
    if residue_descr == 'PRD-3_EHEM':
        s = " %13s: %0.2f " % (residue_descr, pka)
        print s


#############
### PRAa3 ###
#############

combined_results = kbp2.kbp_results.FrameKbpResults()
file = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/Pr/praa3/pra-/results.pkl'
frame_result = pickle.load( open( file, "rb" ) )
combined_results.add_task(frame_result.kbp_results[0])
file = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/Pr/praa3/prah1/results.pkl'
frame_result = pickle.load( open( file, "rb" ) )
combined_results.add_task(frame_result.kbp_results[0])
file = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/Pr/praa3/prah2/results.pkl'
frame_result = pickle.load( open( file, "rb" ) )
combined_results.add_task(frame_result.kbp_results[0])
combined_results.combine_frames_karlsberg(cpus=3)

final_pkas = combined_results.combined_results.pkas
for residue_descr, pka in final_pkas.iteritems():
    if residue_descr == 'PRA-1_EHEM':
        s = " %13s: %0.2f " % (residue_descr, pka)
        print s