
import shutil
import os

for frame in range(202, 501, 2):
    folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/new_models/PRD_H_PRA_a3/titrate/md_Pr_prd-/done/'
    frame_folder = folder + 'frame%i/' % frame
    if os.path.exists(frame_folder):
        shutil.rmtree(frame_folder)