
import shutil

folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_PF_crystal_PRDa-/pls_md_PRDa-/frames_voda_old/'
new_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_PF_crystal_PRDa-/pls_md_PRDa-/frames_voda/'

for i in range(0,500,1):
    shutil.copy2(folder + 'frame%i.pdb' % i, new_folder + 'frame%i.pdb' % (i+1))
