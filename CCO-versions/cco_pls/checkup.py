# coding=utf-8



mds = ['md_pra-', 'md_prah1', 'md_prah2']
folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_pra_a3/10ang_no_cavities_voda_tapbs01/'
folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_pra_a3/10ang_no_cavities_voda_tapbs01/'

# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_xtra/prd/10ang_no_cavities_voda_tapbs01/'
# mds = ['md_F_prd-', 'md_F_prdh1', 'md_F_prdh2']

folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_xtra/prd/10ang_no_cavities_voda_tapbs01/'
mds = ['md_PF_prd-', 'md_PF_prdh1', 'md_PF_prdh2']


for md in mds:
    check_folder = folder + md + '/done/'
    import os
    folds = os.listdir(check_folder)
    for frame in folds:
        f_to_check = check_folder + frame + '/titrate_%s.py' % frame
        f = open(f_to_check)
        checkup = False
        for line in f:
            if 'pka_calculation_tapbs_hr' in line:
                checkup = True
        if not checkup:
            print f_to_check




