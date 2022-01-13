# coding=utf-8

import numpy as np
import shutil
import os

# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/cavities_voda/'
# frames = np.arange(0,601,2)
# for md in ['md_pra-', 'md_prah2', 'md_prah1']:
#     md_folder = folder + md + '/'
#     for frame in frames:
#         if os.path.exists(md_folder + 'done/frame%i/' % frame):
#             if not os.path.exists(md_folder + 'done/frame%i/results.pkl' % frame):
#                 shutil.rmtree(md_folder + 'done/frame%i/' % frame)

folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_prd/cavities_voda/'
frames = np.arange(0,801,2)
for md in ['md_Pr_prd-', 'md_Pr_prdh2', 'md_Pr_prdh1']:
    md_folder = folder + md + '/'
    for frame in frames:
        if os.path.exists(md_folder + 'done/frame%i/' % frame):
            if not os.path.exists(md_folder + 'done/frame%i/results.pkl' % frame):
                shutil.rmtree(md_folder + 'done/frame%i/' % frame)

# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra_a3/cavities_voda/'
# frames = np.arange(0,801,2)
# # for md in ['md_pra-', 'md_prah2', 'md_prah1']:
# for md in ['md_pra-', 'md_prah1']:
#     md_folder = folder + md + '/'
#     for frame in frames:
#         if os.path.exists(md_folder+'done/frame%i/'%frame):
#             if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#                 shutil.rmtree(md_folder+'done/frame%i/'%frame)

# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_md_crystal/cavities_voda_prd_a3/'
# frames = np.arange(0,801,2)
# md_folder = folder
# for frame in frames:
#     if os.path.exists(md_folder+'done/frame%i/'%frame):
#         if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#             shutil.rmtree(md_folder+'done/frame%i/'%frame)

#
#
# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/cavities_voda/'
# frames = np.arange(0,801,2)
# for md in ['md_Pr_prd-', 'md_Pr_prdh2', 'md_Pr_prdh1']:
#     print md
#     md_folder = folder + md + '/'
#     for frame in frames:
#         glu_ok = False
#         asp_ok = False
#         lys_ok = True
#         if os.path.exists(md_folder+'done/frame%i/'%frame):
#             if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#                 shutil.rmtree(md_folder+'done/frame%i/'%frame)
#             else:
#                 init_file = md_folder+'done/frame%i/frame%i_init.pdb'%(frame, frame)
#                 f = open(init_file, 'r')
#                 for line in f:
#                     if 'HE2 GLU A 286' in line:
#                         glu_ok = True
#                     if 'HZ3 LYS A 362' in line:
#                         lys_ok = False
#                     if 'HD2 ASP A 407' in line:
#                         asp_ok = True
#         else:
#             continue
#         if not glu_ok:
#             print 'glu not ok!', md, frame
#         if not asp_ok:
#             print 'asp not ok!', md, frame
#         if not lys_ok:
#             print 'lys not ok!', md, frame
#
# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_voda/'
# frames = np.arange(0,801,2)
# for md in ['md_Pr_prd-', 'md_Pr_prdh2', 'md_Pr_prdh1']:
#     print md
#     md_folder = folder + md + '/'
#     for frame in frames:
#         glu_ok = True
#         asp_ok = True
#         lys_ok = True
#         if os.path.exists(md_folder+'done/frame%i/'%frame):
#             if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#                 shutil.rmtree(md_folder+'done/frame%i/'%frame)
#             else:
#                 init_file = md_folder + 'done/frame%i/frame%i_init.pdb' % (frame, frame)
#                 f = open(init_file, 'r')
#                 for line in f:
#                     if 'GLU A 286' in line:
#                         glu_ok = False
#                     if 'HZ3 LYS A 362' in line:
#                         lys_ok = False
#                     if 'ASP A 407' in line:
#                         asp_ok = False
#         else:
#             continue
#         if not glu_ok:
#             print 'glu not ok!', md, frame
#         if not asp_ok:
#             print 'asp not ok!', md, frame
#         if not lys_ok:
#             print 'lys not ok!', md, frame
#
# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_prd/cavities_voda/'
# frames = np.arange(0,801,2)
# for md in ['md_Pr_prd-', 'md_Pr_prdh2', 'md_Pr_prdh1']:
#     print md
#     md_folder = folder + md + '/'
#     for frame in frames:
#         glu_ok = False
#         asp_ok = False
#         lys_ok = True
#         if os.path.exists(md_folder+'done/frame%i/'%frame):
#             if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#                 shutil.rmtree(md_folder+'done/frame%i/'%frame)
#             else:
#                 init_file = md_folder+'done/frame%i/frame%i_init.pdb'%(frame, frame)
#                 f = open(init_file, 'r')
#                 for line in f:
#                     if 'HE2 GLU A 286' in line:
#                         glu_ok = True
#                     if 'HZ3 LYS A 362' in line:
#                         lys_ok = False
#                     if 'HD2 ASP A 407' in line:
#                         asp_ok = True
#         else:
#             continue
#         if not glu_ok:
#             print 'glu not ok!', md, frame
#         if not asp_ok:
#             print 'asp not ok!', md, frame
#         if not lys_ok:
#             print 'lys not ok!', md, frame
#
#
# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_prd/10ang_voda/'
# frames = np.arange(0,801,2)
# for md in ['md_Pr_prd-', 'md_Pr_prdh2', 'md_Pr_prdh1']:
#     print md
#     md_folder = folder + md + '/'
#     for frame in frames:
#         glu_ok = True
#         asp_ok = True
#         lys_ok = True
#         if os.path.exists(md_folder+'done/frame%i/'%frame):
#             if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#                 shutil.rmtree(md_folder+'done/frame%i/'%frame)
#             else:
#                 init_file = md_folder + 'done/frame%i/frame%i_init.pdb' % (frame, frame)
#                 f = open(init_file, 'r')
#                 for line in f:
#                     if 'GLU A 286' in line:
#                         glu_ok = False
#                     if 'HZ3 LYS A 362' in line:
#                         lys_ok = False
#                     if 'ASP A 407' in line:
#                         asp_ok = False
#         else:
#             continue
#         if not glu_ok:
#             print 'glu not ok!', md, frame
#         if not asp_ok:
#             print 'asp not ok!', md, frame
#         if not lys_ok:
#             print 'lys not ok!', md, frame
#
#
# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/cavities_voda/'
# frames = np.arange(0,601,2)
# for md in ['md_pra-', 'md_prah2', 'md_prah1']:
#     print md
#     md_folder = folder + md + '/'
#     for frame in frames:
#         glu_ok = False
#         asp_ok = False
#         lys_ok = True
#         if os.path.exists(md_folder+'done/frame%i/'%frame):
#             if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#                 shutil.rmtree(md_folder+'done/frame%i/'%frame)
#             else:
#                 init_file = md_folder+'done/frame%i/frame%i_init.pdb'%(frame, frame)
#                 f = open(init_file, 'r')
#                 for line in f:
#                     if 'HE2 GLU A 286' in line:
#                         glu_ok = True
#                     if 'HZ3 LYS A 362' in line:
#                         lys_ok = False
#                     if 'HD2 ASP A 407' in line:
#                         asp_ok = True
#         else:
#             continue
#         if not glu_ok:
#             print 'glu not ok!', md, frame
#         if not asp_ok:
#             print 'asp not ok!', md, frame
#         if not lys_ok:
#             print 'lys not ok!', md, frame
#
#
# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/10ang_voda/'
# frames = np.arange(0,601,2)
# for md in ['md_pra-', 'md_prah2', 'md_prah1']:
#     print md
#     md_folder = folder + md + '/'
#     for frame in frames:
#         glu_ok = True
#         asp_ok = True
#         lys_ok = True
#         if os.path.exists(md_folder+'done/frame%i/'%frame):
#             if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#                 shutil.rmtree(md_folder+'done/frame%i/'%frame)
#             else:
#                 init_file = md_folder + 'done/frame%i/frame%i_init.pdb' % (frame, frame)
#                 f = open(init_file, 'r')
#                 for line in f:
#                     if 'GLU A 286' in line:
#                         glu_ok = False
#                     if 'HZ3 LYS A 362' in line:
#                         lys_ok = False
#                     if 'ASP A 407' in line:
#                         asp_ok = False
#         else:
#             continue
#         if not glu_ok:
#             print 'glu not ok!', md, frame
#         if not asp_ok:
#             print 'asp not ok!', md, frame
#         if not lys_ok:
#             print 'lys not ok!', md, frame
#
#
#
#
#
#
#
# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra_a3/cavities_voda/'
# frames = np.arange(0,601,2)
# for md in ['md_pra-', 'md_prah1', 'md_prah2']:
#     print md
#     if md == 'md_prah2':
#         frames = np.arange(0, 401, 2)
#     md_folder = folder + md + '/'
#     for frame in frames:
#         glu_ok = False
#         asp_ok = True
#         lys_ok = True
#         if os.path.exists(md_folder+'done/frame%i/'%frame):
#             if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#                 shutil.rmtree(md_folder+'done/frame%i/'%frame)
#             else:
#                 init_file = md_folder+'done/frame%i/frame%i_init.pdb'%(frame, frame)
#                 f = open(init_file, 'r')
#                 for line in f:
#                     if 'HE2 GLU A 286' in line:
#                         glu_ok = True
#                     if 'HZ3 LYS A 362' in line:
#                         lys_ok = False
#                     if 'HD2 ASP A 407' in line:
#                         asp_ok = False
#         else:
#             continue
#
#         if not glu_ok:
#             print 'glu not ok!', md, frame
#         if not asp_ok:
#             print 'asp not ok!', md, frame
#         if not lys_ok:
#             print 'lys not ok!', md, frame


# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/10ang_voda/'
# frames = np.arange(0,601,2)
# for md in ['md_pra-', 'md_prah2', 'md_prah1']:
#     md_folder = folder + md + '/'
#     for frame in frames:
#         glu_ok = True
#         asp_ok = True
#         lys_ok = True
#         if os.path.exists(md_folder+'done/frame%i/'%frame):
#             if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#                 shutil.rmtree(md_folder+'done/frame%i/'%frame)
#             else:
#                 init_file = md_folder + 'done/frame%i/frame%i_init.pdb' % (frame, frame)
#                 f = open(init_file, 'r')
#                 for line in f:
#                     if 'GLU A 286' in line:
#                         glu_ok = False
#                     if 'HZ3 LYS A 362' in line:
#                         lys_ok = False
#                     if 'ASP A 407' in line:
#                         asp_ok = False
#         else:
#             continue
#         if not glu_ok:
#             print 'glu not ok!', md, frame
#         if not asp_ok:
#             print 'asp not ok!', md, frame
#         if not lys_ok:
#             print 'lys not ok!', md, frame
#
#
# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_md_crystal/10ang_voda/'
# frames = np.arange(0,801,2)
# md_folder = folder
# md = 'crystal_10ang_voda'
# for frame in frames:
#     glu_ok = True
#     asp_ok = True
#     lys_ok = True
#     if os.path.exists(md_folder+'done/frame%i/'%frame):
#         if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#             shutil.rmtree(md_folder+'done/frame%i/'%frame)
#         else:
#             init_file = md_folder + 'done/frame%i/frame%i_init.pdb' % (frame, frame)
#             f = open(init_file, 'r')
#             for line in f:
#                 if 'GLU A 286' in line:
#                     glu_ok = False
#                 if 'HZ3 LYS A 362' in line:
#                     lys_ok = False
#                 if 'ASP A 407' in line:
#                     asp_ok = False
#     else:
#         continue
#     if not glu_ok:
#         print 'glu not ok!', md, frame
#     if not asp_ok:
#         print 'asp not ok!', md, frame
#     if not lys_ok:
#         print 'lys not ok!', md, frame
#
#
# folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_md_crystal/cavities_voda_prd_a3/'
# frames = np.arange(0,801,2)
# md_folder = folder
# md = 'crystal_cavities_prd_a3'
# for frame in frames:
#     glu_ok = False
#     asp_ok = False
#     lys_ok = True
#     if os.path.exists(md_folder+'done/frame%i/'%frame):
#         if not os.path.exists(md_folder+'done/frame%i/results.pkl'%frame):
#             shutil.rmtree(md_folder+'done/frame%i/'%frame)
#         else:
#             init_file = md_folder+'done/frame%i/frame%i_init.pdb'%(frame, frame)
#             f = open(init_file, 'r')
#             for line in f:
#                 if 'HE2 GLU A 286' in line:
#                     glu_ok = True
#                 if 'HZ3 LYS A 362' in line:
#                     lys_ok = False
#                 if 'HD2 ASP A 407' in line:
#                     asp_ok = True
#     else:
#         continue
#
#     if not glu_ok:
#         print 'glu not ok!', md, frame
#     if not asp_ok:
#         print 'asp not ok!', md, frame
#     if not lys_ok:
#         print 'lys not ok!', md, frame




