# coding=utf-8

import kbp2
import os
import cPickle as pickle

def titrate_structure( pdb_structure, workfolder, state, residues_tt=None, cavities=[], membrane_charge=True, glu_asp_lys_patch=[False, False, True], dolly=False, jobname=''):

    #############
    ### CALCS ###
    #############

    top = []
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw_clean_kb.inp")
    # top.append('/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp')
    if membrane_charge:
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
    else:
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid_no_charge_POPC.rtf")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.rtf")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf")

    par = []
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")


    kbp2_settings = kbp2.pka_calculation.PkaCalcSettings()
    kbp2_settings.top = top
    kbp2_settings.par = par
    kbp2_settings.processes = 1
    kbp2_settings.quiet_mode = True

    decisions = []
    decisions.append(('rename__HOH_TIP3', 'keep'))
    decisions.append(('disu_bridges', 'closed'))
    decisions.append(('cap_termini', 'dont_cap'))

    kbp2_settings.modelling_decisions = decisions
    kbp2_settings.preopt = {}

    protocol = []
    protocol.append((  7, 'h_min'))
    kbp2_settings.protocol = protocol
    if cavities:
        kbp2_settings.cavity_par = cavities

    titratable_yaml = '/scratch/scratch/jdragelj/KB_files/titratable_cco_noglupasp.yaml'
    # titratable_yaml = '/scratch/scratch/jdragelj/KB_files/titratable_cco_extended.yaml'

    kbp2_settings.set_yaml(titratable_yaml)
    titratable_definitions = kbp2.kbp_tools.parse_titratable_yaml(titratable_yaml)
    ##############################################################################

                    #########
                    ### 1 ###        # 1. reduces of oxidized enzime
                    #########


    kbp2_settings.patches = []
    kbp2_settings.patches_no_autogen = []
    charge_patches = []
    bond_patches = []


    if state == 'Pm':
        print 'heme A is reduced'
        charge_patches.append({'AHE2': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                        'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'A3H4': ['HEM-2_EHEM', 'HSD-419_ACHA', 'HAO2-1_FEOH']})
        charge_patches.append({'CBP2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISO': ['HAO2-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                      'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
    if state == 'PF':
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                        'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'A343': ['HEM-2_EHEM', 'HSD-419_ACHA', 'HAO2-1_FEOH']})
        charge_patches.append({'CBP2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISO': ['HAO2-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                      'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
    if state == 'Pr':
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                        'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'A3H4': ['HEM-2_EHEM', 'HSD-419_ACHA', 'HAO2-1_FEOH']})
        charge_patches.append({'CBPN': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISO': ['HAO2-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                      'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
    # if state == 'Pr_deprot_His334':
    #     charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
    #     charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
    #                                     'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
    #     charge_patches.append({'A3H4': ['HEM-2_EHEM', 'HSD-419_ACHA', 'HAO2-1_FEOH']})
    #     charge_patches.append({'CBDD': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA','HSE-284_ACHA']})
    #
    #     bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
    #     bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
    #     bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
    #     bond_patches.append({'EISO': ['HAO2-1_FEOH', 'HEM-2_EHEM']})
    #     bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC']})
    #     bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
    #                                   'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    if state == 'F':
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                        'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'A3H4': ['HEM-2_EHEM', 'HSD-419_ACHA', 'HAO2-1_FEOH']})
        charge_patches.append({'CBF2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISO': ['HAO2-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUBP': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                      'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    #protonation patches
    if glu_asp_lys_patch[0] == True:
        print 'GLU286 patched!'
        charge_patches.append({'GLUP':['GLU-286_ACHA']})
    if glu_asp_lys_patch[1] == True:
        print 'ASP407 patched!'
        charge_patches.append({'ASPP':['ASP-407_ACHA']})
    if glu_asp_lys_patch[2] == True:
        print 'LYS362 patched!'
        charge_patches.append({'LSN':['LYS-362_ACHA']})

    # charge_patches.append({'PAH1': ['PRA-1_EHEM']})


    #CHARGE PACTHES
    for patch in charge_patches:
        # new_patch = cco_config.copy_patch(patch)
        kbp2_settings.patches.append(patch)
        # kbp2_settings.patches_no_autogen.append(patch)

    #BOND PACTHES -> HEME AND COPPER
    for patch in bond_patches:
        # new_patch = cco_config.copy_patch(patch)
        kbp2_settings.patches_no_autogen.append(patch)

    if residues_tt:
        kbp2_settings.set_selected_residues(residues_tt)
    else:
        kbp2_settings.set_excluded_residues(['HSD-102_ACHA', 'HSD-421_ACHA', 'HSD-419_ACHA', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', \
                                         'HSE-260_BCHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'TYR-288_ACHA', 'HSE-284_ACHA'])

    kbp2_settings.preopt = { 'carb_oxi_relax' : False, \
                            'init_die_4' :  False \
                            }

    kbp2_settings.tapbs_bin = '/scratch/scratch/tmeyer/kbplus2/tapbs_1.3_cav_enere'
    # kbp2_settings.remove_folders = 'all'
    kbp2_settings.remove_folders = 'keep'
    kbp2_settings.set_structure(pdb_structure)
    kbp2_settings.apbs_res = 0.3


    if not dolly:
        kbp2_settings.set_workdir(workfolder)
        kbp2.pka_calculation.calc_pkas(kbp2_settings, titratable_definitions)
    else:

        run_folder =  workfolder
        run_folder_tmp = workfolder + jobname + '/'
        done_folder = workfolder
        if not os.path.exists(run_folder_tmp):
            os.mkdir(run_folder_tmp)

        kbp2_set_pkl_file = run_folder_tmp + 'kbp_set_%s.pkl' % jobname
        python_script_file = run_folder_tmp + 'titrate_%s.py' % jobname
        pfile = open(python_script_file, 'w')
        python_script = """
# coding=utf-8

import kbp2
import cPickle as pickle

kbp2_settings_pickle = '%s' """ % kbp2_set_pkl_file + """
kbp2_settings =  pickle.load( open( kbp2_settings_pickle, "rb" ) )
titratable_definitions = kbp2.kbp_tools.parse_titratable_yaml(kbp2_settings.titratable_yaml)
kbp2.pka_calculation.calc_pkas(kbp2_settings, titratable_definitions)
        """
        pfile.write(python_script)
        pfile.close()

        local_dolly_folder = '/public/scratch/jdragelj/kbp2/cco_kbp_%s/' % (jobname)
        kbp2_settings.set_workdir(local_dolly_folder)
        pklfile = open(kbp2_set_pkl_file, 'wb')
        pickle.dump(kbp2_settings, pklfile)
        pklfile.close()
        submitt_kbp2_job(done_folder=done_folder, run_folder=run_folder, local_dolly_folder=local_dolly_folder,
                         jobname=jobname, md_name='cco_kbp')

def submitt_kbp2_job(done_folder, run_folder, local_dolly_folder, jobname, md_name):

    submitt_script_file = run_folder + jobname + '/titr_%s_%s.sh' % (md_name,jobname)
    sfile = open(submitt_script_file, 'w')

    submitt_script = """
#!/bin/csh

echo "Running on `hostname`"

setenv PYTHONPATH /itch/itch/local_python/lib/python2.7/site-packages:/itch/itch/local_python/src

set SDIR=%s%s/"""%(run_folder,jobname) + """

set FDIR=%s/"""%(done_folder) + """

set JDIR=%s"""%local_dolly_folder + """

mkdir -p $JDIR
cp $SDIR/* $JDIR

cd $JDIR

python titrate_%s.py """% jobname + """

cd ..

cp $JDIR/* $FDIR

rm -r $JDIR/*
rm -r $SDIR/*
rmdir $JDIR
rmdir $SDIR
"""


    sfile.write(submitt_script)
    sfile.close()

    import subprocess
    shell = subprocess.Popen('csh\n', \
                             stdin=subprocess.PIPE, \
                             stdout=subprocess.PIPE, \
                             stderr=subprocess.PIPE, \
                             shell=True \
                             )

    shell.stdin.write('qsub -q D61.q,D62.q,D64.q -o %s -e %s %s\n' % (run_folder, run_folder, submitt_script_file) )
    # shell.stdin.write('qsub -q D62.q,D64.q -o %s -e %s %s\n' % (run_folder, run_folder, submitt_script_file) )
    shell.stdin.write('exit\n')



if __name__ == '__main__':

    # ############# CRYSTAL STRUCTURES ##################
    # # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_2gsm/'
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cav_test/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # # pdb_structure = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_2gsm/cco.pdb'
    # pdb_structure = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_PF/md_PF_prd-_xtra6/read_xtra6/cco_prd-_xtra_min2.pdb'
    # jobname = 'crystalPF'
    # state = 'PF'
    # titrate_structure(pdb_structure, workfolder, state, membrane_charge=False,
    #                   cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, False], dolly=False, jobname=jobname)

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_2gsm_cav/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/cco_and_water.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure)
    # pdb.create_struct()
    # new_2 = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
    # pdb_structure_cco = workfolder + '/cco.pdb'
    # new_2.write_pdb(pdb_structure_cco)
    # state = 'Pr'
    # titrate_structure(pdb_structure_cco, workfolder, state, cavities=True)

    # ############## MODELLED OPEN SB PRE CONS MD STRUCTURES ##################
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_model_pull_min/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/move_22_25_min_500_prdh2/cco_moved_prdh2.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb_structure2 = workfolder + '/cco_moved_kb.pdb'
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # new_2 = pdb.copy(segname=['QBH2','PAH2','SWAT','WAT'], exclude=True)
    # new_2.write_pdb(pdb_structure2)
    # state = 'Pr'
    # titrate_structure(pdb_structure2, workfolder, state, membrane_charge=False, cavities=False)

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_model_pull_min_cav/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/move_22_25_min_500_prdh2/cco_moved_mini.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb_structure2 = workfolder + '/cco_moved_kb.pdb'
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # new_2 = pdb.copy(segname=['QBH2','PAH2','SWAT','WAT'], exclude=True)
    # new_2.write_pdb(pdb_structure2)
    # state = 'Pr'
    # titrate_structure(pdb_structure2, workfolder, state, membrane_charge=False, cavities=True)
    #
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_consMD_prdh1/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prdh1/read_cco/cco_prdh1_for_md.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb_structure2 = workfolder + '/cco_kb.pdb'
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # new_2 = pdb.copy(segname=['QBH2','PAH2','SWAT','WAT'], exclude=True)
    # new_2.write_pdb(pdb_structure2)
    # state = 'Pr'
    # titrate_structure(pdb_structure2, workfolder, state, membrane_charge=False, cavities=True)
    #
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_consMD_prdh2/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prdh2/read_cco/cco_prdh2_for_md.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb_structure2 = workfolder + '/cco_kb.pdb'
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # new_2 = pdb.copy(segname=['QBH2','PAH2','SWAT','WAT'], exclude=True)
    # new_2.write_pdb(pdb_structure2)
    # state = 'Pr'
    # titrate_structure(pdb_structure2, workfolder, state, membrane_charge=False, cavities=True)

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_consMD_prd-/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prd-/read_cco/cco_prd-_for_md.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb_structure2 = workfolder + '/cco_kb.pdb'
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # new_2 = pdb.copy(segname=['QBH2','PAH2','SWAT','WAT'], exclude=True)
    # new_2.write_pdb(pdb_structure2)
    # state = 'Pr'
    # titrate_structure(pdb_structure2, workfolder, state, membrane_charge=False, cavities=True)

    ############## DEPROT HIS 334 STRUCTURES ##################
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/frames_sb_open/prdh2/frame0_Pr_His_dep/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/move_22_25_min_500_prdh2/frames_prdh2/frame0.pdb'
    # state = 'Pr_deprot_His334'
    # titrate_structure(pdb_structure, workfolder, state, residues_tt=None)

    ############## MEMBRANE TESTS ##################
    # # esps = [2, 4, 6, 10, 15, 20, 30, 40, 80]
    # esps = [4]
    # for esp in esps:
    #     # workfolder = '/scratch/scratch/jdragelj/tests/mem_test_esp%i_smallmemb/' % esp
    #     # workfolder = '/scratch/scratch/jdragelj/tests/mem_test_esp%i_smallmemb_cav/' % esp
    #     workfolder = '/scratch/scratch/jdragelj/tests/mem_test_esp%i_33tests_cav_3/' % esp
    #     if not os.path.exists(workfolder):
    #         os.mkdir(workfolder)
    #     pdb_structure_cco = '/scratch/scratch/jdragelj/tests/cco.pdb'
    #     state = 'Pr'
    #     residues_tt = ['GLU-286_ACHA']
    #     #todo: read from membrane_thickness.dat
    #     membrane_par = [0.0, 0.0, -16.9, 0.0, 0.0, 33.3, esp] #eps=20.0
    #     # membrane_par = [0.0, 0.0, -10.9, 0.0, 0.0, 21.0, esp] #eps=20.0
    #     # titrate_structure(pdb_structure_cco, workfolder, state, residues_tt=residues_tt, cavities=True, membrane_par=membrane_par)
    #     titrate_structure(pdb_structure_cco, workfolder, state, residues_tt=residues_tt, membrane_par=membrane_par)
    #     f = open(workfolder+'results.dat', 'r')
    #     for line in f:
    #         if line:
    #             print esp, line
    # import time
    # start = time.time()
    # workfolder = '/scratch/scratch/jdragelj/tests/crystal_explicit/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/tests/cco_memb.pdb'
    # state = 'Pr'
    # # residues_tt = ['GLU-286_ACHA', 'LYS-362-ACHA']
    # # titrate_structure(pdb_structure_cco, workfolder, state, residues_tt=residues_tt)
    # titrate_structure(pdb_structure_cco, workfolder, state)
    # end = time.time()
    # print(end - start)/4
    # start = time.time()
    # workfolder = '/scratch/scratch/jdragelj/tests/crystal/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/tests/cco.pdb'
    # state = 'Pr'
    # # residues_tt = ['GLU-286_ACHA', 'LYS-362-ACHA']
    # # titrate_structure(pdb_structure_cco, workfolder, state, residues_tt=residues_tt)
    # titrate_structure(pdb_structure_cco, workfolder, state)
    # end = time.time()
    # print(end - start)/4

    #
    # ############# PRA MINIMIZED STRUCTURES ##################
    # # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/deprot_cav/'
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/deprot/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a/cco_mini_pra-_2.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb_structure2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a/pra-_kb_nomem.pdb'
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # new_2 = pdb.copy(segname=['QBH2','PAH2','SWAT','WAT','MEMB'], exclude=True)
    # new_2.write_pdb(pdb_structure2)
    # state = 'Pr'
    # residues_tt = ['PRA-1_GHEM']
    # # titrate_structure(pdb_structure2, workfolder, state, residues_tt=residues_tt, cavities=True)
    # titrate_structure(pdb_structure2, workfolder, state, residues_tt=residues_tt, cavities=False)
    #
    # # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/prah1_cav/'
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/prah1/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a/cco_mini_prah1_2.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb_structure2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a/prah1_2_kb_nomem.pdb'
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # new_2 = pdb.copy(segname=['QBH2','PAH2','SWAT','WAT','MEMB'], exclude=True)
    # new_2.write_pdb(pdb_structure2)
    # state = 'Pr'
    # residues_tt = ['PRA-1_GHEM']
    # # titrate_structure(pdb_structure2, workfolder, state, residues_tt=residues_tt, cavities=True)
    # titrate_structure(pdb_structure2, workfolder, state, residues_tt=residues_tt, cavities=False)
    #
    # # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/prah2_cav/'
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_pra/prah2/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a/cco_mini_prah2_2.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb_structure2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a/prah2_2_kb_nomem.pdb'
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # new_2 = pdb.copy(segname=['QBH2','PAH2','SWAT','WAT','MEMB'], exclude=True)
    # new_2.write_pdb(pdb_structure2)
    # state = 'Pr'
    # residues_tt = ['PRA-1_GHEM']
    # # titrate_structure(pdb_structure2, workfolder, state, residues_tt=residues_tt, cavities=True)
    # titrate_structure(pdb_structure2, workfolder, state, residues_tt=residues_tt, cavities=False)

    # ############## CONS MD STRUCTURES ##################
    # frames = [0,9]
    # for frame in frames:
    #     workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_frames_cons_md_nt/prd-/frame%i/' % frame
    #     if not os.path.exists(workfolder):
    #         os.mkdir(workfolder)
    #     pdb_structure = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/move_22_25_min_500_prdh2/frames_prd-/frame%i.pdb' % frame
    #     state = 'Pr'
    #     titrate_structure(pdb_structure, workfolder, state, membrane_charge=False)
    #
    #     workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_frames_cons_md_nt/prdh1/frame%i/' % frame
    #     if not os.path.exists(workfolder):
    #         os.mkdir(workfolder)
    #     pdb_structure = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/move_22_25_min_500_prdh2/frames_prdh1/frame%i.pdb' % frame
    #     state = 'Pr'
    #     titrate_structure(pdb_structure, workfolder, state, membrane_charge=False)
    #
    #     workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_frames_cons_md_nt/prdh2/frame%i/' % frame
    #     if not os.path.exists(workfolder):
    #         os.mkdir(workfolder)
    #     pdb_structure = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/move_22_25_min_500_prdh2/frames_prdh2/frame%i.pdb' % frame
    #     state = 'Pr'
    #     titrate_structure(pdb_structure, workfolder, state, membrane_charge=False)



    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_models/Pr_pra_a3/prah_cav/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/move/cco_moved_mini_prah_envir4.pdb'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/move/cco_prah_model_voda.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # for seg in pdb.struct.iter_segments():
    #     for res in seg.iter_residues():
    #         change = False
    #         if res.resname == 'TIP3':
    #             if res.segname == 'QBH2':
    #                 if res.resid in [42, 48, 55]:
    #                     change = True
    #             if change:
    #                 for atm in res.iter_atoms():
    #                     atm['segname'] = 'VODA'
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'Pr'
    # residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #             'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, cavities=True)
    #
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_models/Pr_pra_a3/pra2_cav/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/move/cco_moved_mini_pra2_envir4.pdb'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/move/cco_pra2_model_voda.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # for seg in pdb.struct.iter_segments():
    #     for res in seg.iter_residues():
    #         change = False
    #         if res.resname == 'TIP3':
    #             if res.segname == 'QBH2':
    #                 if res.resid in [42, 48, 55]:
    #                     change = True
    #             if change:
    #                 for atm in res.iter_atoms():
    #                     atm['segname'] = 'VODA'
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'Pr'
    # residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #             'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, cavities=True)

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_models/Pr_pra_a3/pra-_cav/'
    # if not os.path.exists(workfolder):
    #     os.mkdir(workfolder)
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/move/cco_moved_mini_pra-_envir4.pdb'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/move/cco_pra-_model_voda.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # for seg in pdb.struct.iter_segments():
    #     for res in seg.iter_residues():
    #         change = False
    #         if res.resname == 'TIP3':
    #             if res.segname == 'QBH2':
    #                 if res.resid in [42, 48, 55]:
    #                     change = True
    #             if change:
    #                 for atm in res.iter_atoms():
    #                     atm['segname'] = 'VODA'
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'Pr'
    # residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #             'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, cavities=True)

    # main_workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/'
    # states = ['F', 'Pm', 'PF', 'Pr']
    # prots = ['prd-', 'prdh1', 'prdh2']
    # for state in states:
    #     for prot in prots:
    #         workfolder = main_workfolder + '%s/prda3/%s/' % (state, prot)
    #         if not os.path.exists(workfolder):
    #             os.mkdir(workfolder)
    #         pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_test_models/%s/cco_%s.pdb' % (state, prot)
    #         pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_test_models/%s/cco_%s_voda.pdb' % (state, prot)
    #         pdb = kbp2.file_parser.Simple_struct_parser()
    #         pdb.read_pdb(pdb_structure_cco)
    #         pdb.create_struct()
    #         for seg in pdb.struct.iter_segments():
    #             for res in seg.iter_residues():
    #                 change = False
    #                 if res.resname == 'TIP3':
    #                     if res.segname == 'QBH2':
    #                         if res.resid in [42, 48, 55]:
    #                             change = True
    #                     if change:
    #                         for atm in res.iter_atoms():
    #                             atm['segname'] = 'VODA'
    #         pdb.create_struct()
    #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    #         pdb_w.write_pdb(pdb_structure_cco2)
    #         state = state
    #         residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #                     'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    #         jobname = 'cav09_prda3_%s_%s' % (state,prot)
    #         titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                           cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, True], dolly=True, jobname=jobname)
    #
    # main_workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/'
    # states = ['F', 'Pm', 'PF', 'Pr']
    # prots = ['pra-', 'prah1', 'prah2']
    # for state in states:
    #     for prot in prots:
    #         workfolder = main_workfolder + '%s/praa/%s/' % (state, prot)
    #         if not os.path.exists(workfolder):
    #             if not os.path.exists(workfolder):
    #                 os.mkdir(workfolder)
    #         pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a/test_models/%s/cco_%s.pdb' % (state, prot)
    #         pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a/test_models/%s/cco_%s_voda.pdb' % (state, prot)
    #         pdb = kbp2.file_parser.Simple_struct_parser()
    #         pdb.read_pdb(pdb_structure_cco)
    #         pdb.create_struct()
    #         for seg in pdb.struct.iter_segments():
    #             for res in seg.iter_residues():
    #                 change = False
    #                 if res.resname == 'TIP3':
    #                     if res.segname == 'QBH2':
    #                         if res.resid in [42, 48, 55]:
    #                             change = True
    #                     if change:
    #                         for atm in res.iter_atoms():
    #                             atm['segname'] = 'VODA'
    #         pdb.create_struct()
    #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2','XT2', 'XTC'], exclude=True)
    #         pdb_w.write_pdb(pdb_structure_cco2)
    #         state = state
    #         residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #                     'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    #         jobname = 'cav09_praa_%s_%s' % (state,prot)
    #         titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                           cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, True], dolly=True, jobname=jobname)

    # main_workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/'
    # states = ['F', 'Pm', 'PF', 'Pr']
    # prots = ['pra-', 'prah1', 'prah2']
    # for state in states:
    #     for prot in prots:
    #         workfolder = main_workfolder + '%s/praa3/%s/' % (state, prot)
    #         if not os.path.exists(workfolder):
    #             if not os.path.exists(workfolder):
    #                 os.mkdir(workfolder)
    #         pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/test_models/%s/cco_%s.pdb' % (state, prot)
    #         pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/test_models/%s/cco_%s_voda.pdb' % (state, prot)
    #         pdb = kbp2.file_parser.Simple_struct_parser()
    #         pdb.read_pdb(pdb_structure_cco)
    #         pdb.create_struct()
    #         for seg in pdb.struct.iter_segments():
    #             for res in seg.iter_residues():
    #                 change = False
    #                 if res.resname == 'TIP3':
    #                     if res.segname == 'QBH2':
    #                         if res.resid in [42, 48, 55]:
    #                             change = True
    #                     if change:
    #                         for atm in res.iter_atoms():
    #                             atm['segname'] = 'VODA'
    #         pdb.create_struct()
    #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2','XT2', 'XTC'], exclude=True)
    #         pdb_w.write_pdb(pdb_structure_cco2)
    #         state = state
    #         residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #                     'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    #         jobname = 'cav09_praa3_%s_%s' % (state,prot)
    #         titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                           cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, True], dolly=True, jobname=jobname)

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/tests/explicit_wat/'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prd-_xtra6/test_frame240.pdb'
    # state = 'Pr'
    # residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #             'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False, cavities=True, glu_asp_lys_patch=[False, False, True])

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_glu_model/7xtra/'
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pr_prd-_xtra_up/read/7xtra/cco_prd-_xtra_min2.pdb'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pr_prd-_xtra_up/read/7xtra/cco_kbp.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # for seg in pdb.struct.iter_segments():
    #     for res in seg.iter_residues():
    #         change = False
    #         if res.resname == 'TIP3':
    #             if res.segname == 'QBH2':
    #                 if res.resid in [42, 48, 55]:
    #                     change = True
    #             if change:
    #                 for atm in res.iter_atoms():
    #                     atm['segname'] = 'VODA'
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'Pm'
    # residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #             'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # jobname = 'glu_up_open_sb_Pm_7xtra'
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                   cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, True], dolly=True, jobname=jobname)

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_glu_model/7glu_deprot/'
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pr_prd-_xtra_up/read/7xtra/cco_prd-_xtra_min2.pdb'
    # pdb_structure_cco2 = workfolder + '/cco_glu-_kbp.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # for seg in pdb.struct.iter_segments():
    #     for res in seg.iter_residues():
    #         change = False
    #         if res.resname == 'TIP3':
    #             if res.segname == 'QBH2':
    #                 if res.resid in [42, 48, 55]:
    #                     change = True
    #             if change:
    #                 for atm in res.iter_atoms():
    #                     atm['segname'] = 'VODA'
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'Pm'
    # residues = ['PRD-3_EHEM']
    # jobname = 'glu_deprot_up_open_sb_Pm_7xtra'
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                   cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, True, True], dolly=True, jobname=jobname)
    #
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_glu_model/glu_prot_down/'
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prd-_xtra6/read_xtra6/cco_prd-_xtra_min2.pdb'
    # pdb_structure_cco2 = workfolder + 'cco_glu-_kbp.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # for seg in pdb.struct.iter_segments():
    #     for res in seg.iter_residues():
    #         change = False
    #         if res.resname == 'TIP3':
    #             if res.segname == 'QBH2':
    #                 if res.resid in [42, 48, 55]:
    #                     change = True
    #             if change:
    #                 for atm in res.iter_atoms():
    #                     atm['segname'] = 'VODA'
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'Pm'
    # residues = ['PRD-3_EHEM']
    # jobname = 'glu_prot_down_open_sb_Pm_7xtra'
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                   cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[True, True, True], dolly=True, jobname=jobname)

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_glu_model/9xtra/'
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pr_prd-_xtra_up/read/9xtra/cco_prd-_xtra_min2.pdb'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pr_prd-_xtra_up/read/9xtra/cco_kbp.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # for seg in pdb.struct.iter_segments():
    #     for res in seg.iter_residues():
    #         change = False
    #         if res.resname == 'TIP3':
    #             if res.segname == 'QBH2':
    #                 if res.resid in [42, 48, 55]:
    #                     change = True
    #             if change:
    #                 for atm in res.iter_atoms():
    #                     atm['segname'] = 'VODA'
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'Pm'
    # residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #             'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # jobname = 'glu_up_open_sb_Pm_9xtra'
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                   cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, True], dolly=True, jobname=jobname)


    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_glu_model/crystal_PF/'
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_Pr_crystal_xtra4/read_xtra4/cco_and_water_xtra_out.pdb'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/cco_PF_kbp.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # for seg in pdb.struct.iter_segments():
    #     for res in seg.iter_residues():
    #         change = False
    #         if res.resname == 'TIP3':
    #             if res.segname == 'QBH2':
    #                 if res.resid in [42, 48, 55]:
    #                     change = True
    #             if change:
    #                 for atm in res.iter_atoms():
    #                     atm['segname'] = 'VODA'
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'PF'
    # residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #             'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False, cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, True])


    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_glu_model/crystal_Pr/'
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_Pr_crystal_xtra4/read_xtra4/cco_and_water_xtra_out.pdb'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_Pr_crystal_xtra4/cco_kbp.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # for seg in pdb.struct.iter_segments():
    #     for res in seg.iter_residues():
    #         change = False
    #         if res.resname == 'TIP3':
    #             if res.segname == 'QBH2':
    #                 if res.resid in [42, 48, 55]:
    #                     change = True
    #             if change:
    #                 for atm in res.iter_atoms():
    #                     atm['segname'] = 'VODA'
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'Pr'
    # residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #             'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False, cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, True])
    #
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pm_glu_model/crystal_Pm/'
    # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_Pr_crystal_xtra4/read_xtra4/cco_and_water_xtra_out.pdb'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/cco_Pm_kbp.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure_cco)
    # pdb.create_struct()
    # for seg in pdb.struct.iter_segments():
    #     for res in seg.iter_residues():
    #         change = False
    #         if res.resname == 'TIP3':
    #             if res.segname == 'QBH2':
    #                 if res.resid in [42, 48, 55]:
    #                     change = True
    #             if change:
    #                 for atm in res.iter_atoms():
    #                     atm['segname'] = 'VODA'
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC''], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'Pm'
    # residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
    #             'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False, cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, True])

    # for state in ['F','PF','Pm','Pr']:
    #
    #     workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/glu286_influence_pra/gluH_up_%s/' % state
    #     if not os.path.exists(workfolder):
    #         os.mkdir(workfolder)
    #     if state == 'F':
    #         pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pm_prd-_xtra_up/read/7xtra/cco_prd-_xtra_min2_F.pdb'
    #     else:
    #         pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pm_prd-_xtra_up/read/7xtra/cco_prd-_xtra_min2.pdb'
    #     pdb_structure_cco2 = workfolder + '/cco_gluH_kbp.pdb'
    #     pdb = kbp2.file_parser.Simple_struct_parser()
    #     pdb.read_pdb(pdb_structure_cco)
    #     pdb.create_struct()
    #     for seg in pdb.struct.iter_segments():
    #         for res in seg.iter_residues():
    #             change = False
    #             if res.resname == 'TIP3':
    #                 if res.segname == 'QBH2':
    #                     if res.resid in [42, 48, 55]:
    #                         change = True
    #                 if change:
    #                     for atm in res.iter_atoms():
    #                         atm['segname'] = 'VODA'
    #     pdb.create_struct()
    #     pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    #     pdb_w.write_pdb(pdb_structure_cco2)
    #     state = state
    #     # residues = ['PRD-3_EHEM']
    #     residues = ['PRA-1_EHEM']
    #     jobname = 'gluH_up_open_sb_%s' % state
    #     titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                       cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[True, True, True], dolly=False, jobname=jobname)
    #
    #     workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/glu286_influence_pra/glu-_up_%s/' % state
    #     if not os.path.exists(workfolder):
    #         os.mkdir(workfolder)
    #     if state == 'F':
    #         pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pm_prd-_xtra_up/read/7xtra/cco_prd-_xtra_min2_F.pdb'
    #     else:
    #         pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pm_prd-_xtra_up/read/7xtra/cco_prd-_xtra_min2.pdb'
    #     pdb_structure_cco2 = workfolder + '/cco_glu-_kbp.pdb'
    #     pdb = kbp2.file_parser.Simple_struct_parser()
    #     pdb.read_pdb(pdb_structure_cco)
    #     pdb.create_struct()
    #     for seg in pdb.struct.iter_segments():
    #         for res in seg.iter_residues():
    #             change = False
    #             if res.resname == 'TIP3':
    #                 if res.segname == 'QBH2':
    #                     if res.resid in [42, 48, 55]:
    #                         change = True
    #                 if change:
    #                     for atm in res.iter_atoms():
    #                         atm['segname'] = 'VODA'
    #     pdb.create_struct()
    #     pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    #     pdb_w.write_pdb(pdb_structure_cco2)
    #     state = state
    #     # residues = ['PRD-3_EHEM']
    #     residues = ['PRA-1_EHEM']
    #     jobname = 'glu-_up_open_sb_%s' % state
    #     titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                       cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, True, True], dolly=False, jobname=jobname)
    #
    #     workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/glu286_influence_pra/gluH_down_%s/' % state
    #     if not os.path.exists(workfolder):
    #         os.mkdir(workfolder)
    #     pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prd-_xtra6/read_xtra6/cco_prd-_xtra_min2.pdb'
    #     if state == 'F':
    #         pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_F/md_F_prd-_xtra6/read_cco/cco_prd-.pdb'
    #     pdb_structure_cco2 = workfolder + 'cco_gluH_kbp.pdb'
    #     pdb = kbp2.file_parser.Simple_struct_parser()
    #     pdb.read_pdb(pdb_structure_cco)
    #     pdb.create_struct()
    #     for seg in pdb.struct.iter_segments():
    #         for res in seg.iter_residues():
    #             change = False
    #             if res.resname == 'TIP3':
    #                 if res.segname == 'QBH2':
    #                     if res.resid in [42, 48, 55]:
    #                         change = True
    #                 if change:
    #                     for atm in res.iter_atoms():
    #                         atm['segname'] = 'VODA'
    #     pdb.create_struct()
    #     pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    #     pdb_w.write_pdb(pdb_structure_cco2)
    #     state = state
    #     # residues = ['PRD-3_EHEM']
    #     residues = ['PRA-1_EHEM']
    #     jobname = 'gluH_down_open_sb_%s' % state
    #     titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                       cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[True, True, True], dolly=False, jobname=jobname)
    #
    #     workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/glu286_influence_pra/glu-_down_%s/' % state
    #     if not os.path.exists(workfolder):
    #         os.mkdir(workfolder)
    #     pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prd-_xtra6/read_xtra6/cco_prd-_xtra_min2.pdb'
    #     if state == 'F':
    #         pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_F/md_F_prd-_xtra6/read_cco/cco_prd-.pdb'
    #     pdb_structure_cco2 = workfolder + 'cco_glu-_kbp.pdb'
    #     pdb = kbp2.file_parser.Simple_struct_parser()
    #     pdb.read_pdb(pdb_structure_cco)
    #     pdb.create_struct()
    #     for seg in pdb.struct.iter_segments():
    #         for res in seg.iter_residues():
    #             change = False
    #             if res.resname == 'TIP3':
    #                 if res.segname == 'QBH2':
    #                     if res.resid in [42, 48, 55]:
    #                         change = True
    #                 if change:
    #                     for atm in res.iter_atoms():
    #                         atm['segname'] = 'VODA'
    #     pdb.create_struct()
    #     pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
    #     pdb_w.write_pdb(pdb_structure_cco2)
    #     state = state
    #     # residues = ['PRD-3_EHEM']
    #     residues = ['PRA-1_EHEM']
    #     jobname = 'glu-_down_open_sb_%s' % state
    #     titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
    #                       cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, True, True], dolly=False, jobname=jobname)

        # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/glu286_influence/glu_up_%s/' % state
        # if not os.path.exists(workfolder):
        #     os.mkdir(workfolder)
        # if state == 'F':
        #     pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pm_prd-_xtra_up/read/7xtra/cco_prd-_xtra_min2_F.pdb'
        # else:
        #     pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pm_prd-_xtra_up/read/7xtra/cco_prd-_xtra_min2.pdb'
        # pdb_structure_cco2 = workfolder + '/cco_glu_kbp.pdb'
        # pdb = kbp2.file_parser.Simple_struct_parser()
        # pdb.read_pdb(pdb_structure_cco)
        # pdb.create_struct()
        # for seg in pdb.struct.iter_segments():
        #     for res in seg.iter_residues():
        #         change = False
        #         if res.resname == 'TIP3':
        #             if res.segname == 'QBH2':
        #                 if res.resid in [42, 48, 55]:
        #                     change = True
        #             if change:
        #                 for atm in res.iter_atoms():
        #                     atm['segname'] = 'VODA'
        # pdb.create_struct()
        # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
        # pdb_w.write_pdb(pdb_structure_cco2)
        # state = state
        # residues = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM', 'ARG-481_ACHA', 'ARG-482_ACHA',
        #             'GLU-286_ACHA', \
        #             'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA',
        #             'LYS-227_BCHA', 'TYR-336_ACHA', \
        #             'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
        # jobname = 'glu_up_open_sb_%s' % state
        # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
        #                   cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, True], dolly=False, jobname=jobname)


        workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cav_test_no/'
        pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_PF/md_PF_prd-_xtra6/read_xtra6/cco_prd-_xtra_min2.pdb'
        pdb_structure_cco2 = workfolder + '/cco_kbp.pdb'
        pdb = kbp2.file_parser.Simple_struct_parser()
        pdb.read_pdb(pdb_structure_cco)
        pdb.create_struct()
        for seg in pdb.struct.iter_segments():
            for res in seg.iter_residues():
                change = False
                if res.resname == 'TIP3':
                    if res.segname == 'QBH2':
                        if res.resid in [42, 48, 55]:
                            change = True
                    if change:
                        for atm in res.iter_atoms():
                            atm['segname'] = 'VODA'
        pdb.create_struct()
        pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
        pdb_w.write_pdb(pdb_structure_cco2)
        state = "PF"
        residues = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM', 'ARG-481_ACHA', 'ARG-482_ACHA',
                    'GLU-286_ACHA', \
                    'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA',
                    'LYS-227_BCHA', 'TYR-336_ACHA', \
                    'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
        jobname = 'glu_up_open_sb_%s' % state
        titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
                          cavities=[], glu_asp_lys_patch=[False, False, True], dolly=False, jobname=jobname)


        # residues = ['PRD-3_EHEM',
        #             'PRD-3_GHEM',
        #             'PRA-1_GHEM',
        #             'GLU-286_ACHA','ARG-481_ACHA', 'ARG-482_ACHA',\
        #             'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA',
        #             'LYS-227_BCHA', 'TYR-336_ACHA', \
        #             'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
        # pdb_structure_cco = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/md_PF_prah_asp_dep/frames_voda/frame0.pdb'
        # pdb_structure_cco2 = '/scratch/scratch/jdragelj/tests/pls_test/cco.pdb'
        # pdb = kbp2.file_parser.Simple_struct_parser()
        # pdb.read_pdb(pdb_structure_cco)
        # pdb.create_struct()
        # for seg in pdb.struct.iter_segments():
        #     for res in seg.iter_residues():
        #         change = False
        #         if res.resname == 'TIP3':
        #             if res.segname == 'QBH2':
        #                 if res.resid in [42, 48, 55]:
        #                     change = True
        #             if change:
        #                 for atm in res.iter_atoms():
        #                     atm['segname'] = 'VODA'
        # pdb.create_struct()
        # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XT2', 'XTC'], exclude=True)
        # pdb_w.write_pdb(pdb_structure_cco2)
        # workfolder = '/scratch/scratch/jdragelj/tests/pls_test/'
        # if not os.path.exists(workfolder):
        #     os.mkdir(workfolder)
        # state = 'PF'
        # jobname = 'test'
        # titrate_structure(pdb_structure_cco2, workfolder, state, residues_tt=residues, membrane_charge=False,
        #                   cavities=[0.9, 0.2, 0.0], glu_asp_lys_patch=[False, False, True], dolly=False, jobname=jobname)


