# coding=utf-8

import kbp2
import os
import shutil
import numpy as np


def titrate_structure(workfolder, frames_folder, frame_range, state, residues_tt=None, cavities=[], membrane_charge=True, kbp2_eval_titration=False,
                      glu_asp_lys_patch=[True, True, True]):

    print 'Starting work in %s' % workfolder

    run_folder = workfolder + 'run/'
    done_folder = workfolder + 'done/'

    if not os.path.exists(workfolder):
        os.mkdir(workfolder)
    if not os.path.exists(run_folder):
        os.mkdir(run_folder)
    if not os.path.exists(done_folder):
        os.mkdir(done_folder)

    frames_to_do = []
    if os.path.exists(workfolder):
        if os.path.exists(done_folder):
            for frame in frame_range:
                if not os.path.exists(done_folder + 'frame%s/' % frame):
                    if os.path.exists(run_folder + 'frame%s/' % frame):
                        shutil.rmtree(run_folder + 'frame%s/' % frame)
                        frames_to_do.append(frame)
                    else:
                        frames_to_do.append(frame)
                else:
                    continue
        else:
            frames_to_do = frame_range

    top = []
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw_clean_kb.inp")
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
        # kbp2_settings.cavity_par = [0.8, 0.2, 0.0]
        # kbp2_settings.cavity_par = [1.0, 0.2, 0.0]
        kbp2_settings.cavity_par = cavities

    if kbp2_eval_titration:
        kbp2_settings.md_evaluation_mode = True
        titratable_yaml = '/scratch/scratch/jdragelj/KB_files/titratable_cco_extended.yaml'
    else:
        titratable_yaml = '/scratch/scratch/jdragelj/KB_files/titratable_cco_noglupasp.yaml'

    kbp2_settings.set_yaml(titratable_yaml)
    titratable_definitions = kbp2.kbp_tools.parse_titratable_yaml(titratable_yaml)
    ##############################################################################

    kbp2_settings.patches = []
    kbp2_settings.patches_no_autogen = []
    charge_patches = []
    bond_patches = []

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
        bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
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

    #CHARGE PACTHES
    for patch in charge_patches:
        kbp2_settings.patches.append(patch)

    #BOND PACTHES -> HEME AND COPPER
    for patch in bond_patches:
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
    kbp2_settings.remove_folders = 'all'

    frames = []
    for frame in frames_to_do:
        if 'voda' not in frames_folder:
            frames.append(frames_folder+'frame%i.pdb'%frame)
        else:
            if not os.path.exists(frames_folder + 'mgcw/'):
                os.mkdir(frames_folder + 'mgcw/')
            pdb = kbp2.file_parser.Simple_struct_parser()
            pdb.read_pdb(frames_folder+'frame%i.pdb'%frame)
            pdb.create_struct()
            for seg in pdb.struct.iter_segments():
                    for res in seg.iter_residues():
                            change = False
                            if res.resname == 'TIP3':
                                    if res.segname == 'QBH2':
                                            if res.resid in [42,48,55]:
                                                    change = True
                                    if change:
                                            for atm in res.iter_atoms():
                                                    atm['segname'] = 'VODA'
            pdb.create_struct()
            pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2'], exclude=True)
            pdb_w.write_pdb(frames_folder + 'mgcw/frame%i.pdb'%frame)
            frames.append(frames_folder + 'mgcw/frame%i.pdb'%frame)


    # kbp2_settings.remove_folders = 'keep'
    # kbp2_settings.force_conf_ene_calc = True
    # kbp2_settings.set_structure(frames[0])
    # kbp2_settings.set_workdir(workfolder)
    # kbp2.pka_calculation.calc_pkas(kbp2_settings, titratable_definitions)

    from time import sleep
    manager = kbp2.pka_manager.Pka_manager(run_folder, done_folder, kbp2_settings, titratable_definitions, cpus=2)
    manager.submit_structures(frames)
    while True:
        (done, crashed, running, queued) = manager.get_status()
        # print'done=', done, 'crashed=', crashed, 'running=', running, 'queued=', queued
        if running + queued == 0:
            break
        sleep(60)
    print '--------------------------------'


if __name__ == '__main__':

    print'PRD a3'

    chunk80 = np.arange(0, 801, 80)
    chunk40 = np.arange(0, 801, 40)
    chunk20 = np.arange(0, 801, 20)
    chunk10 = np.arange(0, 801, 10)
    chunk4 = np.arange(0, 801, 4)
    chunk2 = np.arange(0, 801, 2)
    frame_range = chunk80

    # frame_range = [240]

    for md  in ['md_Pr_prd-','md_Pr_prdh1','md_Pr_prdh2']:

        titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/no_cavities_voda/'
        frames_suffix = 'frames_voda'
        residues = ['PRD-3_EHEM']
        workfolder = titration_folder + '%s/' % md
        if not os.path.exists(workfolder):
            os.mkdir(workfolder)
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/%s_xtra6/%s/' % (md,frames_suffix)
        state = 'Pr'
        titrate_structure(workfolder, frames_folder, frame_range, state, residues_tt = residues, cavities = [], membrane_charge = False,
                          kbp2_eval_titration = True, glu_asp_lys_patch=[True, True, True])

        titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/cavities_09_voda/'
        frames_suffix = 'frames_voda'
        residues = ['PRD-3_EHEM']
        workfolder = titration_folder + '%s/' % md
        if not os.path.exists(workfolder):
            os.mkdir(workfolder)
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/%s_xtra6/%s/' % (md,frames_suffix)
        state = 'Pr'
        titrate_structure(workfolder, frames_folder, frame_range, state, residues_tt = residues, cavities = [0.9, 0.2, 0.0], membrane_charge = False,
                          kbp2_eval_titration = True, glu_asp_lys_patch=[True, True, True])

        titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/cavities_1_voda/'
        frames_suffix = 'frames_voda'
        residues = ['PRD-3_EHEM']
        workfolder = titration_folder + '%s/' % md
        if not os.path.exists(workfolder):
            os.mkdir(workfolder)
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/%s_xtra6/%s/' % (md,frames_suffix)
        state = 'Pr'
        titrate_structure(workfolder, frames_folder, frame_range, state, residues_tt = residues, cavities = [1.0, 0.2, 0.0], membrane_charge = False,
                          kbp2_eval_titration = True, glu_asp_lys_patch=[True, True, True])




