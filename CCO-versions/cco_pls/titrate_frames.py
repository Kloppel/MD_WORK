# coding=utf-8

import kbp2
import os
import shutil
import numpy as np


def titrate_structure(workfolder, frames_folder, frame_range, state, residues_tt=None, cavities=False, membrane_charge=True, kbp2_eval_titration=False):

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
        kbp2_settings.cavity_par = [0.8, 0.2, 0.0]

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
    if state == 'Pr_deprot_His334':
        print 'histidines inverted on purpose'
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                        'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'A3H4': ['HEM-2_EHEM', 'HSD-419_ACHA', 'HAO2-1_FEOH']})
        charge_patches.append({'CBDD': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA','HSE-284_ACHA']})

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
    print 'LYS362 patched!'
    charge_patches.append({'LSN': ['LYS-362_ACHA']})
    # print 'GLU286 patched!'
    # charge_patches.append({'GLUP': ['GLU-286_ACHA']})
    # print 'ASP407 patched!'
    # charge_patches.append({'ASPP':['ASP-407_ACHA']})

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

    # kbp2_settings.set_structure(pdb_structure)
    # kbp2_settings.set_workdir(workfolder)
    # kbp2.pka_calculation.calc_pkas(kbp2_settings, titratable_definitions)

    frames = []
    for frame in frames_to_do:
        if 'voda' not in frames_folder:
            frames.append(frames_folder+'frame%i.pdb'%frame)
        else:
            if os.path.exists(frames_folder + 'mgcw/frame%i.pdb'%frame):
                frames.append(frames_folder + 'mgcw/frame%i.pdb' % frame)
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

    titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_voda/'

    frames_suffix = 'frames_voda/'
    # frames_suffix = 'frames'

    residues = ['PRD-3_EHEM','PRA-1_EHEM','PRD-3_GHEM','PRA-1_GHEM','ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', 'ASP-407_ACHA', 'HSD-411_ACHA', \
                'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', 'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    print len(residues)

    workfolder = titration_folder + 'md_Pr_prdh1/'
    if not os.path.exists(workfolder):
        os.mkdir(workfolder)
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prdh1_xtra6/%s/' % frames_suffix
    frame_range = np.arange(0, 801, 2)
    state = 'Pr'
    titrate_structure(workfolder, frames_folder, frame_range, state=state, residues_tt=residues, membrane_charge=False, kbp2_eval_titration=True)

    workfolder = titration_folder + 'md_Pr_prd-/'
    if not os.path.exists(workfolder):
        os.mkdir(workfolder)
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prd-_xtra6/%s/' % frames_suffix
    frame_range = np.arange(0, 801, 2)
    state = 'Pr'
    titrate_structure(workfolder, frames_folder, frame_range, state=state, residues_tt=residues, membrane_charge=False, kbp2_eval_titration=True)

    workfolder = titration_folder + 'md_Pr_prdh2/'
    if not os.path.exists(workfolder):
        os.mkdir(workfolder)
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prdh2_xtra6/%s/' % frames_suffix
    frame_range = np.arange(0, 801, 2)
    state = 'Pr'
    titrate_structure(workfolder, frames_folder, frame_range, state=state, residues_tt=residues, membrane_charge=False, kbp2_eval_titration=True)

    ################################
    titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/md_crystal/'
    frames_suffix = 'frames_voda/'

    workfolder = titration_folder + '10ang_voda/'
    if not os.path.exists(workfolder):
        os.mkdir(workfolder)
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_crystal_xtra4/%s/' % frames_suffix
    frame_range = np.arange(0, 801, 2)
    state = 'Pr'
    titrate_structure(workfolder, frames_folder, frame_range, state=state, residues_tt=residues, membrane_charge=False, kbp2_eval_titration=True)



