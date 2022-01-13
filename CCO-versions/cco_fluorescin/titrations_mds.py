# coding=utf-8
# #this works:
# import sys
# sys.path.append('/user/jdragelj/python/karlsbergplus/package/kbp2')

import os
import shutil
import kbp2
from kbp2.workspace_jd import pka_manager_jd
from workspace_jd import cco_config
from time import sleep
import numpy as np

def titrate_frames_md(folder, workfolder, subfolder_list, frames_folder_name, sleep_time, frame_range, residues_tt, md_folder = 'md', special=None):

    print('Set conf energy to zero in kbp_results!')

    for subfolder in subfolder_list:

        ###############
        ### FOLDERS ###
        ###############

        print '--------------------------------'
        print 'Starting work in %s' % subfolder


        work_subfolder = workfolder + subfolder + '/'
        run_folder = work_subfolder + 'run/'
        done_folder = work_subfolder + 'done/'

        if not os.path.exists(work_subfolder):
            os.mkdir(work_subfolder)
        if not os.path.exists(run_folder):
            os.mkdir(run_folder)
        if not os.path.exists(done_folder):
            os.mkdir(done_folder)

        frames_to_do = []
        if os.path.exists(work_subfolder):
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

        # with membrane
        sourcefolder = folder + subfolder + '/' + md_folder + '/%s/' % frames_folder_name
        print sourcefolder


        #############
        ### CALCS ###
        #############

        top = []
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp")
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.rtf")
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf")
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_kbp/top_kbp_fluorescin_scp.rtf")

        par = []
        par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
        par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
        par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
        par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")
        par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_fluorescin.str")
        par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_kbp/par_kbp_fluorescin.prm")

        # kbp2_settings = kbp2.pka_calculation.PkaCalcSettings()
        kbp2_settings = kbp2.pka_cal_jovan.PkaCalcSettings()
        kbp2_settings.top = top
        kbp2_settings.par = par
        kbp2_settings.processes = 1
        kbp2_settings.quiet_mode = True

        decisions = []
        decisions.append(('rename__HOH_TIP3', 'keep'))
        decisions.append(('disu_bridges', 'closed'))

        kbp2_settings.modelling_decisions = decisions

        #keep or not keep MD hydrogens and cap or not cap termini
        kbp2_settings.init_modelling_min = False
        # kbp2_settings.init_modelling_min = True
        # cap or do not cap termini
        # kbp2_settings.tmp_settings['no_cap'] = True

        protocol = []
        protocol.append((  7, 'h_min'))
        kbp2_settings.protocol = protocol

        kbp2_settings.md_evaluation_mode = True
        titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2_np.yaml'
        kbp2_settings.set_yaml(titratable_yaml)

        titratable_definitions = kbp2.kbp_tools.parse_titratable_yaml(titratable_yaml)
        ##############################################################################

                        #########
                        ### 1 ###        # 1. reduces of oxidized enzime
                        #########

        if 'red' in subfolder:
            cco_state = 'R'
            cco_settings = cco_config.COX_settings(specie='A_paraccocus', cycle_state='R')
        else:
            cco_state = 'O'
            cco_settings = cco_config.COX_settings(specie='A_paraccocus', cycle_state='O')


        if special == 'flip_charge':
            if cco_state == 'R':
                cco_settings = cco_config.COX_settings(specie='A_paraccocus', cycle_state='O')
            else:
                cco_settings = cco_config.COX_settings(specie='A_paraccocus', cycle_state='R')

        if special == 'flip_charge_Cua_hemea':
            if cco_state == 'R':
                cco_settings = cco_config.COX_settings(specie='A_paraccocus', cycle_state='Cua_hemea_O_hemea3_Cub_R')
            else:
                cco_settings = cco_config.COX_settings(specie='A_paraccocus', cycle_state='Cua_hemea_R_hemea3_Cub_O')



                        #########
                        ### 2 ###       2. histidines should be as they are
                        #########

        kbp2_settings.set_ignored_residues(['histidines'])

                        #########
                        ### 3 ###       3. lysine protonation state 0/+
                        #########

        if 'lys' in subfolder:
            lysine_neutral = False
            raise AssertionError('In all our calculations, lysine should be neutral! Check your modelling files.')
        else :
            lysine_neutral = True

    #PROTONATION PATCHES
        if lysine_neutral:
            initial_protonation = {('EPP', 278, 'ACHA'):3, ('LYS', 354, 'ACHA'):1, ('DPP', 399, 'ACHA'):3, ('EPP', 481, 'ACHA'):3}
        else:
           initial_protonation = {('EPP', 278, 'ACHA'):3, ('LYS', 354, 'ACHA'):2, ('DPP', 399, 'ACHA'):3, ('EPP', 481, 'ACHA'):3}

                        #########
                        ### 4 ###   4. fluorescin present or not
                        #########

        kbp2_settings.patches = []
        rename = False
        if  'flu' in subfolder:

            resname = "FLU"
            states = []
            state = {'pka': 0.0,
                     'name':  'R',
                     'external_patches':[('FLUN', ['CYS-299_ACHA', 'FLU-1_FLUR'])]}
            states.append(state)
            state = {'pka': -100.0,
                     'name':  '0',
                     'external_patches':[('FLUH', ['CYS-299_ACHA', 'FLU-1_FLUR'])]}
            states.append(state)
            state = {'pka': -93.8,
                     'name': 'D',
                     'external_patches': [('FLUX', ['CYS-299_ACHA', 'FLU-1_FLUR'])]}
            states.append(state)
            state = {'pka': -100.0,
                     'name': '0',
                     'external_patches': [('FLU2', ['CYS-299_ACHA', 'FLU-1_FLUR'])]}
            states.append(state)

            # titratable_definition = kbp2.pka_calculation.get_atom_charges(resname, states, top)
            titratable_definition = kbp2.pka_cal_jovan.get_atom_charges(resname, states, top)
            titratable_definition['FLU'][0]['external_atoms'][0][1]['H6'] = 0.001
            titratable_definitions.update(titratable_definition)

            #GENERAL PATCHES
            kbp2_settings.patches.append({'FCYS': ['CYS-299_ACHA', 'FLU-1_FLUR']})

            rename_flu_deprot = False
            rename_flu_prot = False
            rename = False

            if 'flu-' in subfolder:

                # renaming must be done from FLX to FLU
                print 'Renaming from FLX to FLU'
                rename_flu_deprot = True
                rename = True

                        #########
                        ### 5 ###   5. if present what is the initial protonaton state of fluorescin
                        #########

                initial_protonation.update({('FLU', 1, 'FLUR'):2})

            if 'fluh' in subfolder:
                # print'*'

                initial_protonation.update({('FLU', 1, 'FLUR'):1})

            if 'flu2' in subfolder:
                # renaming must be done from FLX to FLU
                print 'Renaming from FLH to FLU'
                rename_flu_prot = True
                rename = True

                        #########
                        ### 5 ###   5. if present what is the initial protonaton state of fluorescin
                        #########

                initial_protonation.update({('FLU', 1, 'FLUR'):3})

        kbp2_settings.set_initial_protonation(initial_protonation)

        #CHARGE PACTHES
        for patch in cco_settings.charge_patches:
            new_patch = cco_config.copy_patch(patch)
            kbp2_settings.patches.append(patch)

        #BOND PACTHES -> HEME AND COPPER
        kbp2_settings.patches_no_autogen = []
        for patch in cco_settings.bond_patches:
            new_patch = cco_config.copy_patch(patch)
            kbp2_settings.patches_no_autogen.append(patch)

        ### FRAMES ###
        frames = []

        if rename or (special == 'flip_charge'):

            if 'flu' in subfolder:
                renamed_folder = sourcefolder + 'renamed/'

            if special == 'flip_charge':
                renamed_folder = sourcefolder + 'renamed_flip/'

            if not os.path.exists(renamed_folder):
                os.mkdir(renamed_folder)

        #convert frames to be done!
        if len(frames_to_do) != len(frame_range):
            frame_range = list(frames_to_do)

        for i in frame_range:

            if not rename and (special != 'flip_charge'):
                frames.append(sourcefolder + 'frame%i.pdb' % i)
            else:
                ##### PDB renaming ######
                if os.path.exists(renamed_folder+'frame%i.pdb' % i):
                    frames.append(renamed_folder+'frame%i.pdb' % i)
                    print 'Already renamed'
                    continue
                pdb = sourcefolder + 'frame%i.pdb' % i
                pdb_mod = kbp2.file_parser.Simple_struct_parser()
                pdb_mod.read_pdb(pdb)

                if 'flu' in subfolder:
                    if rename_flu_deprot:
                        pdb_mod.rename('FLX', 'FLU')
                    elif rename_flu_prot:
                        pdb_mod.rename('FLH', 'FLU')

                # FLIP charge of BNC calculations!
                if special == 'flip_charge':
                    if 'red' in sourcefolder:
                        for seg in pdb_mod.struct.iter_segments():
                            for res in seg.iter_residues():
                                if res.resname == 'HOH':
                                    for atm in res.iter_atoms():
                                        if atm['segname'] ==  'FEOH':
                                            atm['resname'] = 'OHMI'
                                            if atm['name'] == 'XH2':
                                                atm['name'] = 'OH2'
                                            if atm['name'] == 'H2':
                                                pdb_mod.del_atom(atm, no_restruct=False)

                    if 'red' not in sourcefolder:
                        for seg in pdb_mod.struct.iter_segments():
                            for res in seg.iter_residues():
                                if res.resname == 'OHMI':
                                    for atm in res.iter_atoms():
                                        if atm['segname'] ==  'FEOH':
                                            atm['resname'] = 'HOH'
                                            if atm['name'] == 'OH2':
                                                atm['name'] = 'XH2'

                pdb_mod.create_struct()
                pdb_mod.write_pdb(renamed_folder+'frame%i.pdb' % i)
                frames.append(renamed_folder+'frame%i.pdb' % i)

        kbp2_settings.preopt = { 'carb_oxi_relax' : False, \
                                'init_die_4' :  False \
                                }

        kbp2_settings.set_selected_residues(residues_tt)
        # kbp2_settings.tapbs_bin = 'LD_LIBRARY_PATH=/user/jdragelj/bin /scratch/scratch/tmeyer/kbplus2/tapbs_1.3_cav_enere'
        kbp2_settings.tapbs_bin = '/scratch/scratch/tmeyer/kbplus2/tapbs_1.3_cav_enere'


        # ################################  RUN ####################################
        if special == 'flip_charge':
            kbp2_settings.tmp_settings['bomblev_flip_charge'] = True

        # decision about keeping or deleting folders
        kbp2_settings.remove_folders = 'all'


        # # Test for one strcuture -> comment out manager below
        # kbp2_settings.set_structure(frames[0])
        # kbp2_settings.set_workdir('/scratch/scratch/jdragelj/tests/cco_test/')
        # #kbp2.pka_calculation.calc_pkas(kbp2_settings, titratable_definitions)
        # kbp2.pka_cal_jovan.calc_pkas(kbp2_settings, titratable_definitions)


        manager = pka_manager_jd.Pka_manager(run_folder, done_folder, kbp2_settings, titratable_definitions, cpus=3)
        manager.submit_structures(frames)
        while True:
            (done, crashed, running, queued) = manager.get_status()
            print'done=', done, 'crashed=', crashed, 'running=', running, 'queued=', queued
            if running + queued == 0:
                break
            sleep(sleep_time)
        print '--------------------------------'

if __name__ == '__main__':

    # #### K299C-CcO #####
    # print('K299C-CcO')
    # o_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/his526_mds/'
    # o_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/mutant_no_flu/'
    #
    # subfolders_extended = ['hsp73_hse526_red', 'hsd73_hse526_red']
    #
    # subfolders_new = []
    # subfolders = os.listdir(o_folder)
    #
    # for subfolder in subfolders:
    #     if '526' not in subfolder:
    #         continue
    #     subfolders_new.append(subfolder)
    #     frames_folder_name = 'frames'
    #     if subfolder in subfolders_extended:
    #         continue
    # md_f = 'md'
    #
    # for subfolder in subfolders:
    #     if '526' not in subfolder:
    #         continue
    #     if subfolder not in subfolders_extended:
    #         continue
    #     subfolders_new.append(subfolder)
    #
    # md_f = 'md'
    # frames_folder_name = 'frames'
    # sleep_time = 600
    # frame_range = np.arange(0,103,1)
    # residues_tt = ['HSP-526_ACHA', 'HSP-73_BCHA']
    # titrate_frames_md(o_folder, o_workfolder, subfolders_new, frames_folder_name=frames_folder_name, sleep_time=sleep_time, frame_range=frame_range, residues_tt=residues_tt, md_folder = md_f)


    # # ##### wt-CcO #####
    # print('wt-CcO')
    # o_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/wild_type_mds/'
    # # o_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/'
    # o_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/wild_type/'
    # # subfolders = os.listdir(o_folder)
    # # subfolders_new = []
    # # for subfolder in subfolders:
    # #     if '526' not in subfolder:
    # #         continue
    # #     subfolders_new.append(subfolder)
    #
    # subfolders_new = ['hsp73_hsp526', 'hsp73_hse526', 'hsp73_hsp526_red', 'hsd73_hsp526_red']
    # sleep_time = 600
    # frame_range = np.arange(0,405,2)
    # residues_tt = ['HSP-526_ACHA', 'HSP-73_BCHA']
    # # titrate_frames_md(o_folder, o_workfolder, subfolders_new, frames_folder_name='frames_1', sleep_time=sleep_time, frame_range=frame_range, residues_tt=residues_tt, md_folder = 'md_ex_10')
    # titrate_frames_md(o_folder, o_workfolder, subfolders_new, frames_folder_name='frames_1', sleep_time=sleep_time, frame_range=frame_range,
    #                   residues_tt=residues_tt, md_folder = 'md_ex_10', special='flip_charge')


    ##### Flu-K299C-CcO trajectories#####
    print('Flu-K299C-CcO')
    o_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/'
    o_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/'
    # o_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/'
    # o_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/Cua_hemea/'

    # selected_mds_after_interaction_energy = {'hsp73_hsp526_red_flu-_2': [90], \
    #                                          'hsp73_hsp526_red_fluh': [90], \
    #                                          'hsp73_hsp526_red_flu2': [90], \
    #
    #                                          'hsp73_hsp526_fluh': [90], \
    #                                          'hsp73_hsp526_flu2': [90], \
    #                                          'hsp73_hsp526_flu-_2': [90],\
    #
    #                                          'hsp73_hse526_flu-': [90], \
    #                                          'hsp73_hse526_fluh': [90], \
    #                                          'hsp73_hse526_flu2': [90], \
    #
    #                                          'hsp73_hse526_red_flu-': [90], \
    #                                          'hsp73_hse526_red_flu2': [90], \
    #                                          'hsp73_hse526_red_fluh': [90],\
    #
    #                                          'hsd73_hse526_flu-': [90],\
    #                                          'hsd73_hse526_red_flu-': [90],\
    #
    #                                          }
    #
    # print('Selected_from_Newst_trajcs_July2017_30ns')
    # selected_mds_after_interaction_energy = {'hsp73_hsp526_red_fluh': [90], \
    #                                          'hsp73_hsp526_red_flu2': [90], \
    #                                          'hsp73_hsp526_fluh': [90], \
    #                                          'hsp73_hsp526_flu-_2': [90], \
    #                                          'hsp73_hse526_flu-': [90], \
    #                                          }

    print('Selected_from_Newst_trajcs_July2017_30ns')
    selected_mds_after_interaction_energy = {'hsp73_hsp526_red_flu2': [90], \
                                             'hsp73_hsp526_flu-_2': [90], \
                                             }


    for trajectory, position_list in selected_mds_after_interaction_energy.iteritems():
        for position in position_list:
            folder = o_folder + str(position) + '/'
            workfolder = o_workfolder + str(position) + '/'
            sleep_time = 600
            frame_range = np.arange(344,645,2)
            frames_folder_name = 'frames_1'
            residues_tt = ['HSP-526_ACHA', 'HSP-73_BCHA', 'FLU-1_FLUR']
            subfolders = []
            subfolders.append(trajectory)
            # titrate_frames_md(folder, workfolder, subfolders, frames_folder_name, sleep_time, frame_range, residues_tt, md_folder='md_ex_17_ex_10')
            # titrate_frames_md(folder, workfolder, subfolders, frames_folder_name, sleep_time, frame_range, residues_tt, md_folder='md_ex_10')
            titrate_frames_md(folder, workfolder, subfolders, frames_folder_name, sleep_time, frame_range, residues_tt, md_folder='md_ex_10', special='flip_charge')
            # titrate_frames_md(folder, workfolder, subfolders, frames_folder_name, sleep_time, frame_range, residues_tt, md_folder='md_ex_10', special='flip_charge_Cua_hemea')

