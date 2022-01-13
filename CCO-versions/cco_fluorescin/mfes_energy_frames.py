# coding=utf-8

import os
import kbp2
from workspace_jd import cco_config
import re
import shutil


def prepare_calculations_fluorescein(workfolder, folder, frame_selection, subfolder_selection, charmm_EPS, with_membrane, restart, t=None, meshsize=None, del_prev = False):

    subfolders = os.listdir(folder)
    subfolder_new=[]
    for subfolder in subfolders:
        if subfolder[:3] not in ['hsd','hsp']:
            continue
        else:
            subfolder_new.append(subfolder)

    if not os.path.exists(workfolder):
        os.mkdir(workfolder)
                                                            ####################
                                                            ### FRAME CALCS  ###
                                                            ####################

                                                        # Things to be cafeul about:
                                                        # 1. reduces of oxidized enzime
                                                        # 2. histidines should be as they are
                                                        # 3. lysine protonation state 0/+
                                                        # 4. fluorescin state which side is used!
                                                        # 5. if present what is the initial protonaton state of fluorescin

    for subfolder in subfolder_new:

        calculation_info = re.split('[_]', subfolder)

        if subfolder_selection is not None:
            if subfolder not in subfolder_selection:
                continue

        ####################
        ### JOB spliting ###
        ####################

        # [0] = hsd/hsp
        # [1] = lsn/lys
        # [2] = if red -> reduced; if fluh, flu-, flu2 -> oxidised
        # [3] = if exists it is reduced + fluh, flu-, flu2

        #also little check for right naming
        if len(calculation_info) == 4:
            if calculation_info[2] == 'red':
                dolly_number = subfolder
            else:
                continue
                raise AssertionError('Something is wrong with folders')

        elif len(calculation_info) == 3:
            if calculation_info[2] in ['fluh', 'flu-', 'flu2']:
                dolly_number = subfolder

            else:
                continue
        else:
            continue

        print '--------------------------------'
        ###############
        ### FOLDERS ###
        ###############

        work_subfolder = workfolder + subfolder + '/'
        if with_membrane:
            sourcefolder = folder + subfolder + '/' + 'md/titrate73b/frames_test/'
        else:
            sourcefolder = folder + subfolder + '/' + 'md/titrate_his73b/frames/'


        if not os.path.exists(sourcefolder):
            print sourcefolder
            raise AssertionError('folder_missing!')
        else:
            print 'Starting work in %s' % sourcefolder


        # ONLY WITH FLUORESCIN
        nr_of_atr = len(calculation_info)
        if  (nr_of_atr >= 2 and calculation_info[nr_of_atr-1][:3] != 'flu'):
            print "Nothing to titrate for now!", sourcefolder
            continue

        if not os.path.exists(work_subfolder):
            os.mkdir(work_subfolder)


        #############
        ### CALCS ###
        #############


        top = []
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/top_alw.inp")
        # top.append("/scratch/scratch/awoelke/md_cco/toppar/patches.rtf")
        # top.append("/scratch/scratch/awoelke/md_cco/toppar/top_all36_lipid.rtf")
        # top.append("/scratch/scratch/awoelke/md_cco/toppar/top_all36_cgenff.rtf")

        par = []
        # par.append("/scratch/scratch/awoelke/md_cco/toppar/par_all22_prot_plus_heme_and_Cu_kb.inp")
        # par.append("/scratch/scratch/awoelke/md_cco/toppar/patches.prm")
        # par.append("/scratch/scratch/awoelke/md_cco/toppar/par_all36_lipid.prm")
        # par.append("/scratch/scratch/awoelke/md_cco/toppar/par_all36_cgenff.prm")
        #

        # top.append("/user/jdragelj/projects/cco_fluorescin/top_alw.inp")
        top.append("/scratch/scratch/jdragelj/toppar/patches.rtf")
        top.append("/scratch/scratch/jdragelj/toppar/top_all36_lipid.rtf")
        top.append("/scratch/scratch/jdragelj/toppar/top_all36_cgenff.rtf")

        par.append("/scratch/scratch/jdragelj/toppar/par_all22_prot_plus_heme_and_Cu_kb.inp")
        par.append("/scratch/scratch/jdragelj/toppar/patches.prm")
        par.append("/scratch/scratch/jdragelj/toppar/par_all36_lipid.prm")
        par.append("/scratch/scratch/jdragelj/toppar/par_all36_cgenff.prm")



        if calculation_info[nr_of_atr-1][:4] == 'fluh':
            top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_fluh.rtf")
            par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_fluh.prm")

        elif calculation_info[nr_of_atr-1][:4] == 'flu-':
            top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_flux.rtf")
            par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_flux.prm")

        elif calculation_info[nr_of_atr-1][:4] == 'flu2':
            top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_flu2.rtf")
            par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_flu2.prm")

                        #########
                        ### 1 ###        # 1. reduces of oxidized enzime
                        #########

        cco_state = None
        if nr_of_atr > 2:
            if calculation_info[2] == 'red':
                cco_state = 'red'
            else:
                cco_state = 'oxi'
        else:
            cco_state = 'oxi'


        cco_settings = cco_config.COX_settings(state=cco_state, specie='A_paraccocus')

                        #########
                        ### 3 ###       3. lysine protonation state 0/+
                        #########


        charmm_structures = []
        for i in range(0,101):
            if i not in frame_selection:
                continue

            frame = sourcefolder + 'frame%i.pdb' % i
            workdir =  work_subfolder + 'frame%i/' %i
            if workdir[-1] != '/':
                workdir += '/'

            if del_prev:
                shutil.rmtree(workdir)

            if not os.path.exists(workdir):
                os.mkdir(workdir)

            print frame

            ### CHARMM ###
            charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
            charmm_struct.add_structure(frame)
            charmm_struct.add_decision('rename__CA_CAL', 'keep')
            charmm_struct.add_decision('rename__HOH_TIP3', 'keep')
            charmm_struct.add_decision('disu_bridges', 'closed')
            charmm_struct.charmm_bin = "charmm_c39a2q"

            titr_residue_dict_template = charmm_struct.get_titr_residue_dict()
            ######

            lysine_neutral = False
            if calculation_info[1][:3] == 'lys':
                lysine_neutral = False
            elif  calculation_info[1][:3] == 'lsn':
                lysine_neutral = True

            #PROTONATION PATCHES
            if lysine_neutral:
                charmm_struct.set_prot_residue(('LYS', 354, 'ACHA'), charge=0)
                charmm_struct.set_prot_residue(('GLU', 278, 'ACHA'), patch='GLUP')
                charmm_struct.set_prot_residue(('ASP', 399, 'ACHA'), patch='ASPP')
                charmm_struct.set_prot_residue(('GLU', 481, 'ACHA'), patch='GLUP')
            else:
                charmm_struct.set_prot_residue(('GLU', 278, 'ACHA'), patch='GLUP')
                charmm_struct.set_prot_residue(('ASP', 399, 'ACHA'), patch='ASPP')
                charmm_struct.set_prot_residue(('GLU', 481, 'ACHA'), patch='GLUP')

                            #########
                            ### 4 ###   4. fluorescin present or not
                            #########

            if calculation_info[nr_of_atr-1] == 'flu-':
                #GENERAL PATCHES
                patches = [{'FXCY': ['CYS-301_ACHA', 'FLX-1_FLUR']}]

            if calculation_info[nr_of_atr-1] == 'fluh':
                #GENERAL PATCHES
                patches = [{'FCYS': ['CYS-301_ACHA', 'FLU-1_FLUR']}]

            if calculation_info[nr_of_atr-1] == 'flu2':
                #GENERAL PATCHES
                patches = [{'FCYS': ['CYS-301_ACHA', 'FLH-1_FLUR']}]


            #CHARGE PACTHES
            for patch in cco_settings.charge_patches:
                new_patch = cco_config.copy_patch(patch)
                patches.append(patch)

            for patch in patches:
                patch_name = patch.keys()[0]
                residues = patch[patch_name]
                charmm_struct.add_patch(patch_name, residues)

            #BOND PACTHES -> HEME AND COPPER
            patches_no_autogen = []
            for patch in cco_settings.bond_patches:
                new_patch = cco_config.copy_patch(patch)
                patches_no_autogen.append(patch)

            for patch in patches_no_autogen:
                patch_name = patch.keys()[0]
                residues = patch[patch_name]
                charmm_struct.add_patch(patch_name, residues, no_autogen=True)

            charmm_struct.workdir = workdir

            # copy files needed for mfes
            # mfes_input_files = "/scratch/scratch/jdragelj/mfes_files/"
            # needed_folders = os.listdir(mfes_input_files)
            # for filename in needed_folders:
            #     shutil.copy2(mfes_input_files+'/'+filename, workdir+filename)
            # charmm_struct.add_charmm_command('mfes sele all .and. .not. resname tip3 end', adj_task='hbuild')

            mfes_input_files = "/scratch/scratch/jdragelj/mfes_files/"
            charmm_struct = kbp2.mfes.start_mfes_with_charmm(charmm_struct, mfes_input_files, t, meshsize)

            if charmm_EPS:
                charmm_struct.add_charmm_command('ENERgy EPS 4.0 E14Fac 0.0', adj_task='hbuild')

            charmm_struct.check_structures(quiet=True)
            charmm_struct.charmm_instructions['do_minimize'] = False
            charmm_struct.set_run_charmm()
            charmm_structures.append(charmm_struct)

        if restart:
            dolly_number = dolly_number + '_rst'
        kbp2.charmm.set_dolly_run_charmm(workfolder, dolly_number, charmm_structures)



if __name__ == "__main__":

    print "this is mFES related"

    #calculations with membrane
    workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/mFES_test/memb_fluh/'
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/mFES_test/memb/'
    folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_data/his73b_lys354_flu/'
    frame_selection = [20]
    subfolder_selection = ['hsd73b_lsn354_fluh']
    charmm_EPS = True
    t = 30
    meshsize = 0.5
    restart = False
    with_membrane = True
    # prepare_calculations_fluorescein(workfolder, frame_selection, subfolder_selection, charmm_EPS, t, meshsize, with_membrane, restart)
    prepare_calculations_fluorescein(workfolder, folder, frame_selection, subfolder_selection, charmm_EPS, with_membrane = with_membrane, restart = restart, t=t, meshsize=meshsize)

    # # restart_calculstions
    # workfolder = '/scratch/scratch/jdragelj/cco_alexiev/FLU_energies_memb/membrane_calcs/'
    # folder = '/scratch/scratch/jdragelj/cco_alexiev/fluorescin_data_md/'
    # charmm_EPS = False
    # with_membrane = True
    # restart = True
    # t = 28
    # meshsize = 0.45
    # pickle_calculation = pkl.load(open(workfolder+'restart.pkl', 'rb'))
    # for subfolder, frames in pickle_calculation.iteritems():
    #     subfolder_selection = []
    #     subfolder_selection.append(subfolder)
    #     prepare_calculations_fluorescein(workfolder, frames, subfolder_selection, charmm_EPS, with_membrane, restart, t, meshsize, del_prev = True)












