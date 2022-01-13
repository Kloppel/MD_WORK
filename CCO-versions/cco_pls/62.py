# coding=utf-8

import kbp2
import os
import shutil
import cPickle as pickle
import numpy as np
import time


def titrate_structures_dolly(workfolder, frames_folder, frame_range, md_name, state, residues_tt=None, cavities=[],
                             membrane_charge=True, kbp2_eval_titration=False, glu_asp_lys_patch=[True, True, True], dolly=None):

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
                    if not os.path.exists(done_folder + 'frame%s/results.pkl' % frame):
                        shutil.rmtree(done_folder + 'frame%s/' % frame)
                        frames_to_do.append(frame)
                    else:
                        continue
        else:
            frames_to_do = frame_range

    print md_name, len(frames_to_do), 'left from', len(frame_range)

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
        kbp2_settings.cavity_par = cavities

    if kbp2_eval_titration:
        kbp2_settings.md_evaluation_mode = True
        titratable_yaml = '/scratch/scratch/jdragelj/KB_files/titratable_cco_extended.yaml'
    else:
        titratable_yaml = '/scratch/scratch/jdragelj/KB_files/titratable_cco_noglupasp.yaml'

    kbp2_settings.set_yaml(titratable_yaml)
    ##############################################################################

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

    # #protonation patches
    # print 'LYS362 patched!'
    # charge_patches.append({'LSN': ['LYS-362_ACHA']})
    # print 'GLU286 patched!'
    # charge_patches.append({'GLUP': ['GLU-286_ACHA']})
    # # print 'ASP407 patched!'
    # # charge_patches.append({'ASPP':['ASP-407_ACHA']})

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

    # kbp2_settings.remove_folders = 'keep'
    # kbp2_settings.force_conf_ene_calc = True

    t = 10
    print 'Sleep activated for %i seconds' %t

    frames = []
    for frame in frames_to_do:
        if 'voda' not in frames_folder:
            frames.append(frames_folder+'frame%i.pdb'%frame)
            frame_pdb_final = frames_folder+'frame%i.pdb'%frame
            frames.append(frame_pdb_final)
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
            pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XTC', 'XT2'], exclude=True)
            pdb_w.write_pdb(frames_folder + 'mgcw/frame%i.pdb'%frame)
            frame_pdb_final = frames_folder + 'mgcw/frame%i.pdb'%frame
            frames.append(frame_pdb_final)

        run_folder_frame = run_folder + 'frame%i/' % frame
        if not os.path.exists(run_folder_frame):
            os.mkdir(run_folder_frame)
        kbp2_settings.set_structure(frame_pdb_final)
        jobname = 'frame%i' % frame

        kbp2_set_pkl_file = run_folder_frame + 'kbp_set_frame%i.pkl' % frame
        python_script_file = run_folder_frame + 'titrate_%s.py' % jobname
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
        if 'tapbs01' in md_name:
            python_script = """
# coding=utf-8

import kbp2
import cPickle as pickle

kbp2_settings_pickle = '%s' """ % kbp2_set_pkl_file + """
kbp2_settings =  pickle.load( open( kbp2_settings_pickle, "rb" ) )
titratable_definitions = kbp2.kbp_tools.parse_titratable_yaml(kbp2_settings.titratable_yaml)
kbp2.pka_calculation_tapbs_hr.calc_pkas(kbp2_settings, titratable_definitions)

                    """



        pfile.write(python_script)
        pfile.close()

        local_dolly_folder = '/public/scratch/jdragelj/kbp2/%s_%s/' % (md_name, jobname)
        kbp2_settings.set_workdir(local_dolly_folder)
        pklfile = open(kbp2_set_pkl_file, 'wb')
        pickle.dump(kbp2_settings, pklfile)
        pklfile.close()
        submitt_kbp2_job(done_folder=done_folder, run_folder=run_folder, local_dolly_folder=local_dolly_folder,
                            jobname=jobname, md_name=md_name, dolly=dolly)
        time.sleep(10)


def submitt_kbp2_job(done_folder, run_folder, local_dolly_folder, jobname, md_name, dolly=None):

    submitt_script_file = run_folder + jobname + '/titr_%s_%s.sh' % (md_name,jobname)
    sfile = open(submitt_script_file, 'w')

    submitt_script = """
#!/bin/csh

echo "Running on `hostname`"

setenv PYTHONPATH /itch/itch/local_python/lib/python2.7/site-packages:/itch/itch/local_python/src

set SDIR=%s%s/"""%(run_folder,jobname) + """

set FDIR=%s%s/"""%(done_folder,jobname) + """

set JDIR=%s"""%local_dolly_folder + """

mkdir -p $JDIR
cp $SDIR/* $JDIR

cd $JDIR

python titrate_%s.py """% jobname + """

cd ..

mkdir $FDIR
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

    if dolly:
        shell.stdin.write('qsub -q %s.q -o %s -e %s %s\n' % (dolly, run_folder, run_folder, submitt_script_file) )
    else:
        shell.stdin.write('qsub -q D64.q -o %s -e %s %s\n' % (run_folder, run_folder, submitt_script_file) )
        # shell.stdin.write('qsub -q D62.q -o %s -e %s %s\n' % (run_folder, run_folder, submitt_script_file) )
        # shell.stdin.write('qsub -q D61.q,D62.q,D64.q -o %s -e %s %s\n' % (run_folder, run_folder, submitt_script_file) )
        # shell.stdin.write('qsub -q D62.q,D64.q -o %s -e %s %s\n' % (run_folder, run_folder, submitt_script_file) )
        shell.stdin.write('exit\n')



if __name__ == '__main__':

    titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_pra_a3/10ang_no_cavities_voda_tapbs01/'
    residues = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM', 'ARG-481_ACHA', 'ARG-482_ACHA',
                'GLU-286_ACHA', \
                'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA',
                'LYS-227_BCHA', 'TYR-336_ACHA', \
                'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    frames_suffix = 'frames_voda/'

    frame_range = [20, 100, 180]

    titration_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/F_pra_a3/10ang_no_cavities_voda_tapbs01/'
    workfolder = titration_folder + 'md_pra-/'
    if not os.path.exists(workfolder):
        os.mkdir(workfolder)
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/F_state/md_pra-/%s/' % frames_suffix
    state = 'F'
    md_name = 'tapbs01_10ang_F_pra-_a3'
    titrate_structures_dolly(workfolder, frames_folder, frame_range, md_name, state, residues_tt=residues,
                             cavities=[],
                             membrane_charge=False, kbp2_eval_titration=True, glu_asp_lys_patch=[False, False, True], dolly='D62')
    while True:
        for frame in frame_range:
            comp_done = os.path.exists(workfolder + '/done/frame%s/results.pkl' % frame)
            if comp_done == False:
                break
        if comp_done:
            break
        time.sleep(60)


    frame_range = [20, 100, 180]

    workfolder = titration_folder + 'md_prah1/'
    if not os.path.exists(workfolder):
        os.mkdir(workfolder)
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/F_state/md_prah1/%s/' % frames_suffix
    state = 'F'
    md_name = 'tapbs01_10ang_F_prah1_a3'
    titrate_structures_dolly(workfolder, frames_folder, frame_range, md_name, state, residues_tt=residues,
                             cavities=[],
                             membrane_charge=False, kbp2_eval_titration=True, glu_asp_lys_patch=[False, False, True], dolly='D62')
    while True:
        for frame in frame_range:
            comp_done = os.path.exists(workfolder + '/done/frame%s/results.pkl' % frame)
            if comp_done == False:
                break
        if comp_done:
            break
        time.sleep(60)


    frame_range = [20, 100, 180]

    workfolder = titration_folder + 'md_prah2/'
    if not os.path.exists(workfolder):
        os.mkdir(workfolder)
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/F_state/md_prah2/%s/' % frames_suffix
    state = 'F'
    md_name = 'tapbs01_10ang_F_prah2_a3'
    titrate_structures_dolly(workfolder, frames_folder, frame_range, md_name, state, residues_tt=residues,
                             cavities=[],
                             membrane_charge=False, kbp2_eval_titration=True, glu_asp_lys_patch=[False, False, True], dolly='D62')


