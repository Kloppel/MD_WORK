# coding=utf-8

import kbp2
import os
import numpy as np
import shutil


def conf_energy_computation(pdb_structure_cco, modelling_folder, jobname, cavity_parameter=None, apbs=False, small=None,
                            glu_deprot=False, asp_deprot=False):
    if not os.path.exists(modelling_folder):
        os.mkdir(modelling_folder)

    top = []
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp")
    top.append("/scratch/scratch/tmeyer/karlsbergplus/patches.rtf")
    par = []
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
    charge_patches = []
    bond_patches = []
    charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
    charmm_struct.add_structure(pdb_structure_cco)
    decisions = []
    decisions.append(('rename__HOH_TIP3', 'keep'))
    decisions.append(('disu_bridges', 'closed'))
    decisions.append(('cap_termini', 'dont_cap'))

    ## Pr state ##
    charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
    charge_patches.append(
        {'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                  'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
    charge_patches.append({'A3H4': ['HEM-2_EHEM', 'HSD-419_ACHA', 'HAO2-1_FEOH']})
    charge_patches.append({'CBPN': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC',
                                    'TYR-288_ACHA', 'HSE-284_ACHA']})

    bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
    bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
    bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
    bond_patches.append({'EISO': ['HAO2-1_FEOH', 'HEM-2_EHEM']})
    bond_patches.append(
        {'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC']})
    bond_patches.append(
        {'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                  'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    charge_patches.append({'LSN': ['LYS-362_ACHA']})


    if decisions is not None:
        for decision_name, decision in decisions:
            charmm_struct.add_decision(decision_name, decision)
    for patch in charge_patches:
        patch_name = patch.keys()[0]
        residues = patch[patch_name]
        charmm_struct.add_patch(patch_name, residues)
    for patch in bond_patches:
        patch_name = patch.keys()[0]
        residues = patch[patch_name]
        charmm_struct.add_patch(patch_name, residues, no_autogen=True)

    if not glu_deprot:
        charmm_struct.set_prot_residue(('GLU', 286, 'ACHA'), patch='GLUP')
    if not asp_deprot:
        charmm_struct.set_prot_residue(('ASP', 407, 'ACHA'), patch='ASPP')

    charmm_struct.workdir = modelling_folder

    selection_string = ''

    if small:
        selection_string = small + """
ENERgy EPS 4.0 E14Fac 0.0
INTEraction select small end select all .and. .not. small end
bomblev -5

delete atom sele all .and. .not. small end
!scalar charge set 0.0 atom sele all .and. .not. small end

ENERgy EPS 4.0 E14Fac 0.0
"""
    else:
        selection_string = small + """
ENERgy EPS 4.0 E14Fac 0.0
"""

    charmm_struct.add_charmm_command(selection_string, adj_task='hbuild')


    charmm_struct.check_structures(quiet=True)
    charmm_struct.charmm_instructions['do_minimize'] = False
    charmm_struct.run_charmm()

    if apbs:
        premodelled_structure = charmm_struct.get_modelled_structure()
        charmm_ssp = premodelled_structure
        if not charmm_ssp.par_read:
            for par in charmm_struct.par:
                charmm_ssp.read_par(par)

        # charmm_ssp.write_pqr(charmm_struct.workdir + 'test.pqr')

        apbs_bin = "/scratch/scratch/tmeyer/CHARMM_NAMD/apbs"
        coulomb_bin = "/scratch/scratch/tmeyer/CHARMM_NAMD/coulomb"
        jobname = jobname
        run_folder = modelling_folder + 'apbs/'
        print run_folder
        target_res = 0.5
        print target_res
        kbp2.apbs_manager.prepare_job(run_folder, jobname, charmm_ssp, apbs_bin, coulomb_bin, target_res=target_res,
                                      verbose=True, pqr_filename=None, \
                                      water_folder=None, conc=0.1, queues=None,
                                      dolly_run_folder='/public/scratch/jdragelj/', cavity_parameter=cavity_parameter)
        kbp2.apbs_manager.submit_job(run_folder)


def get_inter_charmm_energy(workdir, KJ=True, k=1):
    '''
    Gets a list of energies from output of CHARMM run

    Example output:

    # ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
    # ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
    # ENER CROSS:           CMAPs        PMF1D        PMF2D        PRIMO
    # ENER EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
    #  ----------       ---------    ---------    ---------    ---------    ---------
    # ENER>        0   1522.52937      0.00000     18.03137
    # ENER INTERN>      544.02813   1124.93499    169.23668    713.96413     66.85994
    # ENER CROSS>      -176.10529      0.00000      0.00000      0.00000
    # ENER EXTERN>     -468.39449   -451.99472      0.00000      0.00000      0.00000
    #  ----------       ---------    ---------    ---------    ---------    ---------
    '''

    if workdir[-1] != '/':
        workdir += '/'

    charmm_energies = {}

    files = os.listdir(workdir)
    for file in files:
        if '_charmm.out' in file:
            charmm_out_filename = workdir + file

    f = open(charmm_out_filename)
    i = 0
    for line in f:
        ener_keyword = 'INTE ENR:'
        if ener_keyword in line:
            i += 1
            if i == k:
                types_and_names = []
                types_and_values = []
                while not '----------' in line:
                    line = line.strip('\n')
                    types_and_names.append(line.split(':'))
                    line = f.next()
                line = f.next()
                while not '----------' in line:
                    types_and_values.append(line.split('>'))
                    line = f.next()
                break

    for (types_descr, names), (types, values) in zip(types_and_names, types_and_values):
        names_list = names.split()
        values_list = values.split()
        for name, value in zip(names_list, values_list):
            value = float(value)
            # Convert to kJ/mol
            if KJ:
                value *= 4.184
            if not name in charmm_energies:
                charmm_energies[name] = [value]
            else:
                charmm_energies[name].append(value)
        f.close()

    return charmm_energies


def get_charmm_energy_many(workdir, KJ=True, k=1):
    '''
    Gets a list of energies from output of CHARMM run

    Example output:

    # ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
    # ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
    # ENER CROSS:           CMAPs        PMF1D        PMF2D        PRIMO
    # ENER EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
    #  ----------       ---------    ---------    ---------    ---------    ---------
    # ENER>        0   1522.52937      0.00000     18.03137
    # ENER INTERN>      544.02813   1124.93499    169.23668    713.96413     66.85994
    # ENER CROSS>      -176.10529      0.00000      0.00000      0.00000
    # ENER EXTERN>     -468.39449   -451.99472      0.00000      0.00000      0.00000
    #  ----------       ---------    ---------    ---------    ---------    ---------
    '''

    if workdir[-1] != '/':
        workdir += '/'

    charmm_energies = {}

    files = os.listdir(workdir)
    for file in files:
        if '_charmm.out' in file:
            charmm_out_filename = workdir + file

    f = open(charmm_out_filename)
    i = 0
    for line in f:
        ener_keyword = 'ENER ENR:'
        if ener_keyword in line:
            i += 1
            if i == k:
                types_and_names = []
                types_and_values = []
                while not '----------' in line:
                    line = line.strip('\n')
                    types_and_names.append(line.split(':'))
                    line = f.next()
                line = f.next()
                while not '----------' in line:
                    types_and_values.append(line.split('>'))
                    line = f.next()
                break

    for (types_descr, names), (types, values) in zip(types_and_names, types_and_values):
        names_list = names.split()
        values_list = values.split()
        for name, value in zip(names_list, values_list):
            value = float(value)
            # Convert to kJ/mol
            if KJ:
                value *= 4.184
            if not name in charmm_energies:
                charmm_energies[name] = [value]
            else:
                charmm_energies[name].append(value)
        f.close()

    return charmm_energies


if __name__ == '__main__':

    ######################################### SMALL ###############################################

    string1 = """
define small sele ((resid 44 .and. segid ACHA) -
.or.(resid 48 .and. segid ACHA) -
.or.(resid 51 .and. segid ACHA) -
.or.(resid 52 .and. segid ACHA) -
.or.(resid 55 .and. segid ACHA) -
.or.(resid 95 .and. segid ACHA) -
.or.(resid 96 .and. segid ACHA) -
.or.(resid 97 .and. segid ACHA) -
.or.(resid 98 .and. segid ACHA) -
.or.(resid 99 .and. segid ACHA) -
.or.(resid 100 .and. segid ACHA) -
.or.(resid 101 .and. segid ACHA) -
.or.(resid 102 .and. segid ACHA) -
.or.(resid 103 .and. segid ACHA) -
.or.(resid 104 .and. segid ACHA) -
.or.(resid 105 .and. segid ACHA) -
.or.(resid 106 .and. segid ACHA) -
.or.(resid 107 .and. segid ACHA) -
.or.(resid 108 .and. segid ACHA) -
.or.(resid 109 .and. segid ACHA) -
.or.(resid 110 .and. segid ACHA) -
.or.(resid 111 .and. segid ACHA) -
.or.(resid 112 .and. segid ACHA) -
.or.(resid 169 .and. segid ACHA) -
.or.(resid 170 .and. segid ACHA) -
.or.(resid 171 .and. segid ACHA) -
.or.(resid 172 .and. segid ACHA) -
.or.(resid 173 .and. segid ACHA) -
.or.(resid 174 .and. segid ACHA) -
.or.(resid 175 .and. segid ACHA) -
.or.(resid 176 .and. segid ACHA) -
.or.(resid 178 .and. segid ACHA) -
.or.(resid 179 .and. segid ACHA) -
.or.(resid 194 .and. segid ACHA) -
.or.(resid 242 .and. segid ACHA) -
.or.(resid 243 .and. segid ACHA) -
.or.(resid 246 .and. segid ACHA) -
.or.(resid 249 .and. segid ACHA) -
.or.(resid 250 .and. segid ACHA) -
.or.(resid 253 .and. segid ACHA) -
.or.(resid 273 .and. segid ACHA) -
.or.(resid 275 .and. segid ACHA) -
.or.(resid 276 .and. segid ACHA) -
.or.(resid 277 .and. segid ACHA) -
.or.(resid 278 .and. segid ACHA) -
.or.(resid 279 .and. segid ACHA) -
.or.(resid 280 .and. segid ACHA) -
.or.(resid 281 .and. segid ACHA) -
.or.(resid 282 .and. segid ACHA) -
.or.(resid 283 .and. segid ACHA) -
.or.(resid 284 .and. segid ACHA) -
.or.(resid 285 .and. segid ACHA) -
.or.(resid 286 .and. segid ACHA) -
.or.(resid 287 .and. segid ACHA) -
.or.(resid 288 .and. segid ACHA) -
.or.(resid 289 .and. segid ACHA) -
.or.(resid 290 .and. segid ACHA) -
.or.(resid 291 .and. segid ACHA) -
.or.(resid 327 .and. segid ACHA) -
.or.(resid 328 .and. segid ACHA) -
.or.(resid 330 .and. segid ACHA) -
.or.(resid 331 .and. segid ACHA) -
.or.(resid 332 .and. segid ACHA) -
.or.(resid 333 .and. segid ACHA) -
.or.(resid 334 .and. segid ACHA) -
.or.(resid 335 .and. segid ACHA) -
.or.(resid 336 .and. segid ACHA) -
.or.(resid 337 .and. segid ACHA) -
.or.(resid 338 .and. segid ACHA) -
.or.(resid 340 .and. segid ACHA) -
.or.(resid 348 .and. segid ACHA) -
.or.(resid 352 .and. segid ACHA) -
.or.(resid 355 .and. segid ACHA) -
.or.(resid 391 .and. segid ACHA) -
.or.(resid 393 .and. segid ACHA) -
.or.(resid 394 .and. segid ACHA) -
.or.(resid 395 .and. segid ACHA) -
.or.(resid 396 .and. segid ACHA) -
.or.(resid 397 .and. segid ACHA) -
.or.(resid 398 .and. segid ACHA) -
.or.(resid 399 .and. segid ACHA) -
.or.(resid 400 .and. segid ACHA) -
.or.(resid 401 .and. segid ACHA) -
.or.(resid 402 .and. segid ACHA) -
.or.(resid 403 .and. segid ACHA) -
.or.(resid 404 .and. segid ACHA) -
.or.(resid 406 .and. segid ACHA) -
.or.(resid 407 .and. segid ACHA) -
.or.(resid 408 .and. segid ACHA) -
.or.(resid 409 .and. segid ACHA) -
.or.(resid 410 .and. segid ACHA) -
.or.(resid 411 .and. segid ACHA) -
.or.(resid 412 .and. segid ACHA) -
.or.(resid 413 .and. segid ACHA) -
.or.(resid 414 .and. segid ACHA) -
.or.(resid 415 .and. segid ACHA) -
.or.(resid 416 .and. segid ACHA) -
.or.(resid 417 .and. segid ACHA) -
.or.(resid 418 .and. segid ACHA) -
.or.(resid 419 .and. segid ACHA) -
.or.(resid 420 .and. segid ACHA) -
.or.(resid 421 .and. segid ACHA) -
.or.(resid 422 .and. segid ACHA) -
.or.(resid 423 .and. segid ACHA) -
.or.(resid 424 .and. segid ACHA) -
.or.(resid 425 .and. segid ACHA) -
.or.(resid 426 .and. segid ACHA) -
.or.(resid 427 .and. segid ACHA) -
.or.(resid 428 .and. segid ACHA) -
.or.(resid 429 .and. segid ACHA) -
.or.(resid 467 .and. segid ACHA) -
.or.(resid 468 .and. segid ACHA) -
.or.(resid 471 .and. segid ACHA) -
.or.(resid 472 .and. segid ACHA) -
.or.(resid 475 .and. segid ACHA) -
.or.(resid 479 .and. segid ACHA) -
.or.(resid 480 .and. segid ACHA) -
.or.(resid 481 .and. segid ACHA) -
.or.(resid 482 .and. segid ACHA) -
.or.(resid 483 .and. segid ACHA) -
.or.(resid 216 .and. segid BCHA) -
.or.(resid 217 .and. segid BCHA) -
.or.(resid 218 .and. segid BCHA) -
.or.(resid 227 .and. segid BCHA) -
.or.(resid 229 .and. segid BCHA) -
.or.(resid 252 .and. segid BCHA) -
.or.(resid 253 .and. segid BCHA) -
.or.(resid 254 .and. segid BCHA) -
.or.(resid 255 .and. segid BCHA) -
.or.(resid 256 .and. segid BCHA) -
.or.(resid 260 .and. segid BCHA) -
.or.(resid 263 .and. segid BCHA) -
.or.(resid 1 .and. segid EHEM) -
.or.(resid 2 .and. segid EHEM) -
.or.(resid 3 .and. segid EHEM) -
.or.(resid 1 .and. segid FEOH) -
.or.(resid 1 .and. segid GHEM) -
.or.(resid 2 .and. segid GHEM) -
.or.(resid 3 .and. segid GHEM) -
.or.(resid 1 .and. segid HOHC) -
.or.(resid 1 .and. segid META) -
.or.(resid 2 .and. segid META) -
.or.(resid 3 .and. segid META) -
.or.(resid 4 .and. segid META) -
.or. segid VODA) end
"""

    # tests = ['']
    # for test in tests:
    #     ##### MOVED AND MINIMIZED MODEL PRD MIN ######
    #     # pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/cco_prd-_protmin2.pdb'
    #     pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/cco_prd-_nomin2.pdb'
    #     test_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/e80_05_test%s/' % (test)
    #
    #     if not os.path.exists(test_folder):
    #         os.mkdir(test_folder)
    #     pdb_structure_w = test_folder + '/m_min_cco_voda%s.pdb' % test
    #     jobname  = 'model_min_w%s' % (test)
    #     modelling_folder = test_folder + '%s/' % jobname
    #     if not os.path.exists(modelling_folder+'apbs/apbs.out'):
    #         pdb = kbp2.file_parser.Simple_struct_parser()
    #         pdb.read_pdb(pdb_structure_og)
    #         pdb.create_struct()
    #         for seg in pdb.struct.iter_segments():
    #                 for res in seg.iter_residues():
    #                         change = False
    #                         if res.resname == 'TIP3':
    #                                 if res.segname == 'QBH2':
    #                                         if res.resid in [42,48,55]:
    #                                                 change = True
    #                                 if change:
    #                                         for atm in res.iter_atoms():
    #                                                 atm['segname'] = 'VODA'
    #         pdb.create_struct()
    #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB', 'XTC', 'XT2'], exclude=True)
    #         pdb_w.write_pdb(pdb_structure_w)
    #         ##############
    #
    #         if test in ['_nw']:
    #             string = string1
    #         if test in ['']:
    #             string = ''
    #
    #         modelling_folder = test_folder + '%s/' % jobname
    #         conf_energy_computation(pdb_structure_w, modelling_folder, jobname, cavity_parameter=0.9, apbs=True,
    #                                 small=string)
    #
    #     ##### MINIMIZED CRYSTAL STRUCTURE ######
    #     # pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/cco_and_water_protmin2.pdb'
    #     pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/cco_and_water_nomin2.pdb'
    #     test_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/e80_05_test%s/' % test
    #     if not os.path.exists(test_folder):
    #         os.mkdir(test_folder)
    #     jobname = 'crystal_min_w%s' % test
    #     modelling_folder = test_folder + '%s/' % jobname
    #     pdb_structure_w = test_folder + '/c_min_cco_voda%s.pdb' % test
    #     if not os.path.exists(modelling_folder + 'apbs/apbs.out'):
    #         pdb = kbp2.file_parser.Simple_struct_parser()
    #         pdb.read_pdb(pdb_structure_og)
    #         pdb.create_struct()
    #         for seg in pdb.struct.iter_segments():
    #                 for res in seg.iter_residues():
    #                         change = False
    #                         if res.resname == 'TIP3':
    #                                 if res.segname == 'QBH2':
    #                                         if res.resid in [42,48,55]:
    #                                                 change = True
    #                                 if change:
    #                                         for atm in res.iter_atoms():
    #                                                 atm['segname'] = 'VODA'
    #         pdb.create_struct()
    #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB', 'XTC', 'XT2'], exclude=True)
    #         pdb_w.write_pdb(pdb_structure_w)
    #         ##############
    #         if test in ['_nw']:
    #             string = string1
    #         if test in ['']:
    #             string = ''
    #
    #         modelling_folder = test_folder + '%s/' % jobname
    #         conf_energy_computation(pdb_structure_w, modelling_folder, jobname, cavity_parameter=0.9, apbs=True,
    #                                 small=string)
    #
    # ###################################################################################################
    # ######################################   FRAMES   #################################################
    # ###################################################################################################
    #
    # frames = np.arange(0,601,10)
    # tests = ['']
    # for test in tests:
    #     for frame in frames:
    #         pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prdh1_xtra6/frames_voda/frame%i.pdb' % (frame)
    #         run_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/xtra/e80_05_frames_test%s/model/f%i/' % (test, frame)
    #         pdb_structure_w = run_folder + '/mod_cco_voda%s.pdb' % test
    #         jobname  = 'model_f%i_w%s' % (frame, test)
    #         modelling_folder = run_folder + '%s/' % jobname
    #         if not os.path.exists(modelling_folder+'apbs/apbs.out'):
    #             if not os.path.exists(run_folder):
    #                 os.mkdir(run_folder)
    #             else:
    #                 shutil.rmtree(run_folder)
    #                 os.mkdir(run_folder)
    #             pdb = kbp2.file_parser.Simple_struct_parser()
    #             pdb.read_pdb(pdb_structure_og)
    #             pdb.create_struct()
    #             for seg in pdb.struct.iter_segments():
    #                     for res in seg.iter_residues():
    #                             change = False
    #                             if res.resname == 'TIP3':
    #                                     if res.segname == 'QBH2':
    #                                             if res.resid in [42,48,55]:
    #                                                     change = True
    #                                     if change:
    #                                             for atm in res.iter_atoms():
    #                                                     atm['segname'] = 'VODA'
    #             pdb.create_struct()
    #             pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
    #             pdb_w.write_pdb(pdb_structure_w)
    #             ##############
    #
    #             if test in ['_nw']:
    #                 string = string1
    #             if test in ['']:
    #                 string = ''
    #
    #             conf_energy_computation(pdb_structure_w, modelling_folder, jobname, cavity_parameter=0.9, apbs=True,
    #                                     small=string)
    #         else:
    #             continue
    #
    #     for frame in frames:
    #         ##### MINIMIZED CRYSTAL STRUCTURE ######
    #         pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_Pr_crystal_xtra4/frames_voda/frame%i.pdb' % (frame)
    #         run_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/xtra/e80_05_frames_test%s/crystal/f%i/' % (test, frame)
    #         pdb_structure_w = run_folder + '/crystal_cco_voda%s.pdb' % test
    #         jobname  = 'crystal_f%i_w%s' % (frame, test)
    #         modelling_folder = run_folder + '%s/' % jobname
    #         if not os.path.exists(modelling_folder + 'apbs/apbs.out'):
    #             if not os.path.exists(run_folder):
    #                 os.mkdir(run_folder)
    #             else:
    #                 shutil.rmtree(run_folder)
    #                 os.mkdir(run_folder)
    #             pdb = kbp2.file_parser.Simple_struct_parser()
    #             pdb.read_pdb(pdb_structure_og)
    #             pdb.create_struct()
    #             for seg in pdb.struct.iter_segments():
    #                     for res in seg.iter_residues():
    #                             change = False
    #                             if res.resname == 'TIP3':
    #                                     if res.segname == 'QBH2':
    #                                             if res.resid in [42,48,55]:
    #                                                     change = True
    #                                     if change:
    #                                             for atm in res.iter_atoms():
    #                                                     atm['segname'] = 'VODA'
    #             pdb.create_struct()
    #             pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
    #             pdb_w.write_pdb(pdb_structure_w)
    #             ##############
    #             modelling_folder = run_folder + '%s/' % jobname
    #             if test in ['_nw']:
    #                 string = string1
    #             if test in ['']:
    #                 string = ''
    #
    #             modelling_folder = run_folder + '%s/' % jobname
    #             conf_energy_computation(pdb_structure_w, modelling_folder, jobname, cavity_parameter=0.9, apbs=True,
    #                                     small=string)
    #         else:
    #             continue





################################################################################################################################################################################################################################


                    ##################################################################################################################
                    ######################################### RESULTS MODELS #########################################################
                    ##################################################################################################################


    frames_tests = ['']
    for ft in frames_tests:
        print 'crystal', ft
        modelling_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/e80_05_test/crystal_min_w%s/' % (ft)
        run_folder = modelling_folder + 'apbs/'
        apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
        charmm_energy = get_charmm_energy_many(modelling_folder, k=1)
        conf_energy_crystal = apbs_energy + charmm_energy['ELEC'][0]
        print "**********Conf energy: %.2f KJ/mol (apbs %.2f, coulomb %.2f)" % (conf_energy_crystal, apbs_energy, charmm_energy['ELEC'][0])

        print 'model', ft
        modelling_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/e80_05_test/model_min_w%s/' % (ft)
        run_folder = modelling_folder + 'apbs/'
        apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
        charmm_energy = get_charmm_energy_many(modelling_folder, k=1)
        conf_energy_model = apbs_energy + charmm_energy['ELEC'][0]
        print "**********Conf energy: %.2f KJ/mol (apbs %.2f, coulomb %.2f)" % (conf_energy_model, apbs_energy, charmm_energy['ELEC'][0])
        print '***********'
        print 'Econf(crystal)-Econf(model)=%.2f KJ/mol' %(conf_energy_crystal-conf_energy_model)
        print '***********'


    #################################################################################################################
    ######################################## RESULTS FRAMES #########################################################
    #################################################################################################################


    frames_tests = ['']
    for ft in frames_tests:
        print 'tests', ft
        # frames = np.arange(0, 601, 10)
        frames = np.arange(0, 501, 10)
        energies = []
        for frame in frames:
            modelling_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/xtra/e80_05_frames_test%s/crystal/f%i/crystal_f%i_w%s/' % (ft, frame, frame,ft)
            run_folder = modelling_folder + 'apbs/'
            apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
            charmm_energy = get_charmm_energy_many(modelling_folder, k=1)
            conf_energy_closed_w = apbs_energy + charmm_energy['ELEC'][0]
            energies.append(conf_energy_closed_w + 53200)
            if (conf_energy_closed_w + 52000) < -5000:
                print frame
        avgEc = sum(energies) / len(energies)

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax1 = fig.add_subplot(211)

        frame_scale = []
        for time in frames:
            time = time / 20.0
            frame_scale.append(time)

        frames = np.arange(0, 501, 10)
        avgC = []
        for frame in frames:
            avgC.append(avgEc)

        ax1.plot(frame_scale, energies, lw=1.0, color='red')
        # ax1.plot(frame_scale, energies, lw=1.0, color='green')
        ax1.plot(frame_scale, avgC, lw=1.0, linestyle='--', color='red')


        energies = []
        for frame in frames:
            modelling_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/xtra/e80_05_frames_test%s/model/f%i/model_f%i_w%s/' % (ft, frame, frame,ft)
            run_folder = modelling_folder + 'apbs/'
            apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
            charmm_energy = get_charmm_energy_many(modelling_folder, k=1)
            conf_energy_open_w = apbs_energy + charmm_energy['ELEC'][0]
            energies.append(conf_energy_open_w + 53200)
        avgEm = sum(energies) / len(energies)
        print 'diff model', avgEc - avgEm


        frames = np.arange(0, 501, 10)
        avgM = []
        for frame in frames:
            avgM.append(avgEm)


        # ax1.plot(frame_scale, energies, lw=1.0, color='blue')
        ax1.plot(frame_scale, energies, lw=1.0, color='green')
        ax1.plot(frame_scale, avgM, lw=1.0, linestyle='--', color='green')


        ax1.tick_params(direction='out', top='off', right='off')
        ax1.spines['top'].set_linewidth(0.0)
        ax1.spines['left'].set_linewidth(1.0)
        ax1.spines['right'].set_linewidth(0.0)
        ax1.spines['bottom'].set_linewidth(1.0)
        # fig.set_size_inches(3.2, 2.7)
        # plt.plot([2.6, 2.15])
        plt.rcParams.update({'font.size': 16})
        plt.ylim(0,700)
        plt.savefig('/home/users/j/jdragelj/Desktop/PRD_ConfE_25ns.png' )



        print 'diff', avgEc - avgEm
        # plt.show()







