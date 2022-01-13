# coding=utf-8

import kbp2
import os
from time import sleep
import numpy as np


def conf_energy_computation(pdb_structure_cco, modelling_folder, jobname,  cavity_parameter=None, apbs=False, small=None, propionic_patch=None):

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

        # propionic
        if propionic_patch == 'PRAH':
                charge_patches.append({'PRAH': ['PRA-1_GHEM']})
        if propionic_patch == 'PRA2':
                charge_patches.append({'PRA2': ['PRA-1_GHEM']})
        if propionic_patch == 'PRDH':
                charge_patches.append({'PRDH': ['PRD-3_EHEM']})
        if propionic_patch == 'PRD2':
                charge_patches.append({'PRD2': ['PRD-3_EHEM']})

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
        charmm_struct.set_prot_residue(('GLU', 286, 'ACHA'), patch='GLUP')
        charmm_struct.set_prot_residue(('ASP', 407, 'ACHA'), patch='ASPP')
        charmm_struct.workdir = modelling_folder

        if small:

                selection_string = small + """
ENERgy EPS 4.0 E14Fac 0.0
INTEraction select small end select all .and. .not. small end
bomblev -5
delete atom sele all .and. .not. small end
ENERgy EPS 4.0 E14Fac 0.0
"""
                charmm_struct.add_charmm_command(selection_string, adj_task='hbuild')

        else:
                command_block = """
ENERgy EPS 4.0 E14Fac 0.0
INTEraction select ((resname ARG .and. resid 481) .or. (resname ARG .and. resid 482)) end select all .and. -
                    .not. ((resname ARG .and. resid 481) .or. (resname ARG .and. resid 482)) end

INTEraction select (resname ARG .and. resid 481) end select all .and. .not. (resname ARG .and. resid 481) end

INTEraction select (resname ARG .and. resid 482) end select all .and. .not. (resname ARG .and. resid 482) end
        """
                charmm_struct.add_charmm_command(command_block, adj_task='hbuild')



        charmm_struct.check_structures(quiet=True)
        charmm_struct.charmm_instructions['do_minimize'] = False
        charmm_struct.run_charmm()

        if apbs:
                charmm_ssp = charmm_struct.structure
                premodelled_structure = charmm_struct.get_modelled_structure()
                charmm_ssp = premodelled_structure
                if not charmm_ssp.par_read:
                        for par in charmm_struct.par:
                            charmm_ssp.read_par(par)

                apbs_bin    = "/scratch/scratch/tmeyer/CHARMM_NAMD/apbs"
                coulomb_bin = "/scratch/scratch/tmeyer/CHARMM_NAMD/coulomb"
                jobname = jobname
                run_folder = modelling_folder + 'apbs/'
                print run_folder
                kbp2.apbs_manager.prepare_job(run_folder, jobname, charmm_ssp, apbs_bin, coulomb_bin, target_res=0.3, verbose=True, pqr_filename=None,\
                        water_folder=None, conc=0.1, queues=None, dolly_run_folder='/public/scratch/jdragelj/', cavity_parameter=cavity_parameter)
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
            if i==k:
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
            if i==k:
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

        #### CRYSTAL PRD-####
        # pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/cco_and_water.pdb'
        # #######
        # pdb_structure_no_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/cl_cco.pdb'
        # pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/cl_cco_mgcw.pdb'
        # pdb = kbp2.file_parser.Simple_struct_parser()
        # pdb.read_pdb(pdb_structure_og)
        # pdb.create_struct()
        # pdb_no_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_no_w.write_pdb(pdb_structure_no_w)
        # for seg in pdb.struct.iter_segments():
        #         for res in seg.iter_residues():
        #                 if res.resname == 'TIP3':
        #                         if res.segname == 'QBH2':
        #                                 if res.resid in [42,48,55]:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'MGCW'
        # pdb.create_struct()
        # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_w.write_pdb(pdb_structure_w)
        ################
        # jobname  = 'closed_no_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/' % jobname
        # conf_energy_computation(pdb_structure_no_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)
        # jobname  = 'closed_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)
        # jobname  = 'closed_no_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/inter/' % jobname
        # conf_energy_computation(pdb_structure_no_w, modelling_folder, jobname,  inter=True)
        # jobname  = 'closed_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/inter/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname)

        # #### CRYSTAL PRD- MIN ####
        # pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/cco_and_water_min.pdb'
        # #######
        # pdb_structure_no_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/min_ene/cl_cco.pdb'
        # pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/min_ene/cl_cco_mgcw.pdb'
        # pdb = kbp2.file_parser.Simple_struct_parser()
        # pdb.read_pdb(pdb_structure_og)
        # pdb.create_struct()
        # pdb_no_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_no_w.write_pdb(pdb_structure_no_w)
        # for seg in pdb.struct.iter_segments():
        #         for res in seg.iter_residues():
        #                 if res.resname == 'TIP3':
        #                         if res.segname == 'QBH2':
        #                                 if res.resid in [42,48,55]:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'MGCW'
        # pdb.create_struct()
        # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_w.write_pdb(pdb_structure_w)
        # ################
        # jobname  = 'closed_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/min_ene/%s/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)


        # #### CRYSTAL PRD- PROT MIN ####
        # pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/cco_and_water_protmin.pdb'
        # #######
        # pdb_structure_no_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/min_prot/cl_cco.pdb'
        # pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/min_prot/cl_cco_mgcw.pdb'
        # pdb = kbp2.file_parser.Simple_struct_parser()
        # pdb.read_pdb(pdb_structure_og)
        # pdb.create_struct()
        # pdb_no_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_no_w.write_pdb(pdb_structure_no_w)
        # for seg in pdb.struct.iter_segments():
        #         for res in seg.iter_residues():
        #                 if res.resname == 'TIP3':
        #                         if res.segname == 'QBH2':
        #                                 if res.resid in [42,48,55]:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'MGCW'
        # pdb.create_struct()
        # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_w.write_pdb(pdb_structure_w)
        # ################
        # jobname  = 'closed_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/min_prot/%s/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)

        # #### CRYSTAL PRD- PROT MIN ####
        # pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/cco_and_water_protmin2.pdb'
        # #######
        # pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/min_prot2/cl_cco_mgcw.pdb'
        # pdb = kbp2.file_parser.Simple_struct_parser()
        # pdb.read_pdb(pdb_structure_og)
        # pdb.create_struct()
        # for seg in pdb.struct.iter_segments():
        #         for res in seg.iter_residues():
        #             change = False
        #             if res.resname == 'TIP3':
        #                     if res.segname == 'QBH2':
        #                             if res.resid in [42,48,55]:
        #                                     change = True
        #                     elif res.segname == 'PAH2':
        #                             if res.resid in [3, 15]:
        #                                     change = True
        #                     if change:
        #                             for atm in res.iter_atoms():
        #                                     atm['segname'] = 'VODA'
        # pdb.create_struct()
        # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_w.write_pdb(pdb_structure_w)
        # ################
        # jobname  = 'closed_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/min_prot2/%s/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)


        ###########################################


        # ##### MOVED AND MINIMIZED MODEL PRD-######
        # pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/move_22_25_min_500_prdh2/cco_moved_mini.pdb'
        # pdb_structure_no_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/o_min_cco.pdb'
        # pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/o_min_cco_mgcw.pdb'
        # pdb = kbp2.file_parser.Simple_struct_parser()
        # pdb.read_pdb(pdb_structure_og)
        # pdb.create_struct()
        # pdb_no_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_no_w.write_pdb(pdb_structure_no_w)
        # for seg in pdb.struct.iter_segments():
        #         for res in seg.iter_residues():
        #                 if res.resname == 'TIP3':
        #                         if res.segname == 'QBH2':
        #                                 if res.resid in [42,48,55]:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'MGCW'
        # pdb.create_struct()
        # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_w.write_pdb(pdb_structure_w)
        ###############
        # jobname  = 'open_min_no_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/' % jobname
        # conf_energy_computation(pdb_structure_no_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)
        # jobname  = 'open_min_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)
        # jobname  = 'open_min_no_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/inter/' % jobname
        # conf_energy_computation(pdb_structure_no_w, modelling_folder, jobname)
        # jobname  = 'open_min_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/inter/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname)


        # ##### MOVED AND MINIMIZED MODEL PRD- more MIN ######
        # pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/cco_prd-.pdb'
        # pdb_structure_no_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/mini/o_min_cco.pdb'
        # pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/mini/o_min_cco_mgcw.pdb'
        # pdb = kbp2.file_parser.Simple_struct_parser()
        # pdb.read_pdb(pdb_structure_og)
        # pdb.create_struct()
        # pdb_no_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_no_w.write_pdb(pdb_structure_no_w)
        # for seg in pdb.struct.iter_segments():
        #         for res in seg.iter_residues():
        #                 if res.resname == 'TIP3':
        #                         if res.segname == 'QBH2':
        #                                 if res.resid in [42,48,55]:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'MGCW'
        # pdb.create_struct()
        # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_w.write_pdb(pdb_structure_w)
        # ##############
        # # jobname  = 'open_min_no_w'
        # # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/' % jobname
        # # conf_energy_computation(pdb_structure_no_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)
        # jobname  = 'open_min_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/mini/%s/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)
        # jobname  = 'open_min_no_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/inter/' % jobname
        # conf_energy_computation(pdb_structure_no_w, modelling_folder, jobname)
        # jobname  = 'open_min_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/%s/inter/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname)
        #
        #
        # ##### MOVED AND MINIMIZED MODEL PRD- PROT MIN ######
        # pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/cco_prd-_protmin1.pdb'
        # pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/min_prot/o_min_cco_mgcw.pdb'
        # pdb = kbp2.file_parser.Simple_struct_parser()
        # pdb.read_pdb(pdb_structure_og)
        # pdb.create_struct()
        # pdb_no_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_no_w.write_pdb(pdb_structure_no_w)
        # for seg in pdb.struct.iter_segments():
        #         for res in seg.iter_residues():
        #                 if res.resname == 'TIP3':
        #                         if res.segname == 'QBH2':
        #                                 if res.resid in [42,48,55]:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'MGCW'
        # pdb.create_struct()
        # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_w.write_pdb(pdb_structure_w)
        # ##############
        # jobname  = 'open_min_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/min_prot/%s/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)

        # ##### MOVED AND MINIMIZED MODEL PRD- PROT MIN ######
        # pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/cco_prd-_protmin2.pdb'
        # pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/min_prot2/o_min_cco_mgcw.pdb'
        # pdb = kbp2.file_parser.Simple_struct_parser()
        # pdb.read_pdb(pdb_structure_og)
        # pdb.create_struct()
        # for seg in pdb.struct.iter_segments():
        #         for res in seg.iter_residues():
        #             change = False
        #             if res.resname == 'TIP3':
        #                     if res.segname == 'QBH2':
        #                             if res.resid in [42,48,55]:
        #                                     change = True
        #                     elif res.segname == 'PAH2':
        #                             if res.resid in [3, 15]:
        #                                     change = True
        #                     if change:
        #                             for atm in res.iter_atoms():
        #                                     atm['segname'] = 'VODA'
        # pdb.create_struct()
        # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        # pdb_w.write_pdb(pdb_structure_w)
        # ##############
        # jobname  = 'open_min_w'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/e20_min_prot2/%s/' % jobname
        # conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True)


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

        string2 = """
define small sele (resname CU1 .or. resname MG .or. segid EHEM .or. segid GHEM .or. segid FEOH -
.or. segid VODA .or. (resname CYS .and. (resid 252 .or. resid 256)) .or. (resname ARG .and. (resid 481 .or. resid 482)) -
.or. (resname MET .and. resid 263) .or. (resname TYR .and. resid 288) .or. (resname GLU .and. resid 254) -
.or. (resname ASP .and. resid 412) .or. (segid BCHA .and. (resid 217 .or. resid 260)) .or. segid HOHC -
.or. (segid ACHA .and. (resid 411 .or. resid 102 .or. resid 419 .or. resid 334 .or. resid 333 .or. resid 284 .or. resid 421))) -
end
"""

        # # tests = ['','2']
        # tests = ['_nw','2_nw']
        # for test in tests:
        #     ##### MOVED AND MINIMIZED MODEL PRD- PROT MIN ######
        #     pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/cco_prd-_protmin2.pdb'
        #     test_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/e20_small_test%s/' % test
        #     if not os.path.exists(test_folder):
        #         os.mkdir(test_folder)
        #     pdb_structure_w = test_folder + '/o_min_cco_voda%s.pdb' % test
        #     jobname  = 'open_min_w%s' % test
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
        #                                 # elif res.segname == 'PAH2':
        #                                 #         if res.resid in [3, 15]:
        #                                 #                 change = True
        #                                 if change:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'VODA'
        #         pdb.create_struct()
        #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        #         pdb_w.write_pdb(pdb_structure_w)
        #         ##############
        #         if test in ['','_nw']:
        #             string = string1
        #         if test in ['2', '2_nw']:
        #             string = string2
        #         conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True, small=string)
        #
        #     pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/cco_and_water_protmin2.pdb'
        #     test_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/e20_small_test%s/' % test
        #     if not os.path.exists(test_folder):
        #         os.mkdir(test_folder)
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
        #                                 # elif res.segname == 'PAH2':
        #                                 #         if res.resid in [3, 15]:
        #                                 #                 change = True
        #                                 if change:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'VODA'
        #         pdb.create_struct()
        #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        #         pdb_w.write_pdb(pdb_structure_w)
        #         ##############
        #         jobname  = 'closed_w%s' % test
        #         modelling_folder = test_folder + '%s/' % jobname
        #         if test in ['','_nw']:
        #             string = string1
        #         if test in ['2', '2_nw']:
        #             string = string2
        #         conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True, small=string)

        # ### FRAMES CRYSTAL - NO WAT ####
        # # frames = np.arange(0,801,2)
        # frames=[0]
        # tests = ['_nw']
        # for test in tests:
        #     print 'Conf energy computation: Crystal MD frames with xtra waters PRD-_a3'
        #     for frame in frames:
        #             frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_Pr_crystal_xtra4/frames_voda/'
        #             pdb_structure_og = frames_folder + 'frame%i.pdb' % frame
        #             #######
        #             run_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/xtra/e20_frames_test%s/crystal/f%i/' % (test, frame)
        #             if not os.path.exists(run_folder):
        #                     os.mkdir(run_folder)
        #             else:
        #                 continue
        #             pdb_structure_w = run_folder + '/cl_cco_voda%s.pdb' % test
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
        #                                             # elif res.segname == 'PAH2':
        #                                             #     if 'nw' not in test:
        #                                             #         if res.resid in [3, 15]:
        #                                             #             change = True
        #                                     if change:
        #                                             for atm in res.iter_atoms():
        #                                                     atm['segname'] = 'VODA'
        #             pdb.create_struct()
        #             pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB', 'XTC', 'XT2'], exclude=True)
        #             pdb_w.write_pdb(pdb_structure_w)
        #             ################
        #             jobname = 'closedf%i' % frame
        #             modelling_folder = run_folder + '%s/' % jobname
        #             if test in ['', '_nw']:
        #                 string = string1
        #             if test in ['2', '2_nw']:
        #                 string = string2
        #             conf_energy_computation(pdb_structure_w, modelling_folder, jobname, cavity_parameter=0.9, apbs=True, small=string)
        #
        #     frames = np.arange(0,801,2)
        #     print 'Conf energy computation: Model open sb MD frames with xtra waters PRD-_a3'
        #     for frame in frames:
        #             ##### FRAMES PRD- NO WAT ######
        #             frames_folder =  '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prd-_xtra6/frames_voda/'
        #             pdb_structure_og = frames_folder + 'frame%i.pdb' % frame
        #             #######
        #             run_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/xtra/e20_frames_test%s/model/f%i/' % (test,frame)
        #             if not os.path.exists(run_folder):
        #                     os.mkdir(run_folder)
        #             else:
        #                 continue
        #             pdb_structure_w = run_folder + '/op_cco_voda%s.pdb' % test
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
        #                                     # elif res.segname == 'PAH2':
        #                                     #     if 'nw' not in test:
        #                                     #         if res.resid in [3, 15]:
        #                                     #             change = True
        #                                     if change:
        #                                             for atm in res.iter_atoms():
        #                                                     atm['segname'] = 'VODA'
        #             pdb.create_struct()
        #             pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB', 'XTC', 'XT2'], exclude=True)
        #             pdb_w.write_pdb(pdb_structure_w)
        #             ################
        #             jobname = 'openf%i' % frame
        #             modelling_folder = run_folder + '%s/' % jobname
        #             if test in ['', '_nw']:
        #                 string = string1
        #             if test in ['2', '2_nw']:
        #                 string = string2
        #             conf_energy_computation(pdb_structure_w, modelling_folder, jobname, cavity_parameter=0.9, apbs=True, small=string)

                                    ################################################################################################
                                    ##################################### protonated structures ####################################
                                    ################################################################################################

        ##### MOVED AND MINIMIZED MODEL PRD- PROT MIN ######

        # tests = ['_nw','2_nw']
        # protons = ['h1', 'h2']
        #
        # for proton in protons:
        #     for test in tests:
        #         pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model%s/cco_and_water%s_protmin2.pdb' % (proton, proton)
        #         workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model%s/e20_small_test%s/'  % (proton, test)
        #         if not os.path.exists(workfolder):
        #             os.mkdir(workfolder)
        #         pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model%s/e20_small_test%s/o_min_cco_voda.pdb'  % (proton, test)
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
        #                                 # elif res.segname == 'PAH2':
        #                                 #         if res.resid in [3, 15]:
        #                                 #                 change = True
        #                                 if change:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'VODA'
        #         pdb.create_struct()
        #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        #         pdb_w.write_pdb(pdb_structure_w)
        #
        #         if proton == 'h1':
        #             propionic_patch = 'PRDH'
        #         if proton == 'h2':
        #             propionic_patch = 'PRD2'
        #
        #         if test in ['_nw']:
        #             string = string1
        #         if test in ['2_nw']:
        #             string = string2
        #         jobname  = 'open_min_%s%s' % (proton, test)
        #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model%s/e20_small_test%s/%s/' % (proton, test, jobname)
        #         conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True, small=string, propionic_patch=propionic_patch)
        #
        #         ######################
        #
        #         pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal%s/cco_and_water%s_protmin2.pdb' % (proton, proton)
        #         workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal%s/e20_small_test%s/'  % (proton, test)
        #         if not os.path.exists(workfolder):
        #             os.mkdir(workfolder)
        #         pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal%s/e20_small_test%s/c_min_cco_voda.pdb' % (proton, test)
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
        #                                 # elif res.segname == 'PAH2':
        #                                 #         if res.resid in [3, 15]:
        #                                 #                 change = True
        #                                 if change:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'VODA'
        #         pdb.create_struct()
        #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        #         pdb_w.write_pdb(pdb_structure_w)
        #
        #         if proton == 'h1':
        #             propionic_patch = 'PRAH'
        #         if proton == 'h2':
        #             propionic_patch = 'PRA2'
        #
        #         if test in ['_nw']:
        #             string = string1
        #         if test in ['2_nw']:
        #             string = string2
        #
        #         jobname  = 'closed_%s%s' % (proton, test)
        #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal%s/e20_small_test%s/%s/' % (proton, test, jobname)
        #         conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True, small=string, propionic_patch=propionic_patch)# tests = ['_nw','2_nw']

        # protons = ['h1', 'h2']
        # tests = ['_nw','2_nw']
        # for proton in protons:
        #     for test in tests:
        #         pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model%s/cco_and_water%s_protmin2.pdb' % (proton, proton)
        #         workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model%s/e20_small_test%s/'  % (proton, test)
        #         if not os.path.exists(workfolder):
        #             os.mkdir(workfolder)
        #         pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model%s/e20_small_test%s/o_min_cco_voda.pdb'  % (proton, test)
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
        #                                 # elif res.segname == 'PAH2':
        #                                 #         if res.resid in [3, 15]:
        #                                 #                 change = True
        #                                 if change:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'VODA'
        #         pdb.create_struct()
        #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        #         pdb_w.write_pdb(pdb_structure_w)
        #
        #         if proton == 'h1':
        #             propionic_patch = 'PRDH'
        #         if proton == 'h2':
        #             propionic_patch = 'PRD2'
        #
        #         if test in ['_nw']:
        #             string = string1
        #         if test in ['2_nw']:
        #             string = string2
        #         jobname  = 'open_min_%s%s' % (proton, test)
        #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model%s/e20_small_test%s/%s/' % (proton, test, jobname)
        #         conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True, small=string, propionic_patch=propionic_patch)
        #
        #         ######################
        #
        #         pdb_structure_og = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal%s/cco_and_water%s_protmin2.pdb' % (proton, proton)
        #         workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal%s/e20_small_test%s/'  % (proton, test)
        #         if not os.path.exists(workfolder):
        #             os.mkdir(workfolder)
        #         pdb_structure_w = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal%s/e20_small_test%s/c_min_cco_voda.pdb' % (proton, test)
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
        #                                 # elif res.segname == 'PAH2':
        #                                 #         if res.resid in [3, 15]:
        #                                 #                 change = True
        #                                 if change:
        #                                         for atm in res.iter_atoms():
        #                                                 atm['segname'] = 'VODA'
        #         pdb.create_struct()
        #         pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
        #         pdb_w.write_pdb(pdb_structure_w)
        #
        #         if proton == 'h1':
        #             propionic_patch = 'PRAH'
        #         if proton == 'h2':
        #             propionic_patch = 'PRA2'
        #
        #         if test in ['_nw']:
        #             string = string1
        #         if test in ['2_nw']:
        #             string = string2
        #
        #         jobname  = 'closed_%s%s' % (proton, test)
        #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal%s/e20_small_test%s/%s/' % (proton, test, jobname)
        #         conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True, small=string, propionic_patch=propionic_patch)

        # ### FRAMES PRA_a models protonated NO WAT ####
        # ### REFERENCE STRUCTURES ######
        # print 'Conf energy computation: reference model MD frames PRA_a protonated'
        # frames = np.arange(0, 601, 2)
        # tests = ['_nw']
        # protons = ['h1', 'h2']
        # for proton in protons:
        #     for test in tests:
        #         for frame in frames:
        #             frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a/md_pra%s/frames_voda/' % (proton)
        #             pdb_structure_og = frames_folder + 'frame%i.pdb' % frame
        #             #######
        #             run_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/e20_frames_test%s/ref_model_%s/f%i/' % (test, proton, frame)
        #             if not os.path.exists(run_folder):
        #                 os.mkdir(run_folder)
        #             else:
        #                 continue
        #             pdb_structure_w = run_folder + '/cl_cco_voda%s.pdb' % test
        #             pdb = kbp2.file_parser.Simple_struct_parser()
        #             pdb.read_pdb(pdb_structure_og)
        #             pdb.create_struct()
        #             for seg in pdb.struct.iter_segments():
        #                 for res in seg.iter_residues():
        #                     change = False
        #                     if res.resname == 'TIP3':
        #                         if res.segname == 'QBH2':
        #                             if res.resid in [42, 48, 55]:
        #                                 change = True
        #                             # elif res.segname == 'PAH2':
        #                             #     if 'nw' not in test:
        #                             #         if res.resid in [3, 15]:
        #                             #             change = True
        #                         if change:
        #                             for atm in res.iter_atoms():
        #                                 atm['segname'] = 'VODA'
        #             pdb.create_struct()
        #             pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB', 'XTC', 'XT2'], exclude=True)
        #             pdb_w.write_pdb(pdb_structure_w)
        #             ################
        #             if proton == 'h1':
        #                 propionic_patch = 'PRAH'
        #             if proton == 'h2':
        #                 propionic_patch = 'PRA2'
        #
        #             if test in ['_nw']:
        #                 string = string1
        #             if test in ['2_nw']:
        #                 string = string2
        #
        #             jobname  = 'closed_%s%s' % (proton, test)
        #             modelling_folder = run_folder + '%s/' % jobname
        #             conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True, small=string, propionic_patch=propionic_patch)
        #
        # ### FRAMES PRD_a3 models protonated NO WAT ####
        # print 'Conf energy computation: open sb model MD frames PRD_a3 protonated'
        # frames = np.arange(0, 801, 2)
        # tests = ['_nw']
        # protons = ['h1', 'h2']
        # for proton in protons:
        #     for test in tests:
        #         for frame in frames:
        #             ##### FRAMES PRD- NO WAT ######
        #             frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prd%s_xtra6/frames_voda/' % (proton)
        #             pdb_structure_og = frames_folder + 'frame%i.pdb' % frame
        #             #######
        #             run_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/xtra/e20_frames_test%s/model%s/f%i/' % (test, proton, frame)
        #             if not os.path.exists(run_folder):
        #                 os.mkdir(run_folder)
        #             else:
        #                 continue
        #             pdb_structure_w = run_folder + '/cl_cco_voda%s.pdb' % test
        #             pdb = kbp2.file_parser.Simple_struct_parser()
        #             pdb.read_pdb(pdb_structure_og)
        #             pdb.create_struct()
        #             for seg in pdb.struct.iter_segments():
        #                 for res in seg.iter_residues():
        #                     change = False
        #                     if res.resname == 'TIP3':
        #                         if res.segname == 'QBH2':
        #                             if res.resid in [42, 48, 55]:
        #                                 change = True
        #                             # elif res.segname == 'PAH2':
        #                             #     if 'nw' not in test:
        #                             #         if res.resid in [3, 15]:
        #                             #             change = True
        #                         if change:
        #                             for atm in res.iter_atoms():
        #                                 atm['segname'] = 'VODA'
        #             pdb.create_struct()
        #             pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB', 'XTC', 'XT2'], exclude=True)
        #             pdb_w.write_pdb(pdb_structure_w)
        #             ################
        #             if proton == 'h1':
        #                 propionic_patch = 'PRDH'
        #             if proton == 'h2':
        #                 propionic_patch = 'PRD2'
        #
        #             if test in ['_nw']:
        #                 string = string1
        #             if test in ['2_nw']:
        #                 string = string2
        #
        #             jobname  = 'open_%s%s' % (proton, test)
        #             modelling_folder = run_folder + '%s/' % jobname
        #             conf_energy_computation(pdb_structure_w, modelling_folder, jobname,  cavity_parameter=0.9, apbs=True, small=string, propionic_patch=propionic_patch)

        ######################################################################################################################################################
        ######################################################################################################################################################
        ######################################################################################################################################################
        ######################################################################################################################################################


        # print 'closed sb with water'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/min_prot2/closed_w/'
        # run_folder = modelling_folder + 'apbs/'
        # apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
        # charmm_energy = kbp2.charmm.get_charmm_energy(modelling_folder)
        # conf_energy_closed_w = apbs_energy + charmm_energy['ELEC'][0]
        # print "**********Conf energy: %.2f KJ/mol (apbs %.2f, elec %.2f)" % (conf_energy_closed_w, apbs_energy, charmm_energy['ELEC'][0])
        # # charmm_energy_inter = get_inter_charmm_energy(modelling_folder, k=1)
        # # print '**********Interaction energy 481/482-all: %.2f KJ/mol' % charmm_energy_inter['ELEC'][0]
        # # charmm_energy_inter = get_inter_charmm_energy(modelling_folder, k=2)
        # # print '**********Interaction energy 481-all: %.2f KJ/mol' % charmm_energy_inter['ELEC'][0]
        # # charmm_energy_inter = get_inter_charmm_energy(modelling_folder, k=3)
        # # print '**********Interaction energy 482-all: %.2f KJ/mol' % charmm_energy_inter['ELEC'][0]
        #
        # print 'open min sb with water'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/min_prot2/open_min_w/'
        # run_folder = modelling_folder + 'apbs/'
        # apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
        # charmm_energy = kbp2.charmm.get_charmm_energy(modelling_folder)
        # conf_energy_open_w = apbs_energy + charmm_energy['ELEC'][0]
        # print "**********Conf energy: %.2f KJ/mol (apbs %.2f, elec %.2f)" % (conf_energy_open_w, apbs_energy, charmm_energy['ELEC'][0])
        # print '***********'
        # print 'Econf(closed)-Econf(open)=%.2f KJ/mol' %(conf_energy_closed_w-conf_energy_open_w)
        # print '***********'
        # # charmm_energy_inter = get_inter_charmm_energy(modelling_folder, k=1)
        # # print '**********Interaction energy 481/482-all: %.2f KJ/mol' % charmm_energy_inter['ELEC'][0]
        # # charmm_energy_inter = get_inter_charmm_energy(modelling_folder, k=2)
        # # print '**********Interaction energy 481-all: %.2f KJ/mol' % charmm_energy_inter['ELEC'][0]
        # # charmm_energy_inter = get_inter_charmm_energy(modelling_folder, k=3)
        # # print '**********Interaction energy 482-all: %.2f KJ/mol' % charmm_energy_inter['ELEC'][0]
        #
        # frames_tests = ['','2']
        # for ft in frames_tests:
        #     ################# small ##############################
        #     print 'closed sb with water no minimization', ft
        #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/no_min_crystal/small_test%s/closed_w%s/' % (ft,ft)
        #     run_folder = modelling_folder + 'apbs/'
        #     apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
        #     inter_energy = kbp2.charmm.get_charmm_energy(modelling_folder, inter=True)
        #     charmm_energy = get_charmm_energy_many(modelling_folder, k=2)
        #     conf_energy_closed_w = apbs_energy + inter_energy['ELEC'][0] + charmm_energy['ELEC'][0]
        #     print "**********Conf energy: %.2f KJ/mol (apbs %.2f, coulomb %.2f, inter %.2f)" % (conf_energy_closed_w, apbs_energy, charmm_energy['ELEC'][0], inter_energy['ELEC'][0])
        #
        # print ''
        # print ''

        frames_tests = ['_nw','2_nw']
        for ft in frames_tests:
            ################# small ##############################
            print 'closed sb with water', ft
            modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystal/e20_small_test%s/closed_w%s/' % (ft,ft)
            run_folder = modelling_folder + 'apbs/'
            apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
            inter_energy = kbp2.charmm.get_charmm_energy(modelling_folder, inter=True)
            charmm_energy = get_charmm_energy_many(modelling_folder, k=2)
            conf_energy_closed_w = apbs_energy + inter_energy['ELEC'][0] + charmm_energy['ELEC'][0]
            print "**********Conf energy: %.2f KJ/mol (apbs %.2f, coulomb %.2f, inter %.2f)" % (conf_energy_closed_w, apbs_energy, charmm_energy['ELEC'][0], inter_energy['ELEC'][0])


            print 'open min sb with water', ft
            modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/model/e20_small_test%s/open_min_w%s/' % (ft,ft)
            run_folder = modelling_folder + 'apbs/'
            apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
            inter_energy = kbp2.charmm.get_charmm_energy(modelling_folder, inter=True)
            charmm_energy = get_charmm_energy_many(modelling_folder, k=2)
            conf_energy_open_w = apbs_energy + inter_energy['ELEC'][0] + charmm_energy['ELEC'][0]
            print "**********Conf energy: %.2f KJ/mol (apbs %.2f, coulomb %.2f, inter %.2f)" % (conf_energy_open_w, apbs_energy, charmm_energy['ELEC'][0], inter_energy['ELEC'][0])
            print '***********'
            print 'Econf(closed)-Econf(open)=%.2f KJ/mol' %(conf_energy_closed_w-conf_energy_open_w)
            print '***********'

        frames_tests = ['_nw']
        for ft in frames_tests:
            print 'tests', ft
            # print 'closed sb MD crystal'
            # frames = np.arange(0,401,2)
            frames = np.arange(0,801,2)
            energies = []
            for frame in frames:
                    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/e20_frames_test%s/crystal/f%i/closedf%i/' % (ft,frame, frame)
                    modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/xtra/e20_frames_test%s/crystal/f%i/closedf%i/' % (ft,frame, frame)
                    run_folder = modelling_folder + 'apbs/'
                    apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
                    inter_energy = kbp2.charmm.get_charmm_energy(modelling_folder, inter=True)
                    charmm_energy = get_charmm_energy_many(modelling_folder, k=2)
                    conf_energy_closed_w = apbs_energy + inter_energy['ELEC'][0] + charmm_energy['ELEC'][0]
                    energies.append(conf_energy_closed_w+10900)
            avgEc = sum(energies)/len(energies)

            import matplotlib.pyplot as plt
            plt.plot(frames, energies, color='red')

            fig = plt.figure()
            ax1 = fig.add_subplot(211)

            frame_scale = []
            for time in frames:
                time = time / 20.0
                frame_scale.append(time)

            ax1.plot(frame_scale, energies, lw=0.75, color='red')

            # print 'open sb MD prd-'
            energies = []
            for frame in frames:
                    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/e20_frames_test%s/model/f%i/openf%i/' % (ft,frame, frame)
                    modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/frames_small/xtra/e20_frames_test%s/model/f%i/openf%i/' % (ft,frame, frame)
                    run_folder = modelling_folder + 'apbs/'
                    apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
                    inter_energy = kbp2.charmm.get_charmm_energy(modelling_folder, inter=True)
                    charmm_energy = get_charmm_energy_many(modelling_folder, k=2)
                    conf_energy_open_w = apbs_energy + inter_energy['ELEC'][0] + charmm_energy['ELEC'][0]
                    energies.append(conf_energy_open_w+10900)
            avgEm = sum(energies)/len(energies)

            ax1.plot(frame_scale, energies, lw=0.75, color='blue')
            ax1.tick_params(direction='out', top='off', right='off')

            ax1.spines['top'].set_linewidth(0.0)
            ax1.spines['left'].set_linewidth(0.5)
            ax1.spines['right'].set_linewidth(0.0)
            ax1.spines['bottom'].set_linewidth(0.5)
            # fig.set_size_inches(3.2, 2.7)
            # plt.plot([2.6, 2.15])
            plt.rcParams.update({'font.size': 12})
            # plt.savefig(out_folder + '%s_pka.png' % (out_filename + '_final'))


            print 'diff',avgEc-avgEm
            plt.show()

        # ################# small ##############################
        # print 'closed sb with water'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystalh1/e20_small_test_nw/closed_h1_nw/'
        # run_folder = modelling_folder + 'apbs/'
        # apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
        # inter_energy = kbp2.charmm.get_charmm_energy(modelling_folder, inter=True)
        # charmm_energy = get_charmm_energy_many(modelling_folder, k=2)
        # conf_energy_closed_w1 = apbs_energy + inter_energy['ELEC'][0] + charmm_energy['ELEC'][0]
        # # print "**********Conf energy: %.2f KJ/mol (apbs %.2f, coulomb %.2f, inter %.2f)" % (conf_energy_closed_w1, apbs_energy, charmm_energy['ELEC'][0], inter_energy['ELEC'][0])
        #
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/crystalh2/e20_small_test_nw/closed_h2_nw/'
        # run_folder = modelling_folder + 'apbs/'
        # apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
        # inter_energy = kbp2.charmm.get_charmm_energy(modelling_folder, inter=True)
        # charmm_energy = get_charmm_energy_many(modelling_folder, k=2)
        # conf_energy_closed_w2 = apbs_energy + inter_energy['ELEC'][0] + charmm_energy['ELEC'][0]
        # # print "**********Conf energy: %.2f KJ/mol (apbs %.2f, coulomb %.2f, inter %.2f)" % (conf_energy_closed_w2, apbs_energy, charmm_energy['ELEC'][0], inter_energy['ELEC'][0])
        #
        # conf_energy_closed_w = (conf_energy_closed_w1 + conf_energy_closed_w2) / 2
        # print "**********Conf energy: %.2f KJ/mol " % conf_energy_closed_w
        #
        #
        # print 'open min sb with water'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/modelh1/e20_small_test_nw/open_min_h1_nw/'
        # run_folder = modelling_folder + 'apbs/'
        # apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
        # inter_energy = kbp2.charmm.get_charmm_energy(modelling_folder, inter=True)
        # charmm_energy = get_charmm_energy_many(modelling_folder, k=2)
        # conf_energy_open_w = apbs_energy + inter_energy['ELEC'][0] + charmm_energy['ELEC'][0]
        # print "**********Conf energy: %.2f KJ/mol (apbs %.2f, coulomb %.2f, inter %.2f)" % (
        # conf_energy_open_w, apbs_energy, charmm_energy['ELEC'][0], inter_energy['ELEC'][0])
        # print '***********'
        # print 'Econf(closed)-Econf(open)=%.2f KJ/mol' % (conf_energy_closed_w - conf_energy_open_w)
        # print '***********'
        #
        # print 'open min sb with water'
        # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/conf_ene_Pr_prd_a3/minimized/modelh2/e20_small_test_nw/open_min_h2_nw/'
        # apbs_energy, coulomb_energy = kbp2.apbs_manager.read_result(run_folder, epsilon=4.0)
        # inter_energy = kbp2.charmm.get_charmm_energy(modelling_folder, inter=True)
        # charmm_energy = get_charmm_energy_many(modelling_folder, k=2)
        # conf_energy_open_w = apbs_energy + inter_energy['ELEC'][0] + charmm_energy['ELEC'][0]
        # print "**********Conf energy: %.2f KJ/mol (apbs %.2f, coulomb %.2f, inter %.2f)" % (
        # conf_energy_open_w, apbs_energy, charmm_energy['ELEC'][0], inter_energy['ELEC'][0])
        # print '***********'
        # print 'Econf(closed)-Econf(open)=%.2f KJ/mol' % (conf_energy_closed_w - conf_energy_open_w)
        print '***********'



