# coding=utf-8

import os
import kbp2


def build_structure(protein, modelling_folder, state=None, propionic_patch=None):

    top = []
    # top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw_jd_oxygen.inp")
    top.append("/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_PF_crystal_PRDa-/top_alw_clean.inp")
    # top.append("/mnt/fu-scratch/tmeyer/karlsbergplus/patches.rtf")
    top.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
    # top.append("/mnt/fu-scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf")
    par = []
    par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
    # par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
    par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
    # par.append("//scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")

    charge_patches = []
    bond_patches = []

    charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
    charmm_struct.add_structure(protein)

    decisions = []
    decisions.append(('rename__HOH_TIP3', 'keep'))
    decisions.append(('disu_bridges', 'closed'))
    decisions.append(('cap_termini', 'dont_cap'))

    if state == '21Cua_3hea_3hea3_2Cub':
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'A3H3': ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
        charge_patches.append({'CBF2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISO': ['OHMI-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA','CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    if state == '11Cua_2hea_3hea3_2Cub':
        charge_patches.append({'AHE2': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA11': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'A3H3': ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
        charge_patches.append({'CBF2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISO': ['OHMI-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA','CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    if state == 'O2_bound':
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'OA3O': ['HEM-2_EHEM', 'HSD-419_ACHA', 'PER-1_PERI']})
        charge_patches.append({'OCBO': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'PER-1_PERI', 'TYR-288_ACHA', 'HSE-284_ACHA']})
        # charge_patches.append({'OXIO': ['HEM-2_EHEM', 'PER-1_PERI', 'CU1-1_META']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA','CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    if state == 'O2_bound_YH':
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'OA3O': ['HEM-2_EHEM', 'HSD-419_ACHA', 'PER-1_PERI']})
        charge_patches.append({'OCB1': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'PER-1_PERI', 'TYR-288_ACHA', 'HSE-284_ACHA']})
        # charge_patches.append({'OXIO': ['HEM-2_EHEM', 'PER-1_PERI', 'CU1-1_META']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA','CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    if state == 'O2_bound_YH_HOHC':
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'OA3O': ['HEM-2_EHEM', 'HSD-419_ACHA', 'PER-1_PERI']})
        charge_patches.append({'CB1T': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})
        charge_patches.append({'OXI2': ['HEM-2_EHEM', 'PER-1_PERI']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA','CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})


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

    if state == 'Pm':
        print 'hemea reduced'
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


    charge_patches.append({'LSN': ['LYS-362_ACHA']})


    # #propionic
    # if propionic_patch == 'PRAHa':
    #     charge_patches.append({'PRA2': ['PRA-1_GHEM']})
    # if propionic_patch == 'PRDHa':
    #     charge_patches.append({'PRDH': ['PRD-3_GHEM']})
    # if propionic_patch == 'PRAHa3':
    #     charge_patches.append({'PRA2': ['PRA-1_EHEM']})
    # if propionic_patch == 'PRDHa3':
    #     charge_patches.append({'PRDH': ['PRD-3_EHEM']})

    # propionic
    if propionic_patch == 'PRAH':
        charge_patches.append({'PRAH': ['PRA-1_EHEM']})
    if propionic_patch == 'PRA2':
        charge_patches.append({'PRA2': ['PRA-1_EHEM']})
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
    # charmm_struct.set_prot_residue(('ASP', 407, 'ACHA'), patch='ASPP')

#     print 'WARNING! ADDED BLOCK OF TEXT!'
#     block = """
# CONS FIX SELECT .not. (resname HOH .and. type H*) END
# MINI SD NSTEP 500 TOLG 0.1
# CONS FIX SELECT NONE END
#         """
#     charmm_struct.add_charmm_command(block, adj_task='hbuild')

    charmm_struct.workdir = modelling_folder

    charmm_struct.check_structures(quiet=True)

    charmm_struct.charmm_instructions['do_minimize'] = False

    charmm_struct.run_charmm()


if __name__ == '__main__':

    # #pre prep
    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/input_2GSM/cco_and_water.pdb'
    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/read_2GSM/OO_jd/'
    # state = '21Cua_3hea_3hea3_2Cub'
    # build_structure(protein, modelling_folder, state)
    #
    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/input_2GSM/cco_and_water.pdb'
    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/read_2GSM/RO_jd/'
    # state = '11Cua_2hea_3hea3_2Cub'
    # build_structure(protein, modelling_folder, state)
    #
    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/model_O2_2GSM/cco_and_water_and_O2.pdb'
    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/read_2GSM/withO2/'
    # state = 'O2_bound'
    # build_structure(protein, modelling_folder, state)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/raman_models/final_model/boundO2_YH/O2H/CuH2O/cco_O2_YH.pdb'
    # # protein = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/raman_models/final_model/boundO2_YH/O2H/cco_withO2_YH_min_HOHC.pdb'
    # # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/raman_models/final_model/boundO2_YH/O2H/CuH2O/'
    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/raman_models/final_model/boundO2_YH/O2H/CuH2O/O2_YH_model_read/'
    # state = 'O2_bound_YH_HOHC'
    # build_structure(protein, modelling_folder, state)

    ###################################################################################################

    # propionic_patches = [None, 'PRAHa', 'PRDHa', 'PRAHa3', 'PRDHa3']
    # propionic_patches = ['None']
    #
    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/original_models/cco_O_O.pdb'
    # state = '21Cua_3hea_3hea3_2Cub'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     if propionic_patch is not None:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/test/O_O_%s/' % propionic_patch
    #     else:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/test/O_O/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/original_models/cco_O_O.pdb'
    # state = '21Cua_3hea_3hea3_2Cub'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     if propionic_patch is not None:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/O_O_%s/' % propionic_patch
    #     else:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/O_O/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/original_models/cco_R_O.pdb'
    # state = '11Cua_2hea_3hea3_2Cub'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     if propionic_patch is not None:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/R_O_%s/' % propionic_patch
    #     else:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/R_O/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)
    #
    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/original_models/cco_withO2.pdb'
    # state = 'O2_bound'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     if propionic_patch is not None:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/withO2_%s/' % propionic_patch
    #     else:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/final_model/with_O2/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/raman_models/final_model/original_models/cco_withO2_YH.pdb'
    # state = 'O2_bound_YH'
    # for propionic_patch in propionic_patches:
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/raman_models/final_model/withO2YH/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)


    #
    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/model_open_sb_2GSM/quick_with_water/md_cons_not_Arg_wat/cco_Pr_nomin.pdb'
    # state = 'Pr'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     if propionic_patch is not None:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/model_open_sb_2GSM/quick_with_water/md_cons_not_Arg_wat/model_Pr_%s_nomin/' % propionic_patch
    #     else:q
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/model_open_sb_2GSM/quick_with_water/md_cons_not_Arg_wat/model_Pr_nomin/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # propionic_patches = ['PRDHa3']
    #
    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/model_open_sb_2GSM/quick_with_water/md_cons_not_Arg_wat/cco_moved_prdh1.pdb'
    # state = 'Pr'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     if propionic_patch is not None:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/model_open_sb_2GSM/quick_with_water/md_cons_not_Arg_wat/model_Pr_%s_nomin/' % propionic_patch
    #     else:
    #         modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/model_open_sb_2GSM/quick_with_water/md_cons_not_Arg_wat/model_Pr_nomin/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls_prd/model_open_sb_2GSM/quick_with_water/modelF/cco_and_water.pdb'
    # state = 'F'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/model_open_sb_2GSM/quick_with_water/modelF/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_Pr/md_Pr_prdh2/read_cco/cco_prdh2_for_md.pdb'
    # state = 'Pr'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_Pr/md_Pr_prdh2/read_cco/read/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_Pr/md_Pr_prd-/read_cco/cco_prd-_for_md.pdb'
    # state = 'Pr'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_Pr/md_Pr_prd-/read_cco/read/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_Pr/md_Pr_crystal_xtra4/read_xtra4/cco_and_water_xtra.pdb'
    # state = 'Pr'
    # propionic_patches = [None]
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_Pr/md_Pr_crystal_xtra4/read_xtra4/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prd-_xtra8_up/read/cco_prd-_for_md_xwat_all.pdb'
    # state = 'Pr'
    # propionic_patches = [None]
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pr/md_Pr_prd-_xtra8_up/read/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pr_prd-_xtra_up/read/9xtra/cco_prd-_xtra9.pdb'
    # state = 'Pm'
    # propionic_patches = [None]
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_Pm/md_Pr_prd-_xtra_up/read/9xtra/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_F/cco_and_water_F.pdb'
    # state = 'F'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_F/model_F/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_F/prep_open_sb/cco_and_water.pdb'
    # state = 'F'
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_F/prep_open_sb/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)


    # protein = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/move/cco_moved_mini_prah_envir4.pdb'
    # state = 'Pr'
    # propionic_patches = [None]
    # for propionic_patch in propionic_patches:
    #     print state, propionic_patch
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/remodel/'
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #     build_structure(protein, modelling_folder, state, propionic_patch)

    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_PF_crystal_PRDa-/frame326.pdb'
    protein = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_PF_crystal_PRDa-/cco_and_water_prep.pdb'
    state = 'PF'
    propionic_patches = [None]
    for propionic_patch in propionic_patches:
        print state, propionic_patch
        modelling_folder = '/mnt/fu-scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_PF_crystal_PRDa-/change_prot/'
        if not os.path.exists(modelling_folder):
            os.mkdir(modelling_folder)
        build_structure(protein, modelling_folder, state, propionic_patch)

    # import shutil
    #
    # for md in ['prd-', 'prdh1', 'prdh2']:
    #     protein = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_F/md_F_%s_xtra6/read_cco/cco_%s_xtra_min2.pdb' % (md, md)
    #     state = 'F'
    #     if md == 'prdh1':
    #         propionic_patch = 'PRDH'
    #     elif md == 'prdh2':
    #         propionic_patch = 'PRD2'
    #     else:
    #         propionic_patch = None
    #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_F/md_F_%s_xtra6/read_cco/' % md
    #     if not os.path.exists(modelling_folder):
    #         os.mkdir(modelling_folder)
    #
    #     build_structure(protein, modelling_folder, state, propionic_patch)
    #     shutil.copy2(modelling_folder + 'cco_%s_xtra_min2_out.xplor.psf' % md,modelling_folder + 'cco_%s.psf.xplor' % md)
    #     shutil.copy2(modelling_folder + 'cco_%s_xtra_min2_out.pdb' % md, modelling_folder + 'cco_%s.pdb' % md)
    #     os.remove(modelling_folder + 'cco_%s_xtra_min2_out.xplor.psf' % md)
    #     os.remove(modelling_folder + 'cco_%s_xtra_min2_out.pdb' % md)
    #
    # # for md in ['pra-', 'prah1', 'prah2']:
    # #     protein = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/F_state/md_%s/read_cco/cco_moved_mini_%s_envir4.pdb' % (md, md)
    # #     state = 'F'
    # #     if md == 'prah1':
    # #         propionic_patch = 'PRAH'
    # #     elif md == 'prah2':
    # #         propionic_patch = 'PRA2'
    # #     else:
    # #         propionic_patch = None
    # #     modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/F_state/md_%s/read_cco/' % md
    # #     if not os.path.exists(modelling_folder):
    # #         os.mkdir(modelling_folder)
    # #
    # #     build_structure(protein, modelling_folder, state, propionic_patch)
    # #     shutil.copy2(modelling_folder + 'cco_moved_mini_%s_envir4_out.xplor.psf' % md, modelling_folder + 'cco_%s.psf.xplor' % md)
    # #     shutil.copy2(modelling_folder + 'cco_moved_mini_%s_envir4_out.pdb' % md, modelling_folder + 'cco_%s.pdb' % md)
    # #     os.remove(modelling_folder + 'cco_moved_mini_%s_envir4_out.xplor.psf' % md)
    # #     os.remove(modelling_folder + 'cco_moved_mini_%s_envir4_out.pdb' % md)
