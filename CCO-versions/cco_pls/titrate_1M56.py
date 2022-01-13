# coding=utf-8

import kbp2
from workspace_jd import cco_config


def titrate_structure(pdb_structure, workfolder, state, residues_tt=None):

    #############
    ### CALCS ###
    #############

    top = []
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw_clean_kb.inp")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.rtf")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
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

    protocol = []
    protocol.append((  7, 'h_min'))
    kbp2_settings.protocol = protocol

    # kbp2_settings.md_evaluation_mode = True
    # titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2_np.yaml'

    # titratable_yaml = '/user/jdragelj/python/karlsberg/kbp2/additional_files/titratable.yaml'
    # titratable_yaml = '/scratch/scratch/jdragelj/KB_files/titratable_cco.yaml'
    titratable_yaml = '/scratch/scratch/jdragelj/KB_files/titratable_cco_noglupasp.yaml'

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


    if state == '21Cua_3hea_3hea3_2Cub':
        charge_patches.append({'AHE3': ['HEM-2_F', 'HSD-102_A', 'HSD-421_A']})
        charge_patches.append({'CA21': ['CU1-2_E', 'CU1-3_E', 'HSE-217_B', 'CYS-252_B', 'GLU-254_B', 'CYS-256_B', 'HSE-260_B', 'MET-263_B']})
        charge_patches.append({'A3H3': ['HEM-2_G', 'HSD-419_A', 'OHMI-1_V']})
        charge_patches.append({'CBF2': ['CU1-1_E', 'HSD-333_A', 'HSD-334_A', 'HOH-1_W', 'TYR-288_A', 'HSE-284_A']})
        bond_patches.append({'PHEM': ['HSD-419_A', 'HEM-2_G']})
        bond_patches.append({'PHEM': ['HSD-102_A', 'HEM-2_F']})
        bond_patches.append({'PHE2': ['HSD-421_A', 'HEM-2_F']})
        bond_patches.append({'EISE': ['OHMI-1_V', 'HEM-2_G']})
        bond_patches.append({'CUB2': ['CU1-1_E', 'HSE-284_A', 'HSD-333_A', 'HSD-334_A', 'HOH-1_W']})
        bond_patches.append({'CUAP': ['CU1-2_E', 'CU1-3_E', 'HSE-217_B', 'CYS-252_B', 'GLU-254_B', 'CYS-256_B', 'HSE-260_B', 'MET-263_B']})

    if state == '11Cua_2hea_3hea3_2Cub':
        charge_patches.append({'AHE2': ['HEM-2_F', 'HSD-102_A', 'HSD-421_A']})
        charge_patches.append({'CA11': ['CU1-2_E', 'CU1-3_E', 'HSE-217_B', 'CYS-252_B', 'GLU-254_B', 'CYS-256_B', 'HSE-260_B', 'MET-263_B']})
        charge_patches.append({'A3H3': ['HEM-2_G', 'HSD-419_A', 'OHMI-1_V']})
        charge_patches.append({'CBF2': ['CU1-1_E', 'HSD-333_A', 'HSD-334_A', 'HOH-1_W', 'TYR-288_A', 'HSE-284_A']})
        bond_patches.append({'PHEM': ['HSD-419_A', 'HEM-2_G']})
        bond_patches.append({'PHEM': ['HSD-102_A', 'HEM-2_F']})
        bond_patches.append({'PHE2': ['HSD-421_A', 'HEM-2_F']})
        bond_patches.append({'EISE': ['OHMI-1_V', 'HEM-2_G']})
        bond_patches.append({'CUB2': ['CU1-1_E', 'HSE-284_A', 'HSD-333_A', 'HSD-334_A', 'HOH-1_W']})
        bond_patches.append({'CUAP': ['CU1-2_E', 'CU1-3_E', 'HSE-217_B', 'CYS-252_B', 'GLU-254_B', 'CYS-256_B', 'HSE-260_B', 'MET-263_B']})

    #CHARGE PACTHES
    for patch in charge_patches:
        # new_patch = cco_config.copy_patch(patch)
        kbp2_settings.patches.append(patch)
        # kbp2_settings.patches_no_autogen.append(patch)

    #BOND PACTHES -> HEME AND COPPER
    for patch in bond_patches:
        # new_patch = cco_config.copy_patch(patch)
        kbp2_settings.patches_no_autogen.append(patch)

    kbp2_settings.set_excluded_residues(['HSD-102_A', 'HSD-421_A', 'HSD-419_A', 'HSE-217_B', 'CYS-252_B', 'GLU-254_B', 'CYS-256_B', \
                                         'HSE-260_B', 'HSD-333_A', 'HSD-334_A', 'TYR-288_A', 'HSE-284_A'])

    kbp2_settings.preopt = { 'carb_oxi_relax' : False, \
                            'init_die_4' :  False \
                            }
    kbp2_settings.tapbs_bin = 'LD_LIBRARY_PATH=/user/jdragelj/bin /scratch/scratch/tmeyer/kbplus2/tapbs_1.3_cav_enere'
    kbp2_settings.remove_folders = 'keep'
    kbp2_settings.set_structure(pdb_structure)
    kbp2_settings.set_workdir(workfolder)
    kbp2_settings.apbs_res = 0.3
    kbp2.pka_calculation.calc_pkas(kbp2_settings, titratable_definitions)


if __name__ == '__main__':

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_model_pls_prd/titrate/21Cua_3hea_3hea3_2Cub/'
    # pdb_structure = '/scratch/scratch/jdragelj/projects/cco_model_pls_prd/titrate/cco_3a3_2Cub_HOH.pdb'
    # state = '21Cua_3hea_3hea3_2Cub'
    # print 'state is O O'
    # titrate_structure(pdb_structure, workfolder, state, residues_tt=None)

    workfolder = '/scratch/scratch/jdragelj/projects/cco_model_pls_prd/titrate/11Cua_2hea_3hea3_2Cub/'
    pdb_structure = '/scratch/scratch/jdragelj/projects/cco_model_pls_prd/titrate/cco_3a3_2Cub_HOH.pdb'
    state = '11Cua_2hea_3hea3_2Cub'
    print 'state is R O'
    titrate_structure(pdb_structure, workfolder, state, residues_tt=None)