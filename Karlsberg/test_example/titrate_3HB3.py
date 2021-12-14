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
    # decisions.append(('rename__HOH_TIP3', 'keep'))
    decisions.append(('rename__HOH_TIP3', 'remove'))
    decisions.append(('disu_bridges', 'closed'))
    decisions.append(('cap_termini', 'dont_cap'))

    kbp2_settings.modelling_decisions = decisions

    protocol = []
    protocol.append((  7, 'h_min'))
    kbp2_settings.protocol = protocol
    # kbp2_settings.cavity_par = [0.8, 0.2, 0.0]

    # kbp2_settings.titr_residue_charges_0 = True
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
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                        'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'A3H3': ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
        charge_patches.append(
            {'CBF2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})
        # charge_patches.append({'CBF2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'TIP3-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISO': ['OHMI-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
        # bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'TIP3-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                      'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    if state == '11Cua_2hea_3hea3_2Cub':
        charge_patches.append({'AHE2': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA11': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                        'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'A3H3': ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
        charge_patches.append(
            {'CBF2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISO': ['OHMI-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                      'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    if state == '11Cua_2hea_2hea3_1Cub':
        charge_patches.append({'AHE2': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA11': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                        'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
        charge_patches.append({'A3H2': ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
        charge_patches.append(
            {'CBT4': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISO': ['OHMI-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                      'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    if state == 'O2_bound':
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                        'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

        charge_patches.append({'OA3O': ['HEM-2_EHEM', 'HSD-419_ACHA', 'PER-1_PERI']})
        charge_patches.append(
            {'OCBO': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'PER-1_PERI', 'TYR-288_ACHA', 'HSE-284_ACHA']})


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



    #CHARGE PACTHES
    for patch in charge_patches:
        # new_patch = cco_config.copy_patch(patch)
        kbp2_settings.patches.append(patch)
        # kbp2_settings.patches_no_autogen.append(patch)

    #BOND PACTHES -> HEME AND COPPER
    for patch in bond_patches:
        # new_patch = cco_config.copy_patch(patch)
        kbp2_settings.patches_no_autogen.append(patch)

    kbp2_settings.set_excluded_residues(['HSD-102_ACHA', 'HSD-421_ACHA', 'HSD-419_ACHA', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', \
                                         'HSE-260_BCHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'TYR-288_ACHA', 'HSE-284_ACHA'])

    kbp2_settings.preopt = { 'carb_oxi_relax' : False, \
                            'init_die_4' :  False \
                            }
    kbp2_settings.tapbs_bin = 'LD_LIBRARY_PATH=/user/jdragelj/bin /scratch/scratch/tmeyer/kbplus2/tapbs_1.3_cav_enere'
    kbp2_settings.remove_folders = 'all'
    kbp2_settings.set_structure(pdb_structure)
    kbp2_settings.set_workdir(workfolder)
    kbp2_settings.apbs_res = 0.3
    kbp2.pka_calculation.calc_pkas(kbp2_settings, titratable_definitions)


if __name__ == '__main__':

    print('for project with "Maria Andrea and Inez Weidinger')

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_model_pls_prd/titrate_2GSM/OO_cavities/'
    # pdb_structure = '/scratch/scratch/jdragelj/projects/cco_model_pls_prd/titrate_2GSM/cco_2gsm.pdb'
    # state = '21Cua_3hea_3hea3_2Cub'
    # print 'state is OO'
    # titrate_structure(pdb_structure, workfolder, state, residues_tt=None)

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_model_pls_prd/titrate_2GSM/RO/'
    # pdb_structure = '/scratch/scratch/jdragelj/projects/cco_model_pls_prd/titrate_2GSM/cco_2gsm.pdb'
    # state = '11Cua_2hea_3hea3_2Cub'
    # print 'state is RO'
    # titrate_structure(pdb_structure, workfolder, state, residues_tt=None)

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_pls_prd/titrate_2GSM/RR/'
    # pdb_structure = '/scratch/scratch/jdragelj/projects/cco_pls_prd/titrate_2GSM/cco_2gsm.pdb'
    # state = '11Cua_2hea_2hea3_1Cub'
    # print 'state is RR'
    # titrate_structure(pdb_structure, workfolder, state, residues_tt=None)


