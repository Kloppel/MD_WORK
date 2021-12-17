# coding=utf-8

import kbp2

#from workspace_jd import cco_config


def titrate_structure(pdb_structure, workfolder, state, residues_tt=None):

    #############
    ### CALCS ###
    #############
    #top = topology
    top = []
    top.append("/toppar27/top_alw_clean.inp")#here, but what is the difference between top_alw_clean.inp und top_alw_clean_kb.inp
    top.append("/toppar27/patches.rtf")#here
    top.append("/toppar27/top_all36_lipid.rtf")#here
    top.append("/toppar27/top_all36_cgenff.rtf")#here
    #par = parameters
    par = []
    par.append("/toppar27/par_all22_prot_plus_heme_and_Cu.inp")#here but what is the diference between regular and _kb file?
    par.append("/toppar27/patches.prm")#here
    par.append("/toppar27/par_all36_lipid.prm")#here
    par.append("/toppar27/par_all36_cgenff.prm")#here


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
    #titratable_yaml = 'toppar27/titratable_cco_noglupasp.yaml'
    #kbp2_settings.set_yaml(titratable_yaml)
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
        charge_patches.append({'CBF2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})
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
        charge_patches.append({'CBF2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

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
        charge_patches.append({'CBT4': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

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
        charge_patches.append({'OCBO': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'PER-1_PERI', 'TYR-288_ACHA', 'HSE-284_ACHA']})


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

    kbp2_settings.preopt = { 'carb_oxi_relax' : False,
                            'init_die_4' :  True
                            }#beide false für Titration aus MD, Koordinaten der schweren Atome nicht ändern für Koordinaten aus MD
    kbp2_settings.tapbs_bin = 'LD_LIBRARY_PATH=/home/pbuser/Desktop/PhD_WORK/TARs/kb2plus_package/lib_fort_gnu /home/pbuser/Desktop/PhD_WORK/TARs/kb2plus_package/bin/tapbs_1.3_cav_enere'
    kbp2_settings.remove_folders = 'all'
    kbp2_settings.set_structure(pdb_structure)
    kbp2_settings.set_workdir(workfolder)
    kbp2_settings.apbs_res = 0.3
    kbp2.pka_calculation.calc_pkas(kbp2_settings, titratable_definitions)


if __name__ == '__main__':
    workfolder = '/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/test/titrate_3HB3/'
    pdb_structure = '/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/test/titrate_3HB3/3hb3.pdb'
    state = '21Cua_3hea_3hea3_2Cub'
    titrate_structure(pdb_structure=pdb_structure, workfolder=workfolder, state=state, residues_tt=None)

    #print('for project with "Maria Andrea and Inez Weidinger')

    #workfolder = '/test/titrate_3HB3/'
    #pdb_structure = '/test/titrate_3HB3/3hb3.pdb'
    #state = '21Cua_3hea_3hea3_2Cub'
    #print('state is OO')
    #titrate_structure(pdb_structure, workfolder, state, residues_tt=None)

    #workfolder = '/test/titrate_3HB3/'
    #pdb_structure = '/test/titrate_3HB3/3hb3.pdb'
    #state = '11Cua_2hea_3hea3_2Cub'
    #print('state is RO')
    #titrate_structure(pdb_structure, workfolder, state, residues_tt=None)

    #workfolder = '/test/titrate_3HB3/'
    #pdb_structure = '/test/titrate_3HB3/3hb3.pdb'
    #state = '11Cua_2hea_2hea3_1Cub'
    #print('state is RR')
    #titrate_structure(pdb_structure, workfolder, state, residues_tt=None)
