# coding=utf-8

import kbp2
from workspace_jd import cco_config


def titrate_structure(pdb_structure, workfolder, state, memb_charge_zero=False):

    top = []
    # top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw_clean_kb.inp")
    if not memb_charge_zero:
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
    else:
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid_no_charge_POPC.rtf")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.rtf")
    # top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_kbp/top_kbp_fluorescin_scp.rtf")

    par = []
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_fluorescin.str")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_kbp/par_kbp_fluorescin.prm")

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

    kbp2_settings.md_evaluation_mode = True
    titratable_yaml = '/scratch/scratch/jdragelj/KB_files/titratable_cco_extended.yaml'

    # titratable_yaml = '/scratch/scratch/jdragelj/KB_files/titratable_cco_noglupasp.yaml'

    kbp2_settings.set_yaml(titratable_yaml)
    titratable_definitions = kbp2.kbp_tools.parse_titratable_yaml(titratable_yaml)


    kbp2_settings.patches = []
    kbp2_settings.patches_no_autogen = []

    if state == 'red':
        cco_settings = cco_config.COX_settings(specie='A_paraccocus', cycle_state='R')
    else:
        cco_settings = cco_config.COX_settings(specie='A_paraccocus', cycle_state='O')

    kbp2_settings.set_excluded_residues(['CYS-66_ACHA', 'CYS-80_ACHA', 'HSD-411_ACHA', 'HSD-94_ACHA', 'HSD-413_ACHA', \
                                        'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HSE-181_BCHA', 'CYS-216_BCHA', \
                                        'GLU-278_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'TYR-280_ACHA', 'CYS-299_ACHA'])

    resname = "FLU"
    states = []
    state = {'pka': 0.0,
             'name': 'R',
             'external_patches': [('FLUN', ['CYS-299_ACHA', 'FLU-1_FLUR'])]}
    states.append(state)
    state = {'pka': -100.0,
             'name': '0',
             'external_patches': [('FLUH', ['CYS-299_ACHA', 'FLU-1_FLUR'])]}
    states.append(state)
    state = {'pka': -93.8,
             'name': 'D',
             'external_patches': [('FLUX', ['CYS-299_ACHA', 'FLU-1_FLUR'])]}
    states.append(state)
    state = {'pka': -100.0,
             'name': '0',
             'external_patches': [('FLU2', ['CYS-299_ACHA', 'FLU-1_FLUR'])]}
    states.append(state)

    titratable_definition = kbp2.pka_calculation.get_atom_charges(resname, states, top)
    titratable_definition['FLU'][0]['external_atoms'][0][1]['H6'] = 0.001
    titratable_definitions.update(titratable_definition)

    # GENERAL PATCHES


    #CHARGE PACTHES
    for patch in cco_settings.charge_patches:
        kbp2_settings.patches.append(patch)
    #BOND PACTHES -> HEME AND COPPER
    for patch in cco_settings.bond_patches:
        kbp2_settings.patches_no_autogen.append(patch)

    kbp2_settings.patches.append({'FCYS': ['CYS-299_ACHA', 'FLU-1_FLUR']})
    # kbp2_settings.init_modelling_min = False


    kbp2_settings.preopt = { 'carb_oxi_relax' : False, \
                            'init_die_4' :  False \
                            }
    kbp2_settings.tapbs_bin = '/scratch/scratch/tmeyer/kbplus2/tapbs_1.3_cav_enere'
    # kbp2_settings.tapbs_bin = '/scratch/scratch/tmeyer/kbplus2/tapbs_cavities'
    kbp2_settings.remove_folders = 'all'
    kbp2_settings.set_structure(pdb_structure)
    kbp2_settings.set_workdir(workfolder)
    kbp2.pka_calculation.calc_pkas(kbp2_settings, titratable_definitions)


if __name__ == '__main__':

    # workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/crystal/oxi_memb_flu/'
    # pdb_structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/model/conformation_tests/mini_confs/cco_flu_90_min.pdb'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/cco_flu_memb.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure)
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'oxi'
    # titrate_structure(pdb_structure_cco2, workfolder, state)
    #
    # workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/crystal/oxi_flu/'
    # pdb_structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/model/conformation_tests/mini_confs/cco_flu_90_min.pdb'
    # pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/cco_flu.pdb'
    # pdb = kbp2.file_parser.Simple_struct_parser()
    # pdb.read_pdb(pdb_structure)
    # pdb.create_struct()
    # pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'MEMB'], exclude=True)
    # pdb_w.write_pdb(pdb_structure_cco2)
    # state = 'oxi'
    # titrate_structure(pdb_structure_cco2, workfolder, state)

    workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/crystal/oxi_memb_noc_flu/'
    pdb_structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/model/conformation_tests/mini_confs/cco_flu_90_min.pdb'
    pdb_structure_cco2 = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/cco_flu_memb_noc.pdb'
    pdb = kbp2.file_parser.Simple_struct_parser()
    pdb.read_pdb(pdb_structure)
    pdb.create_struct()
    pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2'], exclude=True)
    pdb_w.write_pdb(pdb_structure_cco2)
    state = 'oxi'
    titrate_structure(pdb_structure_cco2, workfolder, state, memb_charge_zero=True)
