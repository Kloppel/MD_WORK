# coding=utf-8

import kbp2
from workspace_jd import cco_config
import os


# subfolder = 'hsp73_hsp526_flu-_2'
#
# # work_subfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/oxi/fixedH_O/'
# # # structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_flu-_2/md_ex_10/frames_1/renamed/frame620.pdb'
# # structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/oxi/frame620_kbp.pdb'
# # special =  None
#
# work_subfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/oxi/fixedH_O_flip/'
# structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/oxi/frame620_kbp.pdb'
# special = 'flip_charge'
#
# work_subfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/oxi/fixedH_O_allres/'
# # # structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_flu-_2/md_ex_10/frames_1/renamed/frame620.pdb'
# structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/oxi/frame620_kbp.pdb'
# special =  None

# work_subfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/oxi/fixedH_O_allres_flip/'
# structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/oxi/frame620_kbp.pdb'
# special = 'flip_charge'


subfolder = 'hsp73_hsp526_red_flu2'
#
# # work_subfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/red/fixedH_R/'
# # # structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_red_flu2/md_ex_10/frames_1/renamed/frame548.pdb'
# # structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/red/frame548_kbp.pdb'
# # special =  None
#
# work_subfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/red/fixedH_R_flip/'
# structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/red/frame548_kbp.pdb'
# special = 'flip_charge'#
# work_subfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/red/fixedH_R_allres/'
# # # structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_red_flu2/md_ex_10/frames_1/renamed/frame548.pdb'
# structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/red/frame548_kbp.pdb'
# special =  None

# work_subfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/red/fixedH_R_allres_flip/'
structure = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/flip_charge/with_flu/90/frames_influence/red/frame548_kbp.pdb'
special = 'flip_charge'

if not os.path.exists(work_subfolder):
    os.mkdir(work_subfolder)
# residues_tt = ['FLU-1_FLUR']
residues_tt = ['HSP-526_ACHA', 'HSP-73_BCHA', 'FLU-1_FLUR']

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
# kbp2_settings.quiet_mode = True

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

kbp2_settings.set_ignored_residues(['histidines'])

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
        rename_flu_deprot = True
        rename = True

        initial_protonation.update({('FLU', 1, 'FLUR'):2})

    if 'fluh' in subfolder:

        initial_protonation.update({('FLU', 1, 'FLUR'):1})

    if 'flu2' in subfolder:
        # renaming must be done from FLX to FLU
        rename_flu_prot = True
        rename = True

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



# FLIP charge of BNC calculations!
if special == 'flip_charge':
    pdb_mod = kbp2.file_parser.Simple_struct_parser()
    pdb_mod.read_pdb(structure)
    if 'red' in subfolder:
        for seg in pdb_mod.struct.iter_segments():
            for res in seg.iter_residues():
                if res.resname == 'HOH':
                    for atm in res.iter_atoms():
                        if atm['segname'] == 'FEOH':
                            atm['resname'] = 'OHMI'
                            if atm['name'] == 'XH2':
                                atm['name'] = 'OH2'
                            if atm['name'] == 'H2':
                                pdb_mod.del_atom(atm, no_restruct=False)

    if 'red' not in subfolder:
        for seg in pdb_mod.struct.iter_segments():
            for res in seg.iter_residues():
                if res.resname == 'OHMI':
                    for atm in res.iter_atoms():
                        if atm['segname'] == 'FEOH':
                            atm['resname'] = 'HOH'
                            if atm['name'] == 'OH2':
                                atm['name'] = 'XH2'

    pdb_mod.create_struct()
    pdb_mod.write_pdb(work_subfolder + 'flipped.pdb')


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
kbp2_settings.remove_folders = 'keep'


# Test for one strcuture -> comment out manager below
if special == 'flip_charge':
    kbp2_settings.set_structure(work_subfolder + 'flipped.pdb')
else:
    kbp2_settings.set_structure(structure)
kbp2_settings.set_workdir(work_subfolder)
#kbp2.pka_calculation.calc_pkas(kbp2_settings, titratable_definitions)
kbp2.pka_cal_jovan.calc_pkas(kbp2_settings, titratable_definitions)




