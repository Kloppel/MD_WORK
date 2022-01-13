 # coding=utf-8

import kbp2

protein = '/scratch/scratch/jdragelj/projects/cco_model_pls_prd/titrate/cco_3a3_2Cub_HOH.pdb'
modelling_folder = '/scratch/scratch/jdragelj/projects/cco_model_pls_prd/read_1M56/'
state = '21Cua_3hea_3hea3_2Cub'


top = []
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp")
top.append("/scratch/scratch/tmeyer/karlsbergplus/patches.rtf")
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf")
par = []
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
par.append("//scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")

charge_patches = []
bond_patches = []

charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
charmm_struct.add_structure(protein)

decisions = []
decisions.append(('rename__HOH_TIP3', 'keep'))
decisions.append(('disu_bridges', 'closed'))
decisions.append(('cap_termini', 'dont_cap'))

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

charmm_struct.set_prot_residue(('GLU', 286, 'A'), patch='GLUP')
charmm_struct.set_prot_residue(('LYS', 362, 'A'), charge=0)
charmm_struct.set_prot_residue(('ASP', 407, 'A'), patch='ASPP')

charmm_struct.workdir = modelling_folder

charmm_struct.check_structures(quiet=True)

charmm_struct.charmm_instructions['do_minimize'] = False

charmm_struct.run_charmm()