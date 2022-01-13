# coding=utf-8
import numpy as np
import kbp2

# protein = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_red_flu2/change/cco_and_water.pdb'
protein = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_flu2/change/cco_and_water.pdb'
# modelling_folder = '/scratch/scratch/jdragelj/tests/cco_flu_red/'
modelling_folder = '/scratch/scratch/jdragelj/tests/cco_flu_oxi/'

top = []
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp")
top.append("/scratch/scratch/tmeyer/karlsbergplus/patches.rtf")
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf")
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_flu2.rtf")
par = []
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
par.append("//scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_flu2.prm")

charge_patches = []
bond_patches = []

charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
charmm_struct.add_structure(protein)

decisions = []
decisions.append(('rename__HOH_TIP3', 'keep'))
decisions.append(('disu_bridges', 'closed'))
decisions.append(('cap_termini', 'dont_cap'))

# ###### RED
# ### CHARGE PATCHES
# charge_patches.append({'AHE2': ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
# charge_patches.append({'CA11': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA',
#                                      'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
# charge_patches.append(
#     {'CB1T': ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
# charge_patches.append({'A3W2': ['HEM-2_EHEM', 'HSD-411_ACHA', 'HOH-1_FEOH']})
#
# ### BOND PATCHES
# bond_patches.append({'PHEM': ['HSD-411_ACHA', 'HEM-2_EHEM']})
# bond_patches.append({'PHEM': ['HSD-94_ACHA', 'HEM-2_GHEM']})
# bond_patches.append({'PHE2': ['HSD-413_ACHA', 'HEM-2_GHEM']})
# bond_patches.append({'EISW': ['HOH-1_FEOH', 'HEM-2_EHEM']})
# bond_patches.append({'CUBP': ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
# bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA',
#                                    'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
#### OXI
### CHARGE PATCHES
charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA',
                                     'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
charge_patches.append(
    {'CB4': ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
charge_patches.append({'A3H3': ['HEM-2_EHEM', 'HSD-411_ACHA', 'OHMI-1_FEOH']})

### BOND PATCHES
bond_patches.append({'PHEM': ['HSD-411_ACHA', 'HEM-2_EHEM']})
bond_patches.append({'PHEM': ['HSD-94_ACHA', 'HEM-2_GHEM']})
bond_patches.append({'PHE2': ['HSD-413_ACHA', 'HEM-2_GHEM']})
bond_patches.append({'EISE': ['OHMI-1_FEOH', 'HEM-2_EHEM']})
bond_patches.append({'CUBP': ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA',
                                   'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})

charge_patches.append({'LSN': ['LYS-354_ACHA']})
charmm_struct.set_prot_residue(('GLU', 278, 'ACHA'), patch='GLUP')
charmm_struct.set_prot_residue(('ASP', 399, 'ACHA'), patch='ASPP')
charmm_struct.set_prot_residue(('GLU', 481, 'ACHA'), patch='GLUP')
charge_patches.append({'FCYS': ['CYS-299_ACHA','FLH-1_FLUR']})

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

charmm_struct.workdir = modelling_folder
charmm_struct.check_structures(quiet=True)
charmm_struct.charmm_instructions['do_minimize'] = False
charmm_struct.run_charmm()
charmm_struct.update_charges_with_modelled_structure()

structure = charmm_struct.get_modelled_structure()
structure = charmm_struct.get_modelled_structure()

chargesum_a = 0
chargesum_cub = 0
chargesum_a3 = 0
chargesum_cua = 0

for i, new_atom in enumerate(structure.atoms):
    resid = new_atom['resid']
    segname = new_atom['segname']
    name = new_atom['name']

    q = float(structure.struct[segname][resid][name]['charge'])

    if (segname == 'GHEM') or (str(resid) == '94' and new_atom['resname'] == 'HSD') or ((str(resid) == '413' and new_atom['resname'] == 'HSD')):
        chargesum_a += q

    if ((str(resid) == '3' and new_atom['resname'] == 'CU1' and new_atom['segname'] == 'META') or
        (str(resid) == '2' and new_atom['resname'] == 'CU1' and new_atom['segname'] == 'META') or
        (str(resid) == '181' and new_atom['resname'] == 'HSE' and new_atom['segname'] == 'BCHA') or
        (str(resid) == '216' and new_atom['resname'] == 'CYS' and new_atom['segname'] == 'BCHA') or
        (str(resid) == '218' and new_atom['resname'] == 'GLU' and new_atom['segname'] == 'BCHA') or
        (str(resid) == '220' and new_atom['resname'] == 'CYS' and new_atom['segname'] == 'BCHA') or
        (str(resid) == '224' and new_atom['resname'] == 'HSE' and new_atom['segname'] == 'BCHA') or
        (str(resid) == '227' and new_atom['resname'] == 'MET' and new_atom['segname'] == 'BCHA')):
            chargesum_cua += q

    if ((str(resid) == '1' and new_atom['resname'] == 'CU1' and new_atom['segname'] == 'META') or
        (str(resid) == '1' and new_atom['resname'] == 'HOH' and new_atom['segname'] == 'HOHC') or
        (str(resid) == '325' and new_atom['resname'] == 'HSD' and new_atom['segname'] == 'ACHA') or
        (str(resid) == '326' and new_atom['resname'] == 'HSD' and new_atom['segname'] == 'ACHA') or
        (str(resid) == '280' and new_atom['resname'] == 'TYR' and new_atom['segname'] == 'ACHA') or
        (str(resid) == '276' and new_atom['resname'] == 'HSE' and new_atom['segname'] == 'ACHA')):
            chargesum_cub += q

    # if (segname == 'EHEM') or \
    #     (str(resid) == '1' and new_atom['resname'] == 'HOH' and new_atom['segname'] == 'FEOH') or
    #     ((str(resid) == '411' and new_atom['resname'] == 'HSD')):
    #         chargesum_a3 += q

    if (segname == 'EHEM') or \
        (str(resid) == '1' and new_atom['resname'] == 'OHMI' and new_atom['segname'] == 'FEOH') or \
        ((str(resid) == '411' and new_atom['resname'] == 'HSD'and new_atom['segname'] == 'ACHA')):
            chargesum_a3 += q

print chargesum_a, 'charge of heme a'
print chargesum_a3, 'charge of heme a3'
print chargesum_cua, 'charge of heme cua'
print chargesum_cub, 'charge of heme cub'



# DIRECTLY FROM TOPOLOGY
# print 'oxi'
# patches = ['AHE3', 'CA21', 'CB4', 'A3H3']
# top_file = '/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp'
# for patch_name in patches:
#     found = False
#     sum = 0
#     f = open(top_file, 'r')
#     for line in f:
#         if not found:
#             if patch_name in line:
#                 print patch_name
#                 found = True
#         else:
#             k = line.split()
#             if k and k[0] == 'ATOM':
#                 sum += float(k[3])
#             if 'PRES' in k:
#                 found = False
#     print sum
#
# print ' '
# print 'red'
# patches = ['AHE2', 'CA11', 'CB1T', 'A3W2']
# top_file = '/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp'
# for patch_name in patches:
#     found = False
#     sum = 0
#     f = open(top_file, 'r')
#     for line in f:
#         if not found:
#             if patch_name in line:
#                 print patch_name
#                 found = True
#         else:
#             k = line.split()
#             if k and k[0] == 'ATOM':
#                 sum += float(k[3])
#             if not k:
#                 found = False
#     print sum


