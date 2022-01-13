# coding=utf-8

import os
import kbp2


protein = '/scratch/scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
modelling_folder = '/scratch/scratch/jdragelj/projects/cco_heberle/modelled_files/flipPr/'


if not os.path.exists(modelling_folder):
    os.mkdir(modelling_folder)

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


patches = []

charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
charmm_struct.add_structure(protein)
charmm_struct.add_decision('rename__CA_CAL', 'keep')
charmm_struct.add_decision('rename__HOH_TIP3', 'keep')
charmm_struct.add_decision('disu_bridges', 'closed')
# charmm_struct.add_decision('cap_termini', 'ignore_termini')
charmm_struct.add_decision('cap_termini', 'dont_cap')

state = {'charge': 0,
         'patch': 'GLUE',
         'external_patches': None,
         'rename': None,
         'special': None}
titr_residues = charmm_struct.get_titr_residues()
titr_residues['GLU'].append(state)
charmm_struct.set_titr_residues(titr_residues)

charmm_struct.set_prot_residue(('LYS', 362, 'ACHA'), charge=0)
charmm_struct.set_prot_residue(('ASP', 407, 'ACHA'), patch='ASPP')
charmm_struct.set_prot_residue(('GLU', 286, 'ACHA'), patch='GLUE')

charge_patches= []
bond_patches = []

print 'recheck histdiine deprotonation!'
patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
patches.append({'A3H4': ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
patches.append({'CBPN': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

for patch in patches:
    patch_name = patch.keys()[0]
    residues = patch[patch_name]
    charmm_struct.add_patch(patch_name, residues)

bond_patches = []
bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
bond_patches.append({'EISO': ['OHMI-1_FEOH', 'HEM-2_EHEM']})
bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                              'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})


for patch in charge_patches:
    patch_name = patch.keys()[0]
    residues = patch[patch_name]
    charmm_struct.add_patch(patch_name, residues)

#BOND PACTHES -> HEME AND COPPER
patches_no_autogen = []
for patch in bond_patches:
    patches_no_autogen.append(patch)

for patch in patches_no_autogen:
    patch_name = patch.keys()[0]
    residues = patch[patch_name]
    charmm_struct.add_patch(patch_name, residues, no_autogen=True)

block = """
!coor init sele resid 286 .and. resname glu end
coor init sele resname GLU .and. resid 286 .and. .not. (type N .or. type A .or. type O .or. type C) end

IC EDIT
DIHE ACHA 286 CD ACHA 286 CG ACHA 286 CB ACHA 286 CA -150
DIHE ACHA 286 OE1 ACHA 286 CD ACHA 286 CG ACHA 286 CB 60
END

ic build
hbuild

CONS FIX SELE .not. (resname GLU .and. resid 286 .and. .not. (type N .or. type A .or. type O .or. type C)) end
MINI SD NSTEP 800 INBFRQ 20 TOLG 0.1
"""
charmm_struct.add_charmm_command(block, adj_task='hbuild')


charmm_struct.workdir = modelling_folder

charmm_struct.check_structures(quiet=True)
charmm_struct.charmm_instructions['do_minimize'] = False
charmm_struct.run_charmm(submit=False)
