# coding=utf-8

import os
import kbp2



# protein = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hse526_flu-/protonate/ccoHSP.pdb'
# modelling_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hse526_flu-/protonate/'

# protein = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_o_to_red_flu-/change/cco_o_to_red.pdb'
# modelling_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_o_to_red_flu-/change/'

protein = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_oxi_deprotonate_h_to_flu-/prep/cco_flu-.pdb'
# modelling_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_oxi_deprotonate_h_to_flu-/change/'
modelling_folder = '/scratch/scratch/jdragelj/tests/test_sasa/'


if not os.path.exists(modelling_folder):
    os.mkdir(modelling_folder)

top = []
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp")
top.append("/scratch/scratch/tmeyer/karlsbergplus/patches.rtf")
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf")
top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_flux.rtf")
# top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_flu2.rtf")
# top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_fluh.rtf")
par = []
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
par.append("//scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")
par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_flux.prm")
# par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_flu2.prm")
# par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_fluh.prm")

patches = []

charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
charmm_struct.add_structure(protein)
charmm_struct.add_decision('rename__CA_CAL', 'keep')
charmm_struct.add_decision('rename__HOH_TIP3', 'keep')
charmm_struct.add_decision('disu_bridges', 'closed')

charmm_struct.add_patch('FXCY', ['CYS-299_ACHA', 'FLX-1_FLUR'])
# charmm_struct.add_patch('FCYS', ['CYS-299_ACHA', 'FLH-1_FLUR'])

charmm_struct.set_prot_residue(('LYS', 354, 'ACHA'), charge=0)
charmm_struct.set_prot_residue(('GLU', 278, 'ACHA'), charge=0)
charmm_struct.set_prot_residue(('ASP', 399, 'ACHA'), charge=0)
charmm_struct.set_prot_residue(('GLU', 481, 'ACHA'), charge=0)

charge_patches= []
bond_patches = []


charge_patches.append({'AHE3' : ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
charge_patches.append({'CA21' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
charge_patches.append({'CB4' : ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
charge_patches.append({'A3H3' : ['HEM-2_EHEM', 'HSD-411_ACHA', 'OHMI-1_FEOH']})
### BOND PATCHES
bond_patches.append({'PHEM' : ['HSD-411_ACHA', 'HEM-2_EHEM']})
bond_patches.append({'PHEM' : ['HSD-94_ACHA', 'HEM-2_GHEM']})
bond_patches.append({'PHE2' : ['HSD-413_ACHA', 'HEM-2_GHEM']})
bond_patches.append({'EISE' : ['OHMI-1_FEOH', 'HEM-2_EHEM']})
bond_patches.append({'CUBP' : ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
bond_patches.append({'CUAP' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})

# ### CHARGE PATCHES
# charge_patches.append({'AHE2' : ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
# charge_patches.append({'CA11' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
# charge_patches.append({'CB1T' : ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
# charge_patches.append({'A3W2' : ['HEM-2_EHEM', 'HSD-411_ACHA', 'HOH-1_FEOH']})
# ### BOND PATCHES
# bond_patches.append({'PHEM' : ['HSD-411_ACHA', 'HEM-2_EHEM']})
# bond_patches.append({'PHEM' : ['HSD-94_ACHA', 'HEM-2_GHEM']})
# bond_patches.append({'PHE2' : ['HSD-413_ACHA', 'HEM-2_GHEM']})
# bond_patches.append({'EISW' : ['HOH-1_FEOH', 'HEM-2_EHEM']})
# bond_patches.append({'CUBP' : ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
# bond_patches.append({'CUAP' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})



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

charmm_struct.charmm_instructions['cap_termini'] = False
# charmm_struct.charmm_instructions['ignore_termini'] = True
charmm_struct.workdir = modelling_folder

sasa = 1.4
sasa_block = """
coor surf select all end acce rpro %.1f """ % sasa + """
define backbone select type n .or. type ca .or. type c end
SCALAR WMAIN STAT sele ((segid ACHA .and. resid 70) .and. .not. backbone) end
echo ?ELEC ?stot

"""
charmm_struct.add_charmm_command(sasa_block, adj_task='hbuild')



charmm_struct.check_structures(quiet=True)
charmm_struct.charmm_instructions['do_minimize'] = False
charmm_struct.run_charmm(submit=False)