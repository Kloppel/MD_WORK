# -*- coding: utf-8 -*-
import kbp2
import os



# pdb_filename = '/scratch/scratch/tmeyer/md_pka/mod_bio/2ibx_m.pdb1'
# modelling_dir = '/scratch/scratch/tmeyer/md_pka/mod_bio/2ibx_m_charmm/'
# modelling_dir = '/scratch/scratch/tmeyer/md_pka/mod_bio/2ibx_m_charmm2/'
pdb_filename = '/home/pbuser/Desktop/PhD_WORK/MD_WORK/cco_structures/ParacoccusDenitrificans/1ar1.pdb'
modelling_dir = '/home/pbuser/Desktop/PhD_WORK/MD_WORK/cco_structures/ParacoccusDenitrificans/x-ray/O-state/1ar1/'
# pdb_filename = '/scratch/scratch/tmeyer/md_pka/mod_bio/1hng.pdb1'
# modelling_dir = '/scratch/scratch/tmeyer/md_pka/mod_bio/1hng_charmm/'





if not os.path.exists(modelling_dir):
    os.mkdir(modelling_dir)

##########################################
### Read structure and model hydrogens ###
##########################################
s = kbp2.file_parser.Simple_struct_parser()
s.read_pdb(pdb_filename)

#s = s.copy(residuelist=[('A', 0), ('B', 0)], exclude=True)
#s.shift_segment('A', 1)
#s.shift_segment('B', 1)


#top = []
#top.append('/scratch/scratch/tmeyer/karlsbergplus/top.inp')
#top.append('/scratch/scratch/tmeyer/karlsbergplus/patches.rtf')
#par = []
#par.append('/scratch/scratch/tmeyer/karlsbergplus/par.inp')
#par.append('/scratch/scratch/tmeyer/karlsbergplus/patches.prm')


top = []
top.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/top_alw_clean_kb.inp")
top.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/patches.rtf")
top.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/top_all36_lipid.rtf")
top.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/top_all36_cgenff.rtf")

par = []
par.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/par_all22_prot_plus_heme_and_Cu_kb.inp")
par.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/patches.prm")
par.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/par_all36_lipid.prm")
par.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/par_all36_cgenff.prm")

c = kbp2.charmm.Charmm_manager(workdir=modelling_dir, top=top, par=par)
c.add_structure(s)
c.charmm_bin = "charmm"

#for segname in ['A', 'B']:
#    first_resid = c.structure.struct[segname].resids()[0]
#    first_residue = c.structure.struct[segname][first_resid]
#    last_resid = c.structure.struct[segname].resids()[-1]
#    last_residue = c.structure.struct[segname][last_resid]
#    c.add_sequence(first_residue, 'GLY', is_last_residue=True)
#    c.add_sequence(last_residue, 'GLU')

# 2zta
# The N terminus of the molecule is Gly0. This residue is an artifact of the engineered thrombin cleavage site used
# to purify the protein.26 The C terminus of the molecule is Glu32. The pKa value for Glu32 is for the side-chain.

# 2zta
# REMARK 465     GLU A    32
# REMARK 465     ARG A    33
# REMARK 465     GLU B    32
# REMARK 465     ARG B    33

# 2ibx
# REMARK 465     LEU A     1
# REMARK 465     VAL A     2
# REMARK 465     LYS A     3
# REMARK 465     SER A     4
# REMARK 465     GLN A   326
# REMARK 465     ARG A   327
# REMARK 465     GLU A   328

c.add_decision('rename__HOH_TIP3', 'keep')
c.add_decision('cap_termini', 'dont_cap')
c.check_structures(quiet=False)

done = c.run_charmm()




def titrate_structure(pdb_structure, workfolder, state, residues_tt=None):
    #############
    ### CALCS ###
    #############
    #top = topology
    top = []
    top.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/top_alw_clean_kb.inp")
    top.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/patches.rtf")
    top.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/top_all36_lipid.rtf")
    top.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/top_all36_cgenff.rtf")
    #par = parameters
    par = []
    par.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/par_all22_prot_plus_heme_and_Cu_kb.inp")
    par.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/patches.prm")
    par.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/par_all36_lipid.prm")
    par.append("/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/par_all36_cgenff.prm")
    #Add Settings
    kbp2_settings = kbp2.pka_calculation.PkaCalcSettings()
    kbp2_settings.top = top
    kbp2_settings.par = par
    kbp2_settings.processes = 1
    kbp2_settings.quiet_mode = True
    decisions = []
    #decisions.append(('rename__HOH_TIP3', 'keep'))
    decisions.append(('rename__HOH_TIP3', 'remove'))
    decisions.append(('disu_bridges', 'closed'))
    decisions.append(('cap_termini', 'dont_cap'))
    kbp2_settings.modelling_decisions = decisions
    protocol = []
    protocol.append((  7, 'h_min'))
    kbp2_settings.protocol = protocol
    # kbp2_settings.cavity_par = [0.8, 0.2, 0.0]
    # kbp2_settings.titr_residue_charges_0 = True
    titratable_yaml = '/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/toppar27/titratable_cco_noglupasp.yaml'
    kbp2_settings.set_yaml(titratable_yaml)
    titratable_definitions = kbp2.kbp_tools.parse_titratable_yaml(titratable_yaml)
    kbp2_settings.patches = []
    kbp2_settings.patches_no_autogen = []
    charge_patches = []
    bond_patches = []
    if state == 'O':
        ### CHARGE PATCHES
        charge_patches.append({'AHE3' : ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
        charge_patches.append({'CA21' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
        charge_patches.append({'CB4' : ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
        charge_patches.append({'A3H3' : ['HEM-2_EHEM', 'HSD-411_ACHA', 'OHMI-1_FEOH']})
        ### BOND PATCHES
        bond_patches.append({'PHEM' : ['HSD-411_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM' : ['HSD-94_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2' : ['HSD-413_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISE' : ['OHMI-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUBP' : ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
    if state == 'R':
        ### CHARGE PATCHES
        charge_patches.append({'AHE2' : ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
        charge_patches.append({'CA11' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
        charge_patches.append({'CB1T' : ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
        charge_patches.append({'A3W2' : ['HEM-2_EHEM', 'HSD-411_ACHA', 'HOH-1_FEOH']})
        ### BOND PATCHES
        bond_patches.append({'PHEM' : ['HSD-411_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM' : ['HSD-94_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2' : ['HSD-413_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISW' : ['HOH-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUBP' : ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
        # temporary mix -> flip charge only for CuA and hemeA
    if state == 'Cua_hemea_R_hemea3_Cub_O':
        ### CHARGE PATCHES
        charge_patches.append({'AHE2': ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
        charge_patches.append({'CA11': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA','EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
        charge_patches.append({'CB4' : ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
        charge_patches.append({'A3H3' : ['HEM-2_EHEM', 'HSD-411_ACHA', 'OHMI-1_FEOH']})
        ### BOND PATCHES
        bond_patches.append({'PHEM': ['HSD-411_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-94_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-413_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISW': ['HOH-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUBP': ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA','EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
    if state == 'Cua_hemea_O_hemea3_Cub_R':
        ### CHARGE PATCHES
        charge_patches.append({'AHE3' : ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
        charge_patches.append({'CA21' : ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
        charge_patches.append({'CB1T' : ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
        charge_patches.append({'A3W2' : ['HEM-2_EHEM', 'HSD-411_ACHA', 'HOH-1_FEOH']})
        ### BOND PATCHES
        bond_patches.append({'PHEM': ['HSD-411_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-94_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-413_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISW': ['HOH-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUBP': ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA','EPP-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
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
    return

#if __name__ == '__main__':
#    workfolder = '/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/test/titrate_3HB3/'
#    pdb_structure = '/home/pbuser/Desktop/PhD_WORK/MD_WORK/Karlsberg/test_example/test/titrate_3HB3/3hb3.pdb'
#    #state = '21Cua_3hea_3hea3_2Cub'
#    state= 'O'
#    titrate_structure(pdb_structure=pdb_structure, workfolder=workfolder, state=state, residues_tt=None)
