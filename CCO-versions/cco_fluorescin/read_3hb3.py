# coding=utf-8

import os
import kbp2
from workspace_jd import cco_config


def build_structure(protein, modelling_folder, flu=None, state=None):

    top = []
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf")

    if flu == 'flu2':
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_flu2.rtf")
    elif flu == 'fluh':
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_fluh.rtf")
    elif flu =='flu-':
        top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_flux.rtf")

    par = []
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_fluorescin.str")

    if flu == 'flu2':
        par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_flu2.prm")
    elif flu == 'fluh':
        par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_fluh.prm")
    elif flu =='flu-':
        par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_flux.prm")

    charge_patches = []
    bond_patches = []

    charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
    charmm_struct.add_structure(protein)

    decisions = []
    decisions.append(('rename__HOH_TIP3', 'keep'))
    decisions.append(('disu_bridges', 'closed'))
    decisions.append(('cap_termini', 'dont_cap'))


    if state == 'O2_bound_YH_Fe2_Cu1':
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
        charge_patches.append({'OA3R': ['HEM-2_EHEM', 'HSD-411_ACHA', 'PER-1_PERI']})
        charge_patches.append({'OCB1': [ 'CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'PER-1_PERI', 'TYR-280_ACHA', 'HSE-276_ACHA']})
        charge_patches.append({'OXIR': ['HEM-2_EHEM', 'PER-1_PERI', 'CU1-1_META']})

        bond_patches.append({'PHEM': ['HSD-411_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-94_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-413_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
        bond_patches.append({'CUBO': ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA']})

    if state == 'O':
        ### CHARGE PATCHES
        charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
        charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
        charge_patches.append({'CB4': ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
        charge_patches.append({'A3H3': ['HEM-2_EHEM', 'HSD-411_ACHA', 'OHMI-1_FEOH']})

        ### BOND PATCHES
        bond_patches.append({'PHEM': ['HSD-411_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-94_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-413_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISE': ['OHMI-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUBP': ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})

    if state == 'R':
        ### CHARGE PATCHES
        charge_patches.append({'AHE2': ['HEM-2_GHEM', 'HSD-94_ACHA', 'HSD-413_ACHA']})
        charge_patches.append({'CA11': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})
        charge_patches.append({'CB1T': ['CU1-1_META', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC', 'TYR-280_ACHA', 'HSE-276_ACHA']})
        charge_patches.append({'A3W2': ['HEM-2_EHEM', 'HSD-411_ACHA', 'HOH-1_FEOH']})

        ### BOND PATCHES
        bond_patches.append({'PHEM': ['HSD-411_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM': ['HSD-94_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2': ['HSD-413_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISW': ['HOH-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUBP': ['CU1-1_META', 'HSE-276_ACHA', 'HSD-325_ACHA', 'HSD-326_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-181_BCHA', 'CYS-216_BCHA', 'GLU-218_BCHA', 'CYS-220_BCHA', 'HSE-224_BCHA', 'MET-227_BCHA']})

    charge_patches.append({'LSN': ['LYS-354_ACHA']})

    if flu == 'flu2':
        charge_patches.append({'FCYS': ['CYS-299_ACHA', 'FLH-1_FLUR']})
    elif flu == 'fluh':
        charge_patches.append({'FCYS': ['CYS-299_ACHA', 'FLU-1_FLUR']})
    elif flu =='flu-':
        charge_patches.append({'FXCY': ['CYS-299_ACHA', 'FLX-1_FLUR']})

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

    charmm_struct.set_prot_residue(('GLU', 278, 'ACHA'), patch='GLUP')
    charmm_struct.set_prot_residue(('ASP', 399, 'ACHA'), patch='ASPP')
    charmm_struct.set_prot_residue(('GLU', 481, 'ACHA'), patch='GLUP')

    charmm_struct.set_prot_residue(('HIS', 73, 'BCHA'), rename='HSD')
    charmm_struct.set_prot_residue(('HIS', 526, 'ACHA'), rename='HSE')

    charmm_struct.workdir = modelling_folder

    charmm_struct.check_structures(quiet=True)

    charmm_struct.charmm_instructions['do_minimize'] = False

    charmm_struct.run_charmm()


if __name__ == '__main__':

    # protein = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_flu-_2/switch_to_red/last_oxi_renamed.pdb'
    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_flu-_2/switch_to_red/model/'
    # state = 'R'
    # build_structure(protein, modelling_folder, flu='flu-', state=state)

    protein = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_flu-_2/switch_to_red_deprot/last_oxi_renamed.pdb'
    modelling_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_flu-_2/switch_to_red_deprot/model/'
    state = 'R'
    build_structure(protein, modelling_folder, flu='flu-', state=state)

    # protein = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_red_flu2/switch_to_oxi/last_red_renamed.pdb'
    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/hsp73_hsp526_red_flu2/switch_to_oxi/model/'
    # state = 'O'
    # build_structure(protein, modelling_folder, flu='flu2', state=state)


