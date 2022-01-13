# coding=utf-8

import os
import kbp2
from workspace_jd import cco_config


def build_structure(protein, modelling_folder, state=None):
    top = []
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw_jd_oxygen.inp")
    # top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp")
    # top.append("/scratch/scratch/tmeyer/karlsbergplus/patches.rtf")
    top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
    # top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf")
    par = []
    # par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu.inp")
    # par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
    par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
    # par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")

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

    charge_patches.append({'LSN': ['LYS-354_ACHA']})


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

    charmm_struct.workdir = modelling_folder

    charmm_struct.check_structures(quiet=True)

    charmm_struct.charmm_instructions['do_minimize'] = False

    charmm_struct.run_charmm()


if __name__ == '__main__':

    # protein = '/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/read_first/cco_o2.pdb'
    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/read_first/'

    # protein = '/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/split_new_coord_CB/cco_o2.pdb'
    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/split_new_coord_CB/'

    protein = '/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/minimize_CB/model.pdb'
    modelling_folder = '/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/minimize_CB/read_modelO2/'

    state = 'O2_bound_YH_Fe2_Cu1'
    build_structure(protein, modelling_folder, state)

    # protein = '/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/crystal_o2.pdb'
    # modelling_folder = '/scratch/scratch/jdragelj/projects/cco_raman_models/final_model/boundO2_YH/O2H/model_3HB3/read_crystalO2/'
    #
    # state = 'O2_bound_YH_Fe2_Cu1'
    # build_structure(protein, modelling_folder, state)


