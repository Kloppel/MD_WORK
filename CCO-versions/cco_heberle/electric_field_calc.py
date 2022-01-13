# coding=utf-8

import os
import kbp2
import numpy as np
from math import pi, acos
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm


def calc_electric_filed(structure, point):

    # field = []
    field = np.zeros(shape=(1,3))

    # r = np.zeros(shape=(1,3))
    # r = np.array([0.1, 0.2, 0.3])

    r = point

    chargesum = 0

    for i, new_atom in enumerate(structure.atoms):
        resid = new_atom['resid']
        segname = new_atom['segname']
        name = new_atom['name']

        #REMOVE WATER COORDINATED TO FE
        if segname == 'FEOH':
            continue
        # if segname == 'HOHC':
        #     continue

        # q = float(structure.struct[segname][resid][name]['charge']) * 1.6021765 * 10**-19
        q = float(structure.struct[segname][resid][name]['charge'])
        ri = structure.struct[segname][resid][name]['coord']

        if (segname == 'GHEM') or (str(resid) == '102' and new_atom['resname'] == 'HSD') or ((str(resid) == '421' and new_atom['resname'] == 'HSD')):
            chargesum += q

        # r_up = np.subtract(r,ri)
        # r_down = np.linalg.norm(np.subtract(r,ri))
        # fi = q * r_up / (r_down*r_down*r_down)

        r_up = np.subtract(r,ri)
        # r_up = np.subtract(ri,r)
        r_up_mag = np.linalg.norm(r_up)
        unitvector_r  = r_up/r_up_mag
        fi = q * unitvector_r / (r_up_mag*r_up_mag)

        # print field
        # print type(field)
        # print '***'
        # print fi
        # print type(fi)
        # print '----'

        if i == 0:
            field = fi
        else:
            field += fi

    # print chargesum, 'charge of heme a'


    # const = 1/(4*pi*4.0* 8.854187817620*10**-12*10**-20)
    # const = 1.6021765 * (10**-19) /(4*pi*4.0* 8.854187817620*(10**-12)*(10**-20))

    const = 10**10 * 1.6021765 * (10**-19) /(4*pi*4.0* 8.854187817620*(10**-12))

    field = field*const

    magnitude = np.linalg.norm(field)
    # print magnitude, 'electric field in V/Angstrom -> function returns this value'
    # print magnitude*(10**8), 'electric field in V/cm'

    return field, magnitude


def model_structure(modelling_folder, pdb, state_dict):

    top = []
    top.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp")
    top.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/patches.rtf")
    top.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf")
    top.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf")
    par = []
    par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
    par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
    par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
    par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")

    charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
    charmm_struct.add_structure(pdb)
    charmm_struct.add_decision('rename__CA_CAL', 'keep')
    charmm_struct.add_decision('rename__HOH_TIP3', 'keep')
    charmm_struct.add_decision('disu_bridges', 'closed')
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

    # CHARGE PACTHES
    patches = []

    if state_dict['hema'] == 'oxi':
        patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
    elif state_dict['hema'] == 'red':
        patches.append({'AHE2' : ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})

    if state_dict['Glu286_prot'] == 'glup':
        charmm_struct.set_prot_residue(('GLU', 286, 'ACHA'), patch='GLUP')
    elif state_dict['Glu286_prot'] == 'glue':
        charmm_struct.set_prot_residue(('GLU', 286, 'ACHA'), patch='GLUE')

    if state_dict['PRDa3'] == 'prdh':
        charmm_struct.add_patch('PRDH', ['PRD-3_EHEM'])
    elif state_dict['PRDa3'] == 'prd2':
        charmm_struct.add_patch('PRD2', ['PRD-3_EHEM'])

    if state_dict['PRAa3'] == 'prah':
        charmm_struct.add_patch('PRAH', ['PRD-1_EHEM'])
    elif state_dict['PRAa3'] == 'pra2':
        charmm_struct.add_patch('PRA2', ['PRD-1_EHEM'])

    if state_dict['CuA'] == 'oxi':
        patches.append({'CA21' : ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
    elif state_dict['CuA'] == 'red':
        patches.append({'CA11': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})


    # STATE E
    # E:
    # charge patches: A3H3 + CB1T
    # ligands: OHMI at heme a3, HOH at CuB
    # bond patches: EISE + PHEM + CUB2
    #
    # patches.append({'CA21' : ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
    # patches.append({'CB1T' : ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})
    # patches.append({'A3H3' : ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
    # for patch in patches:
    #     patch_name = patch.keys()[0]
    #     residues = patch[patch_name]
    #     charmm_struct.add_patch(patch_name, residues)
    # bond_patches = []
    # bond_patches.append({'PHEM' : ['HSD-419_ACHA', 'HEM-2_EHEM']})
    # bond_patches.append({'PHEM' : ['HSD-102_ACHA', 'HEM-2_GHEM']})
    # bond_patches.append({'PHE2' : ['HSD-421_ACHA', 'HEM-2_GHEM']})
    # bond_patches.append({'EISE' : ['OHMI-1_FEOH', 'HEM-2_EHEM']})
    # bond_patches.append({'CUB2' : ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
    # bond_patches.append({'CUAP' : ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})


    # PR state
    # CBPN 333 prot, CBDD 333 deprot -> you can flip them to obtain results (approx)
    if state_dict['hema3'] == 'Pr':
        print 'recheck histdiine deprotonation!'
        if state_dict['His333'] == 'deprot':
            patches.append({'CBDD' : ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})
        elif state_dict['His334'] == 'deprot':
            patches.append({'CBDD' : ['CU1-1_META', 'HSD-334_ACHA', 'HSD-333_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})
        else:
            patches.append({'CBPN': ['CU1-1_META', 'HSD-334_ACHA', 'HSD-333_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        if state_dict['hema3'] == 'Pr':
            patches.append({'A3H4': ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
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

    elif state_dict['hema3'] == 'E':
        patches.append({'CB1T' : ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})
        patches.append({'A3H3' : ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
        for patch in patches:
            patch_name = patch.keys()[0]
            residues = patch[patch_name]
            charmm_struct.add_patch(patch_name, residues)
        bond_patches = []
        bond_patches.append({'PHEM' : ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM' : ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2' : ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISE' : ['OHMI-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2' : ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP' : ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
#         charmm_command = """
# CONS FIX SELE .not. (resname TYR .and. resid 288 .and. type HH) end
# MINI SD NSTEP 200 INBFRQ 20 TOLG 0.1
#         """
#         charmm_struct.add_charmm_command(charmm_command, adj_task='hbuild')


    elif state_dict['hema3'] == 'E_Tyr-':
        patches.append({'CBT4' : ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})
        patches.append({'A3H3' : ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
        for patch in patches:
            patch_name = patch.keys()[0]
            residues = patch[patch_name]
            charmm_struct.add_patch(patch_name, residues)
        bond_patches = []
        bond_patches.append({'PHEM' : ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM' : ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2' : ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISE' : ['OHMI-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2' : ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP' : ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

    # BOND PACTHES -> HEME AND COPPER
    patches_no_autogen = []
    for patch in bond_patches:
        patches_no_autogen.append(patch)

    for patch in patches_no_autogen:
        patch_name = patch.keys()[0]
        residues = patch[patch_name]
        charmm_struct.add_patch(patch_name, residues, no_autogen=True)

    if modelling_folder[-1] != '/':
        modelling_folder = modelling_folder + '/'
    if not os.path.exists(modelling_folder):
        os.mkdir(modelling_folder)

    charmm_struct.workdir = modelling_folder
    charmm_struct.check_structures(quiet=True)
    charmm_struct.charmm_instructions['do_minimize'] = False
    charmm_struct.run_charmm(submit=False)

    return charmm_struct, modelling_folder

def electric_field(charmm_struct, modelling_folder):

    charmm_ssp = charmm_struct.structure
    premodelled_structure = charmm_struct.get_modelled_structure()
    charmm_ssp = premodelled_structure

    par = []
    par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu_kb.inp")
    par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/patches.prm")
    par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm")
    par.append("/mnt/fu-scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")

    if not charmm_ssp.par_read:
        for par in charmm_struct.par:
            charmm_ssp.read_par(par)

    # #coordinates of C=O -> old model
    # # a- coords of C
    # a = np.array([-3.125, -31.473, 23.398])
    # # b  - coords of O
    # b = np.array([-2.156, -31.832, 23.856])

    # # coordinates of C=O -> new model
    # # a- coords of C
    # a = np.array([-3.554, -31.375, 23.403])
    # # b  - coords of O
    # b = np.array([-3.100, -30.416, 23.792])

    # coordinates of C=O -> linear model
    # a- coords of C
    a = np.array([-3.698, -31.279,  23.484])
    # b  - coords of O
    b = np.array([-3.152, -30.371,  23.878])

    # a = []
    # b = []
    c = []
    for i, new_atom in enumerate(charmm_ssp.atoms):
        resid = new_atom['resid']
        segname = new_atom['segname']
        name = new_atom['name']

        if segname == 'GHEM':
            if name == 'FE':
                c = charmm_ssp.struct[segname][resid][name]['coord']
                # if segname == 'META':
                #     if name == 'MG':
                #         print 'found'
                #         c = charmm_ssp.struct[segname][resid][name]['coord']

    c_o = np.subtract(b, a)
    c_o = b - a
    c_o_mag = np.linalg.norm(c_o)
    unitvector = c_o / c_o_mag
    # point = a + unitvector*2.7
    # point = a + unitvector * 0.514
    point = a + unitvector * 0.565

    # fehea3_c = a - c
    # print np.degrees(np.arccos(dot(c_o,fehea3_c)/norm(c_o)/norm(fehea3_c))), 'vector angle Fehea3-C-O'
    # scalar_projection_random = np.linalg.norm(fehea3_c)*(dot(c_o,fehea3_c)/norm(fehea3_c)/norm(c_o))
    # vector_projection_random = unitvector * scalar_projection_random
    # print np.linalg.norm(scalar_projection_random), 'magnitude random projection'
    # print np.degrees(np.arccos(dot(c_o,vector_projection_random)/norm(c_o)/norm(vector_projection_random))), 'vector proje angle'

    field_vector, field_magnitude = calc_electric_filed(charmm_ssp, point)
    # proof of point where we look for electric field
    for i, new_atom in enumerate(charmm_ssp.atoms):
        resid = new_atom['resid']
        segname = new_atom['segname']
        name = new_atom['name']
        # if segname == 'FEOH':
        #     if name == 'OH2':
        if segname == 'META':
            if str(resid) == '1':
                charmm_ssp.struct[segname][resid][name]['coord'] = point
    charmm_ssp.write_pdb(modelling_folder + 'point.pdb')

    # print field_magnitude*100, 'in MV/cm'
    # print field_vector, 'field'

    return field_vector, field_magnitude, unitvector, c_o


def analyse_electric_fields(folder, unitvector, c_o, E1_field, E2_field, quiet=False):

    correct_sign = False
    filename = folder + 'result.dat'
    log_file = open(filename, 'w')

    angle_unit_E1 = np.degrees(np.arccos(dot(c_o,E1_field)/norm(c_o)/norm(E1_field)))
    if not quiet:
        print angle_unit_E1, 'E1 angle with CO vector'

    angle_unit_E2 = np.degrees(np.arccos(dot(c_o,E2_field)/norm(c_o)/norm(E2_field)))
    if not quiet:
        print angle_unit_E2, 'E2 angle with CO vector'

    if not quiet:
        if not (-1 <= dot(E2_field,E1_field)/norm(E2_field)/norm(E1_field) <= 1):
            print np.degrees(np.arccos(round(dot(E2_field,E1_field)/norm(E2_field)/norm(E1_field),2))), 'E1_E2_angle'
        else:
            print np.degrees(np.arccos(dot(E2_field,E1_field)/norm(E2_field)/norm(E1_field))), 'E1_E2_angle'

    scalar_projection_E1 = np.linalg.norm(E1_field)*(dot(unitvector,E1_field)/norm(unitvector)/norm(E1_field))
    vector_projection_E1 = unitvector * scalar_projection_E1

    if not quiet:
        print np.linalg.norm(vector_projection_E1)*100, 'magnitude E1 projection in MV/cm'
        if not (-1<= dot(c_o,vector_projection_E1)/norm(c_o)/norm(vector_projection_E1) <= 1):
            print np.degrees(np.arccos(round(dot(c_o,vector_projection_E1)/norm(c_o)/norm(vector_projection_E1), 2))), 'vector proje angle E1'
        else:
            print np.degrees(np.arccos(dot(c_o,vector_projection_E1)/norm(c_o)/norm(vector_projection_E1))), 'vector proje angle E1'


    scalar_projection_E2 = np.linalg.norm(E2_field)*(dot(unitvector,E2_field)/norm(unitvector)/norm(E2_field))
    vector_projection_E2= unitvector * scalar_projection_E2

    # #proof
    # diff_e1e2 = E1_field - E2_field
    # scalar_projection_diff = np.linalg.norm(diff_e1e2) * (dot(unitvector, diff_e1e2) / norm(unitvector) / norm(diff_e1e2))
    # vector_projection_diff = unitvector * scalar_projection_diff
    # if np.degrees(np.arccos(round(dot(vector_projection_diff, c_o) / norm(vector_projection_diff) / norm(c_o)))) == 180.0:
    #     print -1*norm(vector_projection_diff)*100, 'final: MV/cm first diff and the priejction'
    # else:
    #     print norm(vector_projection_diff) * 100, 'final: MV/cm first diff and the priejction'

    if not quiet:
        print np.linalg.norm(vector_projection_E2)*100, 'magnitude E2 projection in MV/cm'
        if not (-1<= dot(c_o,vector_projection_E2)/norm(c_o)/norm(vector_projection_E2) <=1):
            print np.degrees(np.arccos(round(dot(c_o,vector_projection_E2)/norm(c_o)/norm(vector_projection_E2), 2))), 'vector proje angle E2'
        else:
            print np.degrees(np.arccos(dot(c_o,vector_projection_E2)/norm(c_o)/norm(vector_projection_E2))), 'vector proje angle E2'

        if not (-1<= dot(vector_projection_E2,vector_projection_E1)/norm(vector_projection_E2)/norm(vector_projection_E1) <=1):
            print np.degrees(np.arccos(round(dot(vector_projection_E2,vector_projection_E1)/norm(vector_projection_E2)/norm(vector_projection_E1),2))), 'E1_E2_angle_projection'
        else:
            print np.degrees(np.arccos(dot(vector_projection_E2,vector_projection_E1)/norm(vector_projection_E2)/norm(vector_projection_E1))), 'E1_E2_angle_projection'

        print np.linalg.norm(vector_projection_E1) - np.linalg.norm(vector_projection_E2), 'magnitude difference in V/Angstrom'

    if np.linalg.norm(vector_projection_E1) - np.linalg.norm(vector_projection_E2) < 0:
        correct_sign = True

    de_new = np.subtract(vector_projection_E1,vector_projection_E2)

    if not quiet:
        if not (-1<= dot(vector_projection_E2,de_new)/norm(vector_projection_E2)/norm(de_new) <=1):
            print np.degrees(np.arccos(round(dot(vector_projection_E2,de_new)/norm(vector_projection_E2)/norm(de_new),2))), 'E2_difference_angle'
        else:
            np.degrees(np.arccos(dot(vector_projection_E2, de_new) / norm(vector_projection_E2) / norm(de_new))), 'E2_difference_angle'
        if not (-1<= dot(vector_projection_E1,de_new)/norm(vector_projection_E1)/norm(de_new) <=1):
            print np.degrees(np.arccos(round(dot(vector_projection_E1,de_new)/norm(vector_projection_E1)/norm(de_new),2))), 'E1_difference_angle'
        else:
            np.degrees(np.arccos(dot(vector_projection_E1, de_new) / norm(vector_projection_E1) / norm(de_new))), 'E1_difference_angle'

    # de = E2_field - E1_field
    de = E1_field - E2_field
    if not quiet:
        if not (-1<= dot(de,c_o)/norm(de)/norm(c_o) <=1):

            print np.degrees(np.arccos(round(dot(de,c_o)/norm(de)/norm(c_o),2))), 'diff_CO_angle'
        else:
            print np.degrees(np.arccos(dot(de,c_o)/norm(de)/norm(c_o))), 'diff_CO_angle'

    if not (-1<= dot(de_new,c_o)/norm(de_new)/norm(c_o) <=1):
        if not quiet:
            print np.degrees(np.arccos(round(dot(de_new,c_o)/norm(de_new)/norm(c_o),2))), 'projection_diff_and_CO_angle'
        p_diff_CO_angle = np.degrees(np.arccos(round(dot(de_new,c_o)/norm(de_new)/norm(c_o),2)))
    else:
        if not quiet:
            print np.degrees(np.arccos(dot(de_new, c_o) / norm(de_new) / norm(c_o)))
        p_diff_CO_angle = np.degrees(np.arccos(dot(de_new, c_o) / norm(de_new) / norm(c_o)))

    if p_diff_CO_angle > 90.0:
        if not correct_sign:
            raise AssertionError('Check sign of the angle again!')

    # project = np.linalg.norm(de)*(dot(unitvector,de)/norm(unitvector)/norm(de))
    # vect = unitvector * project
    # print norm(vect), 'other way, V/Angstrom'
    # de_new = vector_projection_E2 - vector_projection_E1
    # print norm(de_new), 'final: V/Angstrom'

    diff = norm(de_new)*100
    if correct_sign:
        diff = -1*diff
    if not quiet:
        print diff

    log_file.write( str(norm(de_new)*100) + ' in MV/cm \n')

    log_file.write('\n')
    log_file.close()

    return diff


def log_states(folder, state1, state2):

    filename = folder + 'log.dat'
    log_file = open(filename, 'w')

    log_file.write('* state_1: \n')
    for key, value in state1.iteritems():
        log_file.write('%s.....%s \n' % (key, value))

    log_file.write('\n')

    log_file.write('* state_2: \n')
    for key, value in state2.iteritems():
        log_file.write('%s.....%s \n' % (key, value))

    log_file.write('\n')
    log_file.close()


def get_electric_filed_difference(state1, state2, folder, protein):

    if not os.path.exists(folder):
        os.mkdir(folder)

    log_states(folder, state1, state2)

    structure1, modelling_folder_1 = model_structure(folder + 'state_1', protein, state1)
    field_vector1, field_magnitude1, unitvector, c_o = electric_field(structure1, modelling_folder_1)
    E1 = field_magnitude1
    E1_field = field_vector1

    structure2, modelling_folder_2 = model_structure(folder + 'state_2', protein, state2)
    field_vector2, field_magnitude2, unitvector, c_o = electric_field(structure2, modelling_folder_2)
    E2 = field_magnitude2
    E2_field = field_vector2

    difference = analyse_electric_fields(folder, unitvector, c_o, E1_field, E2_field)
    # difference = analyse_electric_fields(folder, unitvector, c_o, E2_field, E1_field) #this was just a test, ignore

    return difference


if __name__=="__main__":

    # print 'new CO model'
    # main_folder = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/results_3ag1/'
    print 'linear CO model'
    main_folder = '/work/jdragelj/projects/cco_E_field/'

    state1_generic = { 'hema'         :  'oxi',
                       'hema3'        :  'Pr',
                       'Glu286_prot'  :  'glu-',
                       'Glu286_conf'  :  'down' ,
                       'PRDa3'        :  'prd-',
                       'PRAa3'        :  'pra-',
                       'His333'       :  'prot',
                       'His334'       :  'prot',
                       'CuA'          :  'oxi',
                      }

    state2_generic = { 'hema'         :  'oxi',
                       'hema3'        :  'Pr',
                       'Glu286_prot'  :  'glu-',
                       'Glu286_conf'  :  'down',
                       'PRDa3'        :  'prd-',
                       'PRAa3'        :  'pra-',
                       'His333'       :  'prot',
                       'His334'       :  'prot',
                       'CuA'          :  'oxi',
                      }

    #### EXPERIMENTS #####
    # print('hema: oxi->red')
    folder = main_folder + 'hemaOtoR/'
    protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
    state1 = dict(state1_generic)
    state2 = dict(state2_generic)
    state2.update({'hema': 'red'})
    diff =  get_electric_filed_difference(state1, state2, folder, protein)
    print'hema: oxi->red, ', diff

    # # print('CuA: oxi->red')
    # folder = main_folder + 'CuAOtoR/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'CuA': 'red'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'CuA: oxi->red,' , diff
    #
    # # print('E286: glu-->glue')
    # folder = main_folder + 'E286-toe/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'Glu286_prot': 'glue'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'E286: glu-->glue, ', diff
    #
    # # print('E286: glu-->glup')
    # folder = main_folder + 'E286-top/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'Glu286_prot': 'glup'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'E286: glu-->glup, ', diff
    #
    # # print('PRD: prd-->prdh')
    # folder = main_folder + 'PRD-toh/'
    # folder = main_folder + 'test1/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'PRDa3': 'prdh'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'PRD: prd-->prdh, ', diff
    #
    # # print('PRD: prd-->prd2')
    # folder = main_folder + 'PRD-to2/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'PRDa3': 'prd2'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'PRD: prd-->prd2, ', diff
    #
    # # print('PRA: pra-->prdh')
    # folder = main_folder + 'PRA-toh/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'PRAa3': 'prah'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'PRA: pra-->prah, ', diff
    #
    # # print('PRA: prd-->prd2')
    # folder = main_folder + 'PRA-to2/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'PRAa3': 'pra2'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'PRA: pra-->pra2, ', diff
    #
    # # print('His334: prot->deprot')
    # folder = main_folder + 'His334ptod/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'His334': 'deprot'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'His334: prot->deprot, ', diff
    #
    # # print('His333: prot->deprot')
    # folder = main_folder + 'His333ptod/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'His333': 'deprot'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'His333: prot->deprot, ', diff
    #
    # ##### E286 Up ######
    #
    # # print('E286 up: glu-->glue')
    # folder = main_folder + 'E286UP-toe/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/flipPr/cco_and_crystal_water_PR_out.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'Glu286_prot': 'glue'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'E286up: glu-->glue, ', diff
    #
    # # print('E286 up: glu-->glup')
    # folder = main_folder + 'E286UP-top/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/flipPr/cco_and_crystal_water_PR_out.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # state2.update({'Glu286_prot': 'glup'})
    # diff =  get_electric_filed_difference(state1, state2, folder, protein)
    # print'E286up: glu-->glup, ', diff
    #
    #
    # #extra TYR deprot
    # # print 'new CO model'
    # # main_folder = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/results_3ag1/'
    # print 'linear CO model'
    # main_folder = '/work/jdragelj/projects/cco_E_field/projects/cco_heberle/results_linear/'
    #
    # state1_generic = {'hema': 'oxi',
    #                   'hema3': 'E_Tyr-',
    #                   'Glu286_prot': 'glu-',
    #                   'Glu286_conf': 'down',
    #                   'PRDa3': 'prd-',
    #                   'PRAa3': 'pra-',
    #                   'His333': 'prot',
    #                   'His334': 'prot',
    #                   'CuA': 'oxi',
    #                   }
    #
    # state2_generic = {'hema': 'oxi',
    #                   'hema3': 'E',
    #                   'Glu286_prot': 'glu-',
    #                   'Glu286_conf': 'down',
    #                   'PRDa3': 'prd-',
    #                   'PRAa3': 'pra-',
    #                   'His333': 'prot',
    #                   'His334': 'prot',
    #                   'CuA': 'oxi',
    #                   }
    #
    # #### EXPERIMENTS #####
    # # print('hema: oxi->red')
    # folder = main_folder + 'tyr_hack_E/'
    # protein = '/mnt/fu-scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_Glu286_flip.pdb'
    # state1 = dict(state1_generic)
    # state2 = dict(state2_generic)
    # diff = get_electric_filed_difference(state1, state2, folder, protein)
    # print'tyr: deprot->prot, ', diff









