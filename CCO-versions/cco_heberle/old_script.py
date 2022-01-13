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

    print chargesum, 'charge of heme a'


    # const = 1/(4*pi*4.0* 8.854187817620*10**-12*10**-20)
    # const = 1.6021765 * (10**-19) /(4*pi*4.0* 8.854187817620*(10**-12)*(10**-20))

    const = 10**10 * 1.6021765 * (10**-19) /(4*pi*4.0* 8.854187817620*(10**-12))

    field = field*const

    magnitude = np.linalg.norm(field)
    # print magnitude, 'electric field in V/Angstrom -> function returns this value'
    # print magnitude*(10**8), 'electric field in V/cm'

    return field, magnitude



if __name__=="__main__":

    E_oxi = 0
    E_red = 0
    for state_hema in ['oxi_hema', 'red_hema']:

        modelling_folder = '/scratch/scratch/jdragelj/projects/cco_heberle/tests/up_model_cco_rspehroides_E_crystal_%s' % state_hema

        # protein = '/media/jdragelj/88bc6e86-bfdd-41a3-b3f0-9921ae91f5b7/awoelke/awoelke/md_glu101b/kb_crystal/crystal/cco_and_crystal_water.pdb'
        protein = '/scratch/scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_water_PR.pdb'
        # protein = '/scratch/scratch/jdragelj/projects/cco_heberle/modelled_files/cco_and_crystal_Glu286_flip.pdb'


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


        state = {'charge' : 0,
                  'patch' : 'GLUE',
                  'external_patches' : None,
                  'rename' : None,
                  'special' : None}
        titr_residues = charmm_struct.get_titr_residues()
        titr_residues['GLU'].append(state)
        charmm_struct.set_titr_residues(titr_residues)

        charmm_struct.set_prot_residue(('LYS', 362, 'ACHA'), charge=0)
        charmm_struct.set_prot_residue(('ASP', 407, 'ACHA'), patch='ASPP')

        #STATE E
          # E:
            # charge patches: A3H3 + CB1T
            # ligands: OHMI at heme a3, HOH at CuB
            # bond patches: EISE + PHEM + CUB2

        #CHARGE PACTHES
        patches = []

        folder_add = ''

        if state_hema == 'oxi_hema':
            print 'HEME-a -> oxidized state'
            patches.append({'AHE3' : ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
            # print 'HEME-a -> reduced state'
            # patches.append({'AHE2' : ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
            # print 'glup'
            # charmm_struct.set_prot_residue(('GLU', 286, 'ACHA'), patch='GLUP')
            # folder_add += '_glup'
            # print 'glue'
            # charmm_struct.set_prot_residue(('GLU', 286, 'ACHA'), patch='GLUE')
            # folder_add += '_glue'
            # print 'PRDH'
            # charmm_struct.add_patch('PRDH', ['PRD-3_EHEM'])
            # folder_add += '_prdh'
            # print 'PRD2'
            # charmm_struct.add_patch('PRD2', ['PRD-3_EHEM'])
            # folder_add += '_prd2'
            # print 'PRAH'
            # charmm_struct.add_patch('PRAH', ['PRD-1_EHEM'])
            # folder_add += '_prah'
            # print 'PRA2'
            # charmm_struct.add_patch('PRA2', ['PRD-1_EHEM'])
            # folder_add += '_pra2'
            print 'CuA oxi'
            patches.append({'CA21' : ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
            folder_add += '_CuA_oxi'
            # print 'CuA red'
            # patches.append({'CA11' : ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
            # folder_add += '_CuA_red'



        if state_hema == 'red_hema':
            print 'HEME-a -> reduced state'
            patches.append({'AHE2' : ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
            # print 'HEME-a -> oxidized state'
            # patches.append({'AHE3' : ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
            # print 'glup'
            # charmm_struct.set_prot_residue(('GLU', 286, 'ACHA'), patch='GLUP')
            # folder_add += '_glup'
            # print 'glue'
            # charmm_struct.set_prot_residue(('GLU', 286, 'ACHA'), patch='GLUE')
            # folder_add += '_glue'
            # print 'PRDH'
            # charmm_struct.add_patch('PRDH', ['PRD-3_EHEM'])
            # folder_add += '_prdh'
            # print 'PRD2'
            # charmm_struct.add_patch('PRD2', ['PRD-3_EHEM'])
            # folder_add += '_prd2'
            # print 'PRAH'
            # charmm_struct.add_patch('PRAH', ['PRD-1_EHEM'])
            # folder_add += '_prah'
            # print 'PRA2'
            # charmm_struct.add_patch('PRA2', ['PRD-1_EHEM'])
            # folder_add += '_pra2'
            print 'CuA oxi'
            patches.append({'CA21' : ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
            folder_add += 'CuA_oxi'
            # print 'CuA red'
            # patches.append({'CA11' : ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
            # folder_add += '_CuA_red'


        modelling_folder += folder_add
        modelling_folder = modelling_folder + '/'
        if not os.path.exists(modelling_folder):
            os.mkdir(modelling_folder)

        # #E state
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


        #PR state
        #todo: CBPN 333 prot, CBDD 333 deprot -> you can flip them to obtain results (approx)

        # print 'CBPN-333'
        # folder_add += '_His333H'
        # patches.append({'CBPN' : ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})
        # print 'CBDD-333'
        # folder_add += '_His333-'
        # patches.append({'CBDD' : ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        print 'CBPN-334'
        folder_add += '_His334H'
        patches.append({'CBPN' : ['CU1-1_META', 'HSD-334_ACHA', 'HSD-333_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})
        # print 'CBDD-334'
        # folder_add += '_His334-'
        # patches.append({'CBDD' : ['CU1-1_META', 'HSD-334_ACHA', 'HSD-333_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

        patches.append({'A3H4' : ['HEM-2_EHEM', 'HSD-419_ACHA', 'OHMI-1_FEOH']})
        for patch in patches:
            patch_name = patch.keys()[0]
            residues = patch[patch_name]
            charmm_struct.add_patch(patch_name, residues)
        bond_patches = []
        bond_patches.append({'PHEM' : ['HSD-419_ACHA', 'HEM-2_EHEM']})
        bond_patches.append({'PHEM' : ['HSD-102_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'PHE2' : ['HSD-421_ACHA', 'HEM-2_GHEM']})
        bond_patches.append({'EISO' : ['OHMI-1_FEOH', 'HEM-2_EHEM']})
        bond_patches.append({'CUB2' : ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
        bond_patches.append({'CUAP' : ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA', 'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})


        #BOND PACTHES -> HEME AND COPPER
        patches_no_autogen = []
        for patch in bond_patches:
            patches_no_autogen.append(patch)

        for patch in patches_no_autogen:
            patch_name = patch.keys()[0]
            residues = patch[patch_name]
            charmm_struct.add_patch(patch_name, residues, no_autogen=True)


        charmm_struct.workdir = modelling_folder

        charmm_struct.check_structures(quiet=True)
        charmm_struct.charmm_instructions['do_minimize'] = False
        charmm_struct.run_charmm(submit=False)


        charmm_ssp = charmm_struct.structure
        premodelled_structure = charmm_struct.get_modelled_structure()
        charmm_ssp = premodelled_structure

        if not charmm_ssp.par_read:
            for par in charmm_struct.par:
                charmm_ssp.read_par(par)

        #coordinates of C=O

        # a- coords of C
        a = np.array([-3.125, -31.473, 23.398])
        # b  - coords of O
        b = np.array([-2.156, -31.832, 23.856])

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

        c_o = np.subtract(b,a)
        c_o = b - a
        c_o_mag = np.linalg.norm(c_o)
        unitvector  = c_o/c_o_mag
        # point = a + unitvector*2.7
        point = a + unitvector*0.514

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
        charmm_ssp.write_pdb( modelling_folder+'point.pdb')

        if state_hema == 'oxi_hema':
            E_oxi = field_magnitude
            # print E_oxi*100, 'in MV/cm, oxi'
            E_oxi_field = field_vector
            # print E_oxi_field, 'oxi field'
        else:
            E_red = field_magnitude
            E_red_field = field_vector
            # print E_red*100, 'in MV/cm, red'

            # print E_red_field, 'red field'


    angle_unit_Eoxi = np.degrees(np.arccos(dot(c_o,E_oxi_field)/norm(c_o)/norm(E_oxi_field)))
    print angle_unit_Eoxi, 'oxi angle with CO vector'
    angle_unit_Ered = np.degrees(np.arccos(dot(c_o,E_red_field)/norm(c_o)/norm(E_red_field)))
    print angle_unit_Ered, 'red angle with CO vector'

    print np.degrees(np.arccos(round(dot(E_red_field,E_oxi_field)/norm(E_red_field)/norm(E_oxi_field)))), 'oxi_red_angle'

    scalar_projection_oxi = np.linalg.norm(E_oxi_field)*(dot(unitvector,E_oxi_field)/norm(unitvector)/norm(E_oxi_field))
    # print scalar_projection_oxi, 'scalar_projection_oxi'
    vector_projection_oxi = unitvector * scalar_projection_oxi
    # print vector_projection_oxi, 'vec_proj_oxi'
    # print np.linalg.norm(vector_projection_oxi), 'magnitude oxi projection'
    print np.linalg.norm(vector_projection_oxi)*100, 'magnitude oxi projection in MV/cm'
    print np.degrees(np.arccos(round(dot(c_o,vector_projection_oxi)/norm(c_o)/norm(vector_projection_oxi)))), 'vector proje angle oxi'



    scalar_projection_red = np.linalg.norm(E_red_field)*(dot(unitvector,E_red_field)/norm(unitvector)/norm(E_red_field))
    # print scalar_projection_red, 'scalar_projection_oxi'
    vector_projection_red= unitvector * scalar_projection_red
    # print vector_projection_red, 'vec_proj_red'
    # print np.linalg.norm(vector_projection_red), 'magnitude red projection'
    print np.linalg.norm(vector_projection_red)*100, 'magnitude red projection in MV/cm'
    print np.degrees(np.arccos(round(dot(c_o,vector_projection_red)/norm(c_o)/norm(vector_projection_red)))), 'vector proje angle red'


    print np.degrees(np.arccos(round(dot(vector_projection_red,vector_projection_oxi)/norm(vector_projection_red)/norm(vector_projection_oxi)))), 'oxi_red_angle_projection'

    print np.linalg.norm(vector_projection_oxi) - np.linalg.norm(vector_projection_red), 'magnitude difference in V/Angstrom'

    de_new = np.subtract(vector_projection_oxi,vector_projection_red)
    print np.degrees(np.arccos(round(dot(vector_projection_red,de_new)/norm(vector_projection_red)/norm(de_new)))), 'red_difference_angle'
    print np.degrees(np.arccos(round(dot(vector_projection_oxi,de_new)/norm(vector_projection_oxi)/norm(de_new)))), 'oxi_difference_angle'

    # de = E_red_field - E_oxi_field
    de = E_oxi_field - E_red_field
    print np.degrees(np.arccos(round(dot(de,c_o)/norm(de)/norm(c_o)))), 'diff_CO_angle'

    # project = np.linalg.norm(de)*(dot(unitvector,de)/norm(unitvector)/norm(de))
    # vect = unitvector * project
    # print norm(vect), 'other way, V/Angstrom'

    # de_new = vector_projection_red - vector_projection_oxi
    print norm(de_new), 'final: V/Angstrom'
    print norm(de_new)*100, 'final: MV/cm'







