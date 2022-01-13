# coding=utf-8

import os
import kbp2
import numpy as np



def get_pbc_set(x_y_z):

    x = float(x_y_z[0])
    y = float(x_y_z[1])
    z = float(x_y_z[2])

    pbc_block = """
calc XSIZ = %.4f """ %x + """
calc YSIZ = %.4f """ %y + """
calc ZSIZ = %.4f """ %z + """

set 1 XSIZ
set 2 YSIZ
set 3 ZSIZ
define SOLUTE sele all .and. .not. (segid SWAT .or. segid WAT .or. segid PAH2O .or. segid QBH2O) end
!open read unit 18 card name "/scratch/scratch/tmeyer/CHARMM_NAMD/charmm/image-trans.img"
!read image unit 18 card
CRYST DEFINE ORTHO @XSIZ @YSIZ @ZSIZ 90. 90. 90.
CRYST BUILD CUTOFF 15
IMAGe BYREs XCEN 0.0 YCEN 0.0 ZCEN 0.0 SELEct .NOT. SOLUTE END
IMAGe BYSEg XCEN 0.0 YCEN 0.0 ZCEN 0.0 SELEct SOLUTE END
close unit 18

"""
    return pbc_block

def prepare_structure(protein, modelling_folder, state, pra_asp_prot = [False, False], x_y_z = None, gbmv=True, submitt=False):

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
        par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm")

        charge_patches = []
        bond_patches = []

        charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
        charmm_struct.add_structure(protein)
        charmm_struct.add_decision('rename__CA_CAL', 'keep')
        charmm_struct.add_decision('rename__HOH_TIP3', 'keep')
        charmm_struct.add_decision('disu_bridges', 'closed')
        charmm_struct.add_decision('cap_termini', 'dont_cap')

        if state == 'PF':
            charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
            charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                            'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
            charge_patches.append({'A343': ['HEM-2_EHEM', 'HSD-419_ACHA', 'HAO2-1_FEOH']})
            charge_patches.append(
                {'CBP2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

            bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
            bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
            bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
            bond_patches.append({'EISO': ['HAO2-1_FEOH', 'HEM-2_EHEM']})
            bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'OMI2-1_HOHC']})
            bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                          'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

        if state == 'F':
            charge_patches.append({'AHE3': ['HEM-2_GHEM', 'HSD-102_ACHA', 'HSD-421_ACHA']})
            charge_patches.append({'CA21': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                            'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})
            charge_patches.append({'A3H4': ['HEM-2_EHEM', 'HSD-419_ACHA', 'HAO2-1_FEOH']})
            charge_patches.append(
                {'CBF2': ['CU1-1_META', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC', 'TYR-288_ACHA', 'HSE-284_ACHA']})

            bond_patches.append({'PHEM': ['HSD-419_ACHA', 'HEM-2_EHEM']})
            bond_patches.append({'PHEM': ['HSD-102_ACHA', 'HEM-2_GHEM']})
            bond_patches.append({'PHE2': ['HSD-421_ACHA', 'HEM-2_GHEM']})
            bond_patches.append({'EISO': ['HAO2-1_FEOH', 'HEM-2_EHEM']})
            bond_patches.append({'CUB2': ['CU1-1_META', 'HSE-284_ACHA', 'HSD-333_ACHA', 'HSD-334_ACHA', 'HOH-1_HOHC']})
            bond_patches.append({'CUAP': ['CU1-2_META', 'CU1-3_META', 'HSE-217_BCHA', 'CYS-252_BCHA', 'GLU-254_BCHA',
                                          'CYS-256_BCHA', 'HSE-260_BCHA', 'MET-263_BCHA']})

        charge_patches.append({'GLUP':['GLU-286_ACHA']})
        charge_patches.append({'LSN':['LYS-362_ACHA']})

        if pra_asp_prot[1]:
            charge_patches.append({'ASPP': ['ASP-407_ACHA']})
        if pra_asp_prot[0]:
            charge_patches.append({'PRAH': ['PRA-1_EHEM']})

        for patch in charge_patches:
            patch_name = patch.keys()[0]
            residues = patch[patch_name]
            charmm_struct.add_patch(patch_name, residues)

        for patch in bond_patches:
            patch_name = patch.keys()[0]
            residues = patch[patch_name]
            charmm_struct.add_patch(patch_name, residues, no_autogen=True)

        charmm_struct.workdir = modelling_folder


        if gbmv:
            gbsw_commands = []
            gbsw_command = """
update ctonnb 10 ctofnb 12 cutnb 13.5 eps 4 vswi switch

GBMV BETA -20 EPSILON 80 DN 1.0 watr 1.4 GEOM -
TOL 1e-8 BUFR 0.5 Mem 10 CUTA 20 HSX1 -0.125 HSX2 0.25 -
ALFRQ 1 EMP 1.5 P4 0.0 P6 8.0 P3 0.70 ONX 1.9 OFFX 2.1 -
WTYP 2 NPHI 38 SHIFT -0.102 SLOPE 0.9085 CORR 1

ENERgy EPS 4 E14Fac 0.0
"""
            gbsw_commands.append(gbsw_command)

            ### adding GB ###
            command_block = " "
            for command in gbsw_commands:
                command_block += '\n'
                command_block += command

            charmm_struct.add_charmm_command(command_block, adj_task='hbuild')

        if x_y_z is not None:
            pbc_block = get_pbc_set(x_y_z)
            charmm_struct.add_charmm_command(pbc_block, adj_task='hbuild')


        charmm_struct.check_structures(quiet=True)
        charmm_struct.charmm_instructions['do_minimize'] = False
        charmm_struct.run_charmm(submit=submitt)


def get_binding_energy(folder, ordered_components, KJ=False):

    binding_energy = 0
    if folder[-1] != '/':
        folder+='/'
    for i, component in enumerate(ordered_components):
        print folder
        charmm_energies = kbp2.charmm.get_charmm_energy(folder + component + '/', KJ=False)
        if 'IMELec' in charmm_energies.keys():
            comp_energy = charmm_energies['ELEC'][0] + charmm_energies['GBEnr'][0] + charmm_energies['IMELec'][0]
        else:
            comp_energy = charmm_energies['ELEC'][0] + charmm_energies['GBEnr'][0]
        if i == 0:
            binding_energy += comp_energy
        else:
            binding_energy -= comp_energy

    return binding_energy


def check_output(out_file):
    status = True
    q = open(out_file, 'r')
    for line in q:
        if 'FATAL ERROR' in line:
            status = False
        if 'ABNORMAL TERMINATION' in line:
            status = False
    return status


def compare_trajectories(folder, positions, subfolder_list, frame):

    csv_file = folder + 'compared_energies.csv'
    f = open(csv_file, 'w')
    first_line = ' ,'
    for position in positions:
        first_line += '%i' % position
        first_line += ','
    f.write(first_line + '\n')

    for subfolder in subfolder_list:
        csv_line = subfolder
        csv_line += ','
        for position in positions:

            energy_file = folder + str(
                position) + '/' + subfolder + '/' + '%i/' % frame + 'interaction_energy.dat'
            ef = open(energy_file, 'r')
            for line in ef:
                if 'KJ/mol' in line:
                    line = line.split()
                    csv_line += line[0]
                    csv_line += ','
        f.write(csv_line + '\n')
    f.close()

def energy_averages(folder, positions, subfolder_list, frames):
    csv_file = folder + 'compared_energies.csv'
    f = open(csv_file, 'w')
    first_line = ' ,'
    for position in positions:
        first_line += '%i' % position
        first_line += ','
    f.write(first_line + '\n')

    for subfolder in subfolder_list:
        csv_line = subfolder
        csv_line += ','
        for position in positions:
            avg = []
            for i, frame in enumerate(frames):
                energy_folder = folder + str(position) + '/' + subfolder + '/' + '%i/' % frame
                energies = kbp2.charmm.get_charmm_energy(energy_folder, inter=True)
                result_ene = energies['ELEC'][0]
                avg.append(result_ene)
            avg = np.mean(avg)
            csv_line += str(avg)
            csv_line += ','
        f.write(csv_line + '\n')
        f.close()

# def extraction_and_minimization(mds_sourcefolder, main_folder, trajectory, tmp_charmm_extraction_folder, minimization_frame_folder,
#                                 psf_file, crd_file, dcd_file, start, step, end, frame, frame_name, pbc_file, step_md):
#
#     minimized_frame = minimization_frame_folder + frame_name + '_out.pdb'
#     if os.path.exists(minimized_frame):
#         shutil.copy2(minimized_frame, main_folder + 'frame%i_water.pdb' % frame)
#         water_segments = ['WAT','SWAT','QBH2','PAH2']
#         remove_water(main_folder + 'frame%i_water.pdb' % frame, main_folder, 'frame%i' % frame, water_segments)
#         shutil.rmtree(tmp_charmm_extraction_folder)
#         shutil.rmtree(minimization_frame_folder)
#         return
#     else:
#         frame_with_water = tmp_charmm_extraction_folder + frame_name + '.pdb'
#         if os.path.exists(frame_with_water):
#             if os.path.exists(minimization_frame_folder + frame_name + '_charmm.out'):
#                 status = check_output(minimization_frame_folder + frame_name + '_charmm.out')
#                 if status == False:
#                     print('Minimizaion crashed!', frame, trajectory)
#                     return
#                 else:
#                     return
#             else:
#                 if not os.path.exists(minimization_frame_folder):
#                     os.mkdir(minimization_frame_folder)
#                 x_y_z = []
#                 pbc_f = open(pbc_file, 'r')
#                 for line in pbc_f:
#                     line = line.split()
#                     if line[0] == step_md:
#                         x_y_z.append(line[1])
#                         x_y_z.append(line[5])
#                         x_y_z.append(line[9])
#                 prepare_structure(frame_with_water, minimization_frame_folder, trajectory, minimize=True, pqr=False,
#                                   qsub_parameter = '-q D47.q,D48.q,D49.q,D50.q,D51.q,D52.q,D61.q,D62.q,D63.q,D64.q', x_y_z = x_y_z)
#                 return
#         else:
#             if not os.path.exists(tmp_charmm_extraction_folder):
#                 os.mkdir(tmp_charmm_extraction_folder)
#             if os.path.exists(tmp_charmm_extraction_folder + 'charmm_frame_water.out'):
#                 status = check_output(tmp_charmm_extraction_folder + 'charmm_frame_water.out')
#                 if status == False:
#                     print('extraction crashed!', frame, trajectory)
#                     return
#                 else:
#                     return
#             shutil.copy2(psf_file, tmp_charmm_extraction_folder+'cco_and_water.psf')
#             shutil.copy2(crd_file, tmp_charmm_extraction_folder+'cco_and_water.crd')
#             shutil.copy2(dcd_file, tmp_charmm_extraction_folder+'3hb3_md_flu.dcd')
#             extract_frame_with_water_charmm(tmp_charmm_extraction_folder, 'cco_and_water.psf', 'cco_and_water.crd',
#                                             '3hb3_md_flu.dcd', start, step, end, frame_name, submitt= True)



if __name__ == '__main__':

    data_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/md_F_prah_asp_dep/'
    frames_folder = data_folder + '/frames_voda/'
    frames = np.arange(284,360,4)
    for frame in frames:
        state = 'F'
        frame_pdb = frames_folder + 'frame%i.pdb' % frame
        pbc_file = data_folder + 'md/output/2gsm_%s.xst' % state
        x_y_z = []
        pbc_f = open(pbc_file, 'r')
        for line in pbc_f:
            line = line.split()
            step_md = frame*25000
            if line[0] == str(step_md):
                x_y_z.append(line[1])
                x_y_z.append(line[5])
                x_y_z.append(line[9])

        protein = frame_pdb
        if not os.path.exists(frames_folder + 'mgcw/'):
            os.mkdir(frames_folder + 'mgcw/')
        pdb = kbp2.file_parser.Simple_struct_parser()
        pdb.read_pdb(frames_folder + 'frame%i.pdb' % frame)
        pdb.create_struct()
        for seg in pdb.struct.iter_segments():
            for res in seg.iter_residues():
                change = False
                if res.resname == 'TIP3':
                    if res.segname == 'QBH2':
                        if res.resid in [42, 48, 55]:
                            change = True
                    if change:
                        for atm in res.iter_atoms():
                            atm['segname'] = 'VODA'
        pdb.create_struct()
        pdb_w = pdb.copy(segname=['SWAT', 'WAT', 'QBH2', 'PAH2', 'XTC', 'XT2'], exclude=True)
        pdb_w.write_pdb(frames_folder + 'mgcw/frame%i.pdb' % frame)
        frame_pdb_final = frames_folder + 'mgcw/frame%i.pdb' % frame

        modelling_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/total_ene_asp_pra/md_F_prah_asp_dep/frame%i/' % frame
        if not os.path.exists(modelling_folder):
            os.mkdir(modelling_folder)
        prepare_structure(frame_pdb_final, modelling_folder, state, pra_asp_prot=[True, False], x_y_z=x_y_z, gbmv=True, submitt=False)


