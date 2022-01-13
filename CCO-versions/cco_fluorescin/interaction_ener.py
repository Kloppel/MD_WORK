# coding=utf-8

import os
import shutil
import kbp2
from workspace_jd import cco_config
import numpy as np
import stat


def extract_frame_with_water_charmm(folder, psf_file, crd_file, dcd_file, start, step, end, frame_name = '', submitt = False):

    if frame_name == '':
        frame_name = '@{current_frame}'
    else:
        frame_name = frame_name + '.pdb'

    charmm_script = folder + 'charmm_frame_water.inp'
    f = open(charmm_script, 'w')

    charmm_input_text = """* read rtf,para,psf
* read coor file ifile n name my_traj.trj
* (n is the n:th frame in the trajectory)
* write coor pdb name frame_n.pdb
* one frame from my_traj


OPEN READ UNIT 42 CARD NAME "/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_alw.inp"
READ rtf CARD UNIT 42
close unit 42

OPEN READ UNIT 42 CARD NAME "/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_lipid.rtf"
READ rtf CARD UNIT 42 APPEND
close unit 42

OPEN READ UNIT 42 CARD NAME "/scratch/scratch/jdragelj/CHARMM_TOPPAR/top_all36_cgenff.rtf"
READ rtf CARD UNIT 42 APPEND
close unit 42

OPEN READ UNIT 42 CARD NAME "/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all22_prot_plus_heme_and_Cu.inp"
READ PARAMETERS UNIT 42 CARD
close unit 42

OPEN READ UNIT 42 CARD NAME "/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_lipid.prm"
READ PARAMETERS UNIT 42 CARD APPEND
close unit 42

OPEN READ UNIT 42 CARD NAME "/scratch/scratch/jdragelj/CHARMM_TOPPAR/par_all36_cgenff.prm"
READ PARAMETERS UNIT 42 CARD APPEND
close unit 42

stream "/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_fluorescin.str"


!Load the structure that already has waterbox
!Read in Sequence and coordinates of the protein in waterbox

OPEN read UNIT 13 CARD NAME "%s" """ %psf_file+ """
Read psf unit 13 card
clos unit 13

OPEN READ UNIT 14 CARD NAME "%s" """ %crd_file+ """
Read coor unit 14 card
clos unit 14


set start_frame %i""" %start + """
set step_frame %i""" %step + """
set end_frame %i""" %end + """

! Now looping over frames
set current_frame @{start_frame}
label loop

       open unit 15 read file name "%s" """ %dcd_file + """
       read coor file unit 15 ifile @{current_frame}
       write coor pdb name %s sele all end """ %frame_name + """
       close unit 15

      incr current_frame by @{step_frame}
if current_frame le @{end_frame} goto loop


stop
"""
    #
    #            write coor pdb name frame@{current_frame}.pdb


    f.write(charmm_input_text)
    f.close()


    charmms_file = folder + 'charmm_run.sh'
    file = open(charmms_file, 'a')
    if not os.path.exists(charmms_file):
        script = """#!/bin/tcsh
"""
    else:
        script = """

"""
    script += """
cd  %s """ %folder + """
#charmm36b1_64 < charmm_frame_water.inp > charmm_frame_water.out
charmm_c39a2q < charmm_frame_water.inp > charmm_frame_water.out

"""
    script += '\n'
    file.write(script)
    file.close()
    st = os.stat(charmms_file)
    os.chmod(charmms_file, st.st_mode | stat.S_IEXEC)

    if submitt:
        kbp2.charmm.submit_charmm(folder, 'charmm_frame_water.inp', 'charmm_frame_water.out')

def remove_water(protein, folder, output_name, water_segments):
    pdb_mod = kbp2.file_parser.Simple_struct_parser()
    pdb_mod.read_pdb(protein)
    new = pdb_mod.copy(segname=water_segments,exclude=True)
    new.write_pdb(folder + output_name + '.pdb')

def split_pqr(pqr_file, residues_to_exclude, exc_pqr_file):

    f = open(pqr_file, 'r')
    q = open(exc_pqr_file, 'a')

    for line in f:
        line_comps = line.split()
        if line_comps[0] == 'ATOM':
            if line_comps[3] not in residues_to_exclude:
                q.write(line)

    f.close()
    q.close()

def add_potential_calculation(folder, subfolder):

    potentials_file = folder + 'potentials.sh'
    file = open(potentials_file, 'a')
    if not os.path.exists(potentials_file):
        script = """#!/bin/tcsh

"""
    else:
        script = """
"""
    script += subfolder + 'pot.sh'
    script += '\n'
    file.write(script)
    file.close()
    st = os.stat(potentials_file)
    os.chmod(potentials_file, st.st_mode | stat.S_IEXEC)

def add_apbs_calculation(main_folder, run_folder, apbs_bin, coulomb_bin, jobname):

    apbs_cals_file = main_folder + 'abps_calcs.sh'
    file = open(apbs_cals_file, 'a')
    if not os.path.exists(apbs_cals_file):
        tcsh_script = """#!/bin/tcsh

"""
    else:
        tcsh_script = """

"""
    tcsh_script += """
cd %s """ % run_folder + """

    """
    tcsh_script += """
%s apbs.in > apbs.out """ %  apbs_bin               + """
%s -e %s > coulomb.out   """ % (coulomb_bin, jobname + '.pqr')
    tcsh_script += '\n'
    file.write(tcsh_script)
    file.close()
    st = os.stat(apbs_cals_file)
    os.chmod(apbs_cals_file, st.st_mode | stat.S_IEXEC)

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

def prepare_structure(protein, modelling_folder, trajectory, minimize = False, pqr=False,
                      x_y_z = None, inter=False, sasa=1.4, delete_struct='', apbs_calc=False,
                      gen_born=None, submitt=False):

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

        if 'flu2' in trajectory:
            top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_flu2.rtf")
            par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_flu2.prm")
            patches.append({'FCYS': ['CYS-299_ACHA', 'FLU-1_FLUR']})
        if 'fluh' in trajectory:
            top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_fluh.rtf")
            par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_fluh.prm")
            patches.append({'FCYS': ['CYS-299_ACHA', 'FLU-1_FLUR']})
        if 'flu-' in trajectory:
            top.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/top_flux.rtf")
            par.append("/scratch/scratch/jdragelj/CHARMM_TOPPAR/toppar_charmm_FLU/split/par_flux.prm")
            patches.append({'FXCY': ['CYS-299_ACHA', 'FLU-1_FLUR']})

        charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
        charmm_struct.add_structure(protein)
        charmm_struct.add_decision('rename__CA_CAL', 'keep')
        charmm_struct.add_decision('rename__HOH_TIP3', 'keep')
        charmm_struct.add_decision('disu_bridges', 'closed')
        charmm_struct.add_decision('cap_termini', 'dont_cap')

        if 'red' in trajectory:
            cco_state = 'red'
            cco_settings = cco_config.COX_settings(state=cco_state, specie='A_paraccocus', cycle_state='R')
        else:
            cco_state = 'oxi'
            cco_settings = cco_config.COX_settings(state=cco_state, specie='A_paraccocus', cycle_state='O')

        charmm_struct.set_prot_residue(('GLU', 278, 'ACHA'), patch='GLUP')
        charmm_struct.set_prot_residue(('LYS', 354, 'ACHA'), charge=0)
        charmm_struct.set_prot_residue(('ASP', 399, 'ACHA'), patch='ASPP')
        charmm_struct.set_prot_residue(('GLU', 481, 'ACHA'), patch='GLUP')

        #CHARGE PACTHES
        for patch in cco_settings.charge_patches:
            new_patch = cco_config.copy_patch(patch)
            patches.append(new_patch)

        for patch in patches:
            patch_name = patch.keys()[0]
            residues = patch[patch_name]
            charmm_struct.add_patch(patch_name, residues)

        #BOND PACTHES -> HEME AND COPPER
        patches_no_autogen = []
        for patch in cco_settings.bond_patches:
            new_patch = cco_config.copy_patch(patch)
            patches_no_autogen.append(new_patch)

        for patch in patches_no_autogen:
            patch_name = patch.keys()[0]
            residues = patch[patch_name]
            charmm_struct.add_patch(patch_name, residues, no_autogen=True)

        charmm_struct.workdir = modelling_folder

        if minimize:
            charmm_struct.add_charmm_command('minimize abnr nsteps 200 tolg 0.01', adj_task='hbuild')
            charmm_struct.add_charmm_command('minimize sd nsteps 200 tolg 0.1', adj_task='hbuild')
            if x_y_z is not None:
                pbc_block = get_pbc_set(x_y_z)
                charmm_struct.add_charmm_command(pbc_block, adj_task='hbuild')
            charmm_struct.check_structures(quiet=True)
            charmm_struct.charmm_instructions['do_minimize'] = False
            charmm_struct.run_charmm(submit=True)
            return

        if inter:
            flu_inter_block = """
ENERgy E14Fac 0.0
INTEraction select segid FLUR end select all .and. .not. segid FLUR end
"""
            charmm_struct.add_charmm_command(flu_inter_block, adj_task='hbuild')

        facts = False
        gbmv = False
        gbsw = False
        if gen_born is not None:


            if  gen_born == 'gbmv':
                gbmv = True
            # if gen_born == 'facts':
            #     facts = True
            # elif  gen_born == 'gbsw':
            #     gbsw = True

        # if facts:
        #     ### Definition of GB ###
        #     facts_commands = []
        #     facts_commands.append("nbond nbxmod 5 atom cdiel eps 1.0 shift vatom vdistance vshift -")
        #     facts_commands.append("cutnb 14.0 ctofnb 12.0 ctonnb 6.5 e14fac 0.4 wmin 1.5")
        #     facts_commands.append("scalar wmain = radius")
        #     facts_commands.append("scalar wmain set 1.0 select hydrogen end")
        #     facts_commands.append("facts tcps 19 teps 2.0 gamm 0.05 tavw true")
        #     ### adding GB ###
        #     command_block = " "
        #     for command in facts_commands:
        #         command_block += '\n'
        #         command_block += command
        #     charmm_struct.add_charmm_command(command_block, adj_task='hbuild')

        # if gbsw:
        #     ### Definition of GB ###
        #     gbsw_commands = []
        #     gbsw_commands.append("stream \"/scratch/scratch/jdragelj/radius_gbsw.str\"")
        #     # gbsw_commands.append("stream \"/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/gbsw/radius_gbsw.str\"")
        #     gbsw_commands.append("scalar wmain statistics select .not. type H* end")
        #     gbsw_commands.append("define check select (.not type H* ) .and. ( prop wmain .eq. 0.0 ) show end")
        #     gbsw_commands.append("if ?nsel ne 0  stop       !some heavy atom have a zero radius")
        #     gbsw_commands.append("GBSW conc 0.1 sgamma 0.03 GBenergy epsp 4 epsw 80 MOLSURF")

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


    	if sasa:
            sasa_block = """
coor surf select all end acce rpro %.1f """ % sasa + """
!define backbone select type n .or. type ca .or. type c .or. type o end
!SCALAR WMAIN STAT sele ((segid BCHA .and. resid 73) .and. .not. backbone) end
define imidazole select type HD* .or. type HE* .or. type NE2 .or. type ND1 .or. type CG .or. type CD2 .or. type CE1 end
SCALAR WMAIN STAT sele ((segid BCHA .and. resid 73) .and. imidazole) end
echo ?ELEC ?stot

"""
#             sasa_block = """
# coor surf select all end acce rpro %.1f """ % sasa + """
# define backbone select type n .or. type ca .or. type c .or. type o end
# SCALAR WMAIN STAT sele ((segid FLUR .and. (type O5 .or. type O6 .or. type H14)) .and. .not. backbone) end
# echo ?ELEC ?stot
#
# """
            charmm_struct.add_charmm_command(sasa_block, adj_task='hbuild')

        if delete_struct == 'flu':
            command_block ="""
bomblev -4

DELETE ATOM SELE segid FLUR END

"""
            charmm_struct.add_charmm_command(command_block, adj_task='hbuild')

        elif delete_struct == 'protein':
            command_block ="""
bomblev -4

DELETE ATOM SELE all .and. .not. segid FLUR END

"""
            charmm_struct.add_charmm_command(command_block, adj_task='hbuild')


        charmm_struct.check_structures(quiet=True)
        charmm_struct.charmm_instructions['do_minimize'] = False
        charmm_struct.run_charmm(submit=submitt)

        if pqr:
            charmm_ssp = charmm_struct.structure
            premodelled_structure = charmm_struct.get_modelled_structure()
            charmm_ssp = premodelled_structure

            if not charmm_ssp.par_read:
                for par in charmm_struct.par:
                    charmm_ssp.read_par(par)

            charmm_ssp.write_pqr(modelling_folder + 'cco_flu.pqr')
            flu0 = False
            if flu0:
                for new_atom in charmm_ssp.atoms:
                    resid = new_atom['resid']
                    segname = new_atom['segname']
                    name = new_atom['name']
                    if segname != 'FLUR':
                        charmm_ssp.struct[segname][resid][name]['charge'] = 0.0
                charmm_ssp.write_pqr(modelling_folder + 'cco_0_flu.pqr')

            if apbs_calc:
                jobname = 'cco_frame'
                structure = charmm_ssp
                run_folder = modelling_folder + 'apbs/'
                if not os.path.exists(run_folder):
                    os.mkdir(run_folder)
                apbs_bin    = "/scratch/scratch/tmeyer/CHARMM_NAMD/apbs"
                coulomb_bin = "/scratch/scratch/tmeyer/CHARMM_NAMD/coulomb"
                kbp2.apbs.create_apbs_input(jobname, run_folder, structure, apbs_bin, coulomb_bin, target_res=0.3, verbose=True, ion_conc=0.1,
                        tmp_run_folder=None)
                kbp2.apbs.run_local(run_folder)

            return charmm_ssp
        else:
            return

def get_binding_energy(folder, ordered_components, KJ=False):

    binding_energy = 0
    if folder[-1] != '/':
        folder+='/'
    for i, component in enumerate(ordered_components):
        print folder
        charmm_energies = kbp2.charmm.get_charmm_energy(folder + component + '/', KJ=False)
        if 'IMELec' in charmm_energies.keys():
            # comp_energy = charmm_energies['ELEC'][0] + charmm_energies['GBEnr'][0]
            comp_energy = charmm_energies['ELEC'][0] + charmm_energies['GBEnr'][0] + charmm_energies['IMELec'][0]
        else:
            comp_energy = charmm_energies['ELEC'][0] + charmm_energies['GBEnr'][0]
        if i == 0:
            binding_energy += comp_energy
        else:
            binding_energy -= comp_energy

    return binding_energy


def run_apbs_int_energy(jobname, run_folder, charmm_ssp, apbs_bin, coulomb_bin, qsub_parameter):

    if not os.path.exists(run_folder):
        os.mkdir(run_folder)
    if qsub_parameter == 'manual':
        kbp2.apbs.create_interaction_energy_apbs_input(jobname, run_folder, charmm_ssp, apbs_bin, coulomb_bin,
                                         target_res=0.3, verbose=True, ion_conc=0.1, only_prepare_calcs = True)
    else:
        kbp2.apbs.create_interaction_energy_apbs_input(jobname, run_folder, charmm_ssp, apbs_bin, coulomb_bin,
                                         target_res=0.3, verbose=True, ion_conc=0.1, only_prepare_calcs = True)
        kbp2.apbs.submitt_apbs_job(run_folder=run_folder, qsub_parameter=qsub_parameter)

def check_output(out_file):
    status = True
    q = open(out_file, 'r')
    for line in q:
        if 'FATAL ERROR' in line:
            status = False
        if 'ABNORMAL TERMINATION' in line:
            status = False
    return status

def produce_interaction_energy(protein, main_folder, frame_workfolder, residues_to_exclude, trajectory, minimize = False, qsub = ''):

    multivalue_source = '/scratch/scratch/tmeyer/CHARMM_NAMD/apbs-1.4_source/tools/bin/multivalue'
    apbs_bin    = "/scratch/scratch/tmeyer/CHARMM_NAMD/apbs"
    coulomb_bin = "/scratch/scratch/tmeyer/CHARMM_NAMD/coulomb"

    jobname = 'cco_0_flu'
    modelling_folder = frame_workfolder + 'modelling/'
    pqr_file = modelling_folder + 'cco_flu.pqr'
    exc_pqr_file = frame_workfolder + 'cco.pqr'
    csv_file = frame_workfolder + 'cco.csv'
    phi_file = frame_workfolder + jobname + '_pot.phi'
    run_folder = frame_workfolder + 'apbs/'
    dx_file = run_folder + jobname + '.dx'

    if not os.path.exists(csv_file):
        charmm_ssp = prepare_structure(protein, modelling_folder, trajectory, minimize=minimize)
        split_pqr(pqr_file, residues_to_exclude, exc_pqr_file)
        kbp2.apbs.pqr2csv(exc_pqr_file)
        run_apbs_int_energy(jobname, run_folder, charmm_ssp, apbs_bin, coulomb_bin, 'manual')
        add_apbs_calculation(main_folder, run_folder, apbs_bin, coulomb_bin, jobname)
        return

    if not os.path.exists(dx_file):
        status = check_output(run_folder+'apbs.out')
        if not status:
            print "This calculation crashed, removing everything!"
            shutil.rmtree(modelling_folder)
            shutil.rmtree(run_folder)
            shutil.rmtree(csv_file)
            shutil.rmtree(exc_pqr_file)
            return
        else:
            print 'Waiting for apbs calculation to be finished!'
            return
    else:
        if not os.path.exists(phi_file):
            kbp2.apbs.apply_potential(frame_workfolder, jobname, multivalue_source, csv_file, dx_file)
            add_potential_calculation(main_folder, frame_workfolder)
            print 'Combine dx to get potential'
            return
        else:
            if not os.path.exists(frame_workfolder + 'interaction_energy.dat'):
                energy = kbp2.apbs.calculate_interaction_energy(frame_workfolder, exc_pqr_file, phi_file)
                return
            else:
                print 'This calculation has been done, moving on!'
                return

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
