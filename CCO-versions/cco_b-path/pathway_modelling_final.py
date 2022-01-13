# coding=utf-8

# import sys
# sys.path.append('/user/tmeyer/workspace/script/protein_toolbox/packages')

import shutil
import os
import kbp2

############# todo: class - protein specific modeling with comments that are inserted ###########


def pathway_modelling(protein, top, par, pdb_file, folder, pathway_def, decisions=None, patches=None, spec_modelling=None, \
                      consecutive = True, sidechain_relax = False, sidechain_envir_relax = False):

    if protein == 'cco':
        print 'modelling cytochrome c oxidase pathway'

    jobname = protein

    additional_titr_residues = {}
    template_state = {'charge' : 0,
                      'patch' : None,
                      'external_patches' : None,
                      'rename' : None,
                      'special' : None}
    additional_titr_residues = {}
    resname = 'SER'
    additional_titr_residues[resname] = [dict(template_state) for i in range(2)]
    additional_titr_residues[resname][0]['charge'] = 0
    additional_titr_residues[resname][1]['charge'] = -1
    additional_titr_residues[resname][1]['patch'] = 'SEPD'
    resname = 'THR'
    additional_titr_residues[resname] = [dict(template_state) for i in range(2)]
    additional_titr_residues[resname][0]['charge'] = 0
    additional_titr_residues[resname][1]['charge'] = -1
    additional_titr_residues[resname][1]['patch'] = 'THPD'
    resname = 'H2O'
    additional_titr_residues[resname] = [dict(template_state) for i in range(2)]
    additional_titr_residues[resname][0]['charge'] = 0
    additional_titr_residues[resname][1]['charge'] = -1
    additional_titr_residues[resname][1]['patch'] = 'HYH2'



    pathway_list = []
    for residue in pathway_def:
        if type(residue) is not list:
            pathway_target_residue = list(residue)
        else:
            pathway_target_residue = residue
        #creating input list for pathways
        pathway_list.append(pathway_target_residue)


    if sidechain_relax or sidechain_envir_relax:
        resid_selection = []
        for target_residue in pathway_list:
            resid = target_residue[2]
            if resid is not None:
                resid_selection.append(str(resid))


    ##### SELECTION STRING FOR CUSTOM MINIMIZATION ####
    selection = 'cons fix select .not. ( none .or. hydrogens )  end'


    # implement script making function
    if  sidechain_relax:
        selection = '''
define selstr1 select (resid 614 .or. resid 630) .and. segid TAH3 end
define selstr2 select (resid 309 .or. resid 244 .or. resid 312 .or. resid 248 .or. resid 315) .and. segid ACHA end
define selstr3 select resid 15 .and. segid BCHA end
define backbone select type n .or. type ca .or. type c .or. type o end
cons fix select .not. ((selstr1 .or. selstr2 .or. selstr3 .or. hydrogens) .and. .not. backbone) end
        ''' % (resid_selection)
        print 'Side chains are realxed'

    elif sidechain_envir_relax:
        selection = '''
define selstr1 select (resid 614 .or. resid 630) .and. segid TAH3 end
define selstr2 select (resid 309 .or. resid 244 .or. resid 312 .or. resid 248 .or. resid 315) .and. segid ACHA end
define selstr3 select resid 15 .and. segid BCHA end
define selstr4 select .byres. ((selstr1 .or. selstr2 .or. selstr3) .around. 5.0) end
define backbone select type n .or. type ca .or. type c .or. type o end
cons fix select .not. ((selstr4 .or. hydrogens) .and. .not. backbone) end
        ''' % (resid_selection)
        print 'Side chains and envir are relaxed'

    resnames = []
    out_files = []
    model_struct_list = []
    charmm_struct_list = []
    folder_list = []
    for i, target_residue in enumerate(pathway_list):
        target_resname, target_segname, target_resid, target_state = target_residue

        resnames.append(target_resname)

        if consecutive:
            if i > 0:
                previous_output = out_files[i-1]
                pdb = previous_output
            else:
                pdb = pdb_file
        else:
            pdb =  pdb_file

        if target_resname != 'basic':
            target = target_resname + target_segname + str(target_resid)
        else:
            target = target_resname

        project_folder = folder + jobname + target + '/'
        if not os.path.exists(project_folder):
            os.mkdir(project_folder)

        ### CHARMM ###
        charmm_struct = kbp2.charmm.Charmm_manager(top=top, par=par)
        charmm_struct.add_structure(pdb)

        titr_residues = charmm_struct.get_titr_residues()
        titr_residues.update(additional_titr_residues)
        charmm_struct.set_titr_residues(titr_residues)

        # todo: move to class with specific instruction for proteins
        #applying decisions - no decision for now
        if decisions is not None:
            for decision in decisions:
                decision_name = decision[0]
                action = decision[1]
                charmm_struct.add_decision(decision_name, action)
        #patches and spec_modelling - patches can be change a protonation state
        if patches is not None:
            for patch_name, residue_tuple_patch, mod_type in patches:
                if mod_type != 'prot':
                    charmm_struct.add_patch(patch_name, residue_tuple_patch)
                else:
                    charmm_struct.set_prot_residue(residue_tuple_patch, patch=patch_name)


        if spec_modelling is not None:
            for spec in spec_modelling:
                spec_command = spec[0]
                after_task = spec[1]
                charmm_struct.add_charmm_command(spec_command, adj_task=after_task)

#         #adding pathway commands
#         ##################### BASIC STRUCTURE ##################
#         if target_resname == 'basic':
#             general_command = '''!standard minimiz'''
#         ############################# ALL ######################
#         else:
#             target_tuple = (target_resname, target_resid, target_segname)
#             charmm_struct.set_prot_residue(target_tuple, state=target_state)
#             general_command = '''
# %s
# minimize sd nsteps 1000 tolg 0.1
# minimize abnr nsteps 5000 tolg 0.01
#     ''' % (selection)
#         ############ COMMAND ADDING AFTER CHOICE #############
#         charmm_struct.add_charmm_command(general_command, adj_task='hbuild')
#         # charmm_struct.add_charmm_command('ener eps 1.0', adj_task='hbuild')


        # #adding pathway commands
        ##################### BASIC STRUCTURE ##################
        if target_resname == 'basic':
            general_command = '''!standard minimiz'''
        ############################# ALL ######################
        else:
            target_tuple = (target_resname, target_resid, target_segname)
            charmm_struct.set_prot_residue(target_tuple, state=target_state)

        general_command = '''
%s
minimize sd nsteps 1000 tolg 0.1
minimize abnr nsteps 5000 tolg 0.01

cons fix select none end
cons fix select .not. selstr4 end

! stream "/user/jdragelj/python/CcO/MD/dynamics_template_cco_tenth_ns.str"
    ''' % (selection)
        ############ COMMAND ADDING AFTER CHOICE #############
        # print "MD runs now!"
        charmm_struct.add_charmm_command(general_command, adj_task='hbuild')


        ####### FOLDER PREPARATION ########
        if consecutive:
            if sidechain_relax or sidechain_envir_relax:
                task_modelling_folder = project_folder + 'cons_mod_relax/'
            else:
                task_modelling_folder = project_folder + 'cons_mod/'
        else:
            task_modelling_folder = project_folder + 'basic_mod'

        pdb_out = task_modelling_folder + target + '_out.pdb'
        out_files.append(pdb_out)
        charmm_struct.charmm_out_prefix = target + '_out'

        ###### FINAL STEP - MODELLING IN CHARMM #######
        charmm_struct.workdir = task_modelling_folder
        charmm_struct.check_structures()
        ### IMPORTANT - SHUTTING DOWN MINIMIZATION FOR THE PURPOSE OF THIS PROCEDURE.
        # Overwriting the decision from check_structure()

        # if target_resname != 'basic':
        charmm_struct.charmm_instructions['do_minimize'] = False

        if not os.path.exists(task_modelling_folder):
            os.mkdir(task_modelling_folder)
            charmm_struct.run_charmm()

        model_struct = charmm_struct.get_modelled_structure()
        model_struct_list.append(model_struct)

        charmm_struct_list.append(charmm_struct)
        folder_list.append(project_folder)

    return charmm_struct_list, model_struct_list, folder_list


def pathway_tapbs_apbs(ref_charmm_struct, model_struct_list, folder_list, titratable_yaml, pathway, top, par):

    #TAPBS
    tapbs_bin = '/scratch/scratch/tmeyer/kbplus2/tapbs_1.3_cav_enere'

    #APBS
    apbs_bin    = "/scratch/scratch/tmeyer/CHARMM_NAMD/apbs"
    coulomb_bin = "/scratch/scratch/tmeyer/CHARMM_NAMD/coulomb"
    qsub_parameter = '-q D47.q,D48.q,D49.q,D50.q,D51.q,D52.q,D53.q,D54.q,D55.q,D56.q,D57.q,D58.q,D59.q,D60.q,D64.q'
    energy_list = []
    pkint_list = []

    titratable_residues_yaml = kbp2.kbp_tools.parse_titratable_yaml(titratable_yaml)

    ref_charmm_struct.top = top
    ref_charmm_struct.par = par

    titr_residues = ref_charmm_struct.get_titr_residues()
    template_state = {'charge' : 0,
                      'patch' : None,
                      'external_patches' : None,
                      'rename' : None,
                      'special' : None}
    resname = 'DPP'
    titr_residues[resname] = [dict(template_state) for i in range(3)]
    titr_residues[resname][0]['charge'] = -1
    titr_residues[resname][1]['charge'] = 0
    titr_residues[resname][1]['patch'] = 'DPP2'
    titr_residues[resname][2]['charge'] = 0
    titr_residues[resname][2]['patch'] = 'DPP1'
    resname = 'EPP'
    titr_residues[resname] = [dict(template_state) for i in range(3)]
    titr_residues[resname][0]['charge'] = -1
    titr_residues[resname][1]['charge'] = 0
    titr_residues[resname][1]['patch'] = 'EPP2'
    titr_residues[resname][2]['charge'] = 0
    titr_residues[resname][2]['patch'] = 'EPP1'
    ref_charmm_struct.set_titr_residues(titr_residues)

    for resname, segname, resid, state in pathway:
        if resid is not None:
            resname_kbp = kbp2.kbp_tools.get_kbp_resname(resname)
            residue_tuple = (resname, resid, segname)
            ref_charmm_struct.rename_residue(residue_tuple, resname_kbp)

    for model_struct, folder, pathway in zip(model_struct_list, folder_list, pathway):

        print "Starting energy calculations for "
        if folder[-1] != '/':
            folder += '/'

        resname, segname, resid, state = pathway
        resname_kbp = kbp2.kbp_tools.get_kbp_resname(resname)
        print "Energy calculations for %s %s %s of state %s" % (resname_kbp, segname, resid, state)

        tapbs_folder = folder + 'tapbs/'
        tapbs_run_folder = tapbs_folder + 'run/'
        apbs_run_folder = folder + 'apbs/'


        if (not os.path.exists(tapbs_folder)) and (not os.path.exists(apbs_run_folder)):
            os.mkdir(tapbs_folder)

            ref_charmm_struct.transfer_coordinates(model_struct, allow_hydrogen_missmatch=True)
            charmm_folder = tapbs_folder + 'modelling/'
            os.mkdir(charmm_folder)
            ref_charmm_struct.workdir = charmm_folder
            ref_charmm_struct.charmm_instructions['do_minimize'] = False
            ref_charmm_struct.run_charmm()

            premodelled_structure = ref_charmm_struct.get_modelled_structure()
            # modelled_structure = premodelled_structure.copy(segname='MEMB',exclude=True)
            modelled_structure = premodelled_structure


            par = ref_charmm_struct.par
            for par_file in par:
                modelled_structure.read_par(par_file)

            jobname = ref_charmm_struct.title

            if resid is not None:
                residues_to_titrate = [(resname_kbp, resid, segname)]

                print "Starting TAPBS for %s %s %s of state %s" % (resname_kbp, segname, resid, state)
                kbp2.tapbs.run_tapbs(tapbs_run_folder, jobname, residues_to_titrate, modelled_structure, titratable_residues_yaml,
                                     tapbs_bin)

                # to run on dollys - jd
                # run_folder = '/scratch/scratch/jdragelj/tmp/apbs_test/run/testrun/'
                # temp_run_folder = '/public/scratch/jdragelj/kb/apbs_runs/'
                #to run on local - create whenever you want - jd

                # apbs_run_folder = apbs_folder + 'run/'
                # apbs_tmp_folder = tapbs_folder + 'tmp/'

            print "Starting APBS for %s %s %s of state %s" % (resname_kbp, segname, resid, state)
            kbp2.apbs.start_apbs_job(jobname, apbs_run_folder, modelled_structure, apbs_bin, coulomb_bin,
                                     target_res=0.4, verbose=True, qsub_parameter=qsub_parameter)

        # else:
        jobname = ref_charmm_struct.title

        print "Parsing TAPBS results for %s %s %s of state %s" % (resname_kbp, segname, resid, state)
        if resid is not None:
            # Read the results - TAPBS
            pkint_filename = tapbs_run_folder + jobname + '.pkint'
            g_filename = tapbs_run_folder + jobname + '.g'
            pkint, g, residue_list_ext = kbp2.kbp_tools.parse_g_pkint(pkint_filename, g_filename)
        else:
            pkint, g, residue_list_ext = (None, None, None)

        pkint_list.append(pkint)

        #read results - APBS
        print "Reading APBS results for %s %s %s of state %s" % (resname_kbp, segname, resid, state)
        results = kbp2.apbs.read_result(apbs_run_folder)
        if len(results) == 1:
            print("No output files found. Job may not be finished yet, come back later.")
        else:
            if None in results:
                print("Job not finished yet or either the solvation or coulomb calculation failed!")
            else:
                print("Everything is fine! Here are the results:\n")
                apbs_energy, coulomb_energy = results
                print("Solvation energy for %s %s %s of state %s : %5.3f" % (resname_kbp, segname, resid, state, apbs_energy))
                print("Coulomb energy for %s %s %s of state %s : %5.3f" % (resname_kbp, segname, resid, state, coulomb_energy))
                energy_list.append((apbs_energy,coulomb_energy))


    return g, pkint_list, energy_list, residue_list_ext



        ############################
        ####    MAIN SCRIPT   ######
        ############################

if __name__ =="__main__":

    ####TESTING FOR CCO####

    protein = 'cco'

    # charmm files
    top = []
    top.append("/user/jdragelj/charmm_toppar/cco/top_alw_jd.inp")
    top.append("/scratch/scratch/awoelke/md_cco/toppar/top_all36_lipid.rtf")
    par = []
    par.append("/scratch/scratch/awoelke/md_cco/toppar/par_all22_prot_plus_heme_and_Cu.inp")
    par.append("/scratch/scratch/awoelke/md_cco/toppar/par_all36_lipid.prm")



    # TYR shifted - modelled structure
    # pdb = '/user/jdragelj/python/karlsbergplus/pdbs/cco_prepared_ref_new.pdb'
    # folder = '/user/jdragelj/python/CcO/Cco_final/Cco_apbs4_eps4_m/'
    # folder = '/user/jdragelj/python/CcO/Cco_final/Cco_apbs4_eps1_m/'



    # TYR not shifted - crystal structure
    # pdb = '/user/jdragelj/python/karlsbergplus/pdbs/cco_start_no_wbox_no_memb_no_water_final_hsw.pdb'
    # like the file above but shorter name only!!!!
    pdb = '/user/jdragelj/python/karlsbergplus/pdbs/cco_cut_final_hsw.pdb'
    # pdb = '/user/jdragelj/python/karlsbergplus/pdbs/cco_tyrc.pdb'

    # folder = '/user/jdragelj/python/CcO/Cco_final/Cco_apbs4_eps4_c/'
    # folder = '/user/jdragelj/python/CcO/Cco_final/Cco_tyrc/'
    # folder = '/user/jdragelj/python/CcO/MD/'
    # folder = '/user/jdragelj/python/CcO/MD/TYR244_all/'

    # folder = '/user/jdragelj/python/CcO/MD/basic_tyr_1ns/'
    folder = '/user/jdragelj/python/CcO/MD/test'
    # folder = '/user/jdragelj/python/CcO/MD/basic_tyr_half/'
    # folder = '/user/jdragelj/python/CcO/MD/basic_tyr_tenth/'
    #


    #patch block todo: move to specific class info for cco
    patches = []
    patches.append(['GLUP', ('GLU', 15, 'BCHA'), 'prot'])
    patches.append(['GLUP', ('GLU', 203, 'ACHA'), 'prot'])
    patches.append(['GLUP', ('GLU', 131, 'BCHA'), 'prot'])
    patches.append(['ASPP', ('ASP', 372, 'ACHA'), 'prot'])


    #after Auto Gene Angl Dihe - patches as comments todo: move to specific class info for cco
    spec_modelling = []

    command = '''!bonds towards heme and copper
PATCH PHEM ACHA 384 EHEM 2 SETUP
PATCH PHEM ACHA 72 GHEM 2 SETUP
PATCH PHE2 ACHA 386 GHEM 2 SETUP
PATCH EISO FEOH 1 EHEM 2 SETUP

PATCH CUBP META 1 ACHA 233 ACHA 282 ACHA 283 HOHC 1 SETUP
PATCH CUAP META 2 META 3 BCHA 114 BCHA 149 BCHA 151 BCHA 153 BCHA 157 BCHA 160 SETUP

'''

    spec_modelling.append([command,'autogen'])

    command = '''!charge patches
PATCH HEB3 GHEM 2 ACHA 72 ACHA 386 SETUP
PATCH CA21 META 2 META 3 BCHA 114 BCHA 149 BCHA 151 BCHA 153 BCHA 157 BCHA 160 SETUP
PATCH CBPP META 1 ACHA 282 ACHA 283 HOHC 1 ACHA 237 ACHA 233 SETUP
PATCH A3H4 EHEM 2 ACHA 384 FEOH 1 SETUP

'''
    spec_modelling.append([command,'autogen'])


    spec_modelling.append(['!no comment', 'autogen'])

    ######Pathway definition#######
    # pathway = [['basic', None, None, None], \
    #            ('H2O', 'TAH3', 614, 1), \
    #            ['SER', 'ACHA', 309, 1,], \
    #            ['TYR', 'ACHA', 244, 1,], \
    #            ['THR', 'ACHA', 312, 1,], \
    #            ['TYR', 'ACHA', 248, 1,], \
    #            ['H2O', 'TAH3', 630, 1], \
    #            ['THR', 'ACHA', 315, 1], \
    #            ['GLU', 'BCHA', 15, 0]]

    pathway = [['basic', None, None, None], \
               ['TYR', 'ACHA', 244, 1,]]

    ##############################


    ###function that models the pathway#### - todo: give charmm object and not top, par, pdb ... rest to class and ONLY leave the CHOICE of freedom!
    charmm_struct_list, model_struct_list, folder_list = pathway_modelling(protein, top, par, pdb, folder, pathway, patches=patches, spec_modelling=spec_modelling, \
                      consecutive=True, sidechain_envir_relax=True)

    ref_charmm_struct = charmm_struct_list[0]
    # titratable_yaml = '/scratch/scratch/tmeyer/kbplus2/titratable.yaml'
    titratable_yaml = '/user/jdragelj/python/toppar_yaml/titratable_standard_pathway.yaml'

    top = []
    # top.append('/user/jdragelj/python/toppar_yaml/top_alw_plus_lipids_wo_charge_kb.inp')
    top.append("/user/jdragelj/charmm_toppar/cco/top_alw_jd.inp")
    # top.append("/user/jdragelj/charmm_toppar/cco/top_alw_jd_advanced_kb.inp")
    top.append("/scratch/scratch/awoelke/md_cco/toppar/top_all36_lipid.rtf")
    top.append('/user/jdragelj/python/toppar_yaml/patches_hyd.rtf')

    par = []
    par.append('/user/jdragelj/python/toppar_yaml/par_all22_prot_cut.inp')
    # par.append("/scratch/scratch/awoelke/md_cco/toppar/par_all36_lipid.prm")
    par.append('/user/jdragelj/python/toppar_yaml/patches.prm')

    #function that does tapbs
    g_list, pkint_list, energy_list, residue_list_ext= pathway_tapbs_apbs(ref_charmm_struct, model_struct_list, folder_list, titratable_yaml, pathway, top, par)

    # for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:

    ener_print_list = []
    ener_print_list_final = []
    for (resname, segname, resid, state), (solv, coul), pkint_list_entry in zip(pathway, energy_list, pkint_list):
        if resname == 'basic':
            ener_print = resname + ' ' + str(solv) + ' ' + str(coul)
            residue_descr = resname
            ener_print = residue_descr + ' ' + str(solv) + ' ' + str(coul)
        else:
            if resname == 'GLU':
                resname = 'EPP'
            residue_descr = "%s-%i_%s" % (resname, resid, segname)
            ener_print = residue_descr + ' ' + str(solv) + ' ' + str(coul)

            if pkint_list_entry[residue_descr] is not None:
                number_of_pkas = len(pkint_list_entry[residue_descr])

                # ener_len.append(number_of_pkas)

            for j in range(0, number_of_pkas):
                ener_print += (' ' + str(pkint_list_entry[residue_descr][j]))
            ener_print_list.append(ener_print)
        # smart calc - check it!!!!
        if resname == 'EPP':
            ener_print_final = pkint_list_entry[residue_descr][1] - pkint_list_entry[residue_descr][2]
        elif resname != 'basic':
            ener_print_final = pkint_list_entry[residue_descr][1] - pkint_list_entry[residue_descr][0]
        else:
            ener_print_final = 0.0
        # Coulomb + solvation
        ener_print_final += solv + coul
        ener_print_list_final.append(ener_print_final)
    # choice of referent state - one of the energies or 0.0 (abs)
    # e_ref = ener_print_list_final[1]
    e_ref = 0.0
    ener_print_list_final = [x - e_ref for x in ener_print_list_final]


    print "Results: Yay!"
    print "residue solvation coulomb"
    print "crystal"
    for l in range(len(ener_print_list)):
        print ener_print_list[l]

    print '-----------------'
    print "Results: Woohoo"
    print "residue solvation coulomb pkint_final"
    print "crystal"
    for energy in ener_print_list_final:
        print energy

    # #todo: implement smart calculation of pKaints
    # #todo: make pqr compare (only topology) - similar to transfer coords





