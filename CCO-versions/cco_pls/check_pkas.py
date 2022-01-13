# coding=utf-8

import re
import kbp2
import matplotlib.pyplot as plt
import numpy as np


def store_pkas(files):
    pkas = {}
    for file in files:
        f = open(file,'r')
        for i, line in enumerate(f):
            if line:
                line = line.strip()
                components = line.split(':')
                resname, resid, segname = re.split('[-_]', components[0])
                if resname == 'EPP':
                    resname = 'GLU'
                if resname == 'DPP':
                    resname = 'ASP'
                if resname == 'HSP':
                    resname = 'HIS'
                identifier = resname + '-' + resid + '_' + segname
                if '>20' in components[1]:
                    components[1] = 20.0
                elif '<-10' in components[1]:
                    components[1] = -10.0
                pka = round(float(components[1]),1)

                if identifier in pkas.keys():
                    pkas[identifier].append(pka)
                else:
                    pkas[identifier]=[]
                    pkas[identifier].append(pka)
        f.close()
    return pkas

def compare_pkas(pkas_dict):

    pkas_to_check = {}
    for residue, pka_list in pkas_dict.iteritems():
        for pka in pka_list:
            if pka >= 6.00 and pka<=8.00:
                if residue not in pkas_to_check.keys():
                    pkas_to_check[residue] = pka_list
    return  pkas_to_check


def check_protonation(files):

    pkas_dict = store_pkas(files)
    pkas_to_check = compare_pkas(pkas_dict)
    for residue, pka_list in pkas_to_check.iteritems():
        avg = round(sum(pka_list)/len(pka_list),2)
        if round(avg,1) > 5.5 and round(avg,1) <8.50:
            print residue, pka_list, round(sum(pka_list)/len(pka_list),2)
        # print residue, pka_list, round(sum(pka_list)/len(pka_list),2)
    return


if __name__ == '__main__':

    # mds = ['md_Pr_prd-', 'md_Pr_prdh1', 'md_Pr_prdh2']
    # for md in mds:
    #     print 'cav_0.8_0.2', 'cav_0.9_0.2', 'cav_1.0_0.2'
    #     files = []
    #     for i in np.arange(80, 801, 80):
    #         file1 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_cavities_voda/%s/done/frame%i/results.dat' % (md,i)
    #         f1 = open(file1, 'r')
    #         file2 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_cavities_09_voda/%s/done/frame%i/results.dat' % (md,i)
    #         f2 = open(file2, 'r')
    #         file3 = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/10ang_cavities_1_voda/%s/done/frame%i/results.dat' % (md,i)
    #         f3 = open(file3, 'r')
    #
    #
    #         file1s = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/cavities_voda/%s/done/frame%i/results.dat' % (md,i)
    #         f1s = open(file1s, 'r')
    #         file2s = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/cavities_09_voda/%s/done/frame%i/results.dat' % (md,i)
    #         f2s = open(file2s, 'r')
    #         file3s = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/Pr_xtra/prd/cavities_1_voda/%s/done/frame%i/results.dat' % (md,i)
    #         f3s = open(file3s, 'r')
    #
    #         for k, line in enumerate(f1):
    #             if line:
    #                 line = line.strip()
    #                 components = line.split(':')
    #                 resname, resid, segname = re.split('[-_]', components[0])
    #                 if resname == 'PRD' and segname == 'EHEM':
    #                     pka1 = round(float(components[1]), 1)
    #         for k, line in enumerate(f2):
    #             if line:
    #                 line = line.strip()
    #                 components = line.split(':')
    #                 resname, resid, segname = re.split('[-_]', components[0])
    #                 if resname == 'PRD' and segname == 'EHEM':
    #                     pka2 = round(float(components[1]), 1)
    #         for k, line in enumerate(f3):
    #             if line:
    #                 line = line.strip()
    #                 components = line.split(':')
    #                 resname, resid, segname = re.split('[-_]', components[0])
    #                 if resname == 'PRD' and segname == 'EHEM':
    #                     pka3 = round(float(components[1]), 1)
    #
    #         for k, line in enumerate(f1s):
    #             if line:
    #                 line = line.strip()
    #                 components = line.split(':')
    #                 resname, resid, segname = re.split('[-_]', components[0])
    #                 if resname == 'PRD' and segname == 'EHEM':
    #                     pka1s = round(float(components[1]), 1)
    #         for k, line in enumerate(f2s):
    #             if line:
    #                 line = line.strip()
    #                 components = line.split(':')
    #                 resname, resid, segname = re.split('[-_]', components[0])
    #                 if resname == 'PRD' and segname == 'EHEM':
    #                     pka2s = round(float(components[1]), 1)
    #         for k, line in enumerate(f3s):
    #             if line:
    #                 line = line.strip()
    #                 components = line.split(':')
    #                 resname, resid, segname = re.split('[-_]', components[0])
    #                 if resname == 'PRD' and segname == 'EHEM':
    #                     pka3s = round(float(components[1]), 1)
    #
    #         # print i,':',pka1, pka1s, 'd', abs(pka1-pka1s),',',pka2,pka2s,'d',abs(pka3-pka3s),',',pka3,pka3s,'d',abs(pka3-pka3s)
    #         print i,':',pka1, pka1s, ',',pka2,pka2s,',',pka3,pka3s
    #         f1.close()
    #         f2.close()
    #         f3.close()
    #         f1s.close()
    #         f2s.close()
    #         f3s.close()
    #     print '---------------------------------------------------------------------'





    # pdb = '/scratch/scratch/jdragelj/projects/cco_pls/cco_pls_prd/model_open_sb_2GSM/quick_with_water/open_sb_Pr/cco_and_water.pdb'
    # results = get_his_protonation(pdb)
    # for residue, pka_list in pkas.iteritems():
    #     if 'HIS' in residue:
    #         resname, resid, segname = re.split('[-_]', residue)
    #         ident = str(resid) + '_' + segname
    #         for key in results.keys():
    #             if ident == key:
    #                 avg = round(sum(pka_list) / len(pka_list), 2)
    #                 if results[key] == 'HSP' and round(avg,1) < 6.6 :
    #                     print results[key], residue, pka_list, avg
    #                 if results[key] in ['HSE','HSD'] and round(avg, 1) > 7.4 :
    #                     print results[key], residue, pka_list, avg
    #
    #                 # if avg > 7.00 and results[key] != 'HSP':
    #                 #     print results[key], residue, pka_list, avg
    #                 # if avg < 7.00 and results[key] not in ['HSD', 'HSE']:
    #                 #     print results[key], residue, pka_list, avg

    # # file1 = '/scratch/scratch/jdragelj/tests/alw_6565/results.dat'
    # file1 = '/media/jdragelj/88bc6e86-bfdd-41a3-b3f0-9921ae91f5b7/awoelke/md_glu101b/f-state/deprot_a3h4_cbpn/md/output/frames_mem_wo/done/frame0/pka_results.dat'
    # # file2 = '/scratch/scratch/jdragelj/tests/alw_9765/results.dat'
    # # file2 = '/scratch/scratch/jdragelj/tests/alw_9797/results.dat'
    # file2 = '/scratch/scratch/jdragelj/tests/kbp_test/results.dat'
    # files = []
    # files.append(file1)
    # files.append(file2)
    # pkas = store_pkas(files)
    #
    # for residue, pka_values in pkas.iteritems():
    #     if len(pka_values) == 2:
    #         if pka_values[0] != pka_values[1]:
    #             if abs(pka_values[0] - pka_values[1]) >= 0.05:
    #                 print residue, pka_values, abs(pka_values[0] - pka_values[1])
    #     else:
    #         print residue, pka_values




    main_workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models_cav09/'
    states = ['Pm', 'PF', 'Pr', 'F']
    prots = ['prd-', 'prdh1', 'prdh2']
    print 'prd a3'
    for state in states:
        print state
        files = []
        for prot in prots:
            file = main_workfolder + '%s/prda3/%s/results.dat' % (state, prot)
            files.append(file)
        pkas = store_pkas(files)
        print pkas['PRD-3_EHEM']
        print round(sum(pkas['PRD-3_EHEM'])/float(len(pkas['PRD-3_EHEM'])),2)
    print '---------'
    print 'pra a3'
    prots = ['pra-', 'prah1', 'prah2']
    for state in states:
        print state
        files = []
        for prot in prots:
            file = main_workfolder + '%s/praa3/%s/results.dat' % (state, prot)
            files.append(file)
        pkas = store_pkas(files)
        print pkas['PRA-1_EHEM']
        print round(sum(pkas['PRA-1_EHEM'])/float(len(pkas['PRA-1_EHEM'])),2)
    print '-----------'
    print 'pra a'
    prots = ['pra-', 'prah1', 'prah2']
    for state in states:
        print state
        files = []
        for prot in prots:
            file = main_workfolder + '%s/praa/%s/results.dat' % (state, prot)
            files.append(file)
        pkas = store_pkas(files)
        print pkas['PRA-1_GHEM']
        print round(sum(pkas['PRA-1_GHEM'])/float(len(pkas['PRA-1_GHEM'])),2)


    # main_workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models/'
    # states = ['Pm', 'PF', 'Pr', 'F']
    # prots = ['prd-', 'prdh1', 'prdh2']
    # residue_list = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM', 'ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', \
    #             'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', \
    #             'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # import cPickle as pickle
    # for state in states:
    #     combined_results = kbp2.kbp_results.FrameKbpResults()
    #     print state
    #     for prot in prots:
    #         file = main_workfolder + '%s/prda3/%s/results.pkl' % (state, prot)
    #         prot_result = pickle.load(open(file, "rb"))
    #         combined_results.add_task(prot_result.kbp_results[0])
    #     combined_results.combine_frames_karlsberg(cpus=3)
    #     final_pkas = combined_results.combined_results.pkas
    #     for residue_descr, pka in final_pkas.iteritems():
    #         if residue_descr == 'PRD-3_EHEM':
    #             s = " %13s: %0.2f " % (residue_descr, pka)
    #             print s
    # print '--------------'


    # main_workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models/'
    # states = ['Pm', 'PF', 'Pr', 'F']
    # prots = ['pra-', 'prah1', 'prah2']
    # residue_list = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM', 'ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', \
    #             'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', \
    #             'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # import cPickle as pickle
    # for state in states:
    #     combined_results = kbp2.kbp_results.FrameKbpResults()
    #     print state
    #     for prot in prots:
    #         file = main_workfolder + '%s/praa3/%s/results.pkl' % (state, prot)
    #         prot_result = pickle.load(open(file, "rb"))
    #         combined_results.add_task(prot_result.kbp_results[0])
    #     combined_results.combine_frames_karlsberg(cpus=3)
    #     final_pkas = combined_results.combined_results.pkas
    #     for residue_descr, pka in final_pkas.iteritems():
    #         if residue_descr == 'PRA-1_EHEM':
    #             s = " %13s: %0.2f " % (residue_descr, pka)
    #             print s
    # print '--------------'

    # main_workfolder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/cycle_models/'
    # states = ['Pm', 'PF', 'Pr', 'F']
    # prots = ['pra-', 'prah1', 'prah2']
    # residue_list = ['PRD-3_EHEM', 'PRA-1_EHEM', 'PRD-3_GHEM', 'PRA-1_GHEM', 'ARG-481_ACHA', 'ARG-482_ACHA', 'GLU-286_ACHA', \
    #             'ASP-407_ACHA', 'HSD-411_ACHA', 'ASP-412_ACHA', 'ARG-52_ACHA', 'TYR-414_ACHA', 'ASP-229_BCHA', 'LYS-227_BCHA', 'TYR-336_ACHA', \
    #             'TYR-415_ACHA', 'ARG-408_ACHA', 'TYR-175_ACHA']
    # import cPickle as pickle
    # for state in states:
    #     combined_results = kbp2.kbp_results.FrameKbpResults()
    #     print state
    #     for prot in prots:
    #         file = main_workfolder + '%s/praa/%s/results.pkl' % (state, prot)
    #         prot_result = pickle.load(open(file, "rb"))
    #         combined_results.add_task(prot_result.kbp_results[0])
    #     combined_results.combine_frames_karlsberg(cpus=3)
    #     final_pkas = combined_results.combined_results.pkas
    #     for residue_descr, pka in final_pkas.iteritems():
    #         if residue_descr == 'PRA-1_EHEM':
    #             s = " %13s: %0.2f " % (residue_descr, pka)
    #             print s
    # print '--------------'












