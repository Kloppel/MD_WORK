# coding=utf-8

import shutil
import os
import kbp2


if __name__ =="__main__":


    #### PREPARATION OF FILE FROM ALENA

    #removal of waterbox
    pdb = '/user/jdragelj/python/karlsbergplus/pdbs/cco_start.pdb'
    pdb_mod = kbp2.file_parser.Simple_struct_parser()
    pdb_mod.read_pdb(pdb)
    new = pdb_mod.copy(segname=['WAT','SWAT','WCH2','VBH2'],exclude=True)
    new.write_pdb('/user/jdragelj/python/karlsbergplus/pdbs/cco_start_no_wbox.pdb')


    #removal of parts of membrane
    pdb1 = '/user/jdragelj/python/karlsbergplus/pdbs/cco_start_no_wbox.pdb'
    pdb_mod_1 = kbp2.file_parser.Simple_struct_parser()
    pdb_mod_1.read_pdb(pdb1)

    print pdb_mod_1.top_content


    resids = [11,34,40,45,46,69,70,75,76,82,91,96,97,98,99,104,105,106,118]
    reslist = []
    for i in range(1,275):
        if i not in resids:
            reslist.append(('MEMB',i))
    new_1 = pdb_mod_1.copy(residuelist = reslist,exclude=True)
    new_1.write_pdb('/user/jdragelj/python/karlsbergplus/pdbs/cco_start_no_wbox_no_memb.pdb')


    #removal of waters
    pdb2 = '/user/jdragelj/python/karlsbergplus/pdbs/cco_start_no_wbox_no_memb.pdb'
    pdb_mod_2 = kbp2.file_parser.Simple_struct_parser()
    pdb_mod_2.read_pdb(pdb2)
    resids = [614,689,630]
    reslist = []
    for i in range(575,692):
        if i not in resids:
            reslist.append(('TAH2',i))
    new_2 = pdb_mod_2.copy(residuelist = reslist,exclude=True)
    new_2.read_top("/user/jdragelj/python/karlsbergplus/examples_kbp2_module/top.inp")
    print new_2.top_content
    new_3 = new_2.copy(residuelist = reslist,exclude=True)
    print new_3.top_content
    new_2.write_pdb('/user/jdragelj/python/karlsbergplus/pdbs/cco_start_no_wbox_no_memb_no_water.pdb')





