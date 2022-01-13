# coding=utf-8

import os
import shutil
import sys

from time import sleep
import re


# print 'With Fluorescein'
# big_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/'
# subfolders = ['hsp73_hsp526_red_flu-_2', \
#                  'hsp73_hsp526_red_fluh', \
#                  'hsp73_hsp526_red_flu2', \
#                  'hsp73_hse526_red_flu-', \
#                  'hsp73_hse526_red_fluh', \
#                  'hsp73_hse526_red_flu2']
#
# # subfolders = ['hsp73_hsp526_fluh', \
# #                  'hsp73_hsp526_flu2', \
# #                  'hsp73_hsp526_flu-_2', \
# #                  'hsp73_hse526_flu-', \
# #                  'hsp73_hse526_fluh', \
# #                  'hsp73_hse526_flu2']
# index = -1
# for subfolder in subfolders:
#     tcl_script_file = '/user/jdragelj/Desktop/traj_script_with_flu_R.tcl'
#     # tcl_script_file = '/user/jdragelj/Desktop/traj_script_with_flu_O.tcl'
#     index += 1
#     tcl_script = """\n"""
#     tcl_script += """
# menu files on
# display resetview
# mol new {%s%s/md_ex_10/output/3hb3_md_flu_ex.dcd} type {dcd} first 0 last -1 step 2 waitfor 1 """ % (big_folder, subfolder) + """
# animate style Loop
# mol addrep %i """ % index + """
# display resetview
# mol addfile {%s%s/change/cco_and_water.psf} type {psf} first 0 last -1 step 2 waitfor 1 %i """ % (big_folder, subfolder, index) + """
# animate style Loop
# mol rename %i {%s} """ % (index, subfolder) + """
# menu files off
#     """
#     f = open(tcl_script_file, 'a')
#     f.write(tcl_script)
#     f.close()
#
# ##############################################################

print 'wild-type'
main_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/wild_type_mds/'
folder = main_folder
subfolders = os.listdir(folder)
subfolder_new = []
for subfolder in subfolders:
    if '526' not in subfolder:
        continue
    subfolder_new.append(subfolder)
index = -1
for subfolder in subfolder_new:
    tcl_script_file = '/user/jdragelj/Desktop/traj_script_wild-type.tcl'
    index += 1
    tcl_script = """\n"""
    tcl_script += """
menu files on
display resetview
mol new {%s%s/md_ex_10/output/3hb3_md_ex.dcd} type {dcd} first 0 last -1 step 2 waitfor 1 """ % (main_folder, subfolder) + """
animate style Loop
mol addrep %i """ % index + """
display resetview
mol addfile {%s%s/change/intermediate_out.psf} type {psf} first 0 last -1 step 2 waitfor 1 %i """ % (main_folder, subfolder, index) + """
animate style Loop
mol rename %i {%s} """ % (index, subfolder) + """
menu files off
    """
    f = open(tcl_script_file, 'a')
    f.write(tcl_script)
    f.close()