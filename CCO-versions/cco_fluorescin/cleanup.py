# coding=utf-8

import os
import numpy as np

if __name__ == '__main__':

    # folders = ['/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/', '/media/88bc6e86-bfdd-41a3-b3f0-9921ae91f5b7/jovan/backup/scratch/projects/cco_alexiev/MD_299/titrations/titrations_with_flu']
    # for folder in folders:
    #     positions = [0,90,180,270]
    #     for pos in positions:
    #         fold = folder + str(pos) + '/'
    #         subf = os.listdir(fold)
    #
    #         for sub in subf:
    #             foldd = fold + sub + '/done/'
    #             framsubs = os.listdir(foldd)
    #
    #             for frame in framsubs:
    #                 if os.path.exists(foldd + frame + '/pac_ph7_h_min/'):
    #                     shutil.rmtree(foldd + frame + '/pac_ph7_h_min/')
    #                 if os.path.exists(foldd + frame + '/initial_modelling/'):
    #                     shutil.rmtree(foldd + frame + '/initial_modelling/')

    mds_sourcefolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/'
    global_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/energy_calculations/new_charmm_total_notermini/'
    frames = np.arange(0,445,2)
    tbd = []
    positions = [90]


    while True:
        tbd = 0
        for position in positions:
            position_workfolder =  global_workfolder + str(position) + '/'
            if not os.path.exists(position_workfolder):
                continue
            subfolders = os.listdir(position_workfolder)
            chosen_trajectories = []
            for subfolder in subfolders:
                # if 'hsp73' not in subfolder:
                #     continue
                chosen_trajectories.append(subfolder)
            for chosen_trajectory in chosen_trajectories:
                trajectory_workfolder = position_workfolder + chosen_trajectory + '/'
                if not os.path.exists(trajectory_workfolder):
                    continue
                for frame in frames:
                    frame_workfolder = trajectory_workfolder + '%i/' % frame
                    if not os.path.exists(frame_workfolder):
                        continue
                    components = ['total']
                    for component in components:
                        new_frame_workfolder = frame_workfolder + component + '/'
                        # new_frame_workfolder = frame_workfolder
                    if os.path.exists(new_frame_workfolder+'frame%i_charmm.out' % frame):
                        out_file = new_frame_workfolder+'frame%i_charmm.out' % frame
                        f_o = open(out_file, 'r')
                        done = False
                        for line in f_o:
                            if 'NORMAL TERMINATION' in line:
                                done = True
                        if done:
                            file_to_delete = os.listdir(new_frame_workfolder)
                            for file_o in file_to_delete:
                                if '_charmm.out' in file_o:
                                        continue
                                elif '_charmm.inp' in file_o:
                                    continue
                                else:
                                    tbd +=1
                                    os.remove(new_frame_workfolder+file_o)


    # folders = ['/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/hsd73/', \
    #            '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/old_flu000/', \
    #            '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/old_oxi_rmsd_weird/', \
    #            '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/old_red_rmsd_weird/', \
    #            '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/pulling/', \
    #            ]
    #
    # for folder in folders:
    #     patterns = os.listdir(folder)
    #     for patt in patterns:
    #         if os.path.exists(folder+patt+'/done/'):
    #             frame_label_folder_list = os.listdir(folder+patt+'/done/')
    #             for frame_folder_label in frame_label_folder_list:
    #                 frame_folder = folder + patt + '/done/' + frame_folder_label + '/'
    #                 files = os.listdir(frame_folder)
    #                 for fil in files:
    #                     if '_init.pdb' in fil:
    #                         os.remove(frame_folder + fil)
    #                     if '_last_reference.pdb' in fil:
    #                         os.remove(frame_folder + fil)


    # folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/'
    # selected_mds_after_interaction_energy = {'hsp73_hsp526_red_flu-_2': [90], \
    #                                          'hsp73_hsp526_red_fluh': [90], \
    #                                          'hsp73_hsp526_red_flu2': [90], \
    #                                          'hsp73_hsp526_fluh': [90], \
    #                                          'hsp73_hsp526_flu2': [90], \
    #                                          'hsp73_hsp526_flu-_2': [90], \
    #                                          'hsp73_hse526_flu-': [90], \
    #                                          'hsp73_hse526_fluh': [90], \
    #                                          'hsp73_hse526_flu2': [90], \
    #                                          'hsp73_hse526_red_flu-': [90], \
    #                                          'hsp73_hse526_red_fluh': [90], \
    #                                          'hsp73_hse526_red_flu2': [90], \
    #                                          'hsd73_hse526_red_flu-': [90], \
    #                                          'hsd73_hse526_flu-': [90], \
    #                                          }
    # mds_selection_reorganized = {}
    # mds_selection_reorganized[0] = []
    # mds_selection_reorganized[90] = []
    # mds_selection_reorganized[180] = []
    # mds_selection_reorganized[270] = []
    # for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
    #     for position in position_list:
    #         if pattern not in mds_selection_reorganized[position]:
    #             mds_selection_reorganized[position].append(pattern)
    #
    # print mds_selection_reorganized[90]
    # for pattern in mds_selection_reorganized[90]:
    #     frames = np.arange(0,445,2)
    #     for frame in frames:
    #         frame_folder = folder + pattern + '/done/frame%i/' % frame
    #         if os.path.exists(frame_folder):
    #             files = os.listdir(frame_folder)
    #             for fil in files:
    #                 if '_init.pdb' in fil:
    #                     os.remove(frame_folder + fil)
    #                 if '_last_reference.pdb' in fil:
    #                     os.remove(frame_folder + fil)





