# coding=utf-8

import os
from kbp2.workspace_jd.cco_fluorescin import interaction_ener

import numpy as np

# def fmt(x, y):
#     return 'x: {x:0.2f}\ny: {y:0.2f}'.format(x = x, y = y)
#
# class DataCursor(object):
#     # http://stackoverflow.com/a/4674445/190597
#     """A simple data cursor widget that displays the x,y location of a
#     matplotlib artist when it is selected."""
#     def __init__(self, artists, x = [], y = [], tolerance = 5, offsets = (-20, 20),
#                  formatter = fmt, display_all = False):
#         """Create the data cursor and connect it to the relevant figure.
#         "artists" is the matplotlib artist or sequence of artists that will be
#             selected.
#         "tolerance" is the radius (in points) that the mouse click must be
#             within to select the artist.
#         "offsets" is a tuple of (x,y) offsets in points from the selected
#             point to the displayed annotation box
#         "formatter" is a callback function which takes 2 numeric arguments and
#             returns a string
#         "display_all" controls whether more than one annotation box will
#             be shown if there are multiple axes.  Only one will be shown
#             per-axis, regardless.
#         """
#         self._points = np.column_stack((x,y))
#         self.formatter = formatter
#         self.offsets = offsets
#         self.display_all = display_all
#         if not cbook.iterable(artists):
#             artists = [artists]
#         self.artists = artists
#         self.axes = tuple(set(art.axes for art in self.artists))
#         self.figures = tuple(set(ax.figure for ax in self.axes))
#
#         self.annotations = {}
#         for ax in self.axes:
#             self.annotations[ax] = self.annotate(ax)
#
#         for artist in self.artists:
#             artist.set_picker(tolerance)
#         for fig in self.figures:
#             fig.canvas.mpl_connect('pick_event', self)
#
#     def annotate(self, ax):
#         """Draws and hides the annotation box for the given axis "ax"."""
#         annotation = ax.annotate(self.formatter, xy = (0, 0), ha = 'right',
#                 xytext = self.offsets, textcoords = 'offset points', va = 'bottom',
#                 bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
#                 arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0')
#                 )
#         annotation.set_visible(False)
#         return annotation
#
#     def snap(self, x, y):
#         """Return the value in self._points closest to (x, y).
#         """
#         idx = np.nanargmin(((self._points - (x,y))**2).sum(axis = -1))
#         return self._points[idx]
#     def __call__(self, event):
#         """Intended to be called through "mpl_connect"."""
#         # Rather than trying to interpolate, just display the clicked coords
#         # This will only be called if it's within "tolerance", anyway.
#         x, y = event.mouseevent.xdata, event.mouseevent.ydata
#         annotation = self.annotations[event.artist.axes]
#         if x is not None:
#             if not self.display_all:
#                 # Hide any other annotation boxes...
#                 for ann in self.annotations.values():
#                     ann.set_visible(False)
#             # Update the annotation in the current axis..
#             x, y = self.snap(x, y)
#             annotation.xy = x, y
#             annotation.set_text(self.formatter(x, y))
#             annotation.set_visible(True)
#             event.canvas.draw()
#
#

if __name__ == '__main__':


    # #####################
    # ## regular frames ###
    # #####################
    # mds_sourcefolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/'
    # global_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/interaction_ene_mds/'
    # frames = np.arange(24,105,1)
    # frame_folder_name = 'frames_1'
    #
    # # positions = [0,90,180,270]
    # positions = [0]
    # for position in positions:
    #     position_workfolder =  global_workfolder + str(position) + '/'
    #     if not os.path.exists(position_workfolder):
    #         os.mkdir(position_workfolder)
    #     subfolders = os.listdir(mds_sourcefolder + str(position) + '/')
    #     chosen_trajectories = []
    #     for subfolder in subfolders:
    #         if '.sh' in subfolder:
    #             continue
    #         if 'hsd73_hsp526_flu-' not in subfolder:
    #             continue
    #         chosen_trajectories.append(subfolder)
    #     for trajectory in chosen_trajectories:
    #         trajectory_workfolder =  position_workfolder + trajectory + '/'
    #         if not os.path.exists(trajectory_workfolder):
    #             os.mkdir(trajectory_workfolder)
    #         for frame in frames:
    #             frame_workfolder = trajectory_workfolder + str(frame) + '/'
    #             if not os.path.exists(frame_workfolder):
    #                 os.mkdir(frame_workfolder)
    #             frame_pdb = mds_sourcefolder + str(position) + '/' + trajectory + '/md/%s/frame%i.pdb' % (frame_folder_name, frame)
    #             # APBS CALCULATION
    #             calculations_folder = frame_workfolder
    #             modelling_folder = calculations_folder + 'modelling/'
    #             apbs_folder = calculations_folder + 'apbs/'
    #             if os.path.exists(calculations_folder + 'interaction_energy.dat'):
    #                 if os.path.exists(modelling_folder):
    #                     shutil.rmtree(modelling_folder)
    #                 if os.path.exists(apbs_folder):
    #                     shutil.rmtree(apbs_folder)
    #             else:
    #                 residues_to_exclude = ['FLU', 'FLX', 'FLH']
    #                 interaction_energy.produce_interaction_energy(frame_pdb, position_workfolder, calculations_folder, residues_to_exclude,
    #                                                               trajectory, minimize = False)

    # components = ['total', 'flu', 'protein']
    # folder1 = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/charmm_binding_energy/4_80/180/hsd73_hse526_flu-/56/'
    # frame56 = interaction_energy.get_binding_energy(folder1, components)
    # print frame56
    # folder2 = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/charmm_binding_energy/4_80/180/hsd73_hse526_flu-/100/'
    # frame100 = interaction_energy.get_binding_energy(folder2, components)
    # print frame100

    ########################################
    ### charmm TOTAL energy PBC AND GBMV ###
    ########################################

    mds_sourcefolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/'
    global_workfolder_prev = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/energy_calculations/new_charmm_total_notermini/'

    # print 'Data in Supplement are in: /media/jdragelj/88bc6e86-bfdd-41a3-b3f0-9921ae91f5b7/jovan/backup/scratch/projects/cco_alexiev/MD_299/energy_calculations/charmm_binding_energy/4_80_membrane/with_pbc/'

    frames = np.arange(0,445,2)
    frame_folder_name = 'frames_1'
    tbd = []

    positions = [90]
    # stem = 'md_ex_10'
    stem = 'md'

    global_workfolder = global_workfolder_prev
    if not os.path.exists(global_workfolder):
        os.mkdir(global_workfolder)
    for position in positions:
        position_workfolder =  global_workfolder + str(position) + '/'
        if not os.path.exists(position_workfolder):
            os.mkdir(position_workfolder)
        subfolders = os.listdir(mds_sourcefolder + str(position) + '/')
        chosen_trajectories = []
        subfolders = ['hsd73_hse526_red_flu-']
        for subfolder in subfolders:
            chosen_trajectories.append(subfolder)
        for chosen_trajectory in chosen_trajectories:
            trajectory_workfolder = position_workfolder + chosen_trajectory + '/'
            if not os.path.exists(trajectory_workfolder):
                os.mkdir(trajectory_workfolder)
            for frame in frames:
                frame_workfolder = trajectory_workfolder + '%i/' % frame
                if not os.path.exists(frame_workfolder):
                    os.mkdir(frame_workfolder)
                    components = ['total']
                    for component in components:
                        new_frame_workfolder = frame_workfolder + component + '/'
                        if not os.path.exists(new_frame_workfolder):
                            os.mkdir(new_frame_workfolder)
                        frame_pdb = mds_sourcefolder + str(position) + '/' + chosen_trajectory + '/%s/%s/frame%i.pdb' % (stem, frame_folder_name, frame)
                        if os.path.exists(frame_pdb):
                            if not os.path.exists(new_frame_workfolder+'frame%i_charmm.out' % frame):
                                step_to_find = frame*25000
                                x_y_z = []
                                pbc_file = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/' + str(position) + '/' + chosen_trajectory + '/md/output/3hb3_md_flu.xst'
                                if frame < 445:
                                    pbc_file = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/' + str(position) + '/' + chosen_trajectory + '/md/output/3hb3_md_flu.xst'
                                else:
                                    pbc_file = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/' + str(position) + '/' + chosen_trajectory + '/%s/output/3hb3_md_flu.xst' % stem
                                pbc_f = open(pbc_file, 'r')
                                for line in pbc_f:
                                    line = line.split()
                                    if line[0] == str(step_to_find):
                                        x_y_z.append(line[1])
                                        x_y_z.append(line[5])
                                        x_y_z.append(line[9])
                                interaction_ener.prepare_structure(frame_pdb, new_frame_workfolder, chosen_trajectory, minimize = False, pqr=False, x_y_z=x_y_z,
                                                    inter=False, delete_struct=component, apbs_calc=False, gen_born='gbmv', sasa=False, submitt=False)
                            else:
                                continue
                        else:
                            print mds_sourcefolder + str(position) + '/' + chosen_trajectory
                            tbd.append(mds_sourcefolder + str(position) + '/' + chosen_trajectory)


    # #
    # # # ###################
    # # # ### charmm sasa ###
    # # # ###################
    #
    # mds_sourcefolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/'
    # # global_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/sasa_calcs/fluorescein/'
    # # global_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/sasa_calcs/His73/'
    # global_workfolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/sasa_calcs/His73_new/'
    # sasa = 1.4
    #
    # selected_mds_after_interaction_energy = {'hsp73_hsp526_red_protonate_fluh': [90], \
    #                                          'hsp73_hsp526_flu-': [90], \
    #                                          'hsd73_hse526_flu-': [90], \
    #                                          'hsd73_hse526_red_flu-': [90]}
    #
    # # selected_mds_after_interaction_energy = {'hsp73_hsp526_red_protonate_fluh': [90], \
    # #                                          'hsd73_hse526_red_flu-': [90]}
    #
    # frames = np.arange(4,45,2)
    #
    # mds_selection_reorganized = {}
    # mds_selection_reorganized[90] = []
    # for pattern, position_list in selected_mds_after_interaction_energy.iteritems():
    #     for position in position_list:
    #         if pattern not in mds_selection_reorganized[position]:
    #             mds_selection_reorganized[position].append(pattern)
    #
    # for position, pattern_list in mds_selection_reorganized.iteritems():
    #     position_workfolder =  global_workfolder + str(position) + '/'
    #     if not os.path.exists(position_workfolder):
    #         os.mkdir(position_workfolder)
    #     for pattern in pattern_list:
    #
    #         # if pattern == 'hsd73_hse526_flu-':
    #         #     frames = np.arange(44, 685, 2)
    #         # elif pattern == 'hsd73_hse526_red_flu-':
    #         #     frames = np.arange(44, 445, 2)
    #
    #         trajectory_workfolder = position_workfolder + pattern + '/'
    #         if not os.path.exists(trajectory_workfolder):
    #             os.mkdir(trajectory_workfolder)
    #         for fr in frames:
    #             frame_workfolder = trajectory_workfolder + '%i/' % fr
    #             frame_pdb = mds_sourcefolder + str(position) + '/' + pattern + '/md/frames_1/frame%i.pdb' % fr
    #
    #             if pattern == 'hsd73_hse526_flu-':
    #                 frame_pdb = mds_sourcefolder + str(position) + '/' + pattern + '/md_ex_17_ex_12/frames_1/frame%i.pdb' % fr
    #             elif pattern == 'hsd73_hse526_red_flu-':
    #                 frame_pdb = mds_sourcefolder + str(position) + '/' + pattern + '/md_ex_17/frames_1/frame%i.pdb' % fr
    #
    #             if os.path.exists(frame_pdb):
    #                 if not os.path.exists(frame_workfolder+'frame%i_charmm.out' % fr):
    #                     interaction_ener.prepare_structure(protein=frame_pdb, modelling_folder=frame_workfolder, trajectory=pattern,
    #                                                        minimize=False, pqr=False, x_y_z=None, inter=False, sasa=1.4,
    #                                                        delete_struct='', apbs_calc=False, gen_born=None, submitt=True)


    # mds_sourcefolder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/wild_type_mds/'
    #
    # pattern_list = ['hsp73_hsp526', 'hsp73_hsp526_red', 'hsd73_hse526']
    # pattern_list = ['hsd73_hse526']
    # position_workfolder = global_workfolder + 'wt/'
    #
    # for pattern in pattern_list:
    #
    #     if pattern == 'hsp73_hsp526':
    #         frames = np.arange(4, 305, 2)
    #     else:
    #         frames = np.arange(4, 205, 2)
    #
    #     trajectory_workfolder = position_workfolder + pattern + '/'
    #     if not os.path.exists(trajectory_workfolder):
    #         os.mkdir(trajectory_workfolder)
    #     for fr in frames:
    #         frame_workfolder = trajectory_workfolder + '%i/' % fr
    #         frame_pdb = mds_sourcefolder + pattern + '/md/frames_1/frame%i.pdb' % fr
    #         if pattern == 'hsp73_hsp526':
    #             frame_pdb = mds_sourcefolder + pattern + '/md_ex_10/frames_1/frame%i.pdb' % fr
    #         if os.path.exists(frame_pdb):
    #             if not os.path.exists(frame_workfolder + 'frame%i_charmm.out' % fr):
    #                 interaction_ener.prepare_structure(protein=frame_pdb, modelling_folder=frame_workfolder,
    #                                                    trajectory=pattern,
    #                                                    minimize=False, pqr=False, x_y_z=None, inter=False, sasa=1.4,
    #                                                    delete_struct='', apbs_calc=False, gen_born=None, submitt=True)








