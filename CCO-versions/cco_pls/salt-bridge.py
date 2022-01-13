from workspace_jd import salt_bridge_distances
from kbp2.workspace_jd.cco_fluorescin import plot_utilities_frames
import numpy as np

if __name__ == '__main__':


    # for state in ['F', 'PF', 'Pm', 'Pr']:
    #     for prot in ['prdh1', 'prdh2', 'prd-']:
    #
    #         out_folder_general = '/user/jdragelj/projects/cco_pls/PRD/salt-bridges/%s/' % state
    #         frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_%s_xtra6/' % (state, state, prot)
    #
    #         y1_lim = [-10, 20]
    #         y1_ticks = np.arange(-10, 21, 1)
    #         y2_lim = [2, 8]
    #         y2_ticks = np.arange(2, 9, 1)
    #
    #         if state != 'Pr':
    #             x_lim = [0, 25]
    #             x_ticks = np.arange(0, 26, 2)
    #         else:
    #             x_lim = [0, 31]
    #             x_ticks = np.arange(0, 31, 2)
    #
    #         if state != 'Pr':
    #             frames = np.arange(0, 501, 4)
    #             time_at_frames1 = np.arange(0, 501, 4)
    #         else:
    #             frames = np.arange(0, 601, 4)
    #             time_at_frames1 = np.arange(0, 601, 4)
    #
    #         frames_folder1 = frames_folder_general + 'frames_voda/'
    #         frame_scale1 = []
    #         for time in time_at_frames1:
    #             time = time / 20.0
    #             frame_scale1.append(time)
    #
    #         if state != 'Pr':
    #             labelsx = np.arange(0, 26, 2)
    #             string_labels_x = []
    #             for i in np.arange(0, 26, 2):
    #                 if i not in labelsx:
    #                     string_labels_x.append('')
    #                 else:
    #                     string_labels_x.append(i)
    #         else:
    #             labelsx = np.arange(0, 31, 2)
    #             string_labels_x = []
    #             for i in np.arange(0, 31, 2):
    #                 if i not in labelsx:
    #                     string_labels_x.append('')
    #                 else:
    #                     string_labels_x.append(i)
    #
    #
    #         labelsy1 = np.arange(-10, 21, 2)
    #         string_labels_y1 = []
    #         for i in np.arange(-10, 21, 1):
    #             if i not in labelsy1:
    #                 string_labels_y1.append('')
    #             else:
    #                 string_labels_y1.append(i)
    #
    #         labelsy2 = np.arange(2, 9, 1)
    #         string_labels_y2 = []
    #         for i in np.arange(2, 9, 1):
    #             if i not in labelsy2:
    #                 string_labels_y2.append('')
    #             else:
    #                 string_labels_y2.append(i)
    #
    #         pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_xtra/prd/10ang_cavities_09_voda/md_%s_%s/' % (state, state, prot)
    #
    #         pkas = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames, 'frame')
    #
    #         name = 'md_%s_%s' % (state, prot)
    #
    #         sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder1, frames, 'PRD-3_EHEM', 'ARG-481_ACHA', name, out_folder_general)
    #
    #         # sb_dist = sb_dist[2:]
    #
    #         plot_utilities_frames.frames_plot_two_properties_pub('md_%s_%s' % (state, prot), 'pka', 'sb_dist', pkas, sb_dist, frame_scale1,
    #                                                              out_folder_general, x_ticks, x_lim, y1_lim, y1_ticks, y2_lim,
    #                                                              y2_ticks, string_labels_x, string_labels_y1, string_labels_y2)


    # for state in ['F', 'PF', 'Pm', 'Pr']:
    #     for prot in ['prah1', 'prah2', 'pra-']:

    for state in ['F']:
        for prot in ['prah2']:

            out_folder_general = '/user/jdragelj/projects/cco_pls/PRAa3/O-O_distance/%s/' % state
            if state != 'Pr':
                frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/%s_state/md_%s/' % (state,prot)
            else:
                frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/md_%s/' % (prot)

            # y1_lim = [-10, 20]
            # y1_ticks = np.arange(-10, 21, 1)
            y1_lim = [4, 20]
            y1_ticks = np.arange(4, 21, 1)
            y2_lim = [2, 8]
            y2_ticks = np.arange(2, 9, 1)

            if state != 'Pr':
                x_lim = [0, 25]
                x_ticks = np.arange(0, 26, 1)
            else:
                x_lim = [0, 31]
                x_ticks = np.arange(0, 31, 2)

            if state != 'Pr':
                frames = np.arange(0, 501, 4)
                time_at_frames1 = np.arange(0, 501, 4)
            else:
                frames = np.arange(0, 601, 4)
                time_at_frames1 = np.arange(0, 601, 4)

            frames_folder1 = frames_folder_general + 'frames_voda/'
            frame_scale1 = []
            for time in time_at_frames1:
                time = time / 20.0
                frame_scale1.append(time)

            if state != 'Pr':
                labelsx = np.arange(0, 26, 5)
                string_labels_x = []
                for i in np.arange(0, 26, 1):
                    if i not in labelsx:
                        string_labels_x.append('')
                    else:
                        string_labels_x.append(i)
            else:
                labelsx = np.arange(0, 31, 2)
                string_labels_x = []
                for i in np.arange(0, 31, 2):
                    if i not in labelsx:
                        string_labels_x.append('')
                    else:
                        string_labels_x.append(i)


            # labelsy1 = np.arange(-10, 21, 2)
            # string_labels_y1 = []
            # for i in np.arange(-10, 21, 1):
            #     if i not in labelsy1:
            #         string_labels_y1.append('')
            #     else:
            #         string_labels_y1.append(i)

            labelsy1 = np.arange(4, 21, 2)
            string_labels_y1 = []
            for i in np.arange(4, 21, 1):
                if i not in labelsy1:
                    string_labels_y1.append('')
                else:
                    string_labels_y1.append(i)

            labelsy2 = np.arange(2, 9, 1)
            string_labels_y2 = []
            for i in np.arange(2, 9, 1):
                if i not in labelsy2:
                    string_labels_y2.append('')
                else:
                    string_labels_y2.append(i)

            pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_pra_a3/10ang_cavities_09_voda/md_%s/' % (state, prot)

            pkas = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'PRA-1_EHEM', frames, 'frame')

            name = 'md_%s_%s' % (state, prot)

            sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder1, frames, 'PRA-1_EHEM', 'ASP-407_ACHA', name, out_folder_general)

            # sb_dist = sb_dist[2:]

            plot_utilities_frames.frames_plot_two_properties_pub('md_%s_%s' % (state, prot), 'pka', 'sb_dist', pkas, sb_dist, frame_scale1,
                                                                 out_folder_general, x_ticks, x_lim, y1_lim, y1_ticks, y2_lim,
                                                                 y2_ticks, string_labels_x, string_labels_y1, string_labels_y2)