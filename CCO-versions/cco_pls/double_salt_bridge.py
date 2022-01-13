from workspace_jd import salt_bridge_distances
from kbp2.workspace_jd.cco_fluorescin import plot_utilities_frames
import numpy as np

if __name__ == '__main__':


    # for state in ['Pm', 'PF', 'Pr', 'F']:
    for state in ['PF']:
    # for state in ['Pr']:
    #     for prot in ['prdh1', 'prdh2', 'prd-']:
    #     for prot in ['prdh1', 'prdh2']:
        for prot in ['prd-']:

            # out_folder_general = '/user/jdragelj/projects/cco_pls/PRD/salt-bridges_double/%s/' % state
            out_folder_general = '/user/jdragelj/Desktop/'
            frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_%s_xtra6/' % (state, state, prot)
            # frames_folder_general = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/crystal_mds/md_PF_crystal/'
            # y1_lim = [-10, 20]
            # y1_ticks = np.arange(-10, 21, 1)
            y1_lim = [-10, 8]
            y1_ticks = np.arange(-10, 9, 1)
            y2_lim = [2, 8]
            y2_ticks = np.arange(2, 9, 1)

            # if state != 'Pr':
            if state != 'X':
                x_lim = [0, 25]
                x_ticks = np.arange(0, 26, 1)
            else:
                x_lim = [0, 30]
                x_ticks = np.arange(0, 31, 1)

            if state != 'Pr':
                frames = np.arange(0, 501, 10)
                time_at_frames1 = np.arange(0, 501, 10)
            else:
                frames = np.arange(0, 601, 10)
                time_at_frames1 = np.arange(0, 601, 10)

            frames_folder1 = frames_folder_general + 'frames_voda/'
            # frames_folder1 = frames_folder_general + 'frames_voda_1/'
            frame_scale1 = []
            for time in time_at_frames1:
                time = time / 20.0
                frame_scale1.append(time)

            # if state != 'Pr':
            if state != 'X':
                labelsx = np.arange(0, 26, 5)
                string_labels_x = []
                for i in np.arange(0, 26, 1):
                    if i not in labelsx:
                        string_labels_x.append('')
                    else:
                        string_labels_x.append(i)
            else:
                # labelsx = np.arange(0, 31, 2)
                # string_labels_x = []
                # for i in np.arange(0, 31, 2):
                #     if i not in labelsx:
                #         string_labels_x.append('')
                #     else:
                #         string_labels_x.append(i)

                labelsx = np.arange(0, 31, 5)
                string_labels_x = []
                for i in np.arange(0, 31, 1):
                    if i not in labelsx:
                        string_labels_x.append('')
                    else:
                        string_labels_x.append(i)


            labelsy1 = np.arange(-10, 9, 2)
            string_labels_y1 = []
            for i in np.arange(-10, 9, 1):
                if i not in labelsy1:
                    string_labels_y1.append('')
                else:
                    string_labels_y1.append(i)

            # labelsy1 = np.arange(-10, 21, 2)
            # string_labels_y1 = []
            # for i in np.arange(-10, 21, 1):
            #     if i not in labelsy1:
            #         string_labels_y1.append('')
            #     else:
            #         string_labels_y1.append(i)

            labelsy2 = np.arange(2, 9, 1)
            string_labels_y2 = []
            for i in np.arange(2, 9, 1):
                if i not in labelsy2:
                    string_labels_y2.append('')
                else:
                    string_labels_y2.append(i)

            # prot = 'prd-_2'
            # prot = 'prd-'
            pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/%s_xtra/prd/10ang_cavities_09_voda/md_%s_%s/' % (state, state, prot)
            # pka_results_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/titrate/PF_md_crystal/10ang_cavities_09_voda/'

            pkas = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEM', frames, 'frame')
            # pkas = plot_utilities_frames.get_pkas_trajectory_residue(pka_results_folder, 'PRD-3_EHEMEA3', frames, 'frame')

            name = 'md_%s_%s' % (state, prot)

            sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder1, frames, 'PRD-3_EHEM', 'ARG-481_ACHA', name, out_folder_general, double='PRD-3_EHEM')
            # sb_dist = salt_bridge_distances.get_distances_trajectory(frames_folder1, frames, 'PRD-3_EHEMEA3', 'ARG-481_ACHAIN', name, out_folder_general, double='PRD-3_EHEMEA3')

            sb_dist_min = []
            sb_dist_other = []
            for dist_tuple in sb_dist:
                sb_dist_min.append(dist_tuple[0])
                sb_dist_other.append(dist_tuple[1])

            csv_file = '/user/jdragelj/Desktop/data.dat'
            f = open(csv_file, 'w')
            # f.write('time,d1,d2 \n')
            for  time, pka, d1, d2 in zip (frame_scale1,pkas,sb_dist_min, sb_dist_other):
                # print time, pka, d1, d2
                f.write('%s,%s,%s \n' % (str(time), str(d1), str(d2)))

            plot_utilities_frames.frames_plot_three_properties_pub('md_%s_%s' % (state, prot), 'pka', 'sb_dist', pkas, sb_dist_min, sb_dist_other, frame_scale1,
                                                                 out_folder_general, x_ticks, x_lim, y1_lim, y1_ticks, y2_lim,
                                                                 y2_ticks, string_labels_x, string_labels_y1, string_labels_y2)
