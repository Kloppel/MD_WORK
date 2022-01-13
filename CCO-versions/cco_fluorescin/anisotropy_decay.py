# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import os
import cPickle as pickle


def get_fluorescein_coordinates(atom_name, file):
    pdb_file = open(file, 'r')

    for line in pdb_file:
        line = line.split()
        if 'FLUR' in line:
            if line[2] == atom_name:
                if line[4] != '1':
                    x = float(line[6])
                    y = float(line[7])
                    z = float(line[8])
                else:
                    x = float(line[5])
                    y = float(line[6])
                    z = float(line[7])
    coords = []
    coords.append(x)
    coords.append(y)
    coords.append(z)
    return coords


def get_green_vector_coords(pdb_file):
    atoms1 = ['C14', 'C13']
    atoms2 = ['C17', 'C18']

    x1m = (
          get_fluorescein_coordinates(atoms1[0], pdb_file)[0] + get_fluorescein_coordinates(atoms1[1], pdb_file)[0]) / 2
    y1m = (
          get_fluorescein_coordinates(atoms1[0], pdb_file)[1] + get_fluorescein_coordinates(atoms1[1], pdb_file)[1]) / 2
    z1m = (
          get_fluorescein_coordinates(atoms1[0], pdb_file)[2] + get_fluorescein_coordinates(atoms1[1], pdb_file)[2]) / 2

    x2m = (
          get_fluorescein_coordinates(atoms2[0], pdb_file)[0] + get_fluorescein_coordinates(atoms2[1], pdb_file)[0]) / 2
    y2m = (
          get_fluorescein_coordinates(atoms2[0], pdb_file)[1] + get_fluorescein_coordinates(atoms2[1], pdb_file)[1]) / 2
    z2m = (
          get_fluorescein_coordinates(atoms2[0], pdb_file)[2] + get_fluorescein_coordinates(atoms2[1], pdb_file)[2]) / 2

    coords_m1 = (x1m, y1m, z1m)
    coords_m2 = (x2m, y2m, z2m)

    return coords_m1, coords_m2



def normalized_vector(crd1, crd2):
    # This function first creates a vector of direction which is c2 - c1 (careful about interpretation later) and then
    # normalization happens -> becomes unit vector
    # this vector can not be used for analysis

    crd1 = np.array(crd1)
    crd2 = np.array(crd2)

    vector = crd2 - crd1

    length = np.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])
    u_x = vector[0] / length
    u_y = vector[1] / length
    u_z = vector[2] / length
    unit_vector = (u_x, u_y, u_z)
    unit_vector = np.array(unit_vector)

    return unit_vector


if __name__ == '__main__':


#     subfolders = {
#                   'hsp73_hsp526_red_fluh': [90], \
#                  'hsp73_hsp526_red_flu2': [90], \
#                  'hsd73_hse526_red_flu-': [90], \
#  \
#                  'hsp73_hsp526_fluh': [90], \
#                  'hsp73_hsp526_flu-_2': [90], \
#                  'hsp73_hse526_flu-': [90], \
# \
#                  'hsd73_hse526_flu-': [90], \
#                  }

    subfolders = {
                  'hsp73_hsp526_red_flu2': [90], \
                  'hsp73_hsp526_flu-_2': [90], \
                  'switch_O_to_R': [90],\
                  'switch_R_to_O': [90],\
                  }
    plt.figure()

    order = []
    for j, subfolder in enumerate(subfolders):

        main_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/'

        print subfolder
        if 'switch_O_to_R' in subfolder:
            appendix = 'red_flu-(pH=8.0)'
            order.append(appendix)
        elif 'switch_R_to_O' in subfolder:
            appendix = 'oxi_flu2(pH=6.6)'
            order.append(appendix)
        elif 'red_flu2' in subfolder:
            appendix = 'red_flu2(pH=6.6)'
            order.append(appendix)
        elif 'flu-_2' in subfolder:
            appendix = 'oxi_flu-(pH=8.0)'
            order.append(appendix)


        if ('hsd' in subfolder) and ('red' in subfolder):
            md_frames_folder = '%s/md_ex_17_ex_10/frames_1' % subfolder
        elif ('hsd' in subfolder) and ('red' not in subfolder):
            md_frames_folder = '%s/md_ex_17_ex_12/frames_1' % subfolder
        else:
            # md_frames_folder = '%s/md_ex_10/frames_1' % subfolder
            md_frames_folder = '%s/md_ex_10_ex_y/frames_1' % subfolder

        if 'switch' in subfolder:
            main_folder = '/scratch/scratch/jdragelj/'
            md_frames_folder = '%s/frames/' % subfolder

        coords_pickle = frame_file = main_folder + md_frames_folder + '/flu_coords_pickle.pkl'
        # if os.path.exists(coords_pickle):
        #     os.remove(coords_pickle)

        if not os.path.exists(coords_pickle):
            data_dict = {}
            if 'switch_O_to_R' in subfolder:
                frame_range = np.arange(0,1125,1)
            elif 'switch_R_to_O' in subfolder:
                frame_range = np.arange(0,997,1)
            else:
                frame_range = np.arange(0,1021,1)

            data_dict[subfolder] = {}
            for frame in frame_range:
                data_dict[subfolder][frame] = []
                frame_file = main_folder + md_frames_folder + '/frame%i.pdb' % frame
                crda, crdb = get_green_vector_coords(frame_file)
                data_dict[subfolder][frame].append(crda)
                data_dict[subfolder][frame].append(crdb)
            f = open(coords_pickle, 'wb')
            pickle.dump(data_dict, f)
            f.close()
        else:
            data_dict = pickle.load( open(coords_pickle, "rb" ))

        if 'switch' in subfolder:
            interval_end = 35.0
            interval_start = 5.0
        else:
            interval_end = 47.2
            interval_start = 17.2

        # if 'switch_O_to_R' in subfolder:
        #     interval_end = 55.0
        #     interval_start = 25.0

        dt = 0.1
        data_points = []
        data_points.append(1.0*2.0/5.0)
        # data_points.append(1.0)
        taus = np.arange(0.05,10.1,0.05)
        # for tau in [1.0, 2.0, 3.0, 4.0, 5.0]:
        for tau in taus:
            N = (30.0-tau)/dt
            print N
            sum = 0.0
            for n in np.arange(0,N,1.0):

                time1 = interval_start + n*dt
                time2 = interval_start + n*dt + tau

                frame1_num = round((time1*20), 0)
                frame2_num = round((time2*20), 0)


                frame1 = main_folder + md_frames_folder + '/frame%i.pdb' % frame1_num
                frame2 = main_folder + md_frames_folder + '/frame%i.pdb' % frame2_num


                crd1 = data_dict[subfolder][frame1_num][0]
                crd2 = data_dict[subfolder][frame1_num][1]
                unit1 = normalized_vector(crd1, crd2)

                crd3 = data_dict[subfolder][frame2_num][0]
                crd4 = data_dict[subfolder][frame2_num][1]
                unit2 = normalized_vector(crd3, crd4)

                dot_product = np.dot(unit1, unit2)
                # dot_transformed = dot_product
                dot_transformed = (3.0*dot_product*dot_product-1.0)/2.0

                sum += dot_transformed


            point_avg = (2.0/5.0)*(sum/(N+1.0))
            # point_avg = (2.0/5.0)*(sum/(N))
            # point_avg = (sum/(N+1.0))
            data_points.append(point_avg)

        # times = []
        # times.append(interval_start)
        # for tau in [1.0,2.0,3.0,4.0,5.0]:
        # # for tau in taus:
        #     times.append((interval_start+tau))


        times = [0]
        # for tau in [1.0,2.0,3.0,4.0,5.0]:
        for tau in taus:
            times.append(tau)
        print times



        if 'red' in subfolder:
            if 'fluh' in subfolder:
                linestyle = '--'
                color = 'black'
            elif 'flu-' in subfolder:
                linestyle = '-.'
                color = 'black'
            else:
                linestyle = '-'
                color = 'black'
        else:
            if 'hsd73_hse526_flu-' in subfolder:
                linestyle = '-.'
                color = 'red'
            elif 'flu-_2' in subfolder:
                linestyle = '-'
                color = 'red'
            elif 'fluh' in subfolder:
                linestyle = '--'
                color = 'red'
            else:
                linestyle = (0, (5, 1))
                color = 'red'

        if 'O_to_R' in subfolder:
            linestyle  = '--'
            color = 'black'
        if 'R_to_O' in subfolder:
            linestyle  = '--'
            color = 'red'
        if 'deprot' in subfolder:
            linestyle  = '-.'
            color = 'black'

        # plt.plot(times,data_points, linestyle=linestyle, lw=3.0, color=color)
        plt.plot(times,data_points, linestyle=linestyle, lw=1.5, color=color)
        plt.ylim(0.0,0.4)

        f = open('/user/jdragelj/Desktop/%s.dat' % appendix, 'w')
        f.write('time,r\n')
        for time, data in zip(times, data_points):
            string = str(time) + ',' + str(data) + '\n'
            f.write(string)
        f.close



        # # import ascii
        #
        # types = ['times', 'r']
        # times = np.array(times)
        # data_points = np.array(data_points)
        # table = {'type': types, 'time': times, 'r': data_points}
        # ascii.write(table, '/user/jdragelj/Desktpo/%s.dat' % subfolder, formats={'time': '%.1f', 'r': '%.5f'})



    # plt.legend(order, prop={'size': 15})
    plt.rcParams.update({'font.size': 15})
    plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/4_traj_10ns.png', dpi=300)
    plt.show()














