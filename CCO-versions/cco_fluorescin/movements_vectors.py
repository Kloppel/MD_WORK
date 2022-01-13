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

    x1m = (get_fluorescein_coordinates(atoms1[0],pdb_file)[0] + get_fluorescein_coordinates(atoms1[1],pdb_file)[0])/2
    y1m = (get_fluorescein_coordinates(atoms1[0],pdb_file)[1] + get_fluorescein_coordinates(atoms1[1],pdb_file)[1])/2
    z1m = (get_fluorescein_coordinates(atoms1[0],pdb_file)[2] + get_fluorescein_coordinates(atoms1[1],pdb_file)[2])/2

    x2m = (get_fluorescein_coordinates(atoms2[0],pdb_file)[0] + get_fluorescein_coordinates(atoms2[1],pdb_file)[0])/2
    y2m = (get_fluorescein_coordinates(atoms2[0],pdb_file)[1] + get_fluorescein_coordinates(atoms2[1],pdb_file)[1])/2
    z2m = (get_fluorescein_coordinates(atoms2[0],pdb_file)[2] + get_fluorescein_coordinates(atoms2[1],pdb_file)[2])/2


    coords_m1 = (x1m,y1m,z1m)
    coords_m2 = (x2m,y2m,z2m)

    return coords_m1, coords_m2

# def get_magenta_vector_coords(pdb_file):
#
#     atoms = ['C7', 'C10']
#
#     x1 = get_fluorescein_coordinates(atoms[0],pdb_file)[0]
#     y1 = get_fluorescein_coordinates(atoms[0],pdb_file)[1]
#     z1 = get_fluorescein_coordinates(atoms[0],pdb_file)[2]
#
#     x2 = get_fluorescein_coordinates(atoms[1],pdb_file)[0]
#     y2 = get_fluorescein_coordinates(atoms[1],pdb_file)[1]
#     z2 = get_fluorescein_coordinates(atoms[1],pdb_file)[2]
#
#     coords_1 = (x1,y1,z1)
#     coords_2 = (x2,y2,z2)
#
#     return coords_1, coords_2

def get_magenta2_vector_coords(pdb_file):

    atoms_in_ring = ['C4','C5','C6','C7','C8','C9']
    x1, y1, z1 = get_centroid_coords(atoms_in_ring, pdb_file)

    atoms = ['C10']
    x2 = get_fluorescein_coordinates(atoms[0],pdb_file)[0]
    y2 = get_fluorescein_coordinates(atoms[0],pdb_file)[1]
    z2 = get_fluorescein_coordinates(atoms[0],pdb_file)[2]

    coords_1 = (x1,y1,z1)
    coords_2 = (x2,y2,z2)

    return coords_1, coords_2


def get_magenta3_vector_coords(pdb_file):

    atoms = ['C6', 'C8']

    x1 = get_fluorescein_coordinates(atoms[0], pdb_file)[0]
    y1 = get_fluorescein_coordinates(atoms[0], pdb_file)[1]
    z1 = get_fluorescein_coordinates(atoms[0], pdb_file)[2]

    x2 = get_fluorescein_coordinates(atoms[1], pdb_file)[0]
    y2 = get_fluorescein_coordinates(atoms[1], pdb_file)[1]
    z2 = get_fluorescein_coordinates(atoms[1], pdb_file)[2]

    coords_1 = (x1, y1, z1)
    coords_2 = (x2, y2, z2)

    return coords_1, coords_2

def get_centroid_coords(atoms_in_ring, pdb_file):

    x=0
    y=0
    z=0
    for i, atom in enumerate(atoms_in_ring):
        x += get_fluorescein_coordinates(atom,pdb_file)[0]
        y += get_fluorescein_coordinates(atom,pdb_file)[1]
        z += get_fluorescein_coordinates(atom,pdb_file)[2]

    cent_x = x / (i+1)
    cent_y = y / (i+1)
    cent_z = z / (i+1)

    return cent_x, cent_y, cent_z


def get_blue_vector_coords(pdb_file):
    
    atoms = ['C11', 'O4']

    x1 = get_fluorescein_coordinates(atoms[0], pdb_file)[0]
    y1 = get_fluorescein_coordinates(atoms[0], pdb_file)[1]
    z1 = get_fluorescein_coordinates(atoms[0], pdb_file)[2]

    x2 = get_fluorescein_coordinates(atoms[1], pdb_file)[0]
    y2 = get_fluorescein_coordinates(atoms[1], pdb_file)[1]
    z2 = get_fluorescein_coordinates(atoms[1], pdb_file)[2]

    coords_1 = (x1, y1, z1)
    coords_2 = (x2, y2, z2)

    return coords_1, coords_2

def normalized_vector(crd1, crd2):

    # This function first creates a vector of direction which is c2 - c1 (careful about interpretation later) and then
    # normalization happens -> becomes unit vector
    # this vector can not be used for analysis

    crd1 = np.array(crd1)
    crd2 = np.array(crd2)

    vector = crd2 - crd1

    length = np.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2])
    u_x = vector[0]/length
    u_y = vector[1]/length
    u_z = vector[2]/length
    unit_vector = (u_x,u_y,u_z)
    unit_vector = np.array(unit_vector)

    return unit_vector

if __name__ == '__main__':


    subfolders= {'hsp73_hsp526_red_fluh': [90], \
                 'hsp73_hsp526_red_flu2': [90], \
                 'hsp73_hsp526_fluh': [90], \
                 'hsp73_hsp526_flu-_2': [90], \
                 'hsd73_hse526_red_flu-': [90], \
                 'hsd73_hse526_flu-': [90], \
                 'hsp73_hse526_flu-': [90]}


    main_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/'

    for j, subfolder in enumerate(subfolders):

        print subfolder

        fig = plt.figure()
        ax = fig.add_subplot(111)

        plt.rcParams.update({'font.size': 10})
        ax.tick_params(direction='out', top='off', right='off', length=3.0)
        ax.spines['top'].set_linewidth(0.5)
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['right'].set_linewidth(0.5)
        ax.spines['bottom'].set_linewidth(0.5)
        plt.plot([2.0, 2.15])
        fig.set_size_inches(3.2, 2.7)


        plt.yticks(np.arange(-1.0, 1.1, 0.1))
        plt.ylim(-1.0, 1.0)

        time_scale = np.arange(4, 645, 1)
        time_at_frames = np.arange(0, 641, 1)
        md_frames_folder = '/md_ex_10/frames_1'
        if subfolder == 'hsd73_hse526_flu-':
            md_frames_folder = '/md_ex_17_ex_12/frames_1'
        if subfolder == 'hsd73_hse526_red_flu-':
            md_frames_folder = '/md_ex_17_ex_10/frames_1'
        plt.xticks(np.arange(0, 33, 1))
        plt.xlim(0, 32)

        # pdb_file_ref = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/four_orient_new/mds/90/%s%s/frame4.pdb' % (subfolder, md_frames_folder)
        pdb_file_ref = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/model/conformation_tests/mini_confs/cco_flu_90.pdb'


        labels = np.arange(0, 33, 4)
        empty_string_labels = []
        for i in np.arange(0, 33, 1):
            if i not in labels:
                empty_string_labels.append('')
            else:
                empty_string_labels.append(i)

        labelsy = np.arange(-1.0, 1.1, 0.2)
        empty_string_labelsy = []
        for i in np.arange(-1.0, 1.1, 0.1):
            if i not in labelsy:
                empty_string_labelsy.append('')
            else:
                if abs(round(i,1))==0.0:
                    to_append = 0.0
                else:
                    to_append = round(i,1)
                empty_string_labelsy.append(to_append)

        frame_scale = []
        for time in time_at_frames:
            time = time/20.0
            # time = time
            frame_scale.append(time)

        folder = main_folder + subfolder + md_frames_folder + '/'

        x = frame_scale

        pickle_file = '/user/jdragelj/projects/cco_fluorescein/movement_vectors/final_choice/' + subfolder + '.pkl'
        loaded = False

        if os.path.exists(pickle_file):
            data_dict = pickle.load(open(pickle_file, 'rb'))
            if (data_dict['info']['time_scale'].all()== time_scale.all()):
                if (data_dict['info']['ref']== pdb_file_ref):
                    y1 = data_dict['green']
                    y2 = data_dict['blue']
                    y3 = data_dict['magenta']
                    loaded = True

        if not loaded:
            data_dict = {}
            data_dict['info'] = {}
            data_dict['info']['time_scale'] = time_scale
            data_dict['info']['ref'] = pdb_file_ref
            data_dict['green'] = []
            data_dict['blue'] = []
            data_dict['magenta'] = []

            ##### GREEN #####
            dot_products = []
            #frame0
            crd1, crd2 = get_green_vector_coords(pdb_file_ref)
            unit0 = normalized_vector(crd1, crd2)

            for frame in time_scale:
                pdb_file = folder + "frame%i.pdb" % frame
                crd1, crd2 = get_green_vector_coords(pdb_file)
                unit_t = normalized_vector(crd1, crd2)
                dot_product = np.dot(unit0, unit_t)
                dot_products.append(dot_product)
            y1 = np.array(dot_products)
            # plt.plot(x, y1, color='green')
            data_dict['green'] = y1
            ###################

            ##### BLUE #####
            dot_products = []
            #frame0
            crd1, crd2 = get_blue_vector_coords(pdb_file_ref)
            unit0 = normalized_vector(crd1, crd2)

            for frame in time_scale:
                pdb_file = folder + "frame%i.pdb" % frame
                crd1, crd2 = get_blue_vector_coords(pdb_file)
                unit_t = normalized_vector(crd1, crd2)
                dot_product = np.dot(unit0, unit_t)
                dot_products.append(dot_product)
            y2 = np.array(dot_products)
            # plt.plot(x, y2, color='blue')
            data_dict['blue'] = y2
            ###################

            ##### MAGENTA2 #####
            dot_products = []
            #frame0
            crd1, crd2 = get_magenta2_vector_coords(pdb_file_ref)
            unit0 = normalized_vector(crd1, crd2)
            angles = []

            for frame in time_scale:
                pdb_file = folder + "frame%i.pdb" % frame
                crd1, crd2 = get_magenta2_vector_coords(pdb_file)
                unit_t = normalized_vector(crd1, crd2)
                dot_product = np.dot(unit0, unit_t)
                dot_products.append(dot_product)
                # angle = np.degrees(np.arccos(np.clip(np.dot(unit0, unit_t), -1.0, 1.0)))
                # angles.append(angle)
            y3 = np.array(dot_products)
            # plt.plot(x, y3, color='magenta')
            data_dict['magenta'] = y3
            pickle.dump(data_dict, open(pickle_file, 'wb'))


        ax.plot(x, y1, color='green', lw=0.4)
        ax.plot(x, y2, color='blue', lw=0.4)
        ax.plot(x, y3, color='magenta', lw=0.4)

        ax.set_xticklabels(empty_string_labels)
        ax.set_yticklabels(empty_string_labelsy)


        # angle_mean = np.mean(angles)
        # print angle_mean, 'mean'
        # ref = 0
        # for angle in angles:
        #     x = angle-angle_mean
        #     if abs(x) > ref:
        #         ref = abs(x)
        # print ref, 'amplitude'
        ###################

        # fig.savefig('/user/jdragelj/projects/cco_fluorescein/movement_vectors/final_choice/%s_over_time_ref4' % (subfolder))
        fig.savefig('/user/jdragelj/projects/cco_fluorescein/movement_vectors/final_choice/%s_over_time_startref.png' % (subfolder), dpi=300)


















