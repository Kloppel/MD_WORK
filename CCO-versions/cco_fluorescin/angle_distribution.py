# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import os
import cPickle as pickle


def get_vectors(coord_file):

    coordinates_file = open(coord_file, 'r')

    c14_coords = []
    c17_coords = []

    vectors = []
    avg_coords_c14 = [0.0,0.0,0.0]
    avg_coords_c17 = [0.0,0.0,0.0]

    for line in coordinates_file:
        line = line.split()
        if 'C14' in line:
        # if 'O4' in line:
            x = float(line[1])
            y = float(line[2])
            z = float(line[3])
            coords = (x,y,z)
            c14_coords.append(coords)
        if 'C17' in line:
        # if 'SG' in line:
            x = float(line[1])
            y = float(line[2])
            z = float(line[3])

            coords = (x,y,z)
            c17_coords.append(coords)

    for c14_a, c17_a in zip(c14_coords, c17_coords):
        x = c17_a[0] - c14_a[0]
        y = c17_a[1] - c14_a[1]
        z = c17_a[2] - c14_a[2]
        vector = np.array((x,y,z))
        vectors.append(vector)

        avg_coords_c14[0] += c14_a[0]
        avg_coords_c14[1] += c14_a[1]
        avg_coords_c14[2] += c14_a[2]

        avg_coords_c17[0] += c17_a[0]
        avg_coords_c17[1] += c17_a[1]
        avg_coords_c17[2] += c17_a[2]

    avg_coords_c14[0] = avg_coords_c14[0]/len(c14_coords)
    avg_coords_c14[1] = avg_coords_c14[1]/len(c14_coords)
    avg_coords_c14[2] = avg_coords_c14[2]/len(c14_coords)

    avg_coords_c17[0] = avg_coords_c17[0]/len(c17_coords)
    avg_coords_c17[1] = avg_coords_c17[1]/len(c17_coords)
    avg_coords_c17[2] = avg_coords_c17[2]/len(c17_coords)

    c14_avg = np.array(avg_coords_c14)
    c17_avg = np.array(avg_coords_c17)

    avg = c17_avg-c14_avg
    avg_vector = avg
    # avg_vector = np.array(avg)
    # print type(avg)

    # print avg_vector, vectors

    return avg_vector, vectors

def normalized_vector(vector):
    # This function first creates a vector of direction which is c2 - c1 (careful about interpretation later) and then
    # normalization happens -> becomes unit vector
    # this vector can not be used for analysis

    length = np.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])
    u_x = vector[0] / length
    u_y = vector[1] / length
    u_z = vector[2] / length
    unit_vector = (u_x, u_y, u_z)
    unit_vector = np.array(unit_vector)

    return unit_vector

def get_angles(vector_list, frame_range):


    angles_from_avg = []

    avg_vector = np.array([0.0,0.0,0.0])
    for vector in vector_list:
        avg_vector += vector
    avg_vector = avg_vector/len(vector_list)
    unit_avg = normalized_vector(avg_vector)

    for frame in frame_range:
        vector = vector_list[frame]
        unit_v = normalized_vector(vector)
        angle = np.degrees(np.arccos(np.clip(np.dot(unit_v, unit_avg))))
        # angle = np.degrees(np.arccos(np.dot(vector, avg_vector) / (np.linalg.norm(vector) * np.linalg.norm(avg_vector))))
        angles_from_avg.append(angle)
        if angle > 180:
            print angle

    return angles_from_avg


if __name__ == '__main__':



    frame_range = np.arange(344,945,1)



    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/O_flu-_blue_protein.dat')
    # angles = get_angles(vector_list, frame_range)
    # width = 0.0
    # for angle in angles:
    #     width += angle*angle
    # width = width/len(angles)
    # figK = plt.figure()
    # bins=np.arange(0.0,50.1,2.5)
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,50,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_align_protein_blue.png', dpi=300)
    # plt.title('oxi align protein, second moment %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_align_protein_labels_blue.png')
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_oxi_align_protein_blue.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_oxi_align_protein_blue.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()
    #
    #
    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/O_flu-_blue_loop.dat')
    # angles = get_angles(vector_list, frame_range)
    # width = 0.0
    # for angle in angles:
    #     width += angle * angle
    # width = width / len(angles)
    # figB = plt.figure()
    # bins=np.arange(0.0,50.1,2.5)
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,50,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_align_loop_blue.png', dpi=300)
    # plt.title('oxi align loop chain A 299-303, second moment %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_align_blue_loop_labels.png')
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_oxi_align_loop_blue.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_oxi_align_loop_blue.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()





    # frame_range = np.arange(344,945,1)
    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/O_flu-.dat')
    # angles = get_angles(vector_list, frame_range)
    # width = 0.0
    # for angle in angles:
    #     width += angle*angle
    # width = width/len(angles)
    # figA = plt.figure()
    # bins=[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,180,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_align_protein.png', dpi=300)
    # plt.title('oxi align protein, second moment %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_align_protein_labels.png')
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_oxi_align_protein.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_oxi_align_protein.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()
    #
    #
    # # plt.show()
    # # plt.close()
    #
    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/O_flu-_noalign.dat')
    # angles = get_angles(vector_list, frame_range)
    # width = 0.0
    # for angle in angles:
    #     width += angle * angle
    # width = width / len(angles)
    # figA = plt.figure()
    # bins=[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,180,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_no_align.png', dpi=300)
    # plt.title('oxi no align, %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_no_align_labels.png')
    # # plt.show()
    # # plt.close()
    #
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_oxi_no_align.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_oxi_no_align.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()
    #
    #
    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/O_flu-_helix.dat')
    # angles = get_angles(vector_list, frame_range)
    # width = 0.0
    # for angle in angles:
    #     width += angle * angle
    # width = width / len(angles)
    # figB = plt.figure()
    # bins=[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,180,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_align_helix.png', dpi=300)
    # plt.title('oxi align helix chain A 280-299, %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_align_helix_labels.png')
    # # plt.show()
    # # plt.close()
    #
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_oxi_align_helix.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_oxi_align_helix.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()
    #
    avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/O_flu-_loop_correct.dat')
    angles = get_angles(vector_list, frame_range)
    width = 0.0
    for angle in angles:
        width += angle * angle
    width = width / len(angles)

    figB = plt.figure()
    bins=[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
    (n, bins, patches) = plt.hist(angles, bins=bins)
    print n
    print bins
    print patches
    y_limit = max(n) + 0.5*max(n)
    plt.axis([0,180,0,y_limit])
    plt.grid(True)
    plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_align_loop.png', dpi=300)

    plt.title('oxi align loop 299-303, %.2f' % width)
    plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/oxi_align_loop_labels.png')
    f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_oxi_align_loop.dat','w')
    f.write('bin_edge_right,no_angles\n')
    for b, v in zip(bins[1:],n):
        f.write('%i;%i\n'%(int(b),int(v)))
    f.close()
    g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_oxi_align_loop.dat','w')
    g.write('time(ns),angle(deg)\n')
    for frame, ang in enumerate(angles):
        g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    g.close()


    frame_range = np.arange(100,701,1)

    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/R_flu-_blue_protein.dat')
    # angles = get_angles(vector_list, frame_range)
    # width = 0.0
    # for angle in angles:
    #     width += angle*angle
    # width = width/len(angles)
    # figA = plt.figure()
    # bins=np.arange(0.0,50.1,2.5)
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,50,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_protein_blue.png', dpi=300)
    # plt.title('red align protein, second moment second moment %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_protein_labels_blue.png')
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_red_align_protein_blue.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_red_align_protein_blue.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()
    #
    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/R_flu-_blue_loop.dat')
    # angles = get_angles(vector_list, frame_range)
    # width = 0.0
    # for angle in angles:
    #     width += angle*angle
    # width = width/len(angles)
    # figC = plt.figure()
    # bins=np.arange(0.0,50.1,2.5)
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,50,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_loop_blue.png', dpi=300)
    # plt.title('red align loop A 299-303, second moment %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_loop_labels_blue.png')
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_red_align_loop_blue.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_red_align_loop_blue.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()

    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/R_flu-.dat')
    # angles = get_angles(vector_list, frame_range)
    # width = 0.0
    # for angle in angles:
    #     width += angle * angle
    # width = width / len(angles)
    # figC = plt.figure()
    # bins=[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,180,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_protein.png', dpi=300)
    # plt.title('red align protein, %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_protein_labels.png')
    # # plt.show()
    # # plt.close()
    #
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_red_align_protein.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_red_align_protein.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()
    #
    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/R_flu-_noalign.dat')
    # angles = get_angles(vector_list, frame_range)
    # width = 0.0
    # for angle in angles:
    #     width += angle * angle
    # width = width / len(angles)
    # figC = plt.figure()
    # bins=[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,180,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_no_align.png', dpi=300)
    # plt.title('red no align, %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_no_align_labels.png')
    # # plt.show()
    # # plt.close()
    #
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_red_no_align.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_red_no_align.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()
    #
    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/R_flu-_helix.dat')
    # angles = get_angles(vector_list, frame_range)
    # figC = plt.figure()
    # bins=[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # width = 0.0
    # for angle in angles:
    #     width += angle * angle
    # width = width / len(angles)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,180,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_helix.png', dpi=300)
    # plt.title('red align helix, %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_chain_labels.png')
    # # plt.show()
    # # plt.close()
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_red_align_helix.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_red_align_helix.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()

    # avg_vector, vector_list = get_vectors('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/R_flu-_loop_correct.dat')
    # angles = get_angles(vector_list, frame_range)
    # figC = plt.figure()
    # bins=[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
    # (n, bins, patches) = plt.hist(angles, bins=bins)
    # y_limit = max(n) + 0.5*max(n)
    # plt.axis([0,180,0,y_limit])
    # plt.grid(True)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_loop.png', dpi=300)
    # plt.title('red align loop 299-303, %.2f' % width)
    # plt.savefig('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/red_align_loop_labels.png')
    # f = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/hist_data_red_align_loop.dat','w')
    # f.write('bin_edge_right,no_angles\n')
    # for b, v in zip(bins[1:],n):
    #     f.write('%i;%i\n'%(int(b),int(v)))
    # f.close()
    # g = open('/user/jdragelj/projects/cco_fluorescein/anisotropy/angle_distribution/angles_red_align_loop.dat','w')
    # g.write('time(ns),angle(deg)\n')
    # for frame, ang in enumerate(angles):
    #     g.write('%.2f;%.4f\n'%(float(frame)/20.0,ang))
    # g.close()
    #
    # plt.show()
