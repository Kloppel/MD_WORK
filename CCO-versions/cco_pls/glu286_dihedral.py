# coding=utf-8

import kbp2
from workspace_jd import salt_bridge_distances
import re
import numpy as np
import cPickle as pickle
import os

def new_dihedral(p0, p1, p2, p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def measure_dihedral(md_name, results_folder, frames_folder, frame_range, residue, atoms):

    results_data = results_folder + '%s_dihedral.dat' % (md_name)
    results_pickle = results_folder + '%s_dihedral.pkl' % (md_name)

    if not os.path.exists(results_pickle):

        f = open(results_data, 'w')

        f_header = 'frame,angle,up/down\n'
        f.write(f_header)

        dihe_results = {}

        for frame in frame_range:
            print frame
            frame_pdb = frames_folder + 'frame%i.pdb' % frame

            pdb = kbp2.file_parser.Simple_struct_parser()
            pdb.read_pdb(frame_pdb)

            entries1 = re.split(r'[-_]', residue)
            resname1, resid1, segname1 = entries1
            resid1 = int(resid1)


            a1 = salt_bridge_distances.get_atom_coords(pdb, resname1, resid1, segname1, atoms[0])
            a2 = salt_bridge_distances.get_atom_coords(pdb, resname1, resid1, segname1, atoms[1])
            a3 = salt_bridge_distances.get_atom_coords(pdb, resname1, resid1, segname1, atoms[2])
            a4 = salt_bridge_distances.get_atom_coords(pdb, resname1, resid1, segname1, atoms[3])

            dihe_angle = new_dihedral(a1,a2,a3,a4)

            if abs(dihe_angle) > 90:
                down = 0
            if abs(dihe_angle) < 90:
                down = 1

            dihe_results[frame] = (abs(dihe_angle), down)
            f.write('%i,%.2f,%i\n' % (frame, abs(dihe_angle), down))

        f = open(results_pickle, 'wb')
        pickle.dump(dihe_results, f)
        f.close()

    else:

        dihe_results = pickle.load( open(results_pickle, "rb" ) )


    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    times=[]
    for frame in frame_range:
        time =  frame/20.0
        times.append(time)
    dihedrals  = []
    for frame in frame_range:
        dihedrals.append(dihe_results[frame][0])
    ax1.plot(times,dihedrals)
    plt.savefig(results_folder+'%s'%md_name)
    plt.close()
    # plt.show()

    return




if __name__ == '__main__':


    # for state in ['F', 'PF', 'Pm', 'Pr']:
    for state in ['Pr']:
        # for md in ['prd-', 'prdh1', 'prdh2']:
        for md in ['prd-']:
            results_folder = '/user/jdragelj/projects/cco_pls/PRD/glu286_dihe/%s/' % (state)
            frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_%s_xtra6/frames_voda/' % (state, state, md)
            frame_range = np.arange(0, 501, 2)
            residue = 'GLU-286_ACHA'
            atoms = ['CA','CB','CG','CD']
            md_name = 'md_%s_%s' % (state, md)
            measure_dihedral(md_name, results_folder, frames_folder, frame_range, residue, atoms)

    # for state in ['F', 'PF', 'Pm']:
    #     for md in ['md_pra-', 'md_prah1', 'md_prah2']:
    #         results_folder = '/user/jdragelj/projects/cco_pls/PRAa3/glu286_dihe/%s/' % (state)
    #         frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/%s_state/%s/frames_voda/' % (state, md)
    #         frame_range = np.arange(0, 501, 2)
    #         residue = 'GLU-286_ACHA'
    #         atoms = ['CA','CB','CG','CD']
    #         md_name = 'md_%s_%s' % (state, md)
    #         measure_dihedral(md_name, results_folder, frames_folder, frame_range, residue, atoms)
    #
    # for md in ['md_pra-', 'md_prah1', 'md_prah2']:
    #     results_folder = '/user/jdragelj/projects/cco_pls/PRAa3/glu286_dihe/Pr/'
    #     frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/%s/frames_voda/' % (md)
    #     frame_range = np.arange(0, 501, 2)
    #     residue = 'GLU-286_ACHA'
    #     atoms = ['CA', 'CB', 'CG', 'CD']
    #     md_name = 'md_Pr_%s' % (md)
    #     measure_dihedral(md_name, results_folder, frames_folder, frame_range, residue, atoms)

