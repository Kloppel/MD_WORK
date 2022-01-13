# coding=utf-8

import kbp2
from workspace_jd import salt_bridge_distances
import re
import numpy as np
import cPickle as pickle
import os

def angle_between(v1, v2):
    # v1 is your firsr vector
    # v2 is your second vector
    angle = np.degrees(np.arccos(np.clip(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)),-1.0,1.0)))
    return angle

def find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d, atom_names_list_a,
                        hydrogen_atom, name, cutoffd = 4.00, cutoffa = 110.00):

    results_data = results_folder + '%s_%s_%s_hydrogen_bonds.dat' % (name, residue_donor, residue_acceptor)
    results_overview = results_folder + '%s_%s_%s_hydrogen_bonds_stats.dat' % (name, residue_donor, residue_acceptor)
    results_pickle = results_folder + '%s_%s_%s_hydrogen_bonds.pkl' % (name, residue_donor, residue_acceptor)

    if not os.path.exists(results_pickle):

        f = open(results_data, 'w')
        f_header = 'frame,distance,angle,hbond\n'
        f.write(f_header)
        g = open(results_overview, 'w')
        g_header = """
criteria: distance <= %.2f, angle >= %i   """ %(cutoffd, cutoffa)  +  """
donor: %s                                 """ %(residue_donor)  +  """
acceptor: %s                              """ %(residue_acceptor)  +  """
\n
"""
        g.write(g_header)
        distances_heavy = []
        angles_hbond = []
        bond_found = 0

        hbond_results = {}


        for frame in frame_range:
            print frame
            frame_pdb = frames_folder + 'frame%i.pdb' % frame

            pdb = kbp2.file_parser.Simple_struct_parser()
            pdb.read_pdb(frame_pdb)
            dict_2 = {}

            entries1 = re.split(r'[-_]', residue_donor)
            resname1, resid1, segname1 = entries1
            resid1 = int(resid1)

            entries2 = re.split(r'[-_]', residue_acceptor)
            resname2, resid2, segname2 = entries2
            resid2 = int(resid2)

            distance = 5000.0
            array_2_chosen = None

            hydrogen_coords = salt_bridge_distances.get_atom_coords(pdb, resname1, resid1, segname1, hydrogen_atom)
            donor_atom_coords = salt_bridge_distances.get_atom_coords(pdb, resname1, resid1, segname1, atom_name_d)

            for atom_name2 in atom_names_list_a:
                dict_2[atom_name2] = salt_bridge_distances.get_atom_coords(pdb, resname2, resid2, segname2, atom_name2)

            for atom2, array_2 in dict_2.iteritems():
                dist = np.linalg.norm(hydrogen_coords - array_2)
                if dist <= distance:
                    distance = dist
                    array_2_chosen = array_2

            dist_heavy = np.linalg.norm(donor_atom_coords - array_2_chosen)
            distances_heavy.append(dist_heavy)

            vector_donor = donor_atom_coords - hydrogen_coords
            vector_between = array_2_chosen - hydrogen_coords

            angle = angle_between(vector_donor, vector_between)
            angles_hbond.append(angle)

            hbond = 0
            if (dist_heavy <= cutoffd) and (angle >= cutoffa):
                bond_found += 1
                hbond = 1

            data_frame = '%i,%.3f,%.3f,%i\n' % (frame, distance, angle, hbond)
            hbond_results[frame] = (distance, angle, hbond)
            f.write(data_frame)

        occ = bond_found*100/len(frame_range)
        g_stat = '%.2f occupancy ' % (occ)
        g.write(g_stat)

        f = open(results_pickle, 'wb')
        pickle.dump(hbond_results, f)
        f.close()

    else:

        f = open(results_data, 'w')
        f_header = 'frame,distance,angle,hbond\n'
        f.write(f_header)
        g = open(results_overview, 'w')
        g_header = """
criteria: distance <= %.2f, angle >= %i   """ % (cutoffd, cutoffa) + """
donor: %s                                 """ % (residue_donor) + """
acceptor: %s                              """ % (residue_acceptor) + """
\n
        """
        g.write(g_header)
        hbond_results =  pickle.load( open( results_pickle, "rb" ) )
        bond_found = 0
        for frame in frame_range:
            bond_found +=hbond_results[frame][2]
        occ = bond_found * 100 / len(frame_range)
        g_stat = '%.2f occupancy ' % (occ)
        g.write(g_stat)

        print 'traj %s: donor %s acceptor %s %.2f occupancy ' % (name, residue_donor, residue_acceptor, occ)



if __name__ == '__main__':

    for state in ['F', 'PF', 'Pm', 'Pr']:

        print state


        results_folder = '/user/jdragelj/projects/cco_pls/PRD/hbonds/%s/' % (state)
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_prdh1_xtra6/frames_voda/' % (state, state)
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'PRD-3_EHEM'
        residue_acceptor = 'PRA-1_EHEM'
        atom_name_d = 'O1D'
        hydrogen_atom = 'H1D'
        atom_names_list_a = ['O1A','O2A']
        name = 'prdh1'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name = name, cutoffd=4.00, cutoffa=110)

        results_folder = '/user/jdragelj/projects/cco_pls/PRD/hbonds/%s/' % (state)
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_prdh1_xtra6/frames_voda/' % (state, state)
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'PRD-3_EHEM'
        residue_acceptor = 'PRD-3_GHEM'
        atom_name_d = 'O1D'
        hydrogen_atom = 'H1D'
        atom_names_list_a = ['O1D','O2D']
        name = 'prdh1'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name = name, cutoffd=4.00, cutoffa=110)

        results_folder = '/user/jdragelj/projects/cco_pls/PRD/hbonds/%s/' % (state)
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_prdh2_xtra6/frames_voda/' % (state, state)
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'PRD-3_EHEM'
        residue_acceptor = 'PRA-1_EHEM'
        atom_name_d = 'O1D'
        hydrogen_atom = 'H2D'
        atom_names_list_a = ['O1A','O2A']
        name = 'prdh2'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name = name, cutoffd=4.00, cutoffa=110)

        results_folder = '/user/jdragelj/projects/cco_pls/PRD/hbonds/%s/' % (state)
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_prdh2_xtra6/frames_voda/' % (state, state)
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'PRD-3_EHEM'
        residue_acceptor = 'PRD-3_GHEM'
        atom_name_d = 'O1D'
        hydrogen_atom = 'H2D'
        atom_names_list_a = ['O1D','O2D']
        name = 'prdh2'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name = name, cutoffd=4.00, cutoffa=110)


        results_folder = '/user/jdragelj/projects/cco_pls/PRD/hbonds/%s/' % (state)
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_prd-_xtra6/frames_voda/' % (state, state)
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'TRP-172_ACHA'
        residue_acceptor = 'PRD-3_EHEM'
        atom_name_d = 'NE1'
        hydrogen_atom = 'HE1'
        atom_names_list_a = ['O1D', 'O2D']
        name = 'prd-'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name = name, cutoffd=4.00, cutoffa=110)
        results_folder = '/user/jdragelj/projects/cco_pls/PRD/hbonds/%s/' % (state)
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_prdh1_xtra6/frames_voda/' % (state, state)
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'TRP-172_ACHA'
        residue_acceptor = 'PRD-3_EHEM'
        atom_name_d = 'NE1'
        hydrogen_atom = 'HE1'
        atom_names_list_a = ['O1D', 'O2D']
        name = 'prdh1'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name = name, cutoffd=4.00, cutoffa=110)
        results_folder = '/user/jdragelj/projects/cco_pls/PRD/hbonds/%s/' % (state)
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRD_a3/open_sb_%s/md_%s_prdh2_xtra6/frames_voda/' % (state, state)
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'TRP-172_ACHA'
        residue_acceptor = 'PRD-3_EHEM'
        atom_name_d = 'NE1'
        hydrogen_atom = 'HE1'
        atom_names_list_a = ['O1D', 'O2D']
        name = 'prdh2'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name = name, cutoffd=4.00, cutoffa=110)



    ##############################################################################################################################

    print "Pr"

    results_folder = '/user/jdragelj/projects/cco_pls/PRAa3/hbonds/Pr/'
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/md_prah1/frames_voda/'
    frame_range = np.arange(0, 501, 2)
    residue_donor = 'PRA-1_EHEM'
    residue_acceptor = 'PRD-3_EHEM'
    atom_name_d = 'O1A'
    hydrogen_atom = 'H1A'
    atom_names_list_a = ['O1D','O2D']
    name = 'prah1'
    find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                        atom_names_list_a, hydrogen_atom, name=name, cutoffd=4.00, cutoffa=110)

    results_folder = '/user/jdragelj/projects/cco_pls/PRAa3/hbonds/Pr/'
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/md_prah2/frames_voda/'
    frame_range = np.arange(0, 501, 2)
    residue_donor = 'PRA-1_EHEM'
    residue_acceptor = 'PRD-3_EHEM'
    atom_name_d = 'O2A'
    hydrogen_atom = 'H2A'
    atom_names_list_a = ['O1D','O2D']
    name = 'prah2'
    find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                        atom_names_list_a, hydrogen_atom, name=name, cutoffd=4.00, cutoffa=110)

    for state in ['F', 'PF', 'Pm']:

        print state
        # results_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/%s_state/md_prah1/' % state
        results_folder = '/user/jdragelj/projects/cco_pls/PRAa3/hbonds/%s/' % state
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/%s_state/md_prah1/frames_voda/' %state
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'PRA-1_EHEM'
        residue_acceptor = 'PRD-3_EHEM'
        atom_name_d = 'O1A'
        hydrogen_atom = 'H1A'
        atom_names_list_a = ['O1D','O2D']
        name = 'prah1'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name = name, cutoffd=4.00, cutoffa=110)

        # results_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/%s_state/md_prah2/' % state
        results_folder = '/user/jdragelj/projects/cco_pls/PRAa3/hbonds/%s/' % state
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/%s_state/md_prah2/frames_voda/' %state
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'PRA-1_EHEM'
        residue_acceptor = 'PRD-3_EHEM'
        atom_name_d = 'O2A'
        hydrogen_atom = 'H2A'
        atom_names_list_a = ['O1D','O2D']
        name = 'prah2'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name = name, cutoffd=4.00, cutoffa=110)

    #### hydrogen bonds PRA-ASP

    print 'Pr'
    results_folder = '/user/jdragelj/projects/cco_pls/PRAa3/hbonds/Pr/'
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/md_prah1/frames_voda/'
    frame_range = np.arange(0, 501, 2)
    residue_donor = 'PRA-1_EHEM'
    residue_acceptor = 'ASP-407_ACHA'
    atom_name_d = 'O1A'
    hydrogen_atom = 'H1A'
    atom_names_list_a = ['OD1','OD2']
    name = 'prah1'
    find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                        atom_names_list_a, hydrogen_atom, name=name, cutoffd=4.00, cutoffa=110)

    results_folder = '/user/jdragelj/projects/cco_pls/PRAa3/hbonds/Pr/'
    frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/md_prah2/frames_voda/'
    frame_range = np.arange(0, 501, 2)
    residue_donor = 'PRA-1_EHEM'
    residue_acceptor = 'ASP-407_ACHA'
    atom_name_d = 'O2A'
    hydrogen_atom = 'H2A'
    atom_names_list_a = ['OD1','OD2']
    name = 'prah2'
    find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                        atom_names_list_a, hydrogen_atom, name=name, cutoffd=4.00, cutoffa=110)

    for state in ['F', 'PF', 'Pm']:

        print state

        results_folder = '/user/jdragelj/projects/cco_pls/PRAa3/hbonds/%s/' % state
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/%s_state/md_prah1/frames_voda/' %state
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'PRA-1_EHEM'
        residue_acceptor = 'ASP-407_ACHA'
        atom_name_d = 'O1A'
        hydrogen_atom = 'H1A'
        atom_names_list_a = ['OD1', 'OD2']
        name = 'prah1'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name=name, cutoffd=4.00, cutoffa=110)

        results_folder = '/user/jdragelj/projects/cco_pls/PRAa3/hbonds/%s/' % state
        frames_folder = '/scratch/scratch/jdragelj/projects/cco_pls/prd_pra/open_sb_2GSM/PRA_a3/%s_state/md_prah2/frames_voda/' %state
        frame_range = np.arange(0, 501, 2)
        residue_donor = 'PRA-1_EHEM'
        residue_acceptor = 'ASP-407_ACHA'
        atom_name_d = 'O2A'
        hydrogen_atom = 'H2A'
        atom_names_list_a = ['OD1', 'OD2']
        name = 'prah2'
        find_hydrogen_bonds(results_folder, frames_folder, frame_range, residue_donor, residue_acceptor, atom_name_d,
                            atom_names_list_a, hydrogen_atom, name=name, cutoffd=4.00, cutoffa=110)


