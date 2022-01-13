# coding=utf-8



import os
from workspace_jd import tools
import numpy as np



def get_protonation_number(file, residue_tuple_string):

    state_number = None

    f = open(file, 'r')
    for line in f:
        if residue_tuple_string+':' in line:
            components = line.split(':')
            for i, component in enumerate(components):
                if residue_tuple_string in component:
                    data = components[i+1].split(',')
                    state_number = data[0]

    return int(state_number)

def get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict):

    csv = open(csv_filepath, 'w')
    csv.write(tools.string_creation([' '] + list(frame_range)))
    csv.write('\n')

    states_percentages = {}
    for state in residue_dict.keys():
        states_percentages[state] = 0


    list_of_states = []
    for frame in frame_range:
        file_folder = titration_folder + 'frame%i/' % frame
        subdirs = os.listdir(file_folder)
        for sub in subdirs:
            if filename in sub:
                filen = sub
        state_num = get_protonation_number(file_folder+filen, residue_tuple_string)
        states_percentages[int(state_num)] +=1
        list_of_states.append(residue_dict[state_num])

    print residue_tuple_string
    for state, sum in states_percentages.iteritems():
        print residue_dict[state], 100*states_percentages[state]/(len(frame_range)-1)


    csv.write(tools.string_creation([residue_tuple_string] + list_of_states))
    csv.write('\n')



if __name__ == '__main__':

    # titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/hsd73_hse526_flu-/done/'
    # frame_range = np.arange(384,685,2)
    # filename = 'h_min_pac.log'
    # csv_filepath = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/titration_plots/90_hsd73_hse526_flu-.csv'
    # residue_tuple_string = """('HSP', 526, 'ACHA')"""
    # residue_dict = {0:'Z', 1:'D', 2:'P', 3:'E'}
    # get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)
    #
    # titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/hsd73_hse526_red_flu-/done/'
    # frame_range = np.arange(144,445,2)
    # filename = 'h_min_pac.log'
    # csv_filepath = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/titration_plots/90_hsd73_hse526_red_flu-.csv'
    # residue_tuple_string = """('HSP', 526, 'ACHA')"""
    # residue_dict = {0:'Z', 1:'D', 2:'P', 3:'E'}
    # get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)
    #
    # titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/270/hsd73_hsp526_flu2/done/'
    # frame_range = np.arange(144,445,2)
    # filename = 'h_min_pac.log'
    # csv_filepath = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/titration_plots/270_hsd73_hsp526_flu2.csv'
    # residue_tuple_string = """('HSP', 526, 'ACHA')"""
    # residue_dict = {0:'Z', 1:'D', 2:'P', 3:'E'}
    # get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)
    #
    # titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/many_residues/90/hsd73_hse526_flu-/done/'
    # frame_range = np.arange(384,685,2)
    # filename = 'h_min_pac.log'
    # csv_filepath = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/titration_plots/90_hsd73_hse526_flu-_epp.csv'
    # residue_tuple_string = """('EPP', 525, 'ACHA')"""
    # residue_dict = {0:'Z', 1:'P1', 2:'DP', 3:'P2'}
    # get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)
    #
    # titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/many_residues/90/hsd73_hse526_flu-/done/'
    # frame_range = np.arange(144,445,2)
    # filename = 'h_min_pac.log'
    # csv_filepath = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/data_collective/titration_plots/90_hsd73_hse526_red_flu-_epp.csv'
    # residue_tuple_string = """('EPP', 525, 'ACHA')"""
    # residue_dict = {0:'Z', 1:'P1', 2:'DP', 3:'P2'}
    # get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)
    #
    # titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/hsp73_hsp526_red/done/'
    # frame_range = np.arange(4,205,2)
    # filename = 'h_min_pac.log'
    # csv_filepath = '/user/jdragelj/Desktop/526.csv'
    # residue_tuple_string = """('HSP', 526, 'ACHA')"""
    # residue_dict = {0:'Z', 1:'D', 2:'P', 3:'E'}
    # get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)

    # titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/hsp73_hse526_flu-/done/'
    # frame_range = np.arange(4,405,2)
    # filename = 'h_min_pac.log'
    # csv_filepath = '/user/jdragelj/Desktop/526oxi.csv'
    # residue_tuple_string = """('HSP', 526, 'ACHA')"""
    # residue_dict = {0:'Z', 1:'D', 2:'P', 3:'E'}
    # residue_dict = {0:'error', 1:'0', 2:'1', 3:'0'}
    # get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)
    #
    # titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/hsp73_hse526_red_flu-/done/'
    # frame_range = np.arange(4,405,2)
    # filename = 'h_min_pac.log'
    # csv_filepath = '/user/jdragelj/Desktop/526red.csv'
    # residue_tuple_string = """('HSP', 526, 'ACHA')"""
    # residue_dict = {0:'Z', 1:'D', 2:'P', 3:'E'}
    # residue_dict = {0:'error', 1:'0', 2:'1', 3:'0'}
    # get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)

    titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/hsp73_hsp526_flu-/done/'
    frame_range = np.arange(4,405,2)
    filename = 'h_min_pac.log'
    csv_filepath = '/user/jdragelj/Desktop/526fluoxi.csv'
    residue_tuple_string = """('HSP', 526, 'ACHA')"""
    residue_dict = {0:'Z', 1:'D', 2:'P', 3:'E'}
    residue_dict = {0:'error', 1:'0', 2:'1', 3:'0'}
    get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)

    titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/titrations_with_flu/90/hsp73_hsp526_oxi_protonate_fluh/done/'
    frame_range = np.arange(4,405,2)
    filename = 'h_min_pac.log'
    csv_filepath = '/user/jdragelj/Desktop/526flured.csv'
    residue_tuple_string = """('HSP', 526, 'ACHA')"""
    residue_dict = {0:'Z', 1:'D', 2:'P', 3:'E'}
    residue_dict = {0:'error', 1:'0', 2:'1', 3:'0'}
    get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)

    titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/hsp73_hsp526_red/done/'
    frame_range = np.arange(4,205,2)
    filename = 'h_min_pac.log'
    csv_filepath = '/user/jdragelj/Desktop/wt526red.csv'
    residue_tuple_string = """('HSP', 526, 'ACHA')"""
    residue_dict = {0:'Z', 1:'D', 2:'P', 3:'E'}
    residue_dict = {0:'error', 1:'0', 2:'1', 3:'0'}
    get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)

    titration_folder = '/scratch/scratch/jdragelj/projects/cco_alexiev/MD_299/titrations/wild_type/hsp73_hsp526/done/'
    frame_range = np.arange(4,305,2)
    filename = 'h_min_pac.log'
    csv_filepath = '/user/jdragelj/Desktop/wt526oxi.csv'
    residue_tuple_string = """('HSP', 526, 'ACHA')"""
    residue_dict = {0:'Z', 1:'D', 2:'P', 3:'E'}
    residue_dict = {0:'error', 1:'0', 2:'1', 3:'0'}
    get_states_md(titration_folder, frame_range, filename, residue_tuple_string, csv_filepath, residue_dict)
