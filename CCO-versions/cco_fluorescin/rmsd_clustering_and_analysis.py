# coding=utf-8

import os
import cPickle as pickle

from kbp2.workspace_jd import R_clustering
from kbp2.workspace_jd.cco_fluorescin import mutual_rmsds


frames = []
for i in range(20,101):
    i = str(i)
    frames.append(i)


###################
# Regular frames #
##################



vmd_data_folder = '/user/jdragelj/projects/cco_fluorescin/Clustering/vmd_data/'
distance_matrix_data_folder = '/user/jdragelj/projects/cco_fluorescin/Clustering/distance_matrix_data/'

####################################
### creating data for clustering ###
####################################

# subfolders = os.listdir(vmd_data_folder)
# # prepare object for usage
# mut_rmsds = mutual_rmsds.Rmsd()
#
# for subfolder in subfolders:
#
#     if subfolder[-8:] != '_all.dat':
#             continue
#     else:
#         check = True
#     if check:
#
#         mut_rmsds.extract_rmsds(subfolder, vmd_data_folder+subfolder, selection=frames)
#         R_clustering.distance_matrix(mut_rmsds.rmsds[subfolder], frames, distance_matrix_data_folder+subfolder+'.csv')
# mut_rmsds.pickle_results('/user/jdragelj/projects/cco_fluorescin/Clustering', 'flu_rmsds')
####################################


###############################
### clustering preparations ###
###############################


workfolder = '/user/jdragelj/projects/cco_fluorescin/Clustering/all_results/results_20_NOC/'
# workfolder = '/user/jdragelj/projects/cco_fluorescin/Clustering/results_20_cutoff/'

# files_to_examine = os.listdir(distance_matrix_data_folder)
# files_to_examine_filter = []
# for file in files_to_examine:
#     finger_print = re.split('[_]', file)[0]
#     if finger_print not in files_to_examine_filter:
#         files_to_examine_filter.append(finger_print)
#
# # data_set is short MD name example: "pyrh"
#
# for data_set in files_to_examine_filter:
#
#     if not os.path.exists(workfolder + data_set + '/'):
#         os.mkdir(workfolder + data_set + '/')
#         shutil.copy2(distance_matrix_data_folder + data_set + '_all.dat.csv', workfolder + data_set + '/'  + data_set + '_all.dat.csv')
#
#     folderpath = workfolder + data_set + '/'
#
#     csv_file = folderpath + data_set + '_all.dat.csv'
#
#     # R_clustering.make_script(folderpath, csv_file, cutoff_range=(0,20,0.5))
#     R_clustering.make_script(folderpath, csv_file, no_of_clusters_range=(2,80,1))
#
#     run_script = 'Rscript %scluster.R %s \n' % (folderpath, csv_file)
#     print 'set JOBNAME%s = %s' % (data_set, data_set)
#     print 'cd %s' % folderpath
#     print run_script

##########################################
### clustering analysis - inter/intra  ###
##########################################

pickle_file = '/user/jdragelj/projects/cco_fluorescin/Clustering/flu_rmsds.pkl'
mut_rmsds = mutual_rmsds.Rmsd()
f = open( pickle_file, "rb" )
mut_rmsds = pickle.load(f)
f.close()

subfolders = os.listdir(workfolder)
for subfolder in subfolders:
    if 'script.sh' in subfolder:
        continue
    if 'inter_intra' in subfolder:
        continue
    R_clustering.inter_intra_distances(workfolder + subfolder + '/cluster_data', mut_rmsds[subfolder+'_all.dat'], '/user/jdragelj/projects/cco_fluorescin/Clustering/all_results/results_20_NOC/inter_intra/'+subfolder+'_')




####################
# No charge frames #
####################