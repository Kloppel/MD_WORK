#!/bin/bash/


mv cEM/F-state/7au3.pdb 7au3.pdb
python pdbtools/pdb_delchain.py -C,D 7au3.pdb>7au3_.pdb
python pdbtools/pdb_keepcoord.py 7au3_.pdb>7au3__.pdb
rm 7au3_.pdb 
rm 7au3.pdb



mv cEM/O-state/7au6.pdb 7au6.pdb
python pdbtools/pdb_delchain.py -C,D 7au6.pdb>7au6_.pdb
python pdbtools/pdb_keepcoord.py 7au6_.pdb>7au6__.pdb
rm 7au6_.pdb
rm 7au6.pdb



mv cEM/P-state/7ate.pdb 7ate.pdb
python pdbtools/pdb_delchain.py -C,D 7ate.pdb>7ate_.pdb
python pdbtools/pdb_keepcoord.py 7ate_.pdb>7ate__.pdb
rm 7ate_.pdb
rm 7ate.pdb



mv cEM/R-state/7atn.pdb 7atn.pdb
python pdbtools/pdb_delchain.py -C,D 7atn.pdb>7atn_.pdb
python pdbtools/pdb_keepcoord.py 7atn_.pdb>7atn__.pdb
rm 7atn_.pdb
rm 7atn.pdb



mv x-ray/O-state/1ar1.pdb 1ar1.pdb
python pdbtools/pdb_delchain.py -C,D 1ar1.pdb>1ar1_.pdb
python pdbtools/pdb_keepcoord.py 1ar1_.pdb>1ar1__.pdb
rm 1ar1_.pdb
rm 1ar1.pdb



mv x-ray/O-state/3hb3.pdb 3hb3.pdb
python pdbtools/pdb_delchain.py -C,D 3hb3.pdb>3hb3_.pdb
python pdbtools/pdb_keepcoord.py 3hb3_.pdb>3hb3__.pdb
rm 3hb3_.pdb
rm 3hb3.pdb



#compare 3hb3 to all cEM structures
python pdbtools/pdb_intersect.py 3hb3__.pdb 7au3__.pdb > 3hb3_7au3.pdb
python pdbtools/pdb_intersect.py 3hb3__.pdb 7au6__.pdb > 3hb3_7au6.pdb
python pdbtools/pdb_intersect.py 3hb3__.pdb 7ate__.pdb > 3hb3_7ate.pdb
python pdbtools/pdb_intersect.py 3hb3__.pdb 7atn__.pdb > 3hb3_7atn.pdb

#compare 1ar1 to all cEM structures
python pdbtools/pdb_intersect.py 1ar1__.pdb 7au3__.pdb > 1ar1_7au3.pdb
python pdbtools/pdb_intersect.py 1ar1__.pdb 7au6__.pdb > 1ar1_7au6.pdb
python pdbtools/pdb_intersect.py 1ar1__.pdb 7ate__.pdb > 1ar1_7ate.pdb
python pdbtools/pdb_intersect.py 1ar1__.pdb 7atn__.pdb > 1ar1_7atn.pdb

#find littlest shared structure
python pdbtools/pdb_intersect.py 1ar1__.pdb 3hb3__.pdb > conserved_observations_xray_cco.pdb
python pdbtools/pdb_intersect.py 7au3__.pdb 7au6__.pdb 7ate__.pdb 7atn__.pdb > conserved_observations_eCM_cco.pdb
python pdbtools/pdb_intersect.py 1ar1__.pdb 3hb3__.pdb 7au3__.pdb 7au6__.pdb 7ate__.pdb 7atn__.pdb > conserved_observations_all_cco.pdb

#clean up resulting PDBs for further use
mv 1ar1__.pdb x-ray/O-state/1ar1.pdb
mv 3hb3__.pdb x-ray/O-state/3hb3.pdb
mv 7atn__.pdb cEM/R-state/7atn.pdb
mv 7ate__.pdb cEM/P-state/7ate.pdb
mv 7au6__.pdb cEM/O-state/7au6.pdb
mv 7au3__.pdb cEM/F-state/7au3.pdb

#clean up common ground PDBs
mv 3hb3_7au3.pdb comp/3hb3_7au3.pdb
mv 3hb3_7au6.pdb comp/3hb3_7au6.pdb
mv 3hb3_7ate.pdb comp/3hb3_7ate.pdb
mv 3hb3_7atn.pdb comp/3hb3_7atn.pdb
mv 1ar1_7au3.pdb comp/1ar1_7au3.pdb
mv 1ar1_7au6.pdb comp/1ar1_7au6.pdb
mv 1ar1_7ate.pdb comp/1ar1_7ate.pdb
mv 1ar1_7atn.pdb comp/1ar1_7atn.pdb
mv conserved_observations_xray_cco.pdb comp/conserved_observations_xray_cco.pdb
mv conserved_observations_eCM_cco.pdb comp/conserved_observations_eCM_cco.pdb
mv conserved_observations_all_cco.pdb comp/conserved_observations_all_cco.pdb
