#Execute using "source remove_selection" in vmd's active TCL-shell (recognizable by "vmd >" beginning of lines. 
#After startup of VMD, often a single Enter will show this console in the shell.)

#set your workfolder here, leave empty if current directory
set path cEM/O-state/7au6
#set your selection of atoms to delete
set sel [atomselect top "lipid within 2 of protein"]

set atomList [lsort -unique [$sel get {segname resid}]]
$sel delete
package require psfgen
resetpsf

#set your input names here again
readpsf $path/7au6_memb_or_solv.psf
coordpdb $path/7au6_memb_or_solv.pdb

foreach atom $atomList {delatom [lindex $atom 0] [lindex $atom 1]}

#set your output names here
writepsf $path/7au6_memb_or_solv_del.psf
writepdb $path/7au6_memb_or_solv_del.pdb
