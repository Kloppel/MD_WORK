#Execute using "source remove_selection" in vmd's active TCL-shell (recognizable by "vmd >" beginning of lines. 
#After startup of VMD, often a single Enter will show this console in the shell.)

#set your workfolder here, leave empty if current directory
set path cEM/O-state/7au6
#set your selection of atoms to delete
set sel [atomselect top "lipid within 1 of protein"]

set atomList [lsort -unique [$sel get {segname resid}]]
$sel delete
package require psfgen
resetpsf

#set your input names here again
readpsf $path/3hb3_memb_or.psf
coordpdb $path/3hb3_memb_or.pdb

foreach atom $atomList {delatom [lindex $atom 0] [lindex $atom 1]}

#set your output names here
writepsf $path/3hb3_memb_or.psf
writepdb $path/3hb3_memb_or.pdb
