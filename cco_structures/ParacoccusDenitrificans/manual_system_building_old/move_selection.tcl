#Execute using "source remove_selection" in vmd's active TCL-shell (recognizable by "vmd >" beginning of lines. 
#After startup of VMD, often a single Enter will show this console in the shell.)

#set your workfolder here, leave empty if current directory
set path x-ray/O-state/3hb3
#set your selection of atoms to delete
set sel [atomselect top "water within 2 of protein or water within 2 of lipid"]

set atomList [lsort -unique [$sel get {segname resid}]]
$sel delete
package require psfgen
resetpsf

#set your input names here again
readpsf $path/3hb3_memb_or_solv_del.psf
coordpdb $path/3hb3_memb_or_solv_del.pdb

foreach atom $atomList {delatom [lindex $atom 0] [lindex $atom 1]}

#set your output names here
writepsf $path/3hb3_memb_or_solv_del.psf
writepdb $path/3hb3_memb_or_solv_del.pdb
