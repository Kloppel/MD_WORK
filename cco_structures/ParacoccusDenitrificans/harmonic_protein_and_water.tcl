set path x-ray/O-state/3hb3
mol new $path/3hb3_memb_or_solv_del.psf
mol addfile $path/3hb3_memb_or_solv_del_fixlipids.pdb
set all [atomselect top "all"]
$all set beta 0
set prot [atomselect top "protein or water"]
$prot set beta 1
$all writepdb $path/3hb3_memb_or_solv_del_constraints.pdb

