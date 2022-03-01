set path cEM/O-state/7au6
mol new $path/7au6_memb_or_solv_del.psf
mol addfile $path/7au6_memb_or_solv_del_fixlipids.pdb
set all [atomselect top "all"]
$all set beta 0
set prot [atomselect top "protein or water"]
$prot set beta 1
$all writepdb $path/7au6_memb_or_solv_del_constraints.pdb

