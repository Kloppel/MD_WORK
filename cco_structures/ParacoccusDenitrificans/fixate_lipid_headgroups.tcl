set path cEM/O-state/7au6
mol new $path/7au6_memb_or_solv_del.psf
mol addfile $path/7au6_memb_or_solv_del.pdb
set all [atomselect top "all"]
$all set beta 0
set fixed [atomselect top "(resname POPC and name O2 P1 O3 O4 O1 C15 H52 H51 H11 C11 H12 N C14 H14 H42 H43 H41 C12 H22 H23 H21 C13 H33 H31 H32)"]
$fixed set beta 1
$all writepdb $path/7au6_memb_or_solv_del_fixlipids.pdb 


