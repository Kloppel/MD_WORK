set path output_building
mol new $path/gag_ma_popc_ionized_solvated_clean.psf
mol addfile $path/gag_ma_popc_ionized_solvated_clean.pdb
set all [atomselect top "all"]
$all set beta 0
set fixed [atomselect top "(resname POPC and name O2 P1 O3 O4 O1 C15 H52 H51 H11 C11 H12 N C14 H14 H42 H43 H41 C12 H22 H23 H21 C13 H33 H31 H32)"]
$fixed set beta 1
$all writepdb $path/gag_ma_popc_ionized_solvated_clean_fixed_lipid_head_groups.pdb 


