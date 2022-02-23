#removes bad water and then saves the new psf
set path x-ray/O-state/1ar1
mol new $path/1ar1_memb_or_solv.psf
mol addfile $path/1ar1_memb_or_solv.pdb

set all [atomselect top all]
$all set beta 0

# find all waters that overlap protein
set badw1 [atomselect top "water within 3 of protein"]
#set badw2 [atomselect top "water within 3 of resname DLPE DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE^J"]
set badw2 [atomselect top "water within 3 of lipid"]

# remove waters that are in the hydrocarbon region. These can only come from the solvate
# plugin and therefore have segnames that begin with WT (or whatever you used)
#set badw2 [atomselect top "segid TIP3X and same residue as abs(z) < 25"]

$badw1 num
$badw2 num

$badw1 set beta 1
$badw2 set beta 1

#set allbadwater [atomselect top "beta > 0"]
set allbadwater [atomselect top "name OH2 and beta > 0"]
set seglistwater [$allbadwater get segid]
set reslistwater [$allbadwater get resid]

# Now lets build a new psf & pdb w/o badwaters
mol delete all
package require psfgen
resetpsf
readpsf $path/1ar1_memb_or_solv.psf
coordpdb $path/1ar1_memb_or_solv.pdb

foreach segid $seglistwater resid $reslistwater {
   delatom $segid $resid
}

writepsf $path/1ar1_memb_or_solv_clean.psf
writepdb $path/1ar1_memb_or_solv_clean.pdb

mol new $path/1ar1_memb_or_solv_clean.psf
mol addfile $path/1ar1_memb_or_solv_clean.pdb

