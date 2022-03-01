#Use this program after reloading the pdb and psf of the newly fused-into-one-file system of protein and membrane.
set path cEM/O-state/7au6

#set atom selections for membrane, only AA-residues, only moieties, and whole molecule (AA's+moieties)
set membrane [atomselect top "segname POPC"]
set protein [atomselect top "protein"]
set moiety [atomselect top "segname EHEM or segname GHEM or segname META"]
set cco [atomselect top "protein or segname EHEM or segname GHEM or segname META"]

#$membrane moveby [vecinvert [measure center $membrane]]
#$cco moveby [vecinvert [measure center $cco]]
#$moiety moveby [vecinvert [measure center $protein]]
#$protein moveby [vecinvert [measure center $protein]]

resetpsf
writepsf $path/7au6_memb.psf
writepdb $path/7au6_memb.pdb

puts "Please do not forget to check if the program did it's job right. "


#try {
#    $membrane moveby [vecinvert [measure center $membrane]]
#} on error puts "Membrane is probably already at the origin. "
#
#try {
#    $cco moveby [vecinvert [measure center $cco]]
#} on error try {
#    $moiety moveby [vecinvert [measure center $protein]]
#    $protein moveby [vecinvert [measure center $protein]]
#} on error puts "Something went wrong with moving protein and moieties. Execute the program line by line without try clauses to find where which error happened. Often errors here stem from wrong structures being read in, so check previous steps with attention to detail. "