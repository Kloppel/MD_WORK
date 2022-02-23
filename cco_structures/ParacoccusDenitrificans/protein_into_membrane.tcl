#Use this program after loading the pure protein structure with VMD and using the Membrane builder Tool with the VMD Main GUI.
#Extensions -> Modeling -> Membrane Builder 
#Recommended Settings are Lipid: POPC, Membrane X-length: 112, Membrane Y-Length: 112, Output prefix: membrane, Topology CHARMM36 (c36).
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