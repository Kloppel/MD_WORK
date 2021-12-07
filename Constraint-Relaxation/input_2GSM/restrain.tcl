# load into vmd
#vmd scripting data type : .tcl
# your data
mol load psf cco_and_water.psf pdb cco_and_water.pdb

#HEAT:
set all [atomselect top all]
#segname hemes metals festhalten, alle anderen constrained
#set fix [atomselect top "(resname AUR S5CD S5CP S6OH NIC F3S F4S SF4)"]
#$fix set occupancy 1

set res1 [atomselect top "name CA C O N"]
set res2 [atomselect top "all and not (hydrogen or ion or water or (name CA C O N))"]
$all set beta 0
$all set occupancy 0
$res1 set beta 1.00
$res2 set beta 0.50
$all writepdb namd_cons_for_heat.pdb

#quit
#Befehöl zum Ausführen
#vmd -dispdev text -e restrain.tcl >restrain.out