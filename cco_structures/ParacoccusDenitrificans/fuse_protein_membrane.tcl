#Use this program after loading the pure protein structure with VMD and using the Membrane builder Tool with the VMD Main GUI.
#Extensions -> Modeling -> Membrane Builder 
#Recommended Settings are Lipid: POPC, Membrane X-length: 112, Membrane Y-Length: 112, Output prefix: membrane, Topology CHARMM36 (c36).

set path cEM/O-state/7au6
resetpsf

readpsf $path/7au6_out.psf
coordpdb $path/7au6_out.pdb

readpsf $path/membrane.psf
coordpdb $path/membrane.pdb

writepsf $path/7au6_memb.pdb
writepdb $path/7au6_memb.pdb