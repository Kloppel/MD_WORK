set path cEM/O-state/7au6
resetpsf

readpsf $path/7au6_out.psf
coordpdb $path/7au6_out.pdb

readpsf $path/membrane.psf
coordpdb $path/membrane.pdb

writepsf $path/7au6_memb.pdb
writepdb $path/7au6_memb.pdb