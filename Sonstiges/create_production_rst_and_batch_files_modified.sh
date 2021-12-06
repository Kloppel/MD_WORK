# create files for production run

i=1	#i=2
startsteps=5000000
endsteps=10000000
numberrestart=2

while [ $i -le $numberrestart ]; do

a=$(( i - 1 ))
s=$(( a * startsteps ))
n=$(( s + startsteps ))
t=$(( n + startsteps ))
#echo $t

file1=prod$i.namd
file2=prod$a.sh

if [ $i -ne 1 ]; then
	cat < prod1.namd |
		sed -e 's/'"firsttimestep       $startsteps"'/'"firsttimestep       $n"'/g'|
	sed -e 's/'"numsteps         $endsteps"'/'"numsteps         $t"'/g'|
	sed -e 's/prod1/'"prod$i"'/g' |
	sed -e 's/prod1.dcd/'"prod$i.dcd"'/g' |
	sed -e 's/prod1.xst/'"prod$i.xst"'/g' |
	sed -e 's/prod0.xsc/'"prod$a.xsc"'/g' |
	sed -e 's/prod0.coor/'"prod$a.coor"'/g' | 
	sed -e 's/prod0.vel/'"prod$a.vel"'/g' > $file1
	
fi



i=$(( i + 1 ))

done


