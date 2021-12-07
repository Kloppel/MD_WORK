
cd prod/

i=13
numberrestart=22

while [ $i -le $numberrestart ]; 
    do
    file=prod$i
    /home/pbuser/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2 +idlepoll +p20 +devices 0 $file.namd > $file.out  
    while [ ! -s prod$i.coor ] ;   
        do
        sleep 80
    done
    i=$(( i + 1 ))
done
