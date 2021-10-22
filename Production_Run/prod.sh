#!/bin/bash

#set loop variables
i=0
stop=1

#create final folder
mkdir prod_all/
mkdir prod_all/logs/

#execution loop (while loop):
while [ $i -le $stop ]; 
    do
    #make a directory numbered by the loop variable
    mkdir prod$i/
    #change into new directory
    cd prod$i/
    #execute namd calculation from inside this new folder, write output into log.out
    #execution on x CUDA kernels (+p(x))
    #prod.inp = input file
    #log.out = output file
    /home/pbuser/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2 +idlepoll +p20 +devices 0 prod.inp > log.out
    # while condition: ! not 
    #-s following file exists AND has a size larger than 0
    while [ ! -s cco_prod.coor ] ; 
        do
        sleep 80
    done
    cd ../
    mv prod$i/log.out prod_all/logs/log$i.out
    mv prod$i/
    #

    i=$(( i + 1 ))
done


#delete all subfolders, comment this for finding errors within this bash script
