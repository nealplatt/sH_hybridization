#there were several nodes that seemed to have issues....blocking them
for node in $(cat block.list); do
    for i in $(seq 1 15); do
       qsub -V -cwd -S /bin/bash -q all.q -j y  -N block -o ~/block.log -pe mpi 1 -l h=$node ~/sleep.sh 
    done
done

