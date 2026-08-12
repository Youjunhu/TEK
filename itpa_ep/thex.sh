#!/bin/bash
np=512
queue=cp6

cat > t.sh <<endmsg
#!/bin/bash
date
yhrun -n ${np}  -p  ${queue} ./TEK
date
endmsg

# module load  GCC/8.5.0
# module load MPI/mpich/4.0.2-mpi-n-gcc8.5
# module load fftw/3.3.10-gcc8.5-mpi-x

module add Intel_compiler/19.0.4
module load MPI/mpich/4.0.2-mpi-x-icc19.0

module load fftw/3.3.10-icc19.0
module load lapack/3.10.0-icc19.0 

make clean
rm *.txt
make on_thex

if [[ $? -ne 0 ]]; then
    echo "*****error in compiling****"
    exit
fi


name=tae15z38
for i in  400
do
  dir=${name}_${i}
  rsync -avz --delete ./ ../${dir}
  cd ../${dir}
  #k=$((i*10))
  #k=$(echo "4.66*10^(19)*$i" | bc -l) #floating computation using bc
  #sed "/nsegment = xxx/c\nsegment = ${i}" cbc/input_tem.nmlt > input.nmlt
  #sed "/nsegment = xxx/c\nsegment = ${i}" cfetr/input_nscan.nmlt > input.nmlt
  sed "s/tfxxx.txt/tf${i}.txt/gI" itpa_ep/input_tae.nmlt > input.nmlt
  #sed "/density_unit = xxx/c\density_unit = ${i}, ${i}" cbc/kbm_turb.nmlt > input.nmlt
  #cp  itpa_ep/input0.nmlt  input.nmlt
  #sed "/density_unit = xxx/c\density_unit = ${i}, ${i}" cbc/input_kbm.nmlt > input.nmlt
  #sed  -i "/nsegment=/c\nsegment=${k}" input.nmlt
  yhbatch -n ${np} -p ${queue}  --job-name=${dir}  ./t.sh
  cd -
done



