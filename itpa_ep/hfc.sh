#!/bin/bash
#here-document
cat > t.sh << 'endmsg'
#!/bin/bash
#SBATCH -J ITG
#SBATCH -p hfacnormal01
#SBATCH -t 20:29:00
#SBATCH -N 4
#SBATCH -n 512
#SBATCH -o slurm%j.loop
#SBATCH -e error%j.loop
#SBATCH --comment=yjtest
echo "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST}"
echo "SLURM_NODELIST=${SLURM_NODELIST}"
date
mpirun -n 512 ./TEK
date
endmsg
make clean
make on_hfc
if [[ $? -ne 0 ]]; then
    echo "*****error in compiling****"
    exit
fi

#name=cfetr_turb
name=z16
for i in  400
do
  dir=${name}_${i}
  rsync -avz --delete ./ ../${dir}
  cd ../${dir}
  sed "s/tfxxx.txt/tf${i}.txt/gI" itpa_ep/input_tae.nmlt > input.nmlt

  sbatch -J ${dir} t.sh
  cd -
done
