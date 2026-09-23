dir=cbc_nlkbm5

rsync -avz --delete thex:~/${dir}/xgrid.txt .

for i in 1; do
   rsync -avz --delete thex:~/${dir}/heat_flux_ns${i}.txt .
   rsync -avz --delete thex:~/${dir}/ptcl_flux_ns${i}.txt .
done


 python flux.py
