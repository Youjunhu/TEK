
dir=cbc_nlkbm3
 rsync --delete -avz thex:~/${dir}/phi_evolution000000000.txt phi.txt
 rsync --delete -avz thex:~/${dir}/apara_evolution000000000.txt apara.txt
 rsync --delete -avz thex:~/${dir}/profiles_ns2.txt profiles.txt
 rsync --delete -avz thex:~/${dir}/input.nmlt input.nmlt

python p.py &

python p00.py &
python p00b.py
