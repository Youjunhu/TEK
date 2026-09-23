dir=antenna13
  rsync -avz --delete thex:~/${dir}/input.nmlt .
#for i in {13..13}
for i in 5
do
    t=$((4000*i+1))
    t=$(printf "%06d" $t)
    rsync -avz --delete thex:~/${dir}/ms/poloidal_plane_t${t}Apara_neq0 .
    rsync -avz --delete thex:~/${dir}/ms/poloidal_plane_t${t}Phi_neq0 .
    rsync -avz --delete thex:~/${dir}/ms/poloidal_plane_t${t}Apara_nneq0 .
    rsync -avz --delete thex:~/${dir}/ms/poloidal_plane_t${t}Phi_nneq0 .

    #rsync -avz --delete thex:~/${dir}/ms/poloidal_plane_t${t}_nh001Apara .
    #rsync -avz --delete thex:~/${dir}/ms/poloidal_plane_t${t}_nh001Phi .
    
done

python p2.py
