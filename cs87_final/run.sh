make clean
make

# OUTFILE = results.txt

declare -a arr=(
                # "16621981"
                # "198729621539"
                # "2516160469693133"
                "18567078082619935259"
                # "6993666669337710100100501"
                # "812945258035564110179904496619"
                # "31015750616613538167589387786383061"
                # "99887766532235934110673400915598003333"
                # "8740616151424572261591492920109047657941"
                )

for epoch in 1; do
  # echo -e "\n ****[EPOCH $epoch]****" >> results.txt
  for N in ${arr[@]}; do
    ./step1 $N
    echo -e "E:$epoch:N:$N" >> results.txt
    (time ./step2 $N) &>> results.txt
  done
done
