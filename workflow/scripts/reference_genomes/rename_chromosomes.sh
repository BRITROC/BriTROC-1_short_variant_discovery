sed -E "s/>([^ ]+)/>chr\1/" $1 |
sed -E 's/\s.*$//' |
sed 's/chrMT/chrM/' |
sed -E 's/\.1//' |
sed 's/GL/gl/g' |
sed -E 's/chr(gl00021[1-9])/chrUn_\1/' |
sed -E 's/chr(gl0002[2-4][0-9])/chrUn_\1/' |
sed 's/chrgl000191/chr1_gl000191_random/' |
sed 's/chrgl000192/chr1_gl000192_random/' |
sed 's/chrgl000193/chr4_gl000193_random/' |
sed 's/chrgl000194/chr4_gl000194_random/' |
sed 's/chrgl000195/chr7_gl000195_random/' |
sed 's/chrgl000196/chr8_gl000196_random/' |
sed 's/chrgl000197/chr8_gl000197_random/' |
sed 's/chrgl000198/chr9_gl000198_random/' |
sed 's/chrgl000199/chr9_gl000199_random/' |
sed 's/chrgl000200/chr9_gl000200_random/' |
sed 's/chrgl000201/chr9_gl000201_random/' |
sed 's/chrgl000202/chr11_gl000202_random/' |
sed 's/chrgl000203/chr17_gl000203_random/' |
sed 's/chrgl000204/chr17_gl000204_random/' |
sed 's/chrgl000205/chr17_gl000205_random/' |
sed 's/chrgl000206/chr17_gl000206_random/' |
sed 's/chrgl000207/chr18_gl000207_random/' |
sed 's/chrgl000208/chr19_gl000208_random/' |
sed 's/chrgl000209/chr19_gl000209_random/' |
sed 's/chrgl000210/chr21_gl000210_random/' |
sed 's/chrNC_007605/NC_007605/' |
sed 's/chrhs37d5/hs37d5/'
