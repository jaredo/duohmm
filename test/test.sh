../bin/duohmm -M genetic_map_chr20_combined_b37.txt -H switched -O observed
diff observed.haps expected.haps
echo "TEST PASSED"
