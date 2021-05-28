#!/usr/bin/env bash
if [[ ! -d test_data ]]; then
 echo "test_data directory missing!"
 exit 1
fi
cd test_data
# array element format:
# 
arrins=("short_reads" "short_reads_and_superreads" "long_reads" "long_reads")
arrparms=("" "" "-L" "-L -G human-chr19_P.gff")
arrout=("short_reads" "short_reads_and_superreads" "long_reads" "long_reads_guided")
arrmsg=("Short reads"  "Short reads and super-reads" \
   "Long reads" "Long reads with annotation guides")
for i in ${!arrmsg[@]}; do
 fout="${arrout[$i]}.out.gtf"
 /bin/rm -f $fout
 fcmp="${arrout[$i]}.out_expected.gtf"
 if [ ! -f $fcmp ]; then
   echo "Error: file $fcmp does not exist! Re-download test data."
   exit 1
 fi
 n=$i
 ((n++))
 echo "Test ${n}: ${arrmsg[$i]}"
 fin=${arrins[$i]}.bam
 ../stringtie ${arrparms[$i]} -o $fout $fin
 if [ ! -f $fout ]; then
   echo "Error: file $fout not created! Failed running stringtie on $fin"
   exit 1
 fi
 if diff -q -I '^#' $fout $fcmp &>/dev/null; then
    echo "  OK."
 else
   echo "Error: test failed, output $fout different than expected ($fcmp)!"
   #exit 1
 fi
done
