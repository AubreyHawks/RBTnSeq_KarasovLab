#Trim1:  
/uufs/chpc.utah.edu/common/home/u0572090/.local/bin/cutadapt -e 0.0 --overlap 30 -a tctctnnnnnnnnnnnnnnnnnnnncgtac --action=retain -o example_seq_files/trimmed30bp_demo_r1.fq.gz --untrimmed-output example_seq_files/fail_trimmed30bp_demo_r1.fq.gz example_seq_files/novogene_p25.c2_r1.fastq.gz 

#Trim2: 
/uufs/chpc.utah.edu/common/home/u0572090/.local/bin/cutadapt -e 0.0 --overlap 5 -a tctct...cgtac -o example_seq_files/final_trimmed_demo_r1.fq.gz -m 20 --too-long-output example_seq_files/fail_final_trimmed_demo_r1.fq.gz example_seq_files/trimmed30bp_demo_r1.fq.gz 
