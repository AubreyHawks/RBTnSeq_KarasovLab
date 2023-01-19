# This script should:
	# 1. Iterate over the r2 files for the insertion reads (ask me if you have multiple files you'd like to process)
	# 2. Call the indexed reference genome
	# 3. Map the positions of the insertions based on the r2 reads in the ref genome (bwa mem)
	# 4. Convert the sam file to a bam file (samtools view)
	# 5. Sort the reads in the bam file (samtools sort) and index them (samtools index)
	# 6. Associate positions with specific genes (bedtools multicov)

# Outputs I need:
	# 1. A sam file for each of the libraries that contains the read id and the position of the insertion
	# 2. A bedtools output file that associates each position range with a gene


#cd /uufs/chpc.utah.edu/common/home/u0572090/tnseq_pipeline_qc

################################################################
#path to your insertion read file:
################################################################
sample=example_seq_files/novogene_p25.c2_r2.fastq.gz

################################################################
#reference genome for your strain of interest:
################################################################
index=/uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/p25_c2/plate25.C2.pilon.contigs_renamed.fasta 


/uufs/chpc.utah.edu/common/home/u0572090/anaconda3/bin/bwa mem $index $sample > $sample.sam 
		
samtools view -bT $index $sample.sam > $sample.bam
		
samtools flagstat $sample.bam

samtools sort -o $sample.bam $sample.bam

samtools index $sample.bam $sample.bai  

################################################################
#annotation gff goes in this line: 
################################################################
bedtools multicov -bams $sample.bam -bed /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/Talia/08_2021/plate25.c2_fin.gff > $sample.gff

