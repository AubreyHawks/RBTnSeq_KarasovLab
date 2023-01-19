###########             FOR RUNNING THIS SCRIPT DO source $name_of_script      #####################################

############################################################
#path to the folder with reads demultiplexed by p2
#############################################################
cd seq_reads/

##############################################################
#you will need to run this individually for each p2 directory
##############################################################
#1. Make sure your 'index.txt' file has been edited to associate the correct sample names with each p1 index barcode
#2. Change the ITxxx to match the correct p2 directory

for i in IT024
do
	echo $i
	cd $i
	echo $i
	echo *L1_1*
	############################################################################################
	#check that the F1 and F2 below correspond to the end of your sequencing read files
	############################################################################################
	je demultiplex F1= *L5_1.fq.gz F2= *L5_2.fq.gz BF=../../index.txt O=demultiplexed BPOS=READ_1
	cd demultiplexed
	for l in *_N*_2*
	do
		echo $l
#when plates are introduced cut needs to be changed for -c1-6
		foldername=$(echo $l | cut -c1-6)
		mkdir $foldername
		mv $l $foldername
		cd $foldername
#cutadapt trims everything before and after the expected barcode
        cutadapt -a GATGTCCACGAGGTCTCT...CGTACGCTGCAGGTCGAC -o trimmed_"$foldername"_r2.fastq.gz --discard-untrimmed $l
#keep only the sequence of the barcode
        zcat trimmed_"$foldername"_r2.fastq.gz | paste - - - - | cut -f2 > "$foldername"_r2_seq.txt
#count the number of barcodes
        reads=$(cat "$foldername"_r2_seq.txt | wc -l)
#############################################################################################################################################################################################
#makes the file that contains the names of the files with the barcodes and the number of reads. I have not managed to have the file as tab separated so one has to edit it with a text editor
# edit the name of the text file here
#############################################################################################################################################################################################
        echo "$foldername"_r2_seq.txt$reads >> ../../../barseq_221018.txt
		cp "$foldername"_r2_seq.txt ../../../files_counting
		cd ../
	done
	cd ../../
done
