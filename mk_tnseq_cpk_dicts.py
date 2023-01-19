from Bio.Seq import Seq
from Bio import SeqIO
import gzip
import pickle

##########################################################################################################
#this will be used to name the pickled output dictionaries and should reflect the barcode pool being used
##########################################################################################################
pool = 'demo'

##########################################################################################################
#input trimmed barcode file from r1 here:
##########################################################################################################
bc_file = 'example_seq_files/final_trimmed_demo_r1.fq.gz'

#This makes dictionaries of the read names to barcodes 
barcodes = SeqIO.to_dict(SeqIO.parse(gzip.open(bc_file, 'rt') , 'fastq'))  

##########################################################################################################
#input sam file associating read names with genome positions which was generated using the insertions from r2
##########################################################################################################
pos_table_file = 'example_seq_files/novogene_p25.c2_r2.fastq.gz.sam'

#########################################################################################################
#input strain annotation gff here
#########################################################################################################
gff_file = '/uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/aubrey/tnseq_pipeline/p25.c2_ann_fin.gff'

#This is a dictionary of the read id (key) to the start position of the insertion (value) in the genome
read_to_pos={}

with open (pos_table_file, 'r') as pos_table:

	for i in pos_table:
		i = i.split('\t')
		if '@' not in i[0]:
			rn = i[0]
#			print(rn)
			pos = i[3]
			read_to_pos[rn] = pos

#This dictionary uses the read id from the read_to_pos dictionary to match positions (values) to barcodes (key)
bc_to_pos={}

for b in barcodes:
	if barcodes[b].id in read_to_pos:
		bc = barcodes[b].seq
		pos = read_to_pos[barcodes[b].id]
		bc_to_pos[bc] = pos	

pos_gene = {}
id_gff = {}
gene_count = {}

with open (gff_file, 'r') as gff:
	
	for g in gff:
		if g[0] != '#':
#			print(g)
			g = g.split("\t")
			start = int(g[3])
			end = int(g[4])
			g_len_10 = (end - start)/10
			start = start + g_len_10
			end = end - g_len_10
			gene_id = g[8].strip('\n')
			desc = [g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],gene_id]
			gene_count[gene_id] = 0
			id_gff[gene_id] = desc
			for key,value in bc_to_pos.items():
				if int(value) >= start and int(value) <= end:
					pos_gene[key] = gene_id

ptg_filename = pool + '_pos_to_gene.cpk'
id_gff_filename = pool + '_id_gff.cpk'
gc_filename = pool + '_gc_dict.cpk'

pickle.dump(pos_gene, open(ptg_filename, "wb"))
pickle.dump(id_gff, open(id_gff_filename, "wb"))
pickle.dump(gene_count, open(gc_filename, "wb"))
			
