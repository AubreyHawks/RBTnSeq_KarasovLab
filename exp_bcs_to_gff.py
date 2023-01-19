import pickle

#########################################################
#the exact string you used to name your cpk dictionaries
#########################################################
pool = 'demo'

ptg_filename = pool + '_pos_to_gene.cpk'
id_gff_filename = pool + '_id_gff.cpk'
gc_filename = pool + '_gc_dict.cpk'

pos_gene = pickle.load(open(ptg_filename, "rb"))
id_gff = pickle.load(open(id_gff_filename, "rb"))
gene_count = pickle.load(open(gc_filename, "rb"))

exp_bc_files = []
######################################################
#a list of your sample files
######################################################
file_order = '/uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/scripts/rbtnseq_pipeline/example_seq_files/bc_file_list.txt'

######################################################
#the path to the directory with your sample files
######################################################
file_path = '/uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/scripts/rbtnseq_pipeline/example_seq_files/'

with open (file_order, 'r') as order:

	for l in order:
                l = l.strip('\n')
#		l = l.split('\t')
                l = l.split('.')
                full_path = file_path + l[0] + '.txt'
                exp_bc_files.append(full_path)

for exp_file in exp_bc_files:
	for i in gene_count.keys():
		gene_count[i] = 0

	with open(exp_file, 'r') as exp:

		for bc in exp:
			bc = bc.strip('\n')
			if bc in pos_gene.keys():
				gene_count[pos_gene[bc]] += 1

	for i in gene_count.keys():
#		print(id_gff[i])
#		print(gene_count[i])
		gff_line = id_gff[i]
		gff_line.append(gene_count[i])
		id_gff[i] = gff_line
	
for ids in id_gff.values():
	out_line = ids
	print(*out_line, sep='\t')

			
