bc_file = "example_seq_files/unique_bcs.txt"

with open (bc_file,'r' ) as bcs:

	for b in bcs:
		b = b.strip('\n')
		b = b.split(' ')
		print(b[-1])
