***Code Starts here***
#Import all packages
import re

#Input and Output files

file_in = open("hg19.fasta", "r")
file_out = open("hg19_matched_guideRNA.txt", "w")

#Declare a list of patterns to search

patterns = ["AAG", "GGG", "CGG", "TGG", "AAG", "GAG", "TAG", "CAG"]

for lines in file_in:
	lines = lines.strip()
	search = re.search(">",lines)
	if search:
		header = lines
		header = re.sub(">", "", header)
		header = header.strip().split("|")
		Geneid = header[0]
		Gene_name = header[1]
		Exon_rank = header[5]
		Exon_id = header[4]

	else:
		for i in range(0,len(lines)-20):
			sequences = lines[i:i+20]
			for i in range(len(patterns)):
				if patterns[i] == sequences[-3:]:
					pattern0 = sequences[11:17]
					count_gc0 = pattern0.count("GC")
					if count_gc0 >=2:
						print>>file_out, Geneid + "\t" +Exon_id +"\t"+ Gene_name +"\t"+ Exon_rank + "\t" + sequences, "\t", count_gc0
file_in.close()
file_out.close()
