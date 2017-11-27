# CRISPR_Cas9_gRNA
gRNAs were designed for exonic regions of specific genes using certain specific conditions to increase specificity of spCAS9 (NAG and NGG) PAM sequences.
##CODE STARTS HERE##
import re
file_in = open("hg19.fasta", "r") #input file
#Patterns to search
pattern0 ='GGG'
pattern1 = 'AGG'
pattern2 = 'CGG'
pattern3= 'TGG'
pattern4 = 'AAG'
pattern5 = 'TAG'
pattern6 = 'CAG'
pattern7 = 'GAG'
for lines in file_in:
	lines = lines.strip()
	search = re.search(">",lines)
	if search:
		header = lines
#		print header
	else:
		#n = 20
		#sequences= [lines[i:i+n] for i in range(0,len(lines)-20,n)]
		for i in range(0,len(lines)-20):
			sequences = lines[i:i+20]
			#print header,sequences
			if pattern0 == sequences[-3:]:
				pattern0 = sequences[:17]
				count_gc0 = pattern0.count("GC")
				if count_gc0 >=1:
					print header, sequences, count_gc0
			elif pattern1 == sequences[-3:]:
				pattern1 = sequences[:17]
				count_gc1 = pattern1.count("GC")
				if count_gc1 >=2:
					print header, sequences, count_gc1
			elif pattern2 == sequences[-3:]:
				pattern2 = sequences[:17]
				count_gc2 = pattern2.count("GC")
                                if count_gc2 >=2:
					print header, sequences, count_gc2
			elif pattern3 == sequences[-3:]:
				pattern3 = sequences[:17]
				count_gc3 = pattern3.count("GC")
				if count_gc3 >=2:
					print header, sequences, count_gc3
			elif pattern4 == sequences[-3:]:
				pattern4 = sequences[:17]
				count_gc4 = pattern4.count("GC")
				if count_gc4 >=2:
					print header, sequences, count_gc4
			elif pattern5 == sequences[-3:]:
                                pattern = sequences[:17]
                                count_gc5 = pattern5.count("GC")
                                if count_gc5 >=2:
                                        print header, sequences, count_gc5
			elif pattern6 == sequences[-3:]:
                                pattern6 = sequences[:17]
                                count_gc6 = pattern6.count("GC")
                                if count_gc6 >=2:
                                        print header, sequences, count_gc6
			elif pattern7 == sequences[-3:]:
                                pattern7 = sequences[:17]
                                count_gc7 = pattern7.count("GC")
                                if count_gc7 >=3:
                                        print header, sequences, count_gc7
				break
