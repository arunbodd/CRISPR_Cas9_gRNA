***Code Starts here***
#Import all packages

import re

#Input and Output files

file_in = open("hg38_new.fasta", "r")
file_out = open("hg38_matched_guideRNA.txt", "w")

#Declare a list of patterns to search

patterns = ["AAG", "GGG", "CGG", "TGG", "AAG", "GAG", "TAG", "CAG"]

for lines in file_in:
        lines = lines.strip()
        search = re.search(">",lines)
        if search:
                header = lines
                header = re.sub(">", "", header)
                header = header.strip().split("|")
                Gene_name = header[1]

#Logic to match the pattern and pull out all the sequences that meet the requirements

        else:
                for i in range(0,len(lines)-20):
                        sequences = lines[i:i+20]
                        for i in range(len(patterns)):
                                if patterns[i] == sequences[-3:]:
                                        pattern0 = sequences[10:17]
                                        seq_len = len(pattern0)
                                        count_G = pattern0.count("G")
                                        percent_of_G = (float(count_G)/seq_len)*100
                                        count_C = pattern0.count("C")
                                        percent_of_C = (float(count_C)/seq_len)*100
                                        if int(percent_of_G) >= 60 or int(percent_of_C) >= 60:
                                                print header,"\t", Gene_name + "\t" + sequences, "\t", int(percent_of_C), "\t", int(percent_of_G)
file_in.close()
file_out.close()
