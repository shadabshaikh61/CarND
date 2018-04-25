input_file_dir = "./Output05"	#directory in which all source and target file are present	
mapping_file_path = "./data type mapping csv.csv"	# this file should be in csv
output_file_dir = "./mismatch"

import os,sys
mapping = {}
f_mapping = open(mapping_file_path,"r")

line = f_mapping.readline()
line = line.strip("'").strip("\n")
col_line = line.lower().split(",")

mapping_index = col_line.index("movement")
connector_index = col_line.index("connector")
dtype_index = col_line.index("ds sql data type")

for line in f_mapping:
	line = "".join(line.split("\n"))
	words = line.split(",")
	
	if(words[connector_index] == ""):
		continue
	
	if(words[connector_index] not in mapping.keys()):
		mapping[words[connector_index]]={}
		
	if(words[mapping_index] not in mapping[words[connector_index]].keys()):
		mapping[words[connector_index]][words[mapping_index]]=[]

	if(words[dtype_index] not in mapping[words[connector_index]][words[mapping_index]]):
		mapping[words[connector_index]][words[mapping_index]].append(words[dtype_index])

for root,dire,files in os.walk(output_file_dir):
	for single_file in files:
		current_file = open(root+"/"+single_file,'r')
		error_file = open()
		current_file.close()

		
