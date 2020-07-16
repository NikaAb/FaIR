"""
Program to format the IMGT High/v-quest output.
The file to format is the summary file, the output of tha annotation step of 
high-Vquest.
The program gather all sequences with the same V gene and allele, J gene and
Junction AA with an identity percentage value higher than Z into one clonotype.
"""

import sys
from optparse import OptionParser
from Levenshtein import distance as levenshtein_distance
# =============================================================================
#               	Read IMGT high/vquest output 
# =============================================================================	

def read_output_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines

def filter_seq(lines,output_file):
	file_name = output_file + "_unannotated_seq.txt"
	filetowrite=open(file_name,"w")
	filtered_lines = []
	for l in range(0,len(lines)):
		seq= lines[l].split("\t")
		if seq[1].rstrip() == "_" or seq[2].rstrip() == "_" or seq[3].rstrip() == "_" :
			seq_unannotated =  seq[0]+ "\t" + seq[1] + "\t" + seq[2] + "\t" + seq[3] +"\n"
			#print seq
			filetowrite.write(seq_unannotated)
		else :
			filtered_lines.append(lines[l])			
	return filtered_lines

def delete_duplicate(lines):
	uniq_seq_dico = {}
	filtered_lines =[]
	dup_corresp = {}
	#print "lines before filter : ",len(lines)
	for l in range(0,len(lines)):
		#print lines[l]
		seq= lines[l].split("\t")
		#print seq
		dup_id = seq[1].rstrip()+"_"+seq[2].rstrip()+"_"+seq[3].rstrip()
		#print dup_id
		if dup_id in  dup_corresp.keys() :
			#print seq[0]
			dup_corresp[dup_id].append(seq[0])
		else :
			dup_corresp[dup_id] = [seq[0]]
			filtered_lines.append(lines[l])
	#print dup_corresp,"dup_corresp"
	#print "lines after filter : ",len(filtered_lines)
	for key in dup_corresp.keys():
		#print dup_corresp[key]
		uniq_seq_dico[dup_corresp[key][0]] = dup_corresp[key][1:]
	print uniq_seq_dico
	return uniq_seq_dico,filtered_lines


#=============================================================================#

def group_same_VJ(lines,formatted_file):
	#print lines
	dico_same_VJ = {} # same Vgene, same J gene
	dicoSeq = {}
	#print lines
	for l in range(0,len(lines)):
		NumClone= lines[l].split("\t")
		dicoSeq[NumClone[0]] = [NumClone[1].rstrip(),NumClone[2].rstrip(),NumClone[3].rstrip()]
		Clone_identity = ""
		#Clone_identity = str(NumClone[1]+"_"+NumClone[2]) # same V gene and allele + same J gene and allele
		#Clone_identity = str(NumClone[1]+"_"+NumClone[2].split("*")[0]) # same V gene and allele + same J gene
		Clone_identity = str(NumClone[1].split("*")[0]+"_"+NumClone[2].split("*")[0]) # same V gene  + same J gene
		if Clone_identity in dico_same_VJ.keys():
			dico_same_VJ[Clone_identity].append(NumClone[0])
		else:
			dico_same_VJ[Clone_identity] = [NumClone[0]]
	#print "dico_same_VJ",len(dico_same_VJ['IGHV3-74_IGHJ4'])
	#print "dico seq", dicoSeq
	return dico_same_VJ,dicoSeq

# =============================================================================
#				
# =============================================================================	

def hamming_distance(s1, s2):
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))
# =============================================================================	
def write_clone_VJ_cdr3(dico_same_VJ,dicoSeq,uniq_seq_dico,output_file):
	file_name = output_file+"_sameVJ_noallele_CDR3_70.txt"
	filetowrite=open(file_name,"w")
	clone_same_VJ_CDR3 = {}
	dicoclone_vj_cdr3 = {}
	written_seq =[]
	#print uniq_seq_dico
	#print dico_same_VJ
	#print "ham_dist",1-(hamming_distance("CNTIGNAVRAW", "CYTIENTVKVW")/float(len("CYTIENTVKVW")))
	for clone_same_VJ in dico_same_VJ.keys():
		print clone_same_VJ,"clone_same_VJ,clone_same_VJ"
		clone_same_VJ_CDR3[clone_same_VJ] = {}
		print dico_same_VJ,"dico_same_VJ"
		
		for seq in dico_same_VJ[clone_same_VJ]:
			list_cdr3 =[]
			print seq,"here is the seq"
			print "line 104"
			if len(clone_same_VJ_CDR3[clone_same_VJ].keys()) != 0 :
				print "line 106"
				print clone_same_VJ_CDR3[clone_same_VJ].keys(),"keys"

				
				for CDR3 in clone_same_VJ_CDR3[clone_same_VJ].keys() :
					print CDR3,"CDR3"
					if len(dicoSeq[seq.rstrip()][2]) == len(CDR3) :
						print "line 110"
						print 1-(hamming_distance(dicoSeq[seq.rstrip()][2] , CDR3)/float(len(CDR3)))
						if 1-(hamming_distance(dicoSeq[seq.rstrip()][2] , CDR3)/float(len(CDR3))) >= 0.7 :
							list_cdr3.append("+")
							#clone_same_VJ_CDR3[clone_same_VJ][CDR3].append(seq)
							#written_seq.append(seq)
							#print written_seq,"written_seq 1"
						#	print "we are here"
						else:
							print "line 116"
							list_cdr3.append("-")
							#clone_same_VJ_CDR3[clone_same_VJ][dicoSeq[seq.rstrip()][2]]= [seq]
							#written_seq.append(seq)
							#print written_seq,"written_seq 2"
							#print "then we are here"
					elif dicoSeq[seq.rstrip()][1].split("*")[0][-1] == "6" :
						length = max(len(dicoSeq[seq.rstrip()][2]),len(CDR3))
						#print length,'lllll'
						#print 1-(levenshtein_distance(dicoSeq[seq.rstrip()][2] , CDR3)/float(length))
						if 1-(levenshtein_distance(dicoSeq[seq.rstrip()][2] , CDR3)/float(length))>= 0.7  : 
							list_cdr3.append("+")
							#clone_same_VJ_CDR3[clone_same_VJ][CDR3].append(seq)
							#written_seq.append(seq)
							#print written_seq,"written_seq 3"
						else :
							#clone_same_VJ_CDR3[clone_same_VJ][dicoSeq[seq.rstrip()][2]]= [seq]
							#written_seq.append(seq)
							#print written_seq,"written_seq 4"
							list_cdr3.append("-")

					else :
						#print dicoSeq[seq.rstrip()][2],"cdr3"
						#print seq,"seq"
						clone_same_VJ_CDR3[clone_same_VJ][dicoSeq[seq.rstrip()][2]]= [seq]
						written_seq.append(seq)
						print written_seq,"written_seq 5"
					print list_cdr3,"list_cdr3"

			else :
				print "line 131"		
				clone_same_VJ_CDR3[clone_same_VJ][dicoSeq[seq.rstrip()][2]]= [seq]
				#print "we are in else 1"
				#print clone_same_VJ_CDR3
				written_seq.append(seq)
		print written_seq,"written_seq"
		print "wwritten"
		print clone_same_VJ_CDR3,"clone_same_VJ_CDR3"
		potential = clone_same_VJ_CDR3[clone_same_VJ].values()
		seq_list = dico_same_VJ[clone_same_VJ]
		#print seq_list,"premier"
		list_loc = []
		i =0
		while len(seq_list) != 0 :
			print i
	  		if seq_list in potential :
	  			#print seq_list,"here"
	  			list_loc.append(seq_list)
	  			seq_list = []
	  		else :
	  			#print seq_list,"seq_list"
	  			#print "\n"
	  			#print potential,"pot"
	  			#print "\n"
	  			if potential != [] :
		  			list_max_element = max(potential, key=len)
					potential.remove(list_max_element)
					#print potential,"potentiel"
					#list_max_element
					list_loc.append(list_max_element)
					#print list_max_element,"list_max_element"
					#print potential
					#print seq_list,"seq_list"
					for x in list_max_element:
							#print x,"xxxxxx"
							if x in seq_list :
				  				seq_list.remove(str(x))
				  			for p in potential :
				  				if x in p :
				  					for pp in p :
				  						if pp not in list_max_element :
				  							list_max_element.append(pp)
				  					potential.remove(p)
				  	#print seq_list,"seq_list"
				  	#print potential
				  	i+=1
				else :
				  	list_loc.append(seq_list)
				  	seq_list = []
				  	i+=1
  		dicoclone_vj_cdr3[clone_same_VJ] = list_loc
  	print dicoclone_vj_cdr3,"dicoclone_vj_cdr3"
	clone_number = 0
	for VJ in dicoclone_vj_cdr3.keys():
		#print VJ
		for cdr3 in dicoclone_vj_cdr3[VJ]:
			#print cdr3
			#print clone_number
			for seq in cdr3 :
				#print seq,"here",dicoSeq[seq.rstrip()][0] 
				sequence = str(clone_number) + "\t" + str(seq) + "\t" + dicoSeq[seq.rstrip()][0] + "\t" +dicoSeq[seq.rstrip()][1] +"\t" +dicoSeq[seq.rstrip()][2]+ "\n"
				filetowrite.write(sequence)
				if len(uniq_seq_dico[seq]) != 0 :
					for dup in uniq_seq_dico[seq] :
						sequence = str(clone_number) + "\t" + str(dup) + "\t" + dicoSeq[seq.rstrip()][0] + "\t" +dicoSeq[seq.rstrip()][1] +"\t" +dicoSeq[seq.rstrip()][2]+ "\n"
						filetowrite.write(sequence)
				"""
				if seq in uniq_seq_dico.keys() :
					for identical_seq in uniq_seq_dico[seq] :
						seq = str(clone_number) + "\t" + str(identical_seq) + "\t" + dicoSeq[seq.rstrip()][0] + "\t" +dicoSeq[seq.rstrip()][1] +"\t" +dicoSeq[seq.rstrip()][2]+ "\n"
						filetowrite.write(seq)	
				"""
			clone_number += 1
	filetowrite.close()
	return 0


def find_max_list(list):
    list_len = [len(i) for i in list]
    return(max(list_len))
#=============================================================================#

def main():
    usage = "python  groupe_clone.py -i <formated IMGT highvquest statistics output> -o <output file name>\n "
    parser = OptionParser(usage)
    parser.add_option("-i", "--hv_stat_output", dest="hv_stat_output",
          help="formated IMGT highvquest statistics output")
    parser.add_option("-o", "--output_file_name", dest="output_file_name",
          help="the name for the file to write")
    
    (options, args) = parser.parse_args()
    if len(sys.argv) != 5:
        parser.error("incorrect number of arguments")
    IMGT_file = options.hv_stat_output
    output_file_name = options.output_file_name
    CloneID = read_output_file(IMGT_file)
    print "all seq : " , len(CloneID)
    annotated_CloneID = filter_seq(CloneID,output_file_name)
    print "annotated seq : ", len(annotated_CloneID)
    uniq_seq_dico,filtered_lines = delete_duplicate(annotated_CloneID)
    print "uniq seq : ", len(filtered_lines)
    dico_same_VJ,dicoSeq = group_same_VJ(filtered_lines,output_file_name)
    
    #write_clone_VJ(dico_same_VJ,dicoSeq,output_file_name)
    
    write_clone_VJ_cdr3(dico_same_VJ,dicoSeq,uniq_seq_dico,output_file_name)

    print ("Done!")

#=============================================================================#

if __name__ == "__main__":
    main()
