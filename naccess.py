import os
import urllib
import sys
import subprocess

def naccess_func(protein):
	#os.system("./naccess "+protein+".pdb")
	subprocess.call(["./naccess "+protein+".pdb"])	
	return

def calc_asa(naccess_asa):
	#current_chain=''
	#all_chains=[]
	#asa_per_monomer={}
	ASA=0.0
	for atom in naccess_asa:
		columns=atom.split()
		#if not(columns[4]==current_chain):
		#	all_chains.append(columns[4])
		#	current_chain=columns[4]
		#	asa_per_monomer[current_chain]=0
		#asa_per_monomer[current_chain]+=float(columns[10])

		ASA+=float(columns[10])#change if different column
	return ASA #, all_chains, asa_per_monomer

def download_file(protein):
	flag=1
	try:
		urllib.urlretrieve("http://iris.physics.iisc.ernet.in/cgi-bin/pdbid.pl?pdbId="+protein, filename=protein+"/"+protein+"_temp.pdb")
	except:
		try:
			urllib.urlretrieve("http://iris.physics.iisc.ernet.in/cgi-bin/pdbid.pl?pdbId="+protein,filename=protein+"/"+protein+"_temp.pdb")
		except:
			print "Error in downloading protein file ", protein
			flag=0
	return flag
	
def split_protein_file(protein, chains):
	chain_count={}	#number of unique chains in protein
	for c in chains:
		chain_count[c]=0
		
	with open(protein+"/"+protein+"_temp.pdb") as in_file:
		current_chain=''
		for line in in_file:
			items=line.split()
			if items[0]=='ATOM':
				if not(items[4]==current_chain):# increments the chain_count of chain when new chain starts
					if not(current_chain==''):
						F.close()
					if items[4] in chain_count:#if chain isnt in chain_count, then its a duplicate and won't be counted
						chain_count[items[4]]+=1
						current_chain=items[4]
						F=open(protein+"/"+protein+"_"+current_chain+".pdb", 'w')
				if items[4] in chain_count and  chain_count[items[4]]==1:#only writes the atoms of unique chain in file.
					F.write(line)
	return
	
					
def calc_bsa(ASA,chain_ASA ):
	sum_of_monomers=0.0
	for c in chain_ASA:
		sum_of_monomers+=chain_ASA[c]
	return ASA-sum_of_monomers

def delete_files(protein,ind_chains):
	for f in os.listdir(protein):
		os.remove(protein+"/"+f)	
	os.rmdir(protein)
	os.remove(protein+".asa")#IF CWD OF NACCESS OUTPUT IS DIFFERENT CHANGE PATH HERE
	os.remove(protein+".log")#IF CWD OF NACCESS OUTPUT IS DIFFERENT CHANGE PATH HERE
	os.remove(protein+".rsa")#IF CWD OF NACCESS OUTPUT IS DIFFERENT CHANGE PATH HERE
	for c in ind_chains:
		os.remove(protein+"_"+c+".asa")#IF CWD OF NACCESS OUTPUT IS DIFFERENT CHANGE PATH HERE
		os.remove(protein+"_"+c+".log")#IF CWD OF NACCESS OUTPUT IS DIFFERENT CHANGE PATH HERE
		os.remove(protein+"_"+c+".rsa")#IF CWD OF NACCESS OUTPUT IS DIFFERENT CHANGE PATH HERE
	return
##################################################################

script, protein_list=sys.argv
output_file=open(protein_list.split(".")[0]+"_"+"ASA_BSA.txt", 'w')
input_file=open(protein_list)
for line in input_file:
	try:
		items=line.split()
	
		if len(items)>2:
			protein, symm, chain, num=items
		else:
			protein, symm=items
			num=0
		os.mkdir(protein)	
		naccess_func(protein)#IF CWD OF IPAC2 PDB OUTPUT IS DIFFERENT CHANGE PATH HERE
		try:
			naccess_asa=open(protein+".asa")#IF CWD OF NACCESS OUTPUT IS DIFFERENT CHANGE PATH HERE
		except:
			print "Error in naccess program", protein
			continue

		ASA=calc_asa(naccess_asa)#calculates asa for the IPAC2 protein pdb file
	
		if (num==1):
			BSA=0.0
			output_file.write(protein+" "+str(ASA)+" "+str(BSA))
		elif(num==0):
			BSA="-"
			output_file.write(protein+" "+str(ASA)+" "+str(BSA)) #case when pdb file exists but not enough info about chains
		else:	
			chain_ASA={}
			flag=download_file(protein)
			if flag==0:
			
				continue
			ind_chains=[x[0] for x in chain.split("-")]
			split_protein_file(protein, ind_chains)
			for c in ind_chains:
				naccess_func(protein+"_"+c)
				try:
					naccess_asa=open(protein+"_"+c+".asa")#IF CWD OF NACCESS OUTPUT IS DIFFERENT CHANGE PATH HERE
				except:
					print "Error in naccess program for monomer in", protein
					flag=0
				if flag==0:
					continue
				chain_ASA[c]=calc_asa(naccess_asa)
			if not(flag==0):
				BSA=calc_bsa(ASA,chain_ASA )
				output_file.write(protein+" "+str(ASA)+" "+str(BSA))
				delete_files(protein, ind_chains)
	except:
		print "Error in", protein
output_file.close()
