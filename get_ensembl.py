#!/usr/bin/python3 

import sys
import os
import tempfile
import gzip
import argparse
from ftplib import FTP

os.environ['HTTP_PROXY']="http://www-int.dkfz-heidelberg.de:80"
os.environ['HTTPS_PROXY']="http://www-int.dkfz-heidelberg.de:80"

def parse_attributes(kvps):
	''' Parse the 9th column of GTF into a dictionary of attributes.'''
	f_dict = {}
	kvp = kvps.split(';')[:-1]
	for i in kvp:
		ri = i.replace(' "', '="')
		try:
			k, v = ri.split('=')
			rk = k.replace(' ', '')
			f_dict[rk] = v
		except ValueError:
			print('Malformed file.\nThis line isn\'t good. %s' % kvp)
			exit()
	return f_dict
def read_gtf_as_dict(gtf_file):
	# function reads gtf file and returns a dict with scheme
	# {ensembl_gene_id: HGNC symbol}
	output_dict = {}
	temp_dict = {}
	for line in gtf_file:
		if line.startswith('#'):
			continue
		else:
			if line.split('\t')[2] == 'gene':
				if "gene_name" in line.split("\t")[8].rstrip():
					temp_dict=parse_attributes(line.split('\t')[8].rstrip())
					try:
						ensgid = temp_dict['gene_id']
						ensgid_r = ensgid.replace('"','')
						hgnc_symbol = temp_dict['gene_name']
						hgnc_symbol_r = hgnc_symbol.replace('"','')
						output_dict[ensgid_r] = hgnc_symbol_r
					except KeyError:
						print('Malformed GTF file. Report this error on Gist.)')
						exit()
				else:
					continue			
			else:
				continue
	return output_dict
def avail_organism(version):
	'''
	Function returns available organism for provided version
	'''
	ftp = FTP("ftp.ensembl.org")
	ftp.login()
	ftp.cwd('/pub/release-'+str(version)+'/gtf')
	list_organism = ftp.nlst()
	return list_organism
def get_gtf(organism, version, organism_list):
	'''
	This function downloades GTF from the provided organism.
	First it looks for correspondance than downloads the file
	'''
	if organism not in organism_list:
		print('Provided organism is not in list of Ensembl annotated organism.')
		print('Please check the available list with argument: --list Y')
		sys.exit()
	else:
		print('-------------------------------')
		print('Downloading GTF for '+organism)
		print('-------------------------------')
		ftp = FTP("ftp.ensembl.org")
		ftp.login()
		ftp.cwd('/pub/release-'+str(version)+'/gtf/'+organism)
		list_dir = ftp.nlst()
		to_download = [x for x in list_dir if x.endswith(f'{version}.gtf.gz')==True]
		tmp_file = tempfile.TemporaryFile()
		with open("tmp_file.gz", 'wb') as tmp_file:
			ftp.retrbinary(str('RETR '+to_download[0]), tmp_file.write)

parser = argparse.ArgumentParser(description="Retrieve gene ID conversion table using gene annotation from Ensembl")
parser.add_argument("--organism", help="Organism name")
parser.add_argument("--version", type=int, default=106)
parser.add_argument("-l", "--list", "--list-available", dest='list', help="Append any char to show available organisms.")
args=parser.parse_args()

if __name__ == '__main__':
	available_organism = avail_organism(args.version)
	if args.list is not None:
		for x in available_organism:
			print(x)
		exit()
	elif args.organism == None:
		print('No organism provided. Exiting..')
		exit()
	else:
		get_gtf(organism=args.organism, version=args.version, organism_list=available_organism)
		gtf = gzip.open('tmp_file.gz', 'rt')
		genes_dict = read_gtf_as_dict(gtf)
		output_filename = args.organism+'_gene_annotation.tsv'
		with open(output_filename, 'w') as oput:
			oput.write("Ensembl_gene_ID"+'\t'+"HGNC_symbol"+'\n')		
			for ensgid in genes_dict.keys():
				oput.write(ensgid+'\t'+genes_dict[ensgid]+'\n')
		os.remove("tmp_file.gz")
		print('Gene list file saved under this directory at: --> '+output_filename)
		print('-------------------------------')
		print('Report any error on Gist. Enjoy')