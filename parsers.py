from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys,csv

def jgi_csv_reader(csvfile,dialect):
	count = 0
	with open(csvfile,'rb') as f:
		reader = csv.reader(f,dialect)
		reader.next() # skip header
		for line in reader:
			orgString,protID = tuple(line[0].split("|"))
			newline = [protID,orgString] + [i.strip("%") for i in line[1:]]
			yield tuple(newline)
			count +=1
			if count % 10 == 0:
				print count
	print "Read %i blast lines" % count
					
def jgi_fasta_reader(fastaFile):
	count = 0
	records = SeqIO.parse(fastaFile,'fasta')
	for rec in records:
		orgString,protID = tuple(rec.id.split("|")[1:3])
		yield (protID,orgString,str(rec.seq))
		count +=1
		if count % 10 == 0:
			print count
	print "Read %i fasta lines" % count


def blast_hit_gen(infile):
	'''Given blast xml output (--outfmt 5). Assumes queries had JGI 
	description lines.'''
	count = 0
	skip_count = 0
	with open(infile) as f:
		records = NCBIXML.parse(f)
		for rec in records:
			try:
				orgString,protID = tuple(rec.query.split("|")[1:3])
				## remove blast's "gnl|..." descriptor
				hit = ' ' .join(rec.alignments[0].title.split()[1:])
				yield (protID, orgString, hit, rec.descriptions[0].e)
				count +=1
			except IndexError: # no hits below evalue threshold
				skip_count +=1
				continue # skip
			if count % 10 == 0:
				print count
	print "Found %i hits above threshold" % count
	print "%i hits were below threshold" % skip_count
