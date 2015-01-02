import csv
import sqlite3 as sql
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class JGIBlastDB:
	
	def __init__(self,db_name):
		self.con = sql.connect(db_name)
	
	def close(self):
		'''Close database'''
		self.con.close()
		
	def table_info(self):
		'''Show tables and number of rows and columns in each table.'''
		print "\t".join(["Table","Columns","Rows"])
		tables = self.con.execute('''
			SELECT NAME FROM sqlite_master WHERE TYPE="table"''').fetchall()
		table_list = [t[0] for t in tables]
		for table in table_list:
			num_cols = len(self.con.execute('''
				PRAGMA table_info(%s)''' % table).fetchall())
			num_rows = self.con.execute('''
				SELECT Count() FROM %s''' % table).fetchone()[0]
			print "\t".join([table,str(num_cols),str(num_rows)])
	
		
	def _jgi_blast_reader(self,csvfile,dialect):
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
					
	def _jgi_fasta_reader(self,fastaFile):
		count = 0
		records = SeqIO.parse(fastaFile,'fasta')
		for rec in records:
			orgString,protID = tuple(rec.id.split("|")[1:3])
			yield (protID,orgString,str(rec.seq))
			count +=1
			if count % 10 == 0:
				print count
		print "Read %i fasta lines" % count


	def load_jgi_csv(self,csv_file,dialect='excel'):
		'''
		'''
		try:
			self.con.execute('''
				CREATE TABLE
				JGIBlast
				(protID TEXT PRIMARY KEY,
				orgString TEXT,
				hsp TEXT,
				score INTEGER,
				evalue REAL,
				alnLength INTEGER,
				percHitCover REAL,
				hitIdent INTEGER,
				perIdent REAL,
				queryName TEXT,
				queryStart INTEGER,
				queryEnd INTEGER,
				hitName TEXT,
				hitStart INTEGER,
				hitEnd INTEGER,
				organism TEXT,
				dataset TEXT)
				''')
		except sql.Error as e: # catch if table already exists
			return e # Could make this more sophisticated, user input?
		with self.con:
			self.con.executemany('''
				INSERT INTO JGIBlast 
				VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
				''', self._jgi_blast_reader(csv_file,dialect))
			
	def load_jgi_fasta(self,fastaFile):
		'''
		'''
		try:
			self.con.execute('''
				CREATE TABLE
				JGIFasta
				(protID TEXT PRIMARY KEY,
				orgString TEXT,
				protein TEXT,
				FOREIGN KEY(protID) REFERENCES JGIBlast(protID))''')
		except sql.Error as e: # catch if table already exists
			return e # Could make this more sophisticated, user input?
		with self.con:
			self.con.executemany('''
			INSERT INTO JGIFasta VALUES (?,?,?)
			''', self._jgi_fasta_reader(fastaFile))
