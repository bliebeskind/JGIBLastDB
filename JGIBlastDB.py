import csv
import parsers
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
	
	## Loading Functions

	def load_jgi_csv(self,csv_file,dialect='excel'):
		'''
		Load a JGI blast csv file into a SQLite database. Creates a new table
		called JGIBlast.
		
		Columns correspond to the csv except the first csv column, "Hit" 
		is broken into two columns in the database: "protID" and "orgString".
		protID is the primay key and is used as a foreign key in other tables.
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
				''', parsers.jgi_csv_reader(csv_file,dialect))
			
	def load_jgi_fasta(self,fastaFile):
		'''
		Load the corresponding fasta into a table called JGIFasta. Assumption
		is that fasta file will have proteins.
		
		Columns:
		
		protID: 	JGI protein ID (TEXT PRIMARY KEY)
		orgString: 	JGI organism identifier (TEXT)
		protein: 	Protein sequence (TEXT)
		** protID is a FOREIGN KEY referencing JGIBlast(protID)
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
			''', parsers.jgi_fasta_reader(fastaFile))
			
	def load_reciprocal_blast(self,blastfile):
		'''
		Load a table of reciprocal blast hits. Infile must be in xml format
		(-outfmt 5).
		
		Columns:
		
		protID: 	JGI protein ID (TEXT PRIMARY KEY)
		orgString: 	JGI organism identifier (TEXT)
		hit: 		Description line of top hit (TEXT)
		evalue: 	REAL
		** protID is a FOREIGN KEY referencing JGIBlast(protID)
		'''
		try:
			self.con.execute('''
				CREATE TABLE
				rBlast
				(protID TEXT PRIMARY KEY,
				orgString TEXT,
				hit TEXT,
				evalue REAL,
				FOREIGN KEY(protID) REFERENCES JGIBlast(protID))''')
		except sql.Error as e:
			return e
		with self.con:
			self.con.executemany('''
				INSERT INTO rBlast VALUES (?,?,?,?)
				''', parsers.blast_hit_gen(blastfile))
			
	## Search Functions
	
	def protIDs_from_rBlast(self,user_string,user_func=None,evalue=10.0):
		'''
		Filter protein ids using a reciprocal blast. Returns a sqlite3 search.
		
		Input string corresponds to the column rBlast.hit. You can pass in 
		a function that parses this column and returns the value matched against
		<user_string>. You can also provide a maximum evalue (must be float), 
		which is 10.0 by default.
		
		Examples:
		>>> search = a.proteinIDs_from_rBlast("jgi|SacceM3707_1|37034|...
		or
		>>> func = lambda x: x.split("|")[2]
		>>> search = a.proteinIDs_from_rBlast("37034",func,evalue=1e-10)
		>>> search.fetchall()
		[(u'196004',),
		(u'36510',)]
		'''
		if user_func:
			try:
				self.con.create_function("uf",1,user_func)
			except sql.OperationalError: # if user function already exists
				pass
			return self.con.execute('''
				SELECT protID from rBlast 
				WHERE uf(rBlast.hit) == ?
				AND rBlast.evalue <= ?''', (user_string,evalue))
		else:
			return self.con.execute('''
				SELECT protID from rBlast
				WHERE rBlast.hit == ?
				AND rBlast.evalue <= ?''', (user_string,evalue))

	def fastas_from_rBlast(self,user_string,user_func=None,evalue=10.0):
		'''
		Filter sequences using a reciprocal blast. Returns generator of 
		SeqRecord objects of sequences that match user input.
		
		Input string corresponds to the column rBlast.hit. You can pass in 
		a function that parses this column and returns the value matched against
		<user_string>. You can also provide a maximum evalue (must be float), 
		which is 10.0 by default.
		
		Examples:
		>>> seqs = a.fastas_from_rBlast("jgi|SacceM3707_1|37034|...
		or
		>>> func = lambda x: x.split("|")[2]
		>>> seqs = a.fastas_from_rBlast("37034",func,evalue=1e-10)
		>>> SeqIO.write(seqs,"rBlast_filtered_seqs.fas",'fasta')
		'''
		if user_func:
			try:
				self.con.create_function("uf",1,user_func)
			except sql.OperationalError: # if user function already exists
				pass
			search = self.con.execute('''
				SELECT b.hitName, f.protein FROM
				JGIBlast b JOIN JGIFasta f on b.protID=f.protID
				JOIN rBlast r on b.protID=r.protID
				WHERE uf(r.hit) == ?
				AND r.evalue <= ?''', (user_string,evalue))
		else:
			search = self.con.execute('''
				SELECT b.hitName, f.protein FROM
				JGIBlast b JOIN JGIFasta f on b.protID=f.protID
				JOIN rBlast r on b.protID=r.protID
				WHERE rBlast.hit == ?
				AND rBlast.evalue <= ?''', (user_string,evalue))
		return (SeqRecord(Seq(j),id=i,description='') for i,j in search)
