import webview
import subprocess 
import os
import json
import tempfile
import pysam
import regex
import platform
from collections import defaultdict
import argparse
import sys

script_path = os.path.dirname(os.path.realpath(__file__))
catalog = ''
config = {
	"Fixed": {
		'CatalogFile':  os.path.join(script_path, 'assets/catalog.json')
	},
	"Dynamic": {}
}

class Api:
	def readConfigurationFile(self):
		"""Reads config file content

		Returns:
			dict: Configuration
		"""
		config_file = os.path.join(script_path, 'assets/config.json')

		with open(config_file, 'r') as f: content = json.load(f)
		
		return content


	def getDynamicConfiguration(self):
		"""Reads configuration (dynamic part) for showing it in GUI's configuration window

		Returns:
			config (dict): Current configuration
		"""
		return config['Dynamic']


	def updateConfiguration(self, new_conf):
		"""Updates configuration file on disk and in program

		Args:
			new_conf (dict): New configuration settings from front end

		Returns:
			bool: True if settings were saved
		"""
		global config
		config_file = os.path.join(script_path, 'assets/config.json')

		if new_conf:
			with open(config_file, 'w') as conf_file:
				json.dump(new_conf, conf_file, indent = 4)
			
			config['Dynamic'] = self.readConfigurationFile()
			
			return True
		else:
			return False


	def readCatalogue(self):
		"""Read in catalogue file

		Returns:
			dict: Catalogues file as dictionary
		"""
		with open(config['Fixed']['CatalogFile'], mode = 'r') as cat_file:
			catalog = json.load(cat_file)

		return catalog


	def getListOfGenesAndCohorts(self):
		"""Get list of genes and cohorts from catalogue which will be displayed in selection boxes

		Returns:
			dict: List of genes and cohorts with coordinates
		"""
		selection_list = {'Genes': [],
						  'Cohorts': []
						 }

		for entry in catalog:
			selection_list['Genes'].append({'Gene': entry['Gene'],
											'ReferenceCoordinates': entry['ReferenceCoordinates']
											})

			for cohort in entry['Cohorts']:
				selection_list['Cohorts'].append({	'Gene': entry['Gene'],
													'ID': cohort['ID'],
													'Name': cohort['Name'],
													'Source': cohort['Source'],
													})


		selection_list["Genes"].sort(key=lambda e: e['Gene'])

		return selection_list


	def selectSeqFile(self):
		"""Open the dialog box in app to select an aligned sequencing file for analysis

		Returns:
			str: Sample's file path
		"""

		file_types = ('BAM file (*.bam)', 'All files (*.*)')
		file_path = window.create_file_dialog(webview.OPEN_DIALOG, allow_multiple = False, file_types = file_types)

		return file_path


	def extractReads(self, file, gene_id, genome_name, tmpbam, tmpfastq, tmpfastq_r2):
		"""Extract out reads from a BAM file and convert them to FASTQ

		Args:
			file (str): File path
			gene_id (str): Gene
			genome_name (str): Genome_
			tmpbam (str): File path for temprary BAM (extracted reads)
			tmpfastq (str): File path of temporary FASTQ (extracted reads)
			tmpfastq_r2 (str): File path of temporary FASTQ R2 reads if found (extracted reads)


		Returns:
			_type_: _description_
		"""
		input_file = pysam.AlignmentFile(file, 'rb')
		output_file = pysam.AlignmentFile(tmpbam, 'wb', template = input_file)
		
		gene_data = next(item for item in catalog if item["Gene"] == gene_id)
		chrom, start, end = regex.split(':|-', gene_data["ReferenceCoordinates"][genome_name])

		# Extract out reads and save to temporary output file
		fastq_reads = []
		fastq_reads_r2 = []

		# first pass get all reads in region, and if pairs found write to fastq files
		read_dict = defaultdict(lambda: [None, None])


		# if pairs are found in region, write to fastq at same time
		for read in input_file.fetch(chrom, int(start), int(end)):
			output_file.write(read)
			
			if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
				continue
			qname = read.query_name
			if qname not in read_dict:
					if read.is_read1:
						read_dict[qname][0] = read
					else:
						read_dict[qname][1] = read
			else:
					if read.is_read1:
						#ikzf1reads.write(read)
						#ikzf1reads.write(read_dict[qname][1])
						fastq_reads.append("@%s\n%s\n+\n%s\n" % (read.query_name, read.query_sequence, "".join(map(lambda x: chr( x+33 ), read.query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr
						fastq_reads_r2.append("@%s\n%s\n+\n%s\n" % (read_dict[qname][1].query_name, read_dict[qname][1].query_sequence, "".join(map(lambda x: chr( x+33 ), read_dict[qname][1].query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr

					else:
						#ikzf1reads.write(read_dict[qname][0])
						#ikzf1reads.write(read)
						fastq_reads.append("@%s\n%s\n+\n%s\n" % (read_dict[qname][0].query_name, read_dict[qname][0].query_sequence, "".join(map(lambda x: chr( x+33 ), read_dict[qname][0].query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr
						fastq_reads_r2.append("@%s\n%s\n+\n%s\n" % (read.query_name, read.query_sequence, "".join(map(lambda x: chr( x+33 ), read.query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr

					del read_dict[qname]


		unmapped_pairs = 0
		pair_locations = set()
		print("Inspecting dictionary...")

		print("Dict length remaining after first pass is",len(read_dict))
		for key,(r1,r2) in read_dict.items():
						#print("Read key is",key)
						if not r1: 
							#print("Read is R2 located on: ",r2.reference_name)
							#print(r2.reference_id)
							#print("Pair is R1 with id: ",r2.next_reference_id)
							#print("Pair is R1 located on: ",r2.next_reference_name)
							#print(r2.next_reference_id)
							#print(r2.next_reference_name)
							#print(r2.next_reference_start)
							if (r2.mate_is_unmapped):
								unmapped_pairs = unmapped_pairs+1
							else:
								pair_locations.add(r2.next_reference_name)

							#for read_2nd in  pysam.AlignmentFile(file, 'rb').fetch(r2.next_reference_name, r2.next_reference_start):
							#print(r2.mate())
							#	if key == r2.reference:
							#		print(read_2nd) 

						else:
							#print("Read is R1 located on: ",r1.reference_name)
							#print("Pair is R2 located on: ",r1.next_reference_name)
							

							if (r1.mate_is_unmapped):
								unmapped_pairs = unmapped_pairs+1
							else:
								pair_locations.add(r1.next_reference_name)
							
							# print(r1.reference_id)
							# print(r1.next_reference_id)
							# print(r1.next_reference_name)
							# print(r1.next_reference_start)
							# print(r1.mate_is_unmapped)
							#print(r1.mate())
 
		print("Dict items with unmapped pairs",unmapped_pairs)
		print("Dict items pair locations",pair_locations)
		print("Rechecking all of chrom for pairs...")

		# co-ordinate may have been wrong, or different chromsome - check each chr again with mate in original region
		for contig in pair_locations:
			for read in input_file.fetch(contig):
							if read.is_unmapped:
								print ("read is unmapped")
							if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
									continue
							qname = read.query_name
							if qname in read_dict: # could be missing pair
								if read.is_read1 and read_dict[qname][0] == None  :
									fastq_reads.append("@%s\n%s\n+\n%s\n" % (read.query_name, read.query_sequence, "".join(map(lambda x: chr( x+33 ), read.query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr
									fastq_reads_r2.append("@%s\n%s\n+\n%s\n" % (read_dict[qname][1].query_name, read_dict[qname][1].query_sequence, "".join(map(lambda x: chr( x+33 ), read_dict[qname][1].query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr

									del read_dict[qname]
								elif read.is_read2 and read_dict[qname][1] == None :
									fastq_reads.append("@%s\n%s\n+\n%s\n" % (read_dict[qname][0].query_name, read_dict[qname][0].query_sequence, "".join(map(lambda x: chr( x+33 ), read_dict[qname][0].query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr
									fastq_reads_r2.append("@%s\n%s\n+\n%s\n" % (read.query_name, read.query_sequence, "".join(map(lambda x: chr( x+33 ), read.query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr

									del read_dict[qname]


		# TODO can we exit these loops if dic is empty? i.e have found all stray pairs
        # Only if there are unmapped pairs *and* it is enabled, do we carry on until EOF for unmapped
		print("after chr re-check pass, dict length remaining is",len(read_dict))
		if len(read_dict) > 0:
			print("Unmapped reads pairs found",len(read_dict))
			if config['Dynamic']['TwoPass']:

				for read in input_file.fetch(until_eof=True):
					if read.is_unmapped:
						print ("read is unmapped")
					if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
							continue
					qname = read.query_name
					if qname in read_dict: # could be missing pair
						if read.is_read1 and read_dict[qname][0] == None  :
							fastq_reads.append("@%s\n%s\n+\n%s\n" % (read.query_name, read.query_sequence, "".join(map(lambda x: chr( x+33 ), read.query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr
							fastq_reads_r2.append("@%s\n%s\n+\n%s\n" % (read_dict[qname][1].query_name, read_dict[qname][1].query_sequence, "".join(map(lambda x: chr( x+33 ), read_dict[qname][1].query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr

							del read_dict[qname]
						elif read.is_read2 and read_dict[qname][1] == None :
							fastq_reads.append("@%s\n%s\n+\n%s\n" % (read_dict[qname][0].query_name, read_dict[qname][0].query_sequence, "".join(map(lambda x: chr( x+33 ), read_dict[qname][0].query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr
							fastq_reads_r2.append("@%s\n%s\n+\n%s\n" % (read.query_name, read.query_sequence, "".join(map(lambda x: chr( x+33 ), read.query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr

							del read_dict[qname]
				print("3rd pass, Dict length remaining is",len(read_dict))
				if len(read_dict) > 0:
					print("Reads not matched:",len(read_dict))
			else:
				print("Unmapped reads not enabled in config.")


		else:
			print("Skipping unmapped check")

			# if read.is_paired:
			# 	if read.is_read1:
			# 		fastq_reads.append("@%s\n%s\n+\n%s\n" % (read.query_name, read.query_sequence, "".join(map(lambda x: chr( x+33 ), read.query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr
			# 	elif read.is_read2: 
			# 		fastq_reads_r2.append("@%s\n%s\n+\n%s\n" % (read.query_name, read.query_sequence, "".join(map(lambda x: chr( x+33 ), read.query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr

			# else:
			# 	fastq_reads.append("@%s\n%s\n+\n%s\n" % (read.query_name, read.query_sequence, "".join(map(lambda x: chr( x+33 ), read.query_qualities)))) # Reads in FASTQ format, convert numeric qualities to chr

		# Write FASTQ file
		with open(tmpfastq, "w") as fastq_out:
			for r in fastq_reads:
				fastq_out.write(r)

		with open(tmpfastq_r2, "w") as fastq_out_r2:
				for r in fastq_reads_r2:
					fastq_out_r2.write(r)

		input_file.close()
		output_file.close()

		return output_file


	def processSequencingFile(self, gene_id, cohort_id, genome_name, bamfile):
		"""Processes the input BAM file - extract out reads and call with tinyT

		Args:
			gene_id (str): Gene
			cohort_id (str): Cohort
			genome_name (str): Genome
			bamfile (str): Input BAM file path

		Returns:
			bool, str, str: Returns three values: whether there was an error (True/False), output of tinyT and Error message
		"""
		tinytHasError = False

		# Assign file index path, try file.bai and file.bam.bai (or cram and crai)
		file_index = bamfile+'.crai' if bamfile.endswith('cram') else bamfile+'.bai'
		if not os.path.isfile(file_index):
			file_index = os.path.splitext(bamfile)[0]+'.crai' if bamfile.endswith('cram') else os.path.splitext(bamfile)[0]+'.bai'

		# Create temporary BAM and FASTQ files which will contain extracted out reads
		tmp_bam = tempfile.NamedTemporaryFile(suffix = '.bam').name
		tmp_bai = None
		tmp_fastq = tempfile.NamedTemporaryFile(suffix = '.fastq').name
		tmp_fastq_r2 = tempfile.NamedTemporaryFile(suffix = '.fastq').name

		try:
			# Check if the file index (.bai) exists. If not then return an error
			if not os.path.exists(file_index):
				tinytHasError, tinyt_stderr, tinyt_output = True, "Error", "The selected BAM file is not indexed."

			if not tinytHasError:
				self.extractReads(bamfile, gene_id, genome_name, tmp_bam, tmp_fastq, tmp_fastq_r2)

				# Check the generated BAM file. If size less than 10 KB then it does not contain any reads and return an error.
				bam_file_size = os.stat(tmp_bam).st_size / 1024 # Generated BAM file size in kilobytes
				if bam_file_size < 10:
					tinytHasError, tinyt_stderr, tinyt_output = True, "Error", "No reads were extracted out."

			# Sort and index the generated BAM file and forward the content to Server along with some metadata
			if not tinytHasError:
				pysam.sort('-o', tmp_bam, tmp_bam)
				pysam.index(tmp_bam)
				tmp_bai = tmp_bam + '.bai'
				content_bam = open(tmp_bam, 'rb')
				content_bai = open(tmp_bai, 'rb')

				files = [('file', tmp_fastq),
						('file', tmp_fastq_r2),
						 ('file', content_bam), 
						 ('file', content_bai)]
				
				metadata = {'gene': gene_id,
							'cohort': cohort_id,
							'genome': genome_name}
				
				try:
					tinyt_output, tinyt_stderr = self.tinytMakeCalls(files, metadata)
				except RuntimeError as e:
					tinyt_output = e
					tinytHasError = True
				else:
					tinytHasError = False
				finally:
					content_bam.close()
					content_bai.close()			
				
		finally:
			if config['Dynamic']['DeleteTempFiles']:
				os.remove(tmp_fastq)
				os.remove(tmp_fastq_r2)

				os.remove(tmp_bam)
				os.remove(tmp_bai)
				print("Temporary files deleted")
			else:
				print("The following temporary files were not deleted:")
				print(tmp_fastq)
				print(tmp_fastq_r2)
				print(tmp_bam)
				print(tmp_bai)

		return tinytHasError, tinyt_output, tinyt_stderr


	def tinytMakeCalls(self, files, metadata):
		"""Use tinyT to make calls on the files containing extracted out reads

		Args:
			files (dict): Path of input files (BAM/FASTQ)
			metadata (dict): Other data important for making calls (such as gene, cohort, genome)

		Raises:
			RuntimeError: TinyT stderr

		Returns:
			str: Output calls (list)
		"""
		try:

			local_platform = platform.system()
			script_path = os.path.dirname(os.path.realpath(__file__))
			os.path.dirname(os.path.realpath(__file__))

			if local_platform == "Darwin":
				print("Platform: macos")
				local_binary  = os.path.join(script_path,"bin/osx/tinyt_macos")
			elif local_platform == "Windows":
				print("Platform: win")

				local_binary = os.path.join(script_path,"bin/win/tinyt.exe")
			else:
				print("Platform: linux")
				#print(os.path.dirname(os.path.realpath(__file__)))	
				local_binary = os.path.join(script_path, "bin/linux/tinyt_amd64")

			# check if we found paired end by size of r2 file
			# files[0][1] = fastq, files[1][1] = fastq_r2, ,files[2][1] = bam and files[3][1] = bai file

			if (os.stat(files[1][1]).st_size == 0):
				print("Single end reads call to tinyt")

				calls = subprocess.run([local_binary,"map","-i",os.path.join(script_path,"assets/ikzf1.idx"),files[0][1]], capture_output=True, text=True)
				calls.check_returncode()
			else:
				#paired reads
				print("Paired end reads call to tinyt")

				calls = subprocess.run([local_binary,"map","-i",os.path.join(script_path,"assets/ikzf1.idx"),files[0][1],files[1][1]], capture_output=True, text=True)
				calls.check_returncode()
				#print(calls.stdout)
		except subprocess.CalledProcessError as e:
			print("Error calling tinyt")
			print(e)

			raise RuntimeError(calls.stderr)
		else:
			return calls.stdout.replace("\n","\n"),calls.stderr.replace("\n","\n")
	

	def getBackgroundData(self, gene_id, cohort_id):
		"""Read in background data from file for a particular Cohort

		Args:
			gene_id (str): Gene
			cohort_id (str): Cohort

		Returns:
			bool, str: Error (T/F) and file content
		"""
		coordinates_data, cohort_data = self.getInfoFromCatalogue(gene_id, cohort_id)
		csv_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'assets/dev/', cohort_data["Source"])

		try:
			file = open(csv_file, 'r')
		except OSError as e:
			hasError = True
			output = e
		else:
			hasError = False
			with file:
				output = file.read()

		return hasError, output


	def getThresholdData(self, gene_id):
		"""Read in threshold data from file for a particular Gene

		Args:
			gene_id (str): Gene

		Returns:
			bool, str: Error (T/F) and file content
		"""
		csv_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'assets/dev/' + gene_id + '_thresholds.csv')

		try:
			file = open(csv_file, 'r')
		except OSError as e:
			hasError = True
			output = e
		else:
			hasError = False
			with file:
				output = file.read()

		return hasError, output


	def getInfoFromCatalogue(self, gene_id, cohort_id):
		"""Get coordinates and cohort data from a catalogue for a specific gene and cohort

		Args:
			gene_id (str): Gene
			cohort_id (str): Cohort

		Returns:
			dict, dict: Coordinates, Cohort data
		"""
		gene_data = next(item for item in catalog if item["Gene"] == gene_id)
		coordinates_data = gene_data["ReferenceCoordinates"]
		cohort_data = next(item for item in gene_data["Cohorts"] if item["ID"] == cohort_id)

		return coordinates_data, cohort_data


	def analyseSample(self, gene, cohort, refgenome, bamfile):
		"""Runs the analysis for a given sample

		Args:
			gene (str): Gene ID
			cohort (str): Cohort ID
			refgenome (str): Reference genome ID
			bamfile (str): BAM file path

		Returns:
			dict: Results from server in Python dict/JavaScript JSON format
		"""

		tinytHadError, results, stderr_output = self.processSequencingFile(gene, cohort, refgenome, bamfile)
		backgroundHadError, background = self.getBackgroundData(gene, cohort)
		thresholdsHadError, thresholds = self.getThresholdData(gene)
		coordinates_data, cohort_data = self.getInfoFromCatalogue(gene, cohort)

		analysis_info = {
			'SampleFile': os.path.basename(bamfile),
			'Gene': gene,
			'GeneCoordinates': coordinates_data[refgenome],
			'CohortName': cohort_data["Name"]
		}

		response = {
			'background': str(background),
			'thresholds': str(thresholds),
			'analysis_info': analysis_info,
			'errors': {
				"tinyt": tinytHadError if tinytHadError else '',
				"background": backgroundHadError,
				"thresholds": thresholdsHadError
			},
			'tinyt_results': str(results), 
			'tinyt_log': str(stderr_output)
		 }

		return response


if __name__ == '__main__':
	"""Starts the app with GUI unless command line options given
	"""

	if len(sys.argv) > 1:

		parser = argparse.ArgumentParser(description='TobleroneApp command line options.')

		#command line mode - extract reads or open gui as per normal
		parser.add_argument("-e","--extract", action='store_true')
		parser.add_argument('-i','--index', type=argparse.FileType('r'), nargs="?", default="./assets/ikzf1.idx")
		parser.add_argument('-g','--gene', default="IKZF1")
		parser.add_argument('-r','--reference', default="hg38")
		parser.add_argument('-b','--bam', type=argparse.FileType('r'), required=True)
		parser.add_argument('-1','--read1',  default=tempfile.NamedTemporaryFile(suffix = '.fastq').name)
		parser.add_argument('-2','--read2',  default=tempfile.NamedTemporaryFile(suffix = '.fastq').name)
	
		args = parser.parse_args()
		with args.index as index,args.bam as bam:
		# Create temporary BAM and FASTQ files which will contain extracted out reads
			tmp_bam = tempfile.NamedTemporaryFile(suffix = '.bam').name
			tmp_bai = None
			tmp_fastq = args.read1
			tmp_fastq_r2 = args.read2

			print(index.name)
			print(bam.name)
			print(args.gene)
			api = Api()
			catalog = api.readCatalogue() # Load catalogue
			extract_result = api.extractReads(bam.name, args.gene, args.reference, tmp_bam, tmp_fastq, tmp_fastq_r2)
			print(tmp_fastq, tmp_fastq_r2)
			print(extract_result)
	else:
		api = Api()
		catalog = api.readCatalogue() # Load catalogue
		config['Dynamic'] = api.readConfigurationFile() # Load dynamic configuration
		window = webview.create_window('Toblerone', 'assets/gui.html', js_api = api, width = 1000, height = 712, resizable = False) # Create GUI window
		webview.start(debug=True)
