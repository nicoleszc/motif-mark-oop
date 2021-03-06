#!/usr/bin/env python

# Motif Mark Script- OOP version
# Nikki Szczepanski

# Import necessary modules
import argparse
import re
import cairo

# Argparse function
def get_args():
    """Define the arguments needed for the script."""
    parser = argparse.ArgumentParser(description='Given a fasta file of genes and a text file of motifs, outputs an vector image\
		of a gene and its motifs to scale.')
    parser.add_argument("-f", "--fasta", type=str, required=True, help='Absolute path to fasta file of genes, where introns are\
		in lower case and exons are in upper case.')
    parser.add_argument("-m", "--motifs", type=str, required=True, help='Absolute path to text file containing a list of motifs,\
		where each motif is on a new line. Accepts all IUPAC degenerate base symbols.')
    args = parser.parse_args()
    return parser.parse_args()
parse_args = get_args()

# Rename variables from argparse arguments
fasta_file = parse_args.fasta
motifs_file = parse_args.motifs

# Define the ouput image file to be created (with the same path as the input fasta file)
image_file = fasta_file.split(".fa")[0] + ".svg"

def get_motifs():
	"""Read in the motifs file to create a list of motifs named 'motif_list'."""
	motif_list = []
	with open(motifs_file, 'r') as motifs:
		for motif in motifs:
			motif = motif.strip("\n")
			motif_list.append(motif)
		# Sort the motifs from longest to shortest, which will aid in the visualization later on.
		motif_list=sorted(motif_list,key=len,reverse=True)
	return motif_list
motif_list = get_motifs()

def degenerate_bases():
    """
	Creates a dictionary in which the keys are degenerate bases and the values are the possible ATCGU base letters.
	This dictionary 'bases' will be used to create regex expressions for motifs including degenerate bases, according to IUPAC.
	"""
    bases = {
		'A': 'A',
		'C': 'C',
		'G': 'G',
		'T': '[TU]',
		'U': '[TU]',
		'W': '[ATU]',
		'S': '[CG]',
		'M': '[AC]',
		'R': '[AG]',
		'Y': '[CTU]',
		'K': '[GTU]',
		'B': '[CGTU]',
		'D': '[AGTU]',
		'H': '[ACTU]',
		'V': '[ACG]',
		'N': '[ACGTU]'}
    return bases
bases = degenerate_bases()

def convert_degenerates():
	"""
	If degenerate bases are present in a motif, create a regex expression that only uses the letters ACTGU, according to IUPAC, using the pre-defined dictionary 'bases'.
	Saves each motif:regex pair in a dictionary 'motif_convmot'. 
	Also creates a dictionary 'converted_motifs' where the key is the regex expression and the value is the length of the original motif.
	Note: Not case sensitive (i.e., can interpret both lower and upper case motifs).
	"""
	converted_motifs = {}
	# Iterate through each motif in the list of motifs
	for motif in motif_list:
		regex = ''
		# Iterate through each letter in the motif
		for char in motif:
			if char.upper() in bases.keys():
				# Add each letter's corresponding regex expression to a string 
				regex += bases[char.upper()]
			else:
				print("Error: Unrecognized base in motif. Please include only IUPAC symbols.")
		# Once the string of regex expressions is complete, set this as the value in the motif:regex dictionary
		converted_motifs[motif]=regex
	return converted_motifs
converted_motifs = convert_degenerates()

def color_scheme():
	"""
	From a list of pre-defined colors, associates each motif to a color.
	Creates a dictionary in which the keys are the motifs as regex expressions and the values are the colors in RGB notation.
	"""
	# Below is the list of default colors that will represent motifs:
	                #blue       magenta     orange      teal        gray        cyan
	color_list=[(0.1,0.4,0.7),(1,0.2,0.4),(1,0.5,0),(0,0.7,0.6),(0.8,0.8,0.8),(0,0.7,1)]
	color_dict={}
	i=0
	for motif in converted_motifs:
		color_dict[motif]=color_list[i]
		i+=1
	return color_dict
color_dict = color_scheme()

# Define global constants (uppercase)
HEIGHT_GENE_GROUP = 50
Y_OFFSET_GENE = 50
LEFT_MARGIN = 50

# Define dimensions of output image
width=1000
height=1000

surface = cairo.SVGSurface(image_file,width,height)
context = cairo.Context(surface)

# Set a white background (otherwise, background is transparent)
context.save()
context.set_source_rgb(1, 1, 1)
context.paint()
context.restore()

class FastaHeader:
	'''
	Stores fasta header string for each gene in fasta
	'''
	# Set fasta headers to be black
	def __init__(self, gene_number, fasta_header):
		self.gene_number = gene_number
		self.fasta_header = fasta_header

	def draw(self):
		context.set_source_rgb(0,0,0)
		y = HEIGHT_GENE_GROUP * self.gene_number + Y_OFFSET_GENE
		context.move_to(LEFT_MARGIN,y)
		context.show_text(self.fasta_header)

class Gene:
	'''
	Stores start (drawing coords) and length of each gene
	'''
	def __init__(self, gene_number, seq):
		self.gene_number = gene_number
		self.seq = seq

	def find_length(self):
		length = len(self.seq)
		return length

	def draw(self):
		context.set_source_rgb(0,0,0)
		y = HEIGHT_GENE_GROUP * self.gene_number + Y_OFFSET_GENE + 20
		context.set_line_width(2)
		context.move_to(LEFT_MARGIN,y)
		context.line_to(LEFT_MARGIN + self.find_length(),y)
		context.stroke()

class Exon:
	'''
	Stores start (drawing coords) and length of each exon for each gene
	'''
	def __init__(self, gene_number, seq):
		self.gene_number = gene_number
		self.seq = seq
	
	def draw(self):
		context.set_source_rgb(0,0,0)
		y = HEIGHT_GENE_GROUP * self.gene_number + Y_OFFSET_GENE + 20
		context.set_line_width(15)
		coords = [(c.span(0)) for c in re.finditer('[ATCGUN]+', self.seq)]
		if coords != []:
			start = coords[0][0]
			end = coords[0][1]
			context.move_to(LEFT_MARGIN + start,y)
			context.line_to(LEFT_MARGIN + end,y)
			context.stroke()

class Motif:
	'''
	Stores each motif start and length of motif
	'''
	def __init__(self, gene_number, seq):
		self.gene_number = gene_number
		self.seq = seq
		self.match_dict = self.find_coords()

	def find_coords(self):
		match_dict={}
		for motif,regex in converted_motifs.items():
			# Regex expression that allows for overlapping matches. Saves starting index for each match. Case insensitive.
			match = [m.start(0) for m in re.finditer(rf'(?=({regex}))', self.seq, re.IGNORECASE)]
			if match != []:
				match_dict[motif] = match
		return match_dict
		
	def draw(self):
		y = HEIGHT_GENE_GROUP * self.gene_number + Y_OFFSET_GENE + 20
		context.set_line_width(10)
		for motif, match in self.match_dict.items():
			i = 0
			R,G,B = color_dict[motif]
			while i < len(match):
				context.set_source_rgb(R,G,B)
				context.move_to(LEFT_MARGIN + match[i],y)
				context.line_to(LEFT_MARGIN + match[i] + len(motif),y)
				context.stroke()
				i += 1

class GeneGroup:
	'''
	Each GeneGroup object stores the FastaHeader, Gene, Exon and Motif objects.
	Tells each object to draw itself.
	'''
	# Responsible for telling its children (genes, motifs) to go draw themselves
	# Make this class as simple as possible when it comes to calculations
	# There is one GeneGroup object for each gene
	def __init__(self, FastaHeader_object, Gene_object, Exon_object, Motif_object):
		self.FastaHeader_object = FastaHeader_object
		self.Gene_object = Gene_object
		self.Exon_object = Exon_object
		self.Motif_object = Motif_object

	def draw_all(self):
		self.FastaHeader_object.draw()
		self.Gene_object.draw()
		self.Exon_object.draw()
		self.Motif_object.draw()

class Legend:
	'''
	Draws the legend according to the number of the motifs in the motif_list
	'''
	def __init__(self, gene_number):
		self.gene_number = gene_number

	def draw(self):
		# Create box for legend
		y = HEIGHT_GENE_GROUP * self.gene_number + Y_OFFSET_GENE + 50
		context.set_source_rgb(0,0,0)
		context.rectangle(LEFT_MARGIN,y,len(max(motif_list,key=len))+100,len(motif_list)*20)
		context.set_line_width(2)
		context.stroke()
		
		i=0
		for motif in motif_list:
			i+=1
			context.move_to(LEFT_MARGIN+30, y+i*16.5)
			context.set_source_rgb(0,0,0)
			context.show_text(motif)
			R,G,B = color_dict[motif]
			context.set_source_rgb(R,G,B)
			context.set_line_width(10)
			context.move_to(LEFT_MARGIN + 12, y+i*15)
			context.line_to(LEFT_MARGIN + 12 + len(motif), y+i*15)
			context.stroke()


def main_function():
	"""
	Outputs one svg of genes and motifs using Pycairo, where each gene is shown one after the other 
	in the same order they appear in the fasta file. The motifs are overlaid on each gene where the 
	sequences match. Inludes a legend that displays the motif color key.
	"""
	with open(fasta_file, 'r') as fasta:
		seq=''
		fasta_header=''
		gene_number = 1
		fasta_header = fasta.readline().split('>')[1]
		FastaHeader_object = FastaHeader(gene_number, fasta_header)
		# Gor through every line of the fasta file
		for line in fasta:
			line = line.strip('\n')
			# Identify headers (gene names) and remove '>' symbol
			if line.startswith('>'):
				Gene_object = Gene(gene_number, seq)
				Exon_object = Exon(gene_number, seq)
				Motif_object = Motif(gene_number, seq)
				GeneGroup_object = GeneGroup(FastaHeader_object,Gene_object,Exon_object,Motif_object)
				GeneGroup_object.draw_all()
				gene_number+=1
				fasta_header = line.split('>')[1]
				FastaHeader_object = FastaHeader(gene_number, fasta_header)
				seq=''
			else:
				seq+=line
		Gene_object = Gene(gene_number, seq)
		Exon_object = Exon(gene_number, seq)
		Motif_object = Motif(gene_number, seq)
		GeneGroup_object = GeneGroup(FastaHeader_object,Gene_object,Exon_object,Motif_object)
		GeneGroup_object.draw_all()

	Legend_object = Legend(gene_number)
	Legend_object.draw()
	return

main_function()
surface.finish()
