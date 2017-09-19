#!/usr/bin/env python

import argparse
import math
import string
import os


class fasta_input(object):
	def __init__(self, seq):
		super(fasta_input, self).__init__()
		self.buffer_size = os.stat(seq).st_blksize
		self.file = open(seq, "r")
		self.contigs_start_positions = {}
		self.line_length = -1
		self.initialize_fasta_reader()

	def initialize_fasta_reader(self):
		contig_id = None
		buffer_offset = 0
		for line in self.file:
			buffer_offset += len(line)
			if line[0] == ">":
				contig_id = line.rstrip().split(" ")[0][1:]
				self.contigs_start_positions[contig_id] = buffer_offset  # no offset based on contig start index because there are that many N at the begining of the contig to compensate
			elif len(line) > self.line_length + 1:
				self.line_length = len(line.rstrip())

	def finalize(self):
		self.file.close()


def make_supertranscript(seq, contig, translation_table, blocks_out, fasta_out, transcripts_out, current_gene_id, starts, ends, transcript_ids, orientations):
	current_gene_id = current_gene_id.replace("\"", "")
	sstarts = sorted(starts.keys())
	sends = sorted(ends.keys())
	i = j = 0
	current_transcripts = set()
	blocks = []  # (start_position, end_position, [transcript_ids])

	last_position = -1
	while i < len(sstarts) and j < len(sends):
		if sstarts[i] == sends[j]:
			if last_position != sstarts[i]:
				blocks.append((last_position, sstarts[i] - 1, set(current_transcripts)))
			for k in starts[sstarts[i]]:
				current_transcripts.add(transcript_ids[k])
			blocks.append((sstarts[i], sstarts[i], set(current_transcripts)))
			for k in ends[sends[j]]:
				current_transcripts.remove(transcript_ids[k])
			last_position = sstarts[i] + 1
			i += 1
			j += 1
		elif sstarts[i] < sends[j]:
			if last_position > 0 and last_position != sstarts[i]:
				blocks.append((last_position, sstarts[i] - 1, set(current_transcripts)))
			for k in starts[sstarts[i]]:
				current_transcripts.add(transcript_ids[k])
			last_position = sstarts[i]
			i += 1
		elif sstarts[i] > sends[j]:
			blocks.append((last_position, sends[j], set(current_transcripts)))
			for k in ends[sends[j]]:
				current_transcripts.remove(transcript_ids[k])
			last_position = sends[j] + 1
			j += 1

	while j < len(sends):  # closing blocks
		blocks.append((last_position, sends[j], set(current_transcripts)))
		for k in ends[sends[j]]:
			current_transcripts.remove(transcript_ids[k])
		last_position = sends[j] + 1
		j += 1

	ori = set(orientations)
	if len(ori) > 1:
		exit("error, one gene contains exons in both orientations")

	current_position = 1
	ori = ori.pop()
	blocks_buffer = ""
	fasta_buffer = ""
	last_transcripts = None
	transcripts_to_output = {}
	next_opening = 1
	transcript_buffer = ""

	if ori == "+":
		offset = sstarts[0]
		fasta_header_exons = ">" + current_gene_id + " loc:" + contig + "|" + str(blocks[0][0]) + "-" + str(blocks[-1][1]) + "|+ "
		fasta_header_exons += "exons:" + str(blocks[0][0])
		fasta_header_segs = " segs:1"

		for block in blocks:  # (start_position, end_position, [transcript_ids])
			if len(block[2]) > 0:
				end_position = current_position + block[1] - block[0]
				if last_transcripts and block[2] != last_transcripts:
					blocks_buffer += current_gene_id + "\trtracklayer\texon\t" + str(next_opening) + "\t" + str(current_position - 1) + "\t.\t+\t.\tgene_id \"" + current_gene_id + "\"; transcript_id \"" + current_gene_id + "\"; ID \"" + current_gene_id + "\"\n"
					for transcript in last_transcripts - block[2]:
						transcripts_to_output[transcript].append(current_position - 1)
					for transcript in block[2] - last_transcripts:
						if transcript not in transcripts_to_output:
							transcripts_to_output[transcript] = [current_position]
						else:
							transcripts_to_output[transcript].append(current_position)
					last_transcripts = block[2]
					next_opening = current_position
				elif not last_transcripts:
					last_transcripts = block[2]
					for transcript in last_transcripts:
						transcripts_to_output[transcript] = [1]
				tmp_fasta = get_fasta(seq, contig, block[0], block[1])
				if not tmp_fasta:
					print "contig " + contig + " from annotation does not exist in the fasta file."
					return
				fasta_buffer += tmp_fasta
				current_position = end_position + 1
			else:  # gap
				offset += block[1] - block[0] + 1
				fasta_header_exons += "-" + str(block[0] - 1) + "," + str(block[1] + 1)  # -1 and +1 because the exons end and start at positions around the gap
				fasta_header_segs += "-" + str(current_position - 1) + "," + str(current_position)

		for transcript in last_transcripts:
			transcripts_to_output[transcript].append(current_position - 1)
		for transcript in sorted(transcripts_to_output):
			for i in xrange(0, len(transcripts_to_output[transcript]), 2):
				transcript_buffer += current_gene_id + "\tsuperTranscript\texon\t" + str(transcripts_to_output[transcript][i]) + "\t" + str(transcripts_to_output[transcript][i + 1]) + "\t.\t+\t.\tgene_id \"" + current_gene_id + "\"; transcript_id \"" + transcript + "\"\n"
		transcripts_out.write(transcript_buffer)

		blocks_buffer += current_gene_id + "\trtracklayer\texon\t" + str(next_opening) + "\t" + str(end_position) + "\t.\t+\t.\tgene_id \"" + current_gene_id + "\"; transcript_id \"" + current_gene_id + "\"; ID \"" + current_gene_id + "\"\n"
		blocks_out.write(blocks_buffer)

		fasta_header_segs += "-" + str(current_position - 1)
		write_fasta(fasta_out, fasta_header_exons + "-" + str(blocks[-1][1]), fasta_header_segs, fasta_buffer, seq.line_length)
	elif ori == "-":
		offset = sends[-1]
		fasta_header = ">" + current_gene_id + " loc:" + contig + "|" + str(blocks[0][0]) + "-" + str(blocks[-1][1]) + "|- "
		fasta_header_exons = ""
		fasta_header_segs = " segs:1"
		for i in xrange(len(blocks) - 1, -1, -1):
			block = blocks[i]  # (start_position, end_position, [transcript_ids])
			if len(block[2]) > 0:
				end_position = current_position + block[1] - block[0]
				if last_transcripts and block[2] != last_transcripts:
					blocks_buffer += current_gene_id + "\trtracklayer\texon\t" + str(next_opening) + "\t" + str(current_position - 1) + "\t.\t+\t.\tgene_id \"" + current_gene_id + "\"; transcript_id \"" + current_gene_id + "\"; ID \"" + current_gene_id + "\"\n"
					for transcript in last_transcripts - block[2]:
						transcripts_to_output[transcript].append(current_position - 1)
					for transcript in block[2] - last_transcripts:
						if transcript not in transcripts_to_output:
							transcripts_to_output[transcript] = [current_position]
						else:
							transcripts_to_output[transcript].append(current_position)
					last_transcripts = block[2]
					next_opening = current_position
				elif not last_transcripts:
					last_transcripts = block[2]
					for transcript in last_transcripts:
						transcripts_to_output[transcript] = [1]
				tmp_fasta = get_fasta(seq, contig, block[0], block[1])
				if not tmp_fasta:
					print "contig " + contig + " from annotation does not exist in the fasta file."
					return
				fasta_buffer = tmp_fasta + fasta_buffer
				current_position = end_position + 1
			else:  # gap
				offset += block[1] - block[0] + 1
				fasta_header_exons = "-" + str(block[0] - 1) + "," + str(block[1] + 1) + fasta_header_exons
				fasta_header_segs += "-" + str(current_position - 1) + "," + str(current_position)

		for transcript in last_transcripts:
			transcripts_to_output[transcript].append(current_position - 1)
		for transcript in sorted(transcripts_to_output):
			for i in xrange(0, len(transcripts_to_output[transcript]), 2):
				transcript_buffer += current_gene_id + "\tsuperTranscript\texon\t" + str(transcripts_to_output[transcript][i]) + "\t" + str(transcripts_to_output[transcript][i + 1]) + "\t.\t+\t.\tgene_id \"" + current_gene_id + "\"; transcript_id \"" + transcript + "\"\n"
		transcripts_out.write(transcript_buffer)

		blocks_buffer += current_gene_id + "\trtracklayer\texon\t" + str(next_opening) + "\t" + str(end_position) + "\t.\t+\t.\tgene_id \"" + current_gene_id + "\"; transcript_id \"" + current_gene_id + "\"; ID \"" + current_gene_id + "\"\n"
		blocks_out.write(blocks_buffer)

		fasta_header_exons = fasta_header + "exons:" + str(blocks[0][0]) + fasta_header_exons + "-" + str(blocks[-1][1])  # [:fasta_header_exons.rfind(",")]
		fasta_header_segs += "-" + str(current_position - 1)
		write_fasta(fasta_out, fasta_header_exons, fasta_header_segs, string.translate(fasta_buffer[::-1], translation_table), seq.line_length)
	else:
		exit("error, invalid orientation " + ori)


def get_fasta(seq, contig, start, end):
	offset = None
	if contig not in seq.contigs_start_positions:
		if len(seq.contigs_start_positions) > 1:
			print "\"" + contig + "\""
			# exit("Error, contig IDs do not match between annotation and fasta file.")
			return None
		else:
			offset = seq.contigs_start_positions.values()[0]
	else:
		offset = seq.contigs_start_positions[contig]

	offset += start - 1 + int(math.floor((start - 1) / seq.line_length))  # number of characters + new lines
	seq.file.seek(offset, 0)
	line_offset = (start - 1) % seq.line_length
	newlines_offset = int(math.floor((end - start + 1 + line_offset) / seq.line_length))
	sequence = seq.file.read(end - start + 1 + newlines_offset)
	sequence = sequence.replace("\n", "")
	return sequence


def write_fasta(fasta_out, exons, segs, sequence, line_length):
	line_length = 70
	fasta_out.write(exons + segs + "\n")
	for i in xrange(0, len(sequence), line_length):
		fasta_out.write(sequence[i:i + line_length] + "\n")


def main():
	usage = "usage: "
	parser = argparse.ArgumentParser(usage)
	parser.add_argument('-gff', dest="gff", required=True, type=str, help="Path to gff file.")
	parser.add_argument('-seq', dest="fasta", required=True, type=str, help="Path to fasta file.")
	parser.add_argument('-ver', dest="version", default="gff3", type=str, choices=["gff2", "gff3"], help="GFF version of annotation (default = gff3)")  # GFF3 or GTF as 2 different exclusive arguments instead?
	parser.add_argument('-o', dest="output", required=True, type=str, help="Name base and path for output")
	args = parser.parse_args()

	blocks_out = open(args.output + ".supertranscripts.blocks.gtf", "w")
	fasta_out = open(args.output + ".supertranscripts.fasta", "w")
	transcripts_out = open(args.output + ".supertranscripts.transcripts.gtf", "w")

	translation_table = string.maketrans("ATCG", "TAGC")
	seq = fasta_input(args.fasta)

	i = 0
	starts = {}
	ends = {}
	previous_contig = ""
	transcript_ids = []
	orientations = []
	current_transcript = ""
	current_gene_id = ""
	if args.version == "gff2":
		with open(args.gff, "r") as annot:
			for l in annot:
				line = l.rstrip().split("\t")
				if line[2] == "exon":
					fields = line[8].lstrip().split("; ")
					current_transcript = ""
					for f in fields:
						if "gene_id" in f:
							if f.split(" ")[1] != current_gene_id:
								if current_gene_id:
									make_supertranscript(seq, previous_contig, translation_table, blocks_out, fasta_out, transcripts_out, current_gene_id, starts, ends, transcript_ids, orientations)
								current_gene_id = f.split(" ")[1]
								i = 0
								starts = {}
								ends = {}
								previous_contig = line[0]
								transcript_ids = []
								orientations = []
						if "transcript_id" in f:
							current_transcript = f.split(" ")[1]
					if int(line[3]) not in starts:
						starts[int(line[3])] = []
					starts[int(line[3])].append(i)
					if int(line[4]) not in ends:
						ends[int(line[4])] = []
					ends[int(line[4])].append(i)
					i += 1
					transcript_ids.append(current_transcript.replace("\"", ""))
					orientations.append(line[6])
		make_supertranscript(seq, previous_contig, translation_table, blocks_out, fasta_out, transcripts_out, current_gene_id, starts, ends, transcript_ids, orientations)
	seq.finalize()


if __name__ == "__main__":
	main()
