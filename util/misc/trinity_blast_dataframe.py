# Universitat Potsdam
# Author Gaurav Sablok
# date: 2024-2-22
import os
import arguably 
import pyprofilers as pp
import modin.pandas as pd
import gget
@pp.profile(sort_by='cumulative', out_lines=30) 
@pp.profile_by_line(exit=1) 
@pp.simple_timer(num=1)
snap = CodeSnap()
snap = CodeSnap(tracer="python") 
snap.start()
@arguably.Command
def trinityAnnotateDraw(transcripts_path = FALSE)
"""
a parallelly implemented python function for the API based blast analysis of
the transcripts and provides a faster and comprehensive approach
to the blast analysis and filter analysis based on the similarity
score. It takes a transcriptome assembly file and provides a
transcriptome wide dataframe as well as the analysis plus the 
accession and the sequences of the blast analysis and alignment.
It also provides a complete information on the taxID and also the 
blast sequences. It will also write and get the accessions identified
in the blast hists and also the sequences for the alignment. Give the assembled 
genome coding regions, or the transcripts or the genes and it will prepare the 
clean headers, transcript headers, blast analysis, accession analysis, taxonomy
identifiers analysis.
"""
start =[]
end = []
strand = []
length = []
names = []
if transcripts_path:
    transcript_final = os.path.join(os.getcwd(), transcripts_path)
fasta_transcript_dict = {}
read_transcripts = [i.strip() for i in open("transcripts_path", "r").readlines()]
for i in read_transcripts:
    if i.startswith(">"):
        path = i.strip()
        if i not in fasta_transcript_dict:
            fasta_transcript_dict[i] = ""
        continue
            fasta_transcript_dict[path] += i.strip()
fasta_sequences = list(fasta_transcript_dict.values())
fasta_names = list(fasta_transcript_dict.keys())
for i in range(len(fasta_names)):
    start.append(list(filter(lambda n: n !=  "-",re.findall(r'[0-9]{0,3}-[0-9]{0,6}',fasta_names[i])))[0].split("-")[0])
    end.append(list(filter(lambda n: n != "-",re.findall(r'[0-9]{0,3}-[0-9]{0,6}',fasta_names[i])))[0].split("-")[1])
    strand.append(fasta_names[i].split()[-1][-2])
    length.append(fasta_names[i].split()[4].split(":")[-1])
    names.append(fasta_names[i].split()[1].split("_")[1])
fasta_dataframe = {}
for i in range(len(names)):
    fasta_dataframe[names[i]] = [start[i],end[i],length[i],strand[i],fasta_sequences[i]]
genome_dataframe = pd.DataFrame(fasta_dataframe).T
genome_dataframe.rename(columns={0:"start", 1:"end", 2:"length", 3: "strand", 4: "sequences"}, inplace=True)
genome_blast = genome_dataframe["sequences"].to_list()
genome_blast_list = {}
for i in range(len(names)):
    genome_blast_list[names[i]] = gget.blast(genome_blast[i], wrap_text=True)[["Description", "Scientific Name", "Taxid", "E value",\
                                                                          "Accession"]].where(gget.blast(genome_blast[i], wrap_text=True)["Per. Ident"] == "100.00%").dropna()
accession_blast_list = {}
for i in range(len(names)):
    accession_blast_list[names[i]] = gget.blast(genome_blast[i], wrap_text=True)[["Description", "Scientific Name", "Taxid", "E value",\
                                                                   "Accession"]].where(gget.blast(genome_blast[i], wrap_text=True)["Per. Ident"] == "100.00%").dropna()["Accession"].to_list()
tax_blast_list = {}
for i in range(len(names)):
    tax_blast_list[names[i]] = gget.blast(genome_blast[i], wrap_text=True)[["Description", "Scientific Name", "Taxid", "E value",\
                                                                    "Accession"]].where(gget.blast(genome_blast[i], wrap_text=True)["Per. Ident"] == "100.00%").dropna()["Taxid"].to_list()                                                                    
blast_accession_ids = []
for i in range(len(list(genome_blast_list.values()))):
    blast_accession_ids.append(list(genome_blast_list.values())[i]["Accession"].to_list())
final_blast_accession_ids = [j for i in blast_accession_ids for j in i]

blast_accession_sequences = {}
for i in range(len(final_blast_accession_ids)):
    blast_accession_sequences[final_blast_accession_ids[i]] = gget.seq(final_blast_accession_ids[i], translate=True, isoforms=True)

with open(os.path.join(os.getcwd(),"taxinformation.txt"), "w") as taxonomy:
    taxonomy.write(f"the identified taxonomy identifiers are:")
    taxonomy.write("\n")
    taxonomy.write(f"{tax_blast_list}")
    taxonomy.close()
with open(os.path.join(os.getcwd(),"accession.txt"), "w") as accession:
    accession.write(f"the identified accession identifiers are:")
    accession.write("\n")
    accession.write(f"{accession_blast_list}")
    accession.close()
with open(os.path.join(os.getcwd(),"write_accession_sequences.txt"), "w") as accession_sequences:
    accession.write(f"the identified accession identifiers sequences are:")
    accession.write("\n")
    accession.write(f"{blast_accession_sequences.keys()}\t{blast_accession_sequences.values()}")
    accession.close()
snap.stop()
snap.save()
if __name__ == "__main__":
    arguably.run()
