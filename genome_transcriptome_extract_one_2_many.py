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
def genomeExtraction(alignmentgenome = FALSE, 
                          reference_genome = FALSE, 
                                 target_genome = FALSE,
                                    percent_match = FALSE):
    """
    Function: genomeExtraction
    Summary: this will take the genome alignment file in the format given below from 
    the lastz or the blast alignments and then will extract the reference and the target
    genome regions. It will extract from both the positive and the negative strand and in 
    the case of the negative strand it reverses the sequences. It also provides the option 
    for the filtering according to the given percentage match and also writes all the coordinates 
    and the sequences as tab delimited files. This functions takes the one target genome and as many
    as the reference genomes. A sample fasta and the alignment file is present for the same.  
    Examples:  
  	genome	query_size	aligned_start	aligned_end	matches	mismatches	%_aligned	%_matched	chromosome	strand	start	end
0	0	taeGut2	8000	1	8000	7788	60	100.00	97.35	1	+	56	100
1	1	taeGut2	8000	6295	7425	1095	22	14.14	96.82	1	+	5	10
3	3	taeGut2	8000	6351	6387	34	1	0.46	91.89	14	+	2	20
4	4	taeGut2	8000	6351	6383	31	2	0.41	93.94	2	+	7	72
5	5	taeGut2	8000	6357	6383	26	1	0.34	96.30	10	-	3	32
    Attributes: 
        @param (multilevelgenome) default=FALSE: 
        @param (splitgenome) default=FALSE: 
    """""""""
if alignmentgenome and reference_genome and target_genome:
    reference_genome = os.path.join(os.getcwd(), reference_genome)
    target_genome = os.path.join(os.getcwd(), target_genome)
    reference_genome_dict = {}
    read_reference_genome = [i.strip() for i in open("reference_genome", "r").readlines()]
    for i in read_reference_genome:
        if i.startswith(">"):
            path = i.strip()
            if i not in reference_genome_dict:
                reference_genome_dict[i] = ""
                continue
        reference_genome_dict[path] += i.strip()
reference_genome_sequences = list(reference_genome_dict.values())
reference_genome_dict_names = list(reference_genome_dict.keys())
target_genome_dict = {}
read_target_genome = [i.strip() for i in open("target_genome", "r").readlines()]
for i in read_target_genome:
    if i.startswith(">"):
        path = i.strip()
        if i not in target_genome_dict:
            target_genome_dict[i] = ""
        continue
    target_genome_dict[path] += i.strip()
target_genome_sequences = list(target_genome_dict.values())
target_genome_dict_names = list(target_genome_dict.keys())
alignment_read = os.path.join(os.getcwd(), alignmentgenome)
alignment_dataframe = pd.read_csv(alignment_read, sep = ",")
aligned_reference_genome = alignment_dataframe[["genome","aligned_start", "aligned_end"]]
reference_aligned_genome = aligned_dataframe[["chromosome","strand", "aligned_start", "aligned_end"]]
reference_aligned_positive_strand = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "+").dropna()
reference_aligned_negative_strand = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "-").dropna()
target_genome_start = aligned_dataframe[["genome","aligned_start", "aligned_end"]]["aligned_start"].to_list()
target_genome_end = aligned_dataframe[["genome","aligned_start", "aligned_end"]]["aligned_end"].to_list()
target_genome_extract_sequences = []
for i in range(len(target_genome_start)):
    target_genome_extract_sequences.append(target_genome_sequences[target_genome_start[i]:target_genome_end[i]])
base_reference_fasta_dict = {}
for k,v in reference_genome_dict.items():
    base_reference_fasta_dict[k.replace(">", "")] = v
seqnames,seq  = list(base_reference_fasta_dict.keys()),list(base_reference_fasta_dict.values())
reference_chrom = []
chrom_names = list(map(str,alignment_dataframe["chromosome"].to_list()))
chrom_start = alignment_dataframe["start"].to_list()
chrom_end = alignment_dataframe["end"].to_list()
for i in range(len(chrom_names)):
        reference_chrom.append([chrom_names[i],chrom_start[i],chrom_end[i]])
seq_tuple = [(i,j)for i,j in zip(seqnames,seq)]
seq_dict = []
for i in range(len(chrom_names)):
    for j in range(len(seq_tuple)):
        if chrom_names[i] == seq_tuple[j][0]:
            seq_dict.append([chrom_names[i],seq_tuple[j][1][chrom_start[i]:chrom_end[i]]]) 
if alignmentgenome and reference_genome and target_genome and percent_match:
    percent_filtering = percent_match
    reference_genome = os.path.join(os.getcwd(), reference_genome)
    target_genome = os.path.join(os.getcwd(), target_genome)
    reference_genome_dict = {}
    read_reference_genome = [i.strip() for i in open("reference_genome", "r").readlines()]
    for i in read_reference_genome:
        if i.startswith(">"):
            path = i.strip()
            if i not in reference_genome_dict:
                reference_genome_dict[i] = ""
                continue
        reference_genome_dict[path] += i.strip()
reference_genome_sequences = list(reference_genome_dict.values())
reference_genome_dict_names = list(reference_genome_dict.keys())
target_genome_dict = {}
read_target_genome = [i.strip() for i in open("target_genome", "r").readlines()]
for i in read_target_genome:
    if i.startswith(">"):
        path = i.strip()
        if i not in target_genome_dict:
            target_genome_dict[i] = ""
        continue
    target_genome_dict[path] += i.strip()
target_genome_sequences = list(reference_genome_dict.values())
target_genome_dict_names = list(reference_genome_dict.keys())
alignment_read = os.path.join(os.getcwd(), alignmentgenome)
alignment_dataframe = pd.read_csv(alignment_read, sep = ",")[pd.read_csv(alignment_read, sep = ",")["%_matched"] > percent_filtering]
aligned_reference_genome = alignment_dataframe[["genome","aligned_start", "aligned_end"]]
reference_aligned_genome = aligned_dataframe[["chromosome","strand", "aligned_start", "aligned_end"]]
reference_aligned_positive_strand = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "+").dropna()
reference_aligned_negative_strand = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "-").dropna()
target_genome_start = aligned_dataframe[["genome","aligned_start", "aligned_end"]]["aligned_start"].to_list()
target_genome_end = aligned_dataframe[["genome","aligned_start", "aligned_end"]]["aligned_end"].to_list()
target_genome_extract_sequences = []
for i in range(len(target_genome_start)):
    target_genome_extract_sequences.append(target_genome_sequences[target_genome_start[i]:target_genome_end[i]])
base_reference_fasta_dict = {}
for k,v in reference_genome_dict.items():
    base_reference_fasta_dict[k.replace(">", "")] = v
seqnames,seq  = list(base_reference_fasta_dict.keys()),list(base_reference_fasta_dict.values())
reference_chrom = []
chrom_names = list(map(str,alignment_dataframe["chromosome"].to_list()))
chrom_start = alignment_dataframe["start"].to_list()
chrom_end = alignment_dataframe["end"].to_list()
for i in range(len(chrom_names)):
        reference_chrom.append([chrom_names[i],chrom_start[i],chrom_end[i]])
seq_tuple = [(i,j)for i,j in zip(seqnames,seq)]
seq_dict = []
for i in range(len(chrom_names)):
    for j in range(len(seq_tuple)):
        if chrom_names[i] == seq_tuple[j][0]:
            seq_dict.append([chrom_names[i],seq_tuple[j][1][chrom_start[i]:chrom_end[i]]]) 
with open("genome_aligned_target.txt", "w") as target:
    target.write(f"{target_genome_start}\t{target_genome_end}\t{target_genome_extract_sequences}")
with open("target_aligned_multiple_reference.txt", "w") as reference:
    reference.write("The multiple reference sequences along with the sequence names are")
    reference.write("\n")
    reference.write(f"seq_dict")
snap.stop()
snap.save()
if __name__ == "__main__":
    arguably.run()
