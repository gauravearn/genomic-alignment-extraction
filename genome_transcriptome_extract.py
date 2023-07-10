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
    and the sequences as tab delimited files. 
    Examples:  
        genome  query_size  aligned_start   aligned_end matches mismatches  %_aligned   %_matched   chromosome  strand  start   end
0   taeGut2 8000    1   8000    7788    60  100 97.35   chr1    +   5637454 5646796
1   taeGut2 8000    6295    7425    1095    22  14.14   96.82   chr1    +   5642077 5643357
2   taeGut2 8000    6332    6383    45  4   0.65    86.54   chr2    -   1365309 1365357
3   taeGut2 8000    6351    6387    34  1   0.46    91.89   chr14   +   2096236 2096628
4   taeGut2 8000    6351    6383    31  2   0.41    93.94   chr2    +   7208240 7208272
5   taeGut2 8000    6357    6383    26  1   0.34    96.3    chr10   -   3293331 3293357
6   taeGut2 8000    3498    3526    25  1   0.36    86.21   chr1A   -   71676009    71676034
7   taeGut2 8000    3602    3624    22  1   0.29    95.65   chr14   +   12695042    12695064
    Attributes: 
        @param (multilevelgenome) default=FALSE: InsertHere
        @param (splitgenome) default=FALSE: InsertHere
    Returns: InsertHere
    """""""""
if alignmentgenome and reference_genome and target_genome:
    reference_genome = os.path.join(os.getcwd(), multilevelgenome)
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
read_target_genome = [i.strip() for i in open("reference_genome", "r").readlines()]
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
reference_aligned_positive_strand_sequences_start = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "+").dropna()["start"].to_list()
reference_aligned_positive_strand_sequences_end = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "+").dropna()["end"].to_list()
reference_genome_positive_strand_extract_sequences = []
for i in range(len(reference_aligned_positive_strand_sequences_start)):
    reference_genome_positive_strand_extract_sequences.append(reference_genome_sequence[reference_genome_positive_strand_extract_sequences[i]]:
                                                              reference_genome_sequence[reference_genome_positive_strand_extract_sequences[i]])
reference_aligned_negative_strand_sequences_start = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "-").dropna()["start"].to_list()
reference_aligned_negative_strand_sequences_end = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "-").dropna()["end"].to_list()
reference_genome_negative_strand_extract_sequences = []
for i in range(len(reference_aligned_negative_strand_sequences_start)):
    reference_genome_negative_strand_extract_sequences.append(reference_genome_sequence[reference_genome_negative_strand_extract_sequences[i]]:
                                                              reference_genome_sequence[reference_genome_negative_strand_extract_sequences[i]])
reference_genome_negative_strand_extract_sequences_reverese = [reversed(reference_genome_negative_strand_extract_sequences[i]) for i in range(len(reference_genome_negative_strand_extract_sequences))]


if alignmentgenome and reference_genome and target_genome and percent_match:
    reference_genome = os.path.join(os.getcwd(), multilevelgenome)
    percent_filtering = percent_match
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
read_target_genome = [i.strip() for i in open("reference_genome", "r").readlines()]
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
reference_aligned_positive_strand_sequences_start = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "+").dropna()["start"].to_list()
reference_aligned_positive_strand_sequences_end = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "+").dropna()["end"].to_list()
reference_genome_positive_strand_extract_sequences = []
for i in range(len(reference_aligned_positive_strand_sequences_start)):
    reference_genome_positive_strand_extract_sequences.append(reference_genome_sequence[reference_genome_positive_strand_extract_sequences[i]]:
                                                              reference_genome_sequence[reference_genome_positive_strand_extract_sequences[i]])
reference_aligned_negative_strand_sequences_start = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "-").dropna()["start"].to_list()
reference_aligned_negative_strand_sequences_end = aligned_dataframe[["chromosome","strand", "start", "end"]].where(align["strand"] == "-").dropna()["end"].to_list()
reference_genome_negative_strand_extract_sequences = []
for i in range(len(reference_aligned_negative_strand_sequences_start)):
    reference_genome_negative_strand_extract_sequences.append(reference_genome_sequence[reference_genome_negative_strand_extract_sequences[i]]:
                                                              reference_genome_sequence[reference_genome_negative_strand_extract_sequences[i]])
reference_genome_negative_strand_extract_sequences_reverese = [reversed(reference_genome_negative_strand_extract_sequences[i]) for i in range(len(reference_genome_negative_strand_extract_sequences))]

with open("genome_aligned_target.txt", "w") as target:
    target.write(f"{target_genome_start}\t{target_genome_end}\t{target_genome_extract_sequences}")
    target.close()
with open("genome_aligned_reference_positive_strand,txt", "w") as positive:
    positive.write(f"{reference_aligned_positive_strand_sequences_start}\t{reference_aligned_positive_strand_sequences_start}\t{reference_aligned_positive_strand_sequences_start}")
    positive.close()
with open("genome_aligned_reference_negative_strand,txt", "w") as negative:
    negative.write(f"{reference_aligned_negative_strand_sequences_start}\t{reference_aligned_negative_strand_sequences_end}\t{reference_genome_negative_strand_extract_sequences}")
    negative.close()
snap.stop()
snap.save()
if __name__ == "__main__":
    arguably.run()
