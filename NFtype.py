import sys
from Bio import SeqIO

# Define the motifs
NF_motif = "TGAAAACCTTTTTTCCATTTACA"
r2_motif = "ATGGTAAAAAAGGTGAAAACCTTTTTTT"
r1_motif = "ATGGTAAAAAAGGTGT"
r1_motif_2 = "ATGGTAAAAAGTGT"
m2_motif = "ATGGTAAAAAAGGTGAAAACCTTTTTT"
m1_motif = "CCATTTACAAAAAAT"
s8_prefix = "AACTCTGTGCAATGGAGCACACGGCTCGTG"

# Define a function to determine the type
def determine_type(sequence):
    # Combine multiple lines into one line and remove any trailing whitespace or tabs
    sequence = sequence.replace("\t", "").replace("\n", "").strip()
    
    # First, count and remove the R2 motifs
    r2_count = sequence.count(r2_motif)
    sequence_without_r2 = sequence.replace(r2_motif, "")
    
    # Then, count the M2 motifs in the sequence without R2
    m2_count = sequence_without_r2.count(m2_motif)
    m1_count = sequence.count(m1_motif)
    r1_count = sequence.count(r1_motif)
    r1_count_2 = sequence.count(r1_motif_2)
    
    s8_index = sequence.find(s8_prefix)
    s8_31st_base = sequence[s8_index + len(s8_prefix)] if s8_index != -1 else None
    
    # First check for Naegleria fowleri
    if (m1_count >= 1 and m2_count >= 1) or (NF_motif in sequence):
    #if m1_count >= 1 and m2_count >= 1:
        # Continue to type determination
        if r2_count == 0 and r1_count == 0 and m2_count == 1 and m1_count == 1 and s8_31st_base == "C":
            return "Type 1"
        elif r2_count == 0 and r1_count == 0 and m2_count == 1 and m1_count == 1 and s8_31st_base == "T":
            return "Type 2"
        elif r2_count == 1 and r1_count == 1 and m2_count == 1 and m1_count == 1 and s8_31st_base == "T":
            return "Type 3"
        elif r2_count == 1 and r1_count == 1 and m2_count == 1 and m1_count == 1 and s8_31st_base == "C":
            return "Type 4"
        elif r2_count == 1 and r1_count == 0 and r1_count_2 == 1 and m2_count == 1 and m1_count == 1 and s8_31st_base == "C":
            return "Type 5"
        elif r2_count == 2 and r1_count == 1 and m2_count == 1 and m1_count == 1 and s8_31st_base == "T":
            return "Type 6"
        elif r2_count == 3 and r1_count == 1 and m2_count == 1 and m1_count == 1 and s8_31st_base == "T":
            return "Type 7"
        elif r2_count == 2 and r1_count == 2 and m2_count == 1 and m1_count == 1 and s8_31st_base == "T":
            return "Type 8"
        else:
            # Ambiguous cases or partial matches
            possible_types = []
            if r1_count_2 == 1:
                possible_types.append("Type 5")
            if r2_count == 3:
                possible_types.append("Type 7")
            if r2_count == 2 and r1_count == 2:
                possible_types.append("Type 8")
            if s8_31st_base == "C":
                possible_types.append("Type 1 or Type 4 or Type 5")
            if s8_31st_base == "T":
                possible_types.append("Type 2 or Type 3 or Type 6 or Type 7 or Type 8")
            if possible_types:
                return f"Naegleria fowleri, Possible type: {' and '.join(possible_types)}"
            else:
                return "Naegleria fowleri, Type not determined"
    else:
        return "Not Naegleria fowleri"

# Define the main function
def main(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        result = determine_type(sequence)
        print(f"{record.id}\t{result}")

        
# Script entry point
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    main(fasta_file)
