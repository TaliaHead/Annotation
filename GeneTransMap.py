from pathlib import Path

def parse_fasta_to_tsv(fasta_file, output_tsv):
    with open(fasta_file, 'r') as fasta, open(output_tsv, 'w') as tsv:
        # Write the TSV file
        for line in fasta:
            if line.startswith(">"):  # Identify header lines
                header = line.strip()
                transcript_id = header.split()[0][1:]  # Remove '>' and extract ID
                gene_id = transcript_id.rsplit('t', 1)[0]  # Remove 't' and trailing number
                tsv.write(f"{gene_id}\t{transcript_id}\n")

# Example usage:
fasta_file = "/path/to.fasta.fasta"  # Replace with the path to your FASTA file
output_tsv = "/path/to/output.tsv"  # Output file path
parse_fasta_to_tsv(fasta_file, output_tsv)

print(f"Tabular data has been written to {output_tsv}")
