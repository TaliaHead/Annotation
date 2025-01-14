import pandas as pd
import argparse
from Bio import SeqIO
import xml.etree.ElementTree as ET
import glob

# Define file paths
fasta_file = "path/to/fasta.fasta"  # Replace with the actual path
blast_file = "path/to/UniProt/Blast/results.tsv"  # Replace with the actual path
interproscan_files = glob.glob("path/to/InterProScan/results*.tabular")  # Update with your directory 
genetransmap_file = "path/to/GeneTransMap.tsv"  # Replace with the actual path
output_file = "path/to/output/annotation.tsv"

# Read GeneTransMap file
genetransmap = pd.read_csv(genetransmap_file, sep="\t", header=None, names=["Gene ID", "Transcript ID"])

# Parse the FASTA file
fasta_data = []
for record in SeqIO.parse(fasta_file, "fasta"):
    header = record.description.split()
    transcript_id = header[0]
    gene_id = transcript_id.split("_")[0]  # Adjust parsing as per header structure
    fasta_data.append({"Transcript ID": transcript_id, "Gene ID": gene_id})
fasta_df = pd.DataFrame(fasta_data)

# Merge FASTA and GeneTransMap data
annotation_df = pd.merge(fasta_df.drop(columns=["Gene ID"]), genetransmap, how="outer", on="Transcript ID")
annotation_df['Transcript ID'] = annotation_df['Transcript ID'].astype(str)
print("...Fasta and Gene map parsed")

# Load the BLAST results into a DataFrame
blastp_cols = ["Transcript ID", "sseqid", "pident", "evalue", "bitscore", "length", "sgi", "Accession", "Description", "unk"]
blast_init = pd.read_csv(blast_file, sep="\t", names=blastp_cols)

# Filter hits with uncharacterized proteins and sort by e-value and bitscore
filtered_hits = (
        blast_init[~blast_init["Description"].str.contains("Uncharacterized protein")]
        .sort_values(by=["Transcript ID", "evalue", "bitscore"], ascending=[True, True, False])
    )
# Select the best hit for each query
best_hits = {}
for transcript_id, group in filtered_hits.groupby("Transcript ID"):
    best_hit = group.iloc[0]
    if best_hit["pident"] >= 50:
    # Extract concise description (only text before "OS=")
        concise_description = best_hit["Description"].split("OS=")[0].strip()
        
    #Populate best_hits dictionary
        best_hits[transcript_id] = {
            "Transcript ID": transcript_id,
            "Accession": best_hit["Accession"],
            "Gene name": None,  # Default to None
            "Description": concise_description,
            "Species": None  # Default to None
        }

	# Parse Gene Name (GN) and Species (OS)
        description = best_hit["Description"]
        if "OS=" in description:
            species = " ".join(description.split("OS=")[1].split()[:2])
            best_hits[transcript_id]["Species"] = species
        if "GN=" in description:
            gene_name = description.split("GN=")[1].split()[0]
            best_hits[transcript_id]["Gene name"] = gene_name

# Convert best_hits dictionary to DataFrame
blastp_data = pd.DataFrame.from_dict(best_hits, orient="index")
print("...Blast hits parsed")

# Merge BLAST results with annotation_df
annotation_df = pd.merge(annotation_df, blastp_data, how="left", on="Transcript ID")
print("...Blast annotation merged")

# Parse multiple InterProScan tabular files
interpro_cols = [
    "Transcript ID", "Ignore1", "Ignore2", "Database Type", "Database ID", "Database Description",
    "Ignore3", "Ignore4", "Ignore5", "Ignore6", "Ignore7", "InterPro Entry", "InterPro Description", "GO terms", "Ignore8"
]

# Read all InterProScan files and concatenate them
interpro_list = []
for file in interproscan_files:
    interpro_df = pd.read_csv(file, sep="\t", names=interpro_cols, usecols=[
        "Transcript ID", "Database Type", "Database ID", "Database Description", "InterPro Entry", "InterPro Description", "GO terms"
    ])
    interpro_list.append(interpro_df)

# Combine all files into a single DataFrame
interpro_df = pd.concat(interpro_list, ignore_index=True)
print("...InterPro files combined")

# Replace '-' with NaN in relevant columns
interpro_df[["InterPro Entry", "InterPro Description", "GO terms"]] = (
    interpro_df[["InterPro Entry", "InterPro Description", "GO terms"]]
    .replace("-", pd.NA)
    .fillna("")
)
interpro_df["GO terms"] = interpro_df["GO terms"].str.replace("|", ",")
print("...Missing values replaced and GO terms formatted")

# Filter relevant Database Types
relevant_types = ["Pfam", "SMART", "CDD", "SUPERFAMILY", "ProSiteProfiles", "PANTHER"]
interpro_df = interpro_df[interpro_df["Database Type"].isin(relevant_types)]
print("...IPS files filtered based on database")

# Group InterProScan data by Transcript ID and Database Type for pivoting
interpro_grouped = (
    interpro_df.groupby(["Transcript ID", "Database Type"])
    .agg({
        "Database ID": lambda x: ",".join(sorted(set(x.dropna()))),
        "Database Description": lambda x: ",".join(sorted(set(x.dropna()))),
    })
    .reset_index()
)

# Pivot InterProScan data for easier merging and rename columns for clarity
interpro_pivot = (
    interpro_grouped.pivot(index="Transcript ID", columns="Database Type", values=["Database ID", "Database Description"])
)
interpro_pivot.columns = [
    col.replace("Database ID_", "").replace("Database Description_", " description") if "Database ID_" in col or "Database Description_" in col else col
    for col in ['_'.join(col).strip() for col in interpro_pivot.columns.values]
]
interpro_pivot.reset_index(inplace=True)
print("...InterPro files pivoted and columns renamed for clarity")

# Group InterPro Entry, InterPro Description, and GO terms separately
interpro_additional = (
    interpro_df.groupby("Transcript ID")
    .agg({
        "InterPro Entry": lambda x: ",".join(sorted(set(x.dropna()))),
        "InterPro Description": lambda x: ",".join(sorted(set(x.dropna()))),
        "GO terms": lambda x: ",".join(sorted(set(x.dropna()))),
    })
    .reset_index()
)
print("...InterPro and GO terms data grouped")

# Merge pivoted data with additional columns
final_interpro_df = pd.merge(interpro_pivot, interpro_additional, on="Transcript ID", how="outer")
columns_to_clean = ["InterPro Entry", "InterPro Description", "GO terms"]
for col in columns_to_clean:
    final_interpro_df[col] = final_interpro_df[col].str.lstrip(",")
print("...InterPro data cleaned and merged back together")

# Merge InterProScan data into annotation_df
annotation_df = pd.merge(annotation_df, final_interpro_df, how="left", on="Transcript ID")
print("...InterPro data merged with annotation")

# Save final output
annotation_df.to_csv(output_file, sep="\t", index=False)
print("ANNOTATION FILE SAVED :)")