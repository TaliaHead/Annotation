# Annotation
Scripts for "quick and dirty" functional annotation of assembled RNAseq databases

For transcriptomic data without a reference genome, this is used to get a "quick" functional annotation of your fasta file(s) based on UniProt BLAST results and InterProScan results.

GeneTransMap.py:
  Purpose - Create "Gene ID" for transcript/contig IDs in your fasta
  Required Input- a fasta file with your transcriptomic data named by the Evidential Gene Pipeline with contig names that are in the format 'Name_EVm#####t#" where 't#' denotes an isoform of the gene 'EVm#####'.
  Function- This uses your fasta file to create a "gene to transcript map", like the one used by Trinotate, without having a reference genome. The first column in the output corresponds to the transcript ID and the second column is the "gene ID"
  Usage- Replace path to fasta and path to output in the script before running


  Gene_Annotation.py:
  Purpose- Create a " functional annotation" for RNAseq fasta files without a reference
  Required Inputs
    1. Fasta file (nucleotide or AA)
    2. Blast results in tabular format (.tsv or .tabular) from UniProt
    3. GeneTransMap output (could edit the code to just create that here if you don't need a separate GeneTransMap- I was just trying a work around for Trinotate)
    4. InterProScan results in tabular format (.tsv or .tabular) 
  Usage-
    1. Check that your InterProScan output format matches that of interpro_cols in the script
      *Right now the code only uses the data from Pfam, SMART, CDD, SUPERFAMILY, ProSiteProfiles, PANTHER, InterPro, and GO terms. Other columns are ignored, change the name of the column and corresponding output lines if you wish to add/delete 
    2. Set your file paths & run
  Output- A tabular sheet with the following columns:
  Transcript ID/Gene ID/Accession/Gene name/Description/Species/CDD/PANTHER/Pfam/ProSiteProfiles/SMART/SUPERFAMILY/descriptionCDD	/descriptionPANTHER/descriptionPfam/descriptionProSiteProfiles/descriptionSMART/descriptionSUPERFAMILY/InterPro Entry/InterPro Description/GO terms
