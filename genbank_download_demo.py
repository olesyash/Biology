from Bio import Entrez, SeqIO

Entrez.email = "jcecomputationalbiology@gmail.com"  # Provide an email address

# First, find out the available search fields in the database you're searching
# Note: the nucleotide database contains DNA/RNA, this is what I want here.
with Entrez.einfo(db="nucleotide") as handle:
    record = Entrez.read(handle)

for field in record["DbInfo"]["FieldList"]:
    print("%(Name)s, %(FullName)s, %(Description)s" % field)

# Returns a printout of search fields like this:
# ACCN, Accession, Accession number of sequence
# PACC, Primary Accession, Does not include retired secondary accessions
# GENE, Gene Name, Name of gene associated with sequence
# PROT, Protein Name, Name of protein associated with sequence
# ECNO, EC/RN Number, EC number for enzyme or CAS registry number
# PDAT, Publication Date, Date sequence added to GenBank

# Fetch the required accessions:
# orangutan: 'AF451972'
# chimp: 'AF176731'
# human: 'X90314'
# rettype="fasta","gb",
output_fasta_file_orangutan = 'AF451972.fasta'
with Entrez.efetch(db="nucleotide", id="AF451972",
                   rettype="fasta", retmode="text") as in_handle:

    # creating an output fasta file:
    with open(output_fasta_file_orangutan, "w") as out_handle:
        out_handle.write(in_handle.read())

    print("Saved file {}".format(output_fasta_file_orangutan))

# reading the fasta
print("Parsing...")
record = SeqIO.read(output_fasta_file_orangutan, "fasta")
print(record)

# id can be a list of accessions:
output_fasta_all = 'all.fasta'
with Entrez.efetch(db="nucleotide", id=["AF451972", "AF176731", "X90314"],
                   rettype="fasta", retmode="text") as in_handle:
    # creating an output fasta file:
    with open(output_fasta_all, "w") as out_handle:
        out_handle.write(in_handle.read())

print("Saved file {}".format(output_fasta_all))


# or, instead of creating a fasta file:
with Entrez.efetch(db="nucleotide", id="AF451972", retmode="xml") as handle:
    features = Entrez.read(handle)[0]

sequence = features['GBSeq_sequence']  # this is your sequence!
print(sequence)

# compare the two sequences fetched via the two ways
assert(str(record.seq).lower() == str(sequence).lower())

# some more info:
print(features['GBSeq_length'])
print(features['GBSeq_update-date'])
print(features['GBSeq_organism'])

# exercise for the class: use it with a list of ids

# Corona demo: get all of the accession numbers for the available
#SARS-CoV-2 sequences using the ORGN field (organism) to limit the search term
with Entrez.esearch(db="nucleotide", term="SARS-CoV-2[ORGN]",
                    idtype="acc", retmax="15") as handle:  # ACCESSION (ID)
    results = Entrez.read(handle)
accs = results["IdList"]  # save the list of accession numbers

print(len(accs), accs)

import os
#Use efetch to download records and save to disk
filename = "SARS-CoV-2.gb"
if not os.path.isfile(filename):
    # Downloading...
    with Entrez.efetch(db="nucleotide", id=accs, rettype="gb", retmode="text") \
            as net_handle:
        with open(filename, "w") as out_handle:
            out_handle.write(net_handle.read())
    print("Saved {}".format(filename))

exit()















with Entrez.efetch(db="nucleotide", id=["AF451972", "AF176731", "X90314"], retmode="xml") as handle:
    all = Entrez.read(handle)
    assert(len(all) == 3)

for i in range(len(all)):
    features = all[i]
    print(i, features['GBSeq_sequence'])
