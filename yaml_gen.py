import yaml
from Bio import SeqIO

def read_fasta_sequence(fasta_file):
    """
    Reads the first sequence from a FASTA file.
    """
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return str(record.seq)

# Read sequences from the FASTA files
atrx_seq = read_fasta_sequence("proteins/atrx.fasta")
idh1_seq = read_fasta_sequence("proteins/idh1.fasta")
p53_seq = read_fasta_sequence("proteins/p53.fasta")

# Define the Vorasidenib SMILES string (from your provided SDF conversion)
vorasidenib_smiles = "Clc1nc(c2nc(N[C@@H](C(F)(F)F)C)nc(N[C@@H](C(F)(F)F)C)n2)ccc1"

# Create the YAML data structure
data = {
    "version": 1,
    "sequences": [
        {
            "protein": {
                "id": ["ATRX"],
                "sequence": atrx_seq,
                "msa": "empty"  # Use auto-generated MSA via --use_msa_server flag
            }
        },
        {
            "protein": {
                "id": ["IDH1"],
                "sequence": idh1_seq,
                "msa": "empty"
            }
        },
        {
            "protein": {
                "id": ["P53"],
                "sequence": p53_seq,
                "msa": "empty"
            }
        },
        {
            "ligand": {
                "id": ["VORASIDENIB"],
                "smiles": vorasidenib_smiles
            }
        }
    ]
}

# Write the YAML data to a file named input.yaml
output_yaml = "input.yaml"
with open(output_yaml, "w") as outfile:
    yaml.dump(data, outfile, sort_keys=False)

print("YAML file created:", output_yaml)