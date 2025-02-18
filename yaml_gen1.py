import yaml
from Bio import SeqIO

def read_fasta_sequence(fasta_file):
    """
    Reads the first sequence from a FASTA file.
    """
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return str(record.seq)

# Read sequences from your FASTA files
atrx_seq = read_fasta_sequence("proteins/atrx.fasta")
idh1_seq = read_fasta_sequence("proteins/idh1.fasta")
p53_seq = read_fasta_sequence("proteins/p53.fasta")

# Vorasidenib SMILES string from your provided conversion
vorasidenib_smiles = "Clc1nc(c2nc(N[C@@H](C(F)(F)F)C)nc(N[C@@H](C(F)(F)F)C)n2)ccc1"

# Create the YAML data structure
data = {
    "version": 1,
    "sequences": [
        {
            "protein": {
                "id": ["A"],  # ATRX chain ID
                "sequence": atrx_seq,
                "msa": "empty"  # Will trigger auto-generation of the MSA
            }
        },
        {
            "protein": {
                "id": ["B"],  # IDH1 chain ID
                "sequence": idh1_seq,
                "msa": "empty"
            }
        },
        {
            "protein": {
                "id": ["C"],  # TP53 chain ID
                "sequence": p53_seq,
                "msa": "empty"
            }
        },
        {
            "ligand": {
                "id": ["L1"],  # Use a ligand-specific ID that won't trigger MSA lookup
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