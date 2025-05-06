#!/usr/bin/env python


import os.path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio.PDB import MMCIFParser, PDBParser
import sys
import argparse
import os


class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Output a sequence from a structure file (mmcif/pdb).")


        # Positional argument for the input JSON file
        self.parser.add_argument("input_structure", help="Path to the input structure file (mmcif or pdb). Note only mmcif will work with alphafold.")

        # Optional positional argument for the output JSON file
        self.parser.add_argument("output_sequence", nargs='?',
                                 help=f"Path to the output sequence file (defaults to <input_structure>{append_string}.fasta)")
        if len(sys.argv) == 1:
            self.usage()
            sys.exit()


    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print("Error: " + '\n'.join(msgs))
        print(" ")
        sys.exit(2)

append_string = '_seq'



def extract_sequence(structure_file):
    if structure_file.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    elif structure_file.endswith(".pdb"):
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Use PDB or mmCIF.")
    
    structure = parser.get_structure("protein", structure_file)

    models = {}

    three_to_one = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
    }

    for model in structure:
        model_id = model.id
        chains = models[model_id] = {}
        for chain in model:
            chain_id = chain.id
            chains[chain_id] = []
            seq = []
            res_ids = []
            for residue in chain:
                if residue.id[0] == " ":  # Ignore heteroatoms
                    seq.append(residue.resname)
                    res_ids.append(residue.id[1])
            sequence = "".join(three_to_one.get(res, "X") for res in seq)
            chains[chain_id].append([sequence, res_ids])
    return models

def main(structure_file, fasta_file):
    filename = os.path.splitext(os.path.basename(structure_file))[0]



    models = extract_sequence(structure_file)
    sequences = []

    for model_id, model in models.items():
        for chain_id, chain in model.items():
            name = filename+'_'+str(model_id)+'_'+str(chain_id)
            sequence = chain[0][0]
            if len(sequence) > 0:
                sequences.append(SeqRecord(Seq(sequence), id=name, description=""))

    with open(fasta_file, "w") as fasta_file_handle:
        SeqIO.write(sequences, fasta_file_handle, "fasta")

    print("FASTA file (%s) written successfully!" % fasta_file)


if __name__ == "__main__":
    config = ArgumentParserConfig()
    args = config.parser.parse_args()

    input_structure = args.input_structure
    output_sequence = args.output_sequence

    # Determine the output filename if not provided
    if output_sequence is None:
        base, ext = os.path.splitext(input_structure)
        output_sequence = f"{base}{append_string}.fasta"

    main(input_structure, output_sequence)
