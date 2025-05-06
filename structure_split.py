import os.path
import copy 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio.PDB import MMCIFParser, PDBParser, PDBIO, MMCIFIO
import sys



import argparse
import os


class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Split the chains of a mmcif structure into multiple files.")


        # Positional argument for the input JSON file
        self.parser.add_argument("input_structure", help="Path to the input structure file (mmcif).")

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




def main(structure_file):
    file_root = os.path.splitext(structure_file)[0]
    if structure_file.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    elif structure_file.endswith(".pdb"):
        raise ValueError('PDB files not currently supported. Use cif.')
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Use PDB or mmCIF.")

    structure = parser.get_structure("protein", structure_file)
    mmcif_dict = parser._mmcif_dict
    models = {}

    three_to_one = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
    }
    out_ext = "cif"

    # Available model numbers
    model_nums = set(mmcif_dict.get("_atom_site.pdbx_PDB_model_num", ["1"]))

    for model in structure:
        model_id = str(model.id + 1)  # model.id is 0-based; mmCIF uses 1-based

        for chain in model:
            chain_id = chain.id

            # Deep copy of the full dict
            chain_dict = copy.deepcopy(mmcif_dict)

            # Filter indices matching both model num and chain ID
            label_asym_ids = chain_dict["_atom_site.label_asym_id"]
            model_nums = chain_dict.get("_atom_site.pdbx_PDB_model_num", ["1"] * len(label_asym_ids))

            keep_indices = [
                i for i, (m, c) in enumerate(zip(model_nums, label_asym_ids))
                if m == model_id and c == chain_id
            ]

            for key in chain_dict:
                if key.startswith("_atom_site.") and isinstance(chain_dict[key], list):
                    chain_dict[key] = [chain_dict[key][i] for i in keep_indices]

            out_file = f"{os.path.splitext(structure_file)[0]}_model_{model_id}_chain_{chain_id}.cif"
            io = MMCIFIO()
            io.set_dict(chain_dict)
            io.save(out_file)
            print(f"Wrote model {model_id}, chain {chain_id} to {out_file}")



if __name__ == "__main__":
    config = ArgumentParserConfig()
    args = config.parser.parse_args()

    input_structure = args.input_structure

    main(input_structure)

