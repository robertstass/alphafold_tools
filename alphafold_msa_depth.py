#!/usr/bin/env python


from Bio import SeqIO
from io import StringIO
from Bio.SeqIO.FastaIO import FastaWriter
import sys

import argparse
import json
import os

class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Set a specific msa depth for an alphafold3 input json file.")

        # Optional depth parameter with a default value of 10
        self.parser.add_argument("--depth", type=int, default=default_depth,
                                 help="The depth for the msa (default: %d)" % default_depth)

        # Positional argument for the input JSON file
        self.parser.add_argument("input_json", help="Path to the input JSON file")

        # Optional positional argument for the output JSON file
        self.parser.add_argument("output_json", nargs='?',
                                 help="Path to the output JSON file (defaults to <input_json>_msa_depth_<depth>.json)")
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


default_depth = 30







####################################################################################################################



class SingleLineFastaWriter(FastaWriter):
    def __init__(self, handle):
        super().__init__(handle, wrap=None)  # Disables line wrapping


def msa_depth(msa, depth):
    msa = msa.replace(r"\n", "\n")
    msa_handle = StringIO(msa)
    alignment = list(SeqIO.parse(msa_handle, "fasta"))
    l = len(alignment)
    if l > depth:
        depth = depth
    else:
        depth = l
        print("MSA input depth is %d. Leaving unchanged." % l)
    subset_alignment = alignment[:depth]
    fasta_string = StringIO()
    writer = SingleLineFastaWriter(fasta_string)
    writer.write_file(subset_alignment)
    escaped_fasta_string = fasta_string.getvalue()
    return escaped_fasta_string




def main(json_file, output_filepath, depth):
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)

        if 'sequences' in data and isinstance(data['sequences'], list):
            for item in data['sequences']:
                if 'protein' in item and isinstance(item['protein'], dict):
                    if 'unpairedMsa' in item['protein']:
                        item['protein']['unpairedMsa'] = msa_depth(item['protein']['unpairedMsa'], depth)
                    else:
                        print(f"Warning: 'unpairedMsa' key not found in a protein entry.")

                    if 'pairedMsa' in item['protein']:
                        item['protein']['pairedMsa'] = msa_depth(item['protein']['pairedMsa'], depth)
                    else:
                        print(f"Warning: 'pairedMsa' key not found in a protein entry.")
                else:
                    print(f"Warning: 'protein' key not found or is not a dictionary in a sequence entry.")
        else:
            print(f"Warning: 'sequences' key not found or is not a list in the input file.")

        with open(output_filepath, 'w') as f:
            json.dump(data, f, indent=4)  # Use indent=4 for pretty formatting

        print(f"Successfully read '{json_file}', modified 'unpairedMsa' and 'pairedMsa' to a depth of {depth}, and wrote to '{output_filepath}'.")

    except FileNotFoundError:
        print(f"Error: Input file '{json_file}' not found.")
    except json.JSONDecodeError as e:
        print(e)
        print(f"Error: Could not decode JSON from '{json_file}'. Please ensure it's a valid JSON file.")



if __name__ == "__main__":
    config = ArgumentParserConfig()
    args = config.parser.parse_args()

    input_json = args.input_json
    output_json = args.output_json
    depth = args.depth

    # Determine the output filename if not provided
    if output_json is None:
        base, ext = os.path.splitext(input_json)
        output_json = f"{base}_msa_depth_{depth}{ext}"

    main(input_json, output_json, depth)




