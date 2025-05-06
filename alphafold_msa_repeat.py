#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
from Bio.SeqIO.FastaIO import FastaWriter
import sys
import argparse
import json
import os

class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Repeat the msa with a linker between for an alphafold3 input json file. (This avoids having to rerun an msa that has already been done.")

        # Optional depth parameter with a default value of 10
        self.parser.add_argument("--linker", type=str, default=default_linker,
                                 help="Amino acid sequence to fit between copies (default: %s)" % default_linker)
        self.parser.add_argument("--copies", type=int, default=default_copies,
                                 help="Number of copies (default: %s)" % default_copies)

        # Positional argument for the input JSON file
        self.parser.add_argument("input_json", help="Path to the input JSON file")

        # Optional positional argument for the output JSON file
        self.parser.add_argument("output_json", nargs='?',
                                 help="Path to the output JSON file (defaults to <input_json>_msa_repeated.json)")
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


default_linker = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" # 32xG
default_copies = 2






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

def msa_repeat(msa, depth, copies):
    msa = msa.replace(r"\n", "\n")
    msa_handle = StringIO(msa)
    alignment = list(SeqIO.parse(msa_handle, "fasta"))
    query = alignment[0]
    query_len = len(query)
    linker_len = len(linker)
    new_query_len = (query_len*copies)+(linker_len*(copies-1))
    new_query_str = ((str(query.seq)+linker)*(copies-1))+str(query.seq)
    query.seq = Seq(new_query_str)
    linker_gap = '-'*linker_len
    subset_alignment = [query]
    for ali in alignment[1:]:
        ali_str = ((str(ali.seq)+linker_gap)*(copies-1))+str(ali.seq)
        ali.seq = Seq(ali_str)
        subset_alignment.append(ali)
    fasta_string = StringIO()
    writer = SingleLineFastaWriter(fasta_string)
    writer.write_file(subset_alignment)
    escaped_fasta_string = fasta_string.getvalue()
    return str(query.seq), escaped_fasta_string



def main(json_file, output_filepath, linker, copies):
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)

        if 'sequences' in data and isinstance(data['sequences'], list):
            for item in data['sequences']:
                if 'protein' in item and isinstance(item['protein'], dict):
                    if 'unpairedMsa' in item['protein']:
                        query, fasta = msa_repeat(item['protein']['unpairedMsa'], linker, copies)
                        item['protein']['unpairedMsa'] = fasta
                        item['protein']['sequence'] = query
                        ids = item['protein']['id']
                        if len(ids) > 1:
                            item['protein']['id'] = ids[0:int(len(ids)/copies)]
                        print("New query sequence:")
                        print(query)
                    else:
                        print(f"Warning: 'unpairedMsa' key not found in a protein entry.")

                    if 'pairedMsa' in item['protein']:
                        query, fasta = msa_repeat(item['protein']['pairedMsa'], linker, copies)
                        item['protein']['pairedMsa'] = fasta
                        #item['protein']['sequence'] = query
                    else:
                        print(f"Warning: 'pairedMsa' key not found in a protein entry.")
                else:
                    print(f"Warning: 'protein' key not found or is not a dictionary in a sequence entry.")
        else:
            print(f"Warning: 'sequences' key not found or is not a list in the input file.")

        with open(output_filepath, 'w') as f:
            json.dump(data, f, indent=4)  # Use indent=4 for pretty formatting

        print(f"Successfully read '{json_file}', modified 'unpairedMsa' and 'pairedMsa', and wrote to '{output_filepath}'.")

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
    linker = args.linker
    copies = args.copies

    # Determine the output filename if not provided
    if output_json is None:
        base, ext = os.path.splitext(input_json)
        output_json = f"{base}_msa_repeated{ext}"

    main(input_json, output_json, linker, copies)




