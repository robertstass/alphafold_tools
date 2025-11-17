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
        self.parser = argparse.ArgumentParser(description="Combine the MSAs (with a linker between) for an alphafold3 input json file. (This avoids having to rerun an msa that has already been done).")

        # Optional depth parameter with a default value of 10
        self.parser.add_argument("--linker", type=str, default=default_linker,
                                 help="Amino acid sequence to fit between copies (default: %s)" % default_linker)
        self.parser.add_argument("--sequence_numbers", "-ns", type=str,
                                 help="Which protein sequences to combine .Zero based. (default: None = all)")

        # Positional argument for the input JSON file
        self.parser.add_argument("input_json", help="Path to the input JSON file")

        # Optional positional argument for the output JSON file
        self.parser.add_argument("output_json", nargs='?',
                                 help=f"Path to the output JSON file (defaults to <input_json>{append_string}.json)")
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


append_string="_msa_combined"



####################################################################################################################


def csv_to_int_list(csv_string, delim=","):
    try:
        return [int(item.strip()) for item in csv_string.split(delim)]
    except ValueError:
        raise ValueError("Input contains non-integer values.")

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

def msa_combine(msas, verbose=True):
    combined_seq = []
    combined_msa = []
    print(f'New starting residue numbers:') if verbose else None
    for i, msa in enumerate(msas):
        msa = msa.replace(r"\n", "\n")
        msa_handle = StringIO(msa)
        alignment = list(SeqIO.parse(msa_handle, "fasta"))
        query = alignment[0]
        query_len = len(query)
        linker_len = len(linker)
        combined_seq.append(str(query.seq))
        print(sum([len(seq) for seq in combined_seq[:i*2]])+1) if verbose else None
        combined_msa.append(alignment[1:])
        if i != len(msas)-1:
            combined_seq.append(linker)
            linker_gap = '-'*linker_len
            #combined_msa.append([Seq(linker_gap)])
            combined_msa.append(None)

    sequence = "".join(combined_seq)
    query.seq = Seq(sequence)
    new_alignment = [query]

    filled_msa = []
    for i,msa in enumerate(combined_msa):
        if msa != None:
            pre_fill_len = sum([len(seq) for seq in combined_seq[:i]])
            post_fill_len = sum([len(seq) for seq in combined_seq[i + 1:]])
            for seq in msa:
                seq.seq = '-'*pre_fill_len+seq.seq+'-'*post_fill_len
                filled_msa.append(seq)

    new_alignment = new_alignment+filled_msa

    fasta_string = StringIO()
    writer = SingleLineFastaWriter(fasta_string)
    writer.write_file(new_alignment)
    escaped_fasta_string = fasta_string.getvalue()
    return sequence, escaped_fasta_string



def main(json_file, output_filepath, linker, sequence_numbers):
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)


        sequence_count = 0
        sequence_index = 0
        sequences = []
        unpaired_msas = []
        paired_msas = []

        if 'sequences' in data and isinstance(data['sequences'], list):
            for item in data['sequences']:
                if 'protein' in item and isinstance(item['protein'], dict):
                    if sequence_numbers is None or sequence_count == sequence_numbers[sequence_index]:
                        if 'unpairedMsa' in item['protein']:
                            sequences.append(item['protein']['sequence'])
                            unpaired_msas.append(item['protein']['unpairedMsa'])
                        else:
                            print(f"Warning: 'unpairedMsa' key not found in a protein entry.")
                        if 'pairedMsa' in item['protein']:
                            paired_msas.append(item['protein']['pairedMsa'])
                        else:
                            print(f"Warning: 'pairedMsa' key not found in a protein entry.")
                        sequence_index+=1
                    sequence_count+=1

                else:
                    print(f"Warning: 'protein' key not found or is not a dictionary in a sequence entry.")
        else:
            print(f"Warning: 'sequences' key not found or is not a list in the input file.")

        sequence, combined_unpaired_msas = msa_combine(unpaired_msas, verbose=True)
        seq2, combined_paired_msas = msa_combine(paired_msas, verbose=False)

        sequence_count = 0
        sequence_index = 0
        sequences_to_remove = []

        if 'sequences' in data and isinstance(data['sequences'], list):
            for i,item in enumerate(data['sequences']):
                if 'protein' in item and isinstance(item['protein'], dict):
                    if sequence_numbers is None or sequence_count == sequence_numbers[sequence_index]:
                        if sequence_index == 0:
                            if 'unpairedMsa' in item['protein']:
                                item['protein']['sequence'] = sequence
                                item['protein']['unpairedMsa'] = combined_unpaired_msas
                            else:
                                print(f"Warning: 'unpairedMsa' key not found in a protein entry.")
                            if 'pairedMsa' in item['protein']:
                                item['protein']['pairedMsa'] = combined_paired_msas
                            else:
                                print(f"Warning: 'pairedMsa' key not found in a protein entry.")
                        else:
                            sequences_to_remove.append(i)
                        sequence_index += 1

                    sequence_count += 1

                else:
                    print(f"Warning: 'protein' key not found or is not a dictionary in a sequence entry.")
        else:
            print(f"Warning: 'sequences' key not found or is not a list in the input file.")

        if 'sequences' in data and isinstance(data['sequences'], list):
            data['sequences'] = [item for i, item in enumerate(data['sequences']) if i not in sequences_to_remove]

        if 'name' in data and isinstance(data['name'], str):
            data['name'] = data['name']+append_string

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
    sequence_numbers = args.sequence_numbers


    delim=","
    if sequence_numbers != None:
        sequence_numbers = csv_to_int_list(sequence_numbers, delim=delim)

    # Determine the output filename if not provided
    if output_json is None:
        base, ext = os.path.splitext(input_json)
        output_json = f"{base}{append_string}{ext}"

    main(input_json, output_json, linker, sequence_numbers)




