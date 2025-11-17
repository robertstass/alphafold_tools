#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import argparse
import os
import json

default_sequence_num = 0
append_string = '_mut'
default_mutation = "X"

class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Mutate a sequence.")

        self.parser.add_argument("--sequence_number", "-n", type=int, default=default_sequence_num,
                                 help="Which number sequence to use from an input file. (zero-based. default: %d)" % default_sequence_num)

        # accept as string so we can allow comma-separated lists like "5,10,20"
        self.parser.add_argument("--index", "-i", type=str, default="1",
                                 help="Index or comma-separated list of indices (1-based). Example: 5 or 5,10,20. Indices refer to original (pre-shift) coordinates.")

        # accept as string so we can allow comma-separated lists like "A,V,GG"
        self.parser.add_argument("--mutation", "-m", type=str, default=default_mutation,
                                 help=f"The new residue(s) or comma-separated list of residues. Default: {default_mutation}")

        self.parser.add_argument("--insert", action="store_true",
                                 help="Insert the mutated residues before the index rather than replacing values. (applies to all provided mutations)")

        # Positional argument for the input JSON or FASTA file
        self.parser.add_argument("input_sequence", help="Path to the input sequence file (fasta or an alphafold3 json).")

        # Optional positional argument for the output file
        self.parser.add_argument("output_sequence", nargs='?',
                                 help=f"Path to the output file (defaults to <input_sequence>{append_string}.ext)")

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


def parse_json(json_file):
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        print(f"Error: Input file '{json_file}' not found.")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(e)
        print(f"Error: Could not decode JSON from '{json_file}'. Please ensure it's a valid JSON file.")
        sys.exit(1)


def sequence_from_json(data, sequence_number):
    sequences = []
    if 'sequences' in data and isinstance(data['sequences'], list):
        for item in data['sequences']:
            if 'protein' in item and isinstance(item['protein'], dict):
                if 'sequence' in item['protein']:
                    sequences.append(item['protein']['sequence'])
                else:
                    print(f"Warning: 'sequence' key not found in a protein entry.")
            else:
                print(f"Warning: 'protein' key not found or is not a dictionary in a sequence entry.")
    else:
        print(f"Warning: 'sequences' key not found or is not a list in the input file.")
    sequences_len = len(sequences)
    if sequence_number < sequences_len:
        return sequences[sequence_number]
    else:
        raise ValueError(f"Only {sequences_len} sequences in input json. Cannot select sequence {sequence_number}.")


def replace_sequence_in_json(data, sequence_number, sequence):
    i = 0
    replaced = False
    if 'sequences' in data and isinstance(data['sequences'], list):
        for item in data['sequences']:
            if 'protein' in item and isinstance(item['protein'], dict):
                if 'sequence' in item['protein']:
                    if i == sequence_number:
                        item['protein']['sequence'] = sequence
                        replaced = True
                    i += 1
                else:
                    print(f"Warning: 'sequence' key not found in a protein entry.")
            else:
                print(f"Warning: 'protein' key not found or is not a dictionary in a sequence entry.")
    else:
        print(f"Warning: 'sequences' key not found or is not a list in the input file.")
    if replaced:
        return data
    else:
        raise ValueError(f"Only {i} sequences in input json. Cannot select sequence {sequence_number}.")


def mutate_sequence(sequence, index, mutation, insert):
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")

    if index < 1:
        raise ValueError("index must be 1 or greater (1-based)")

    seq_len = len(sequence)
    pos = index - 1  # convert to 0-based

    if insert:
        # can insert at pos == seq_len (append)
        if pos > seq_len:
            raise ValueError(f"Insert position {index} is beyond sequence length {seq_len}.")
        new_seq = sequence[:pos] + mutation + sequence[pos:]
    else:
        # replacement: pos must point to an existing residue
        if pos >= seq_len:
            raise ValueError(f"Replacement index {index} is beyond sequence length {seq_len}.")
        end = pos + len(mutation)
        new_seq = sequence[:pos] + mutation + sequence[end:]
    return new_seq


def _parse_index_list(index_arg):
    """
    Parse a string like "5" or "5,10,20" into a list of ints.
    Strips whitespace and validates integer conversion.
    """
    parts = [p.strip() for p in index_arg.split(',') if p.strip() != ""]
    if len(parts) == 0:
        raise ValueError("No indices supplied.")
    try:
        ints = [int(p) for p in parts]
    except ValueError:
        raise ValueError("Indices must be integers (comma-separated).")
    return ints


def _parse_mutation_list(mutation_arg):
    """
    Parse a string like "A" or "A,V,GG" into a list of mutation strings.
    Strips whitespace from each element; empty elements are not allowed.
    """
    parts = [p.strip() for p in mutation_arg.split(',')]
    if len(parts) == 0 or any(p == "" for p in parts):
        raise ValueError("Invalid mutation list provided (empty element found).")
    return parts


def main(sequence_number, index_arg, mutation_arg, insert, input_sequence, output_sequence):

    index_list = _parse_index_list(index_arg)
    mutation_list = _parse_mutation_list(mutation_arg)

    if len(index_list) != len(mutation_list):
        print("Error: --index and --mutation must contain the same number of items.")
        print(f"Indices provided: {len(index_list)}; Mutations provided: {len(mutation_list)}")
        sys.exit(2)

    input_base, input_ext = os.path.splitext(input_sequence)

    if input_ext.lower() == ".json":
        json_input = True
        data = parse_json(input_sequence)
        sequence = sequence_from_json(data, sequence_number)
    else:
        json_input = False
        try:
            records = list(SeqIO.parse(input_sequence, "fasta"))
        except FileNotFoundError:
            print(f"Error: Input file '{input_sequence}' not found.")
            sys.exit(1)
        if len(records) == 0:
            print(f"Error: No FASTA records found in '{input_sequence}'.")
            sys.exit(1)
        if sequence_number < 0 or sequence_number >= len(records):
            raise ValueError(f"Only {len(records)} sequences in input fasta. Cannot select sequence {sequence_number}.")
        record = records[sequence_number]
        sequence = str(record.seq)

    print(f'Applying mutation(s) to {input_sequence}')

    # Build list of mutations with original order index so we can sort deterministically.
    # Each item: (index (int), mutation (str), orig_pos (int))
    edits = [(idx, mut, orig_i) for orig_i, (idx, mut) in enumerate(zip(index_list, mutation_list))]

    # Sort edits by index descending so earlier edits (higher indices) don't shift later ones.
    # For equal indices we sort by orig_pos descending so that the final order in the sequence
    # matches the order the user supplied on the command line.
    edits_sorted = sorted(edits, key=lambda t: (-t[0], -t[2]))

    current_seq = sequence
    for idx, mut, _orig in edits_sorted:
        if idx < 1:
            print(f"Error: index values must be 1 or greater. Got: {idx}")
            sys.exit(2)
        try:
            current_seq = mutate_sequence(current_seq, idx, mut, insert)
        except Exception as e:
            print(f"Error applying mutation at index {idx} with mutation '{mut}': {e}")
            sys.exit(1)

    mutated_sequence = current_seq

    if json_input:
        # replace in JSON structure and write JSON
        data = replace_sequence_in_json(data, sequence_number, mutated_sequence)
        try:
            with open(output_sequence, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            print(f"Error writing JSON to '{output_sequence}': {e}")
            sys.exit(1)
    else:
        # Modify fasta record and write all records back out
        try:
            record.seq = Seq(mutated_sequence)
        except Exception:
            record.seq = mutated_sequence
        try:
            SeqIO.write(records, output_sequence, "fasta")
        except Exception as e:
            print(f"Error writing FASTA to '{output_sequence}': {e}")
            sys.exit(1)

    print(f"Successfully written to '{output_sequence}'.")


if __name__ == "__main__":
    config = ArgumentParserConfig()
    args = config.parser.parse_args()

    sequence_number = args.sequence_number
    index = args.index
    mutation = args.mutation
    insert = args.insert
    input_sequence = args.input_sequence
    output_sequence = args.output_sequence

    # Determine the output filename if not provided
    if output_sequence is None:
        base, ext = os.path.splitext(input_sequence)
        output_sequence = f"{base}{append_string}{ext}"

    main(sequence_number, index, mutation, insert, input_sequence, output_sequence)
