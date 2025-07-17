#!/usr/bin/env python

import sys
import string
import itertools

import argparse
import os
import json

class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Add N-linked glycans to a alphafold input json file.")

        # Optional depth parameter with a default value of 10
        self.parser.add_argument("--glycans", "-g", required=True, help="Numbered N-link glycan sites as a csv (1-based).")

        self.parser.add_argument("--sequence_number", "-n", type=int, default=default_sequence_num,
                                 help="Which sequence to use from an input json file. (zero-based. default: %d)" % default_sequence_num)


        self.parser.add_argument("--chain_ids", default=None,
                                 help=f"The chain ids to use for each glycan as a csv. Defaults to AA,BB,CC...")

        # Positional argument for the input JSON file
        self.parser.add_argument("input_json", help="Input alphafold json file.")

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



default_sequence_num = 0
append_string = '_glyc'


def csv_to_int_list(csv_string, delim=","):
    try:
        return [int(item.strip()) for item in csv_string.split(delim)]
    except ValueError:
        raise ValueError("Input contains non-integer values.")

def letter_iterator():
    # Define the base sequences
    uppercase = list(string.ascii_uppercase)  # A-Z
    uppercase = [x+x for x in uppercase] # AA, BB, CC
    lowercase = list(string.ascii_lowercase)  # a-z
    lowercase = [x + x for x in lowercase] # aa, bb, cc
    # Combine primary sequences
    yield from uppercase
    yield from lowercase
    # Optional: fallback pattern (A1, B1, ..., Z1, A2, ...)
    for suffix in itertools.count(1):
        for ch in uppercase:
            yield f"{ch}{suffix}"

def parse_json(json_file):
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)

    except FileNotFoundError:
        print(f"Error: Input file '{json_file}' not found.")
    except json.JSONDecodeError as e:
        print(e)
        print(f"Error: Could not decode JSON from '{json_file}'. Please ensure it's a valid JSON file.")
    return data


def sequence_from_json(data, sequence_number):
    sequences = []
    ids = []
    if 'name' in data and isinstance(data['name'], str):
        data['name'] = data['name']+append_string
    if 'sequences' in data and isinstance(data['sequences'], list):
        for item in data['sequences']:
            if 'protein' in item and isinstance(item['protein'], dict):
                if 'sequence' in item['protein']:
                    sequences.append(item['protein']['sequence'])
                    ids.append(item['protein']['id'])
                else:
                    print(f"Warning: 'sequence' key not found in a protein entry.")
            else:
                print(f"Warning: 'protein' key not found or is not a dictionary in a sequence entry.")
    else:
        print(f"Warning: 'sequences' key not found or is not a list in the input file.")
    sequences_len = len(sequences)
    if sequence_number < sequences_len:
        return sequences[sequence_number], ids[sequence_number]
    else:
        raise ValueError(f"Only {sequences_len} sequences in input json. Cannot select sequence {sequence_number}.")


def add_ligands_json_data(data, ligands_data, bonds_data):
    sequence_count = 0
    if 'sequences' in data and isinstance(data['sequences'], list):
        data['sequences'] += ligands_data
    else:
        print(f"Warning: 'sequences' key not found or is not a list in the input file.")
    if "bondedAtomPairs" in data and isinstance(data['bondedAtomPairs'], list):
        data["bondedAtomPairs"] += bonds_data
    else:
        data["bondedAtomPairs"] = bonds_data
    return data


def main(glycans, sequence_number, chain_ids, input_json, output_json):
    json_data = parse_json(input_json)
    sequence,ids = sequence_from_json(json_data, sequence_number)

    if chain_ids != None:
        chain_ids = [item.strip() for item in chain_ids.split(delim)]
        if len(chain_ids) != len(glycans)*len(ids):
            raise ValueError("Length of chain_ids must be the same as glycans (times number of symmetry copies). (Or set to None)")
    else:
        letters = letter_iterator()
        chain_ids = [next(letters) for i in range(0,len(glycans)*len(ids))]

    ligands = []
    bonds = []
    for i, glycan in enumerate(glycans):
        try:
            res = sequence[glycan - 1]
        except:
            raise ValueError("Glycan number outside of sequence residue range")
        if res != "N":
            raise ValueError(f"Residue number {glycan} is not N. (it is {res}).")
        for j, id in enumerate(ids):
            chain_id = chain_ids[i*len(ids)+j]
            ligands.append({"ligand": {"id": [chain_id], "ccdCodes": ["NAG"]}})
            bonds.append([[id, glycan, "ND2"],[chain_id, 1, "C1"]])

    output_json_data = add_ligands_json_data(json_data, ligands, bonds)
    with open(output_json, 'w') as f:
        json.dump(output_json_data, f, indent=4)  # Use indent=4 for pretty formatting
    print(f"Successfully written to '{output_json}'.")


if __name__ == "__main__":
    config = ArgumentParserConfig()
    args = config.parser.parse_args()

    glycans = args.glycans
    sequence_number = args.sequence_number
    chain_ids = args.chain_ids
    input_json = args.input_json
    output_json = args.output_json

    delim=","
    glycans = csv_to_int_list(glycans, delim=delim)


    # Determine the output filename if not provided
    if output_json is None:
        base, ext = os.path.splitext(input_json)
        output_json = f"{base}{append_string}.json"

    main(glycans, sequence_number, chain_ids, input_json, output_json)