#!/usr/bin/env python

import sys

import argparse
import os
import json

class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Copy an msa from one alphafold json file to another.")

        # Optional depth parameter with a default value of 10
        self.parser.add_argument("--from_sequence_number", "-fn", type=int, default=0,
                                 help="Which sequence to use from an input json file. (zero-based. default: %d)" % default_from_sequence_num)
        self.parser.add_argument("--to_sequence_number", "-tn", type=int, default=0,
                                 help="Which sequence to use from an input json file. (zero-based. default: %d)" % default_to_sequence_num)

        # Positional argument for the input JSON file
        self.parser.add_argument("copy_from", help="Path to the input json file to copy msa from.")
        self.parser.add_argument("copy_to", help="Path to the json file to copy msa to.")
        # Optional positional argument for the output JSON file
        self.parser.add_argument("output_json", nargs='?',
                                 help=f"Path to the output JSON file (defaults to <copy_to>{append_string}.json)")
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



default_from_sequence_num = 0
default_to_sequence_num = 0
append_string = '_msa_copy'

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


def msa_from_json(data, sequence_number):
    msas = []
    if 'sequences' in data and isinstance(data['sequences'], list):
        for item in data['sequences']:
            if 'protein' in item and isinstance(item['protein'], dict):
                if 'unpairedMsa' in item['protein']:
                    msa = (item['protein']['unpairedMsa'], [], "")
                    if 'pairedMsa' in item['protein']:
                        msa = (msa[0], item['protein']['unpairedMsa'], "")
                    if 'sequence' in item['protein']:
                        msa = (msa[0], msa[1], item['protein']['sequence'])
                    msas.append(msa)
                else:
                    print(f"Warning: 'unpairedMsa' key not found in a protein entry.")
            else:
                print(f"Warning: 'protein' key not found or is not a dictionary in a sequence entry.")
    else:
        print(f"Warning: 'sequences' key not found or is not a list in the input file.")
    msas_len = len(msas)
    if sequence_number < msas_len:
        return msas[sequence_number]
    else:
        raise ValueError(f"Only {msas_len} msas in input json. Cannot select sequence {sequence_number}.")




def msa_to_json(data, msa, sequence_number):
    sequence_count = 0
    msa_set = False
    if 'sequences' in data and isinstance(data['sequences'], list):
        for item in data['sequences']:
            if 'protein' in item and isinstance(item['protein'], dict):
                if sequence_count == sequence_number:
                    item['protein']['unpairedMsa'] = msa[0]
                    item['protein']['pairedMsa'] = msa[1]
                    if msa[2] != "":
                        item['protein']['sequence'] = msa[2]
                    msa_set = True
                sequence_count += 1
            else:
                print(f"Warning: 'protein' key not found or is not a dictionary in a sequence entry.")
    else:
        print(f"Warning: 'sequences' key not found or is not a list in the input file.")
    if not msa_set:
        raise ValueError('Could not set msa to sequence number %d in copy_to json.' % sequence_number)
    return data



def main(copy_from, copy_to, output_json, from_sequence_number, to_sequence_number):

    from_json = parse_json(copy_from)
    to_json = parse_json(copy_to)

    msa = msa_from_json(from_json, from_sequence_number)

    output_json_data = msa_to_json(to_json, msa, to_sequence_number)

    with open(output_json, 'w') as f:
        json.dump(output_json_data, f, indent=4)  # Use indent=4 for pretty formatting
    print(f"Successfully written to '{output_json}'.")


if __name__ == "__main__":
    config = ArgumentParserConfig()
    args = config.parser.parse_args()

    copy_from = args.copy_from
    copy_to = args.copy_to
    output_json = args.output_json
    from_sequence_number = args.from_sequence_number
    to_sequence_number = args.to_sequence_number

    # Determine the output filename if not provided
    if output_json is None:
        base, ext = os.path.splitext(copy_to)
        output_json = f"{base}{append_string}.json"

    main(copy_from, copy_to, output_json, from_sequence_number, to_sequence_number)