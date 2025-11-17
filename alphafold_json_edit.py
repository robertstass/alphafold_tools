#!/usr/bin/env python
import copy
import sys
import itertools

import argparse
import os
import json

class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Edit top level json key values (for alphafold3 json files).")

        # Optional depth parameter with a default value of 10
        self.parser.add_argument("--key", required=True, help="JSON key to edit.")
        self.parser.add_argument("--value", required=True, help=f"Value to replace.")
        self.parser.add_argument("--type", default=default_type, choices=types.keys(), help="Apply a type to the value.")
        self.parser.add_argument("--secondary_type", default=default_type, choices=types.keys(), help="Secondary type. Eg to apply to values in a list or values of dict.")
        self.parser.add_argument("--append", action="store_true", help="Append the value to the string/list.")


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

types = {"str": str, "list": list, "dict": dict, "int": int, "float": float}
default_type = "str"

append_string = '_edit'



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




def main(key, value, data_type, secondary_type, append, input_json, output_json):
    json_data = parse_json(input_json)

    delim = ','
    if data_type == "list":
        value = value.split(delim)
        value = [types[secondary_type](x) for x in value]
    elif data_type == "dict":
        value = dict(itertools.zip_longest(*[iter(value.split(delim))] * 2, fillvalue=""))
        value2 = copy.deepcopy(value)
        for k, v in value.items():
            value2[k] = types[secondary_type](v)
        value = value2
    elif data_type == "str":
        pass
    else:
        value = types[data_type](value)

    if key in json_data:
        if append:
            assert type(json_data[key]) == types[data_type]
            if data_type == "list" or data_type == "str":
                json_data[key] += value
            elif data_type == "dict":
                json_data[key] = json_data[key] | value
        else:
            json_data[key] = value
    else:
        raise ValueError(f"Key f{key} not in json file.")
    print(f'New value for key {key}:')
    print(json_data[key])

    with open(output_json, 'w') as f:
        json.dump(json_data, f, indent=4)  # Use indent=4 for pretty formatting
    print(f"Successfully written to '{output_json}'.")


if __name__ == "__main__":
    config = ArgumentParserConfig()
    args = config.parser.parse_args()

    key = args.key
    value =args.value
    data_type = args.type
    secondary_type = args.secondary_type
    append = args.append

    input_json = args.input_json
    output_json = args.output_json

    # Determine the output filename if not provided
    if output_json is None:
        base, ext = os.path.splitext(input_json)
        output_json = f"{base}{append_string}.json"

    main(key, value, data_type, secondary_type, append, input_json, output_json)

