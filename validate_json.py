#!/usr/bin/env python
import sys
import argparse
import json
import os

class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Validate that a file is readable json.")

        # Positional argument for the input JSON file
        self.parser.add_argument("input_json", help="Path to the input JSON file")

        # Optional positional argument for the output JSON file
        self.parser.add_argument("output_json", nargs='?',
                                 help="(Optional) Path to the output JSON file if ok (it will be unchanged).")
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










####################################################################################################################




def main(json_file, output_filepath):
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)

        if output_filepath is not None:
            with open(output_filepath, 'w') as f2:
                json.dump(data, f2, indent=4)  # Use indent=4 for pretty formatting
                print(f"Validated '{json_file}' and rewritten to '{output_filepath}'.")
        else:
            print(f"Validated '{json_file}'. It is valid json.")

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


    main(input_json, output_json)




