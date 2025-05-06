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
        self.parser = argparse.ArgumentParser(description="Copy the release date entry from one mmcif file to another.")

        self.parser.add_argument("--force_overwrite", action='store_true', help="Force overwite if release date info already exists.")
        # Positional argument for the input JSON file
        self.parser.add_argument("copy_from", help="mmcif file to copy release date from.")
        self.parser.add_argument("copy_to", help="mmcif file to copy release date to.")

        # Optional positional argument for the output JSON file
        self.parser.add_argument("output_to", nargs='?', help=f"Path to output to. Defaults to <input_sequence>{append_string}.cif)")
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


append_string = '_rd'


def read_structure(structure_file):
    if structure_file.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    elif structure_file.endswith(".pdb"):
        raise ValueError('PDB files not currently supported. Use mmcif.')
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Use mmCIF.")

    structure = parser.get_structure("protein", structure_file)
    mmcif_dict = parser._mmcif_dict
    return structure, mmcif_dict


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == "":
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")

def main(from_structure_file, to_structure_file, output_structure_file, force_overwrite=False):
    from_structure, from_mmcif_dict = read_structure(from_structure_file)
    to_structure, to_mmcif_dict = read_structure(to_structure_file)

    release_date_entry = '_pdbx_audit_revision_history.revision_date'
    revision_history_category_name = '_pdbx_audit_revision_history'

    revision_history_keys = [key for key in from_mmcif_dict.keys() if key.startswith(revision_history_category_name)]
    if not release_date_entry in from_mmcif_dict.keys():
        print('mmCIF must have release date. (%s)' % release_date_entry)
    else:
        print("Release date: %s" % from_mmcif_dict[release_date_entry])
    revision_history_dict = {k: from_mmcif_dict[k] for k in revision_history_keys if k in from_mmcif_dict}


    to_revision_history_keys = [key for key in to_mmcif_dict.keys() if key.startswith(revision_history_category_name)]

    overwrite = False
    if to_revision_history_keys != []:
        if not force_overwrite:
            overwrite = query_yes_no('Warning, structure already has revision history entries. Do you want to overwrite?', default='yes')
            if not overwrite:
                print('Quitting without copying release date.')
                return
            else:
                print('Overwriting existing revision history entries...')
        else:
            print('Warning: Overwriting existing revision history entries...')
            overwrite = True
        if overwrite:
            [to_mmcif_dict.pop(k, None) for k in to_revision_history_keys]
    else:
        overwrite = True
    if overwrite:
        to_mmcif_dict = {**to_mmcif_dict, **revision_history_dict}



        out_file = output_structure_file
        io = MMCIFIO()
        io.set_dict(to_mmcif_dict)
        io.save(out_file)
        print(f"Wrote structure to {out_file}")

if __name__ == "__main__":
    config = ArgumentParserConfig()
    args = config.parser.parse_args()

    copy_from = args.copy_from
    copy_to = args.copy_to
    output_to = args.output_to
    force_overwrite = args.force_overwrite

    # Determine the output filename if not provided
    if output_to is None:
        base, ext = os.path.splitext(copy_to)
        output_to = f"{base}{append_string}{ext}"

    main(copy_from, copy_to, output_to, force_overwrite=force_overwrite)