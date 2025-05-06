from Bio import SeqIO
from Bio import Align, AlignIO
from Bio.PDB import MMCIFParser, PDBParser
import sys


import argparse
import os
import json

class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Align a structural template to a sequence to supply to alphafold.")

        # Optional depth parameter with a default value of 10
        self.parser.add_argument("--sequence_number", "-n", type=int, default=default_sequence_num,
                                 help="Which sequence to use from an input json file. (ignored for fasta. zero-based. default: %d)" % default_sequence_num)

        self.parser.add_argument("--copies", "-c", type=int, default=1,
                                 help="Number of copies of the structure sequence to align to the query sequence. Separate templates will be added to the output json. default: 1")
        self.parser.add_argument("--alignment_scoring", default=default_alignment_scoring,
                                 help=f"The scoring type used by biopython to align the sequences. Default: {default_alignment_scoring}")
        self.parser.add_argument("--alignment_type", default=default_alignment_type,
                                 help=f"Instead of a query sequence, an alignment file can be supplied, with the query sequence aligned to the structure sequence. This allows easier customisation if the script's built in alignment doesn't work correctly. Enter the type of alignment file input here (if supplied). Default: {default_alignment_type}")

        # Positional argument for the input JSON file
        self.parser.add_argument("input_structure", help="Path to the input structure file (mmcif or pdb). Note only mmcif will work with alphafold")
        self.parser.add_argument("input_sequence", help="Path to the input query sequence file (fasta or an alphafold3 json). Alternatively can input a sequence alignment file here (eg from emboss needle).")

        # Optional positional argument for the output JSON file
        self.parser.add_argument("output_json", nargs='?',
                                 help=f"Path to the output JSON file (defaults to <input_sequence>{append_string}.json)")
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
append_string = '_template'
overwrite_template = False #Keep/discard templates that are already present in the input json file.  
default_alignment_type = 'emboss'
default_alignment_scoring = "blastp"

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


def add_template_json_data(data, template_data, sequence_number):
    sequence_count = 0
    if 'sequences' in data and isinstance(data['sequences'], list):
        for item in data['sequences']:
            if 'protein' in item and isinstance(item['protein'], dict):
                if 'sequence' in item['protein']:
                    if sequence_count == sequence_number:
                        if 'templates' in item['protein'] and isinstance(item['protein']['templates'], list) and not overwrite_template:
                            item['protein']['templates'] = item['protein']['templates']+template_data['templates']
                        else:
                            item['protein']['templates'] = template_data['templates']
                    sequence_count += 1
                else:
                    print(f"Warning: 'sequence' key not found in a protein entry.")
            else:
                print(f"Warning: 'protein' key not found or is not a dictionary in a sequence entry.")
    else:
        print(f"Warning: 'sequences' key not found or is not a list in the input file.")
    return data

def extract_sequence(structure_file):
    if structure_file.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    elif structure_file.endswith(".pdb"):
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Use PDB or mmCIF.")
    
    structure = parser.get_structure("protein", structure_file)
    seq = []
    res_ids = []

    
    if hasattr(parser, "_mmcif_dict"): #doesn't work
        release_date_entry = '_pdbx_audit_revision_history.revision_date'
        if not release_date_entry in parser._mmcif_dict.keys():
            print('mmCIF must have release date. (%s)' % release_date_entry)
        else:
            print("Release date: %s" % parser._mmcif_dict[release_date_entry])

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":  # Ignore heteroatoms
                    seq.append(residue.resname)
                    res_ids.append(residue.id[1])
            break  # Only take the first chain
        break  # Only take the first model
    
    three_to_one = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
    }
    
    sequence = "".join(three_to_one.get(res, "X") for res in seq)
    return sequence, res_ids

def align_sequences(seq1, seq2, alignment_scoring):
    aligner = Align.PairwiseAligner(scoring=alignment_scoring)
    aligner.mode = "global"
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    alignment = aligner.align(seq1, seq2)[0]
    return alignment

def compute_index_mapping(alignment, seq1, seq2, copies):
    seq1_len = len(seq1)
    if copies > 1:
        seq1_len = seq1_len/copies
    seq1_indices, seq2_indices = [[] for i in range(copies)], [[] for i in range(copies)]
    idx1, idx2 = 0, 0
    copy = 0
    for a, b in zip(alignment[0], alignment[1]):
        if a != "-":
            if b != "-":
                seq1_indices[copy].append(idx1)
                seq2_indices[copy].append(idx2)
            idx1 += 1
            if idx1 >= seq1_len:
                idx1 = 0
                copy += 1
        if b != "-":
            idx2 += 1
    return seq1_indices, seq2_indices

def main(input_structure, input_sequence, output_json, sequence_number, copies, alignment_type, alignment_scoring):
    structure_seq, structure_res_ids = extract_sequence(input_structure)

    input_ali = False
    if os.path.splitext(input_sequence)[1] == '.json':
        json_input = True
        json_data = parse_json(input_sequence)
        fasta_seq = sequence_from_json(json_data, sequence_number)
    else:
        json_input = False
        try:
            input_alignment = AlignIO.read(input_sequence, alignment_type)
            s1 = str(input_alignment[0].seq).strip('-')
            fasta_seq = str(input_alignment[1].seq).strip('-')
            input_ali = True
        except:
            pass

        if input_ali:
            if s1 != structure_seq:
                raise IOError('Alignment sequence 1 must match structure sequence')
        else:
            fasta_seq = str(next(SeqIO.parse(input_sequence, "fasta")).seq)

    if copies > 1:
        structure_seq = structure_seq*copies
    if not input_ali:
        alignment = align_sequences(structure_seq, fasta_seq, alignment_scoring)
    else:
        alignment = input_alignment

    print("Alignment:")
    print(alignment)
    
    mapping = compute_index_mapping(alignment, structure_seq, fasta_seq, copies)
    indent = "    "
    print("Index Mapping:")
    print("")
    template_json_start = f'''
"templates": ['''
    template_json_end = '''
    ]'''
    template_json_middle = ""
    for map1, map2 in zip(mapping[1], mapping[0]):
        middle = f'''
{indent}{{
{indent}"mmcifPath": "{input_structure}",
{indent}"queryIndices": {map1},
{indent}"templateIndices": {map2}
{indent}}},'''
        template_json_middle = template_json_middle+middle
    template_json_middle = template_json_middle[0:-1] #remove trailing comma
    template_json = template_json_start+template_json_middle+template_json_end
    print(template_json)
    if json_input:
        template_json = json.loads(f"{{\n{template_json}\n}}")
        output_json_data = add_template_json_data(json_data, template_json, sequence_number)
        if 'name' in output_json_data and isinstance(output_json_data['name'], str):
            output_json_data['name'] = output_json_data['name']+append_string
        with open(output_json, 'w') as f:
            json.dump(output_json_data, f, indent=4)  # Use indent=4 for pretty formatting
    else:
        name = os.path.basename(os.path.splitext(input_sequence)[0])+append_string
        template_json = template_json.replace('\n', f'\n{indent}{indent}')
        output_json_data = f'''{{
"dialect": "alphafold3",
"version": 2,
"name": "{name}",
"sequences": [
{indent}{{
{indent}"protein": {{
{indent}{indent}"id": ["A"],
{indent}{indent}"sequence": "{fasta_seq}",
{indent}{indent}"modifications": [],
{indent}{indent}"unpairedMsa": "",
{indent}{indent}"pairedMsa": "",{template_json}
{indent}{indent}}}
{indent}}}
],
"modelSeeds": [1],
"bondedAtomPairs": null,
"userCCD": null
}}'''
        with open(output_json, 'w') as f:
            f.writelines(output_json_data)
        print("Warning: Output json is msa-free.")
    print(f"Successfully written to '{output_json}'.")


if __name__ == "__main__":
    config = ArgumentParserConfig()
    args = config.parser.parse_args()

    input_structure = args.input_structure
    input_sequence = args.input_sequence
    output_json = args.output_json
    sequence_number = args.sequence_number
    copies = args.copies
    alignment_type = args.alignment_type
    alignment_scoring = args.alignment_scoring

    # Determine the output filename if not provided
    if output_json is None:
        base, ext = os.path.splitext(input_sequence)
        output_json = f"{base}{append_string}.json"

    main(input_structure, input_sequence, output_json, sequence_number, copies, alignment_type, alignment_scoring)
