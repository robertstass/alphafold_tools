#!/usr/bin/env python

import os.path
import copy
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, MMCIFIO
import sys


import argparse
import os


class ArgumentParserConfig:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="Combine the chains of a mmcif structure into a single chain.")


        # Positional argument for the input JSON file
        self.parser.add_argument("input_structure", help="Path to the input structure file (mmcif).")

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



three_to_one = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
}
out_ext = "cif"


def orig_patch_additional_cif_categories(new_dict, new_label_asym_id="A"):
    """
    Update other mmCIF categories beyond _atom_site to reflect new unified chain.
    This includes things like _pdbx_poly_seq_scheme, _struct_asym, etc.
    """
    def patch_asym_id(category, field="asym_id"):
        key = f"{category}.{field}"
        if key in new_dict:
            new_dict[key] = [new_label_asym_id] * len(new_dict[key])

    def patch_seq_id(category, field="seq_id"):
        key = f"{category}.{field}"
        if key in new_dict:
            new_dict[key] = [str(i + 1) for i in range(len(new_dict[key]))]

    # Update polymer sequence schemes
    patch_asym_id("_pdbx_poly_seq_scheme")
    patch_seq_id("_pdbx_poly_seq_scheme")

    # Update entity polymer sequence
    patch_seq_id("_pdbx_entity_poly_seq")

    # Update branch scheme if it exists
    patch_asym_id("_pdbx_branch_scheme")
    patch_seq_id("_pdbx_branch_scheme")

    # Update structural connection (bonding) info
    for prefix in ["ptnr1", "ptnr2"]:
        for suffix in ["label_asym_id", "auth_asym_id"]:
            key = f"_struct_conn.{prefix}_{suffix}"
            if key in new_dict:
                new_dict[key] = [new_label_asym_id] * len(new_dict[key])

        for suffix in ["label_seq_id", "auth_seq_id"]:
            key = f"_struct_conn.{prefix}_{suffix}"
            if key in new_dict:
                try:
                    new_dict[key] = [str(int(x) + 1) for x in new_dict[key]]  # shift to new numbering
                except ValueError:
                    pass  # sometimes these contain '?' or '.' placeholders

    # Collapse _struct_asym to a single chain
    if "_struct_asym.id" in new_dict:
        count = len(new_dict["_struct_asym.id"])
        new_dict["_struct_asym.id"] = [new_label_asym_id]
        if "_struct_asym.entity_id" in new_dict:
            first_entity = new_dict["_struct_asym.entity_id"][0]
            new_dict["_struct_asym.entity_id"] = [first_entity]
        for k in list(new_dict.keys()):
            if k.startswith("_struct_asym.") and k != "_struct_asym.id" and k != "_struct_asym.entity_id":
                new_dict[k] = [new_dict[k][0]]  # take the first value

    # --- Additional category patches ---

    # _struct_sheet_range
    if "_struct_sheet_range.sheet_id" in new_dict:
        patch_asym_id("_struct_sheet_range", "asym_id")
        for field in ["beg_seq_id", "end_seq_id"]:
            key = f"_struct_sheet_range.{field}"
            if key in new_dict:
                new_dict[key] = [str(int(x) + 1) if x.isdigit() else x for x in new_dict[key]]

    # _pdbx_struct_sheet_hbond
    for direction in ["range_1", "range_2"]:
        patch_asym_id("_pdbx_struct_sheet_hbond", f"{direction}_asym_id")
        key = f"_pdbx_struct_sheet_hbond.{direction}_seq_id"
        if key in new_dict:
            new_dict[key] = [str(int(x) + 1) if x.isdigit() else x for x in new_dict[key]]

    # _struct_site / _struct_site_gen
    patch_asym_id("_struct_site_gen", "auth_asym_id")
    if "_struct_site_gen.auth_seq_id" in new_dict:
        new_dict["_struct_site_gen.auth_seq_id"] = [
            str(int(x) + 1) if x.isdigit() else x for x in new_dict["_struct_site_gen.auth_seq_id"]
        ]

    # _pdbx_validate_torsion
    patch_asym_id("_pdbx_validate_torsion", "auth_asym_id")
    if "_pdbx_validate_torsion.auth_seq_id" in new_dict:
        new_dict["_pdbx_validate_torsion.auth_seq_id"] = [
            str(int(x) + 1) if x.isdigit() else x for x in new_dict["_pdbx_validate_torsion.auth_seq_id"]
        ]

    # _pdbx_unobs_or_zero_occ_residues
    patch_asym_id("_pdbx_unobs_or_zero_occ_residues", "label_asym_id")
    patch_seq_id("_pdbx_unobs_or_zero_occ_residues", "label_seq_id")


    return new_dict



def patch_additional_cif_categories(new_dict, new_label_asym_id="A"):
    """
    Update other mmCIF categories beyond _atom_site to reflect new unified chain.
    This includes things like _pdbx_poly_seq_scheme, _struct_asym, etc.
    """
    def patch_asym_id(category, field="asym_id"):
        key = f"{category}.{field}"
        if key in new_dict:
            new_dict[key] = [new_label_asym_id] * len(new_dict[key])

    def patch_seq_id(category, field="seq_id"):
        key = f"{category}.{field}"
        if key in new_dict:
            new_dict[key] = [str(i + 1) for i in range(len(new_dict[key]))]

    # Update polymer sequence schemes
    patch_asym_id("_pdbx_poly_seq_scheme")
    patch_asym_id("_pdbx_poly_seq_scheme", field='pdb_strand_id')
    patch_seq_id("_pdbx_poly_seq_scheme")

    # Update entity polymer sequence
    patch_seq_id("_pdbx_entity_poly_seq")

    # Update branch scheme if it exists
    patch_asym_id("_pdbx_branch_scheme")
    patch_seq_id("_pdbx_branch_scheme")

    patch_asym_id("_pdbx_nonpoly_scheme")
    patch_asym_id("_pdbx_nonpoly_scheme", field='pdb_strand_id')

    return new_dict


def remove_keys(dict, keys_to_remove):
    dict_keys = dict.copy().keys()
    for key_start in keys_to_remove:
        for key in dict_keys:
            if key.startswith(key_start):
                del dict[key]
    return dict

keys_to_remove = [
'_pdbx_poly_seq_scheme.',
'_pdbx_nonpoly_scheme.',
'_entity_poly_seq',
'_struct_ref_seq.',
'_struct_ref_seq_dif.',
'_struct_asym.',
'_pdbx_struct_assembly.',
'_pdbx_struct_assembly_prop.',
'_pdbx_struct_assembly_gen.',
'_struct_conf.',
'_struct_conf_type.',
'_struct_conn.',
'_struct_conn_type.',
'_pdbx_struct_conn_angle.',
'_pdbx_modification_feature.',
'_struct_mon_prot_cis.',
'_struct_sheet.',
'_struct_sheet_order.',
'_struct_sheet_range.',
'_pdbx_struct_sheet_hbond.',
'_struct_site.',
'_struct_site_gen.',
'_pdbx_validate_torsion.',
'_pdbx_unobs_or_zero_occ_residues',
'_pdbx_poly_scheme.',
'_pdbx_branch_scheme.'
]
'''
keys_to_remove = [
'_entity_poly_seq.',
]
'''

def main(structure_file, keep_hetatms=False):
    file_root = os.path.splitext(structure_file)[0]
    if structure_file.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    elif structure_file.endswith(".pdb"):
        raise ValueError('PDB files not currently supported. Use mmcif.')
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Use mmCIF.")

    structure = parser.get_structure("protein", structure_file)
    mmcif_dict = parser._mmcif_dict
    models = {}

    label_asym_ids = mmcif_dict["_atom_site.label_asym_id"]
    label_seq_ids = mmcif_dict["_atom_site.auth_seq_id"]  #mmcif_dict["_atom_site.label_seq_id"]
    auth_seq_ids = mmcif_dict.get("_atom_site.auth_seq_id", label_seq_ids)  # fallback if not present
    model_nums = mmcif_dict.get("_atom_site.pdbx_PDB_model_num", ["1"] * len(label_asym_ids))
    label_comp_ids = mmcif_dict["_atom_site.label_comp_id"]

    # We'll collect indices of atoms to keep, reindexing as we go
    new_label_asym_id = "A"
    new_indices = []
    new_label_seq_ids = []
    new_auth_seq_ids = []

    current_residue_number = 1  # Start numbering residues from 1

    seen_residues = set()
    for i, (model, chain, resnum, comp) in enumerate(zip(model_nums, label_asym_ids, label_seq_ids, label_comp_ids)):
        key = (model, chain, resnum)
        if keep_hetatms or comp in three_to_one.keys():
            if key not in seen_residues:
                seen_residues.add(key)
                current_residue_number += 1
            new_indices.append(i)
            new_label_seq_ids.append(str(current_residue_number))
            new_auth_seq_ids.append(str(current_residue_number))  # update if needed

    # Make a deep copy of the mmCIF dictionary
    new_dict = copy.deepcopy(mmcif_dict)

    for key in new_dict:
        if key.startswith("_atom_site.") and isinstance(new_dict[key], list):
            new_dict[key] = [new_dict[key][i] for i in new_indices]

    # Overwrite label_asym_id and sequence IDs
    new_dict["_atom_site.label_asym_id"] = [new_label_asym_id] * len(new_indices)
    new_dict["_atom_site.auth_asym_id"] = [new_label_asym_id] * len(new_indices)
    new_dict["_atom_site.label_seq_id"] = new_label_seq_ids
    new_dict["_atom_site.auth_seq_id"] = new_auth_seq_ids

    #new_dict = patch_additional_cif_categories(new_dict)
    new_dict = remove_keys(new_dict, keys_to_remove)

    # Save the new combined structure
    out_file = f"{os.path.splitext(structure_file)[0]}_combined_chain.cif"
    io = MMCIFIO()
    io.set_dict(new_dict)
    io.save(out_file)

    print(f"Wrote combined chain to {out_file}")

if __name__ == "__main__":
    config = ArgumentParserConfig()
    args = config.parser.parse_args()

    input_structure = args.input_structure

    main(input_structure)
