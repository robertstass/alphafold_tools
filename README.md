# Alphafold tools 
A collection of useful scripts for running Alphafold3 with templates to enforce particular conformations.

Installation instructions are at the end of this file. 

# structure_sequence_alignment.py
This script takes a structure (mmcif/pdb) and a sequence as input. It aligns the structure's sequence to the query sequence and outputs json formatted information that defines the mapping of one to the other. (eg "queryIndices" and "templateIndices"). If an untemplated run has already been performed the output .json containing the msa can be used to provide the input sequence to this script. The output json file will then contain both template information and an msa (Otherwise it is msa-free). 

```commandline
structure_sequence_alignment.py <structure.cif> <sequence.fasta>
```
or
```commandline
structure_sequence_alignment.py <structure.cif> <sequence.json>
```

# copy_structure_release_date.py
Alphafold3 requires templates in mmcif format that contain release date information (otherwise it crashes). This metadata is often present in a cif file downloaded from the pdb but is often lost if the file has to be manipulated in some way. If this is the case, this script can be used to copy the release date information from the original cif and put it in the modified file.
```commandline
copy_structure_release_date.py <copy_from.cif> <copy_to.cif> 
```

# alphafold_msa_depth.py
Sometimes, even with a template supplied, alphafold will produce an output with an undesired conformation. In these cases the msa is overpowering the template. To help mitigate this you can weaken the effect of the msa by reducing the number of entries. The number to use will vary but typically reducing to ~50 will force adoption of the template's conformation. Be mindful that the lower the number, the less accurate alphafold's predictions will be. 

```commandline
alphafold_msa_depth.py --depth 50 <input.json>
```
# structure_split.py
Alphafold3 accepts only mmcif files with a single chain as templates. This script takes a multichain mmcif file and splits it into multiple files for each chain. 
```commandline
structure_split.py <structure.cif>
```

# structure_combine.py
Occassionally it is useful to treat a multichain template as if it is a single chain structure. Eg to run a prediction with multiple proteins connected by a flexible linker. Normally, separate proteins can be templated individually but a linked structure is required to provide a template that defines how the components fit together. This script will take a multichain mmcif file and combine all the chains into one chain. (Note: This messes up the format of the mmcif file a bit so it's not recommended to use these for anything other than creating alphafold templates!).
```commandline
structure_combine.py <structure.cif>
```

# alphafold_msa_repeat.py
This script creates a repeated msa linked by a linker. If you have a multimeric protein, you can generate an msa for the monomer by running a simple job and taking the .json output file. This script takes the msa contained within and repeats it linked by a glycine linker. eg KVFGR becomes KVFGRGGGGGGKVFGRGGGGGGKVFGRG with 3 copies and a 6 glycine linker. The resulting msa obviously doesn't cover the linker but we don't need it to. This is a much faster way to generate the msa (particularly for large sequences).
```commandline
alphafold_msa_repeat.py <input.json>
```
# structure_to_sequence.py
This script outputs a fasta sequence file from the coordinates of a structure file.
```commandline
structure_to_sequence.py <structure.cif>
```

# validate_json.py
Running alphafold often involves making small changes to json files and it can be easy to mess up the formatting. This script validates the json input file to minimize accidental job failures.
```commandline
validate_json.py <input.json>
```

# Installation
This script has the following requirements:  
-biopython (tested on v1.81)   
The easiest way to get access to these is with an [Anaconda](https://www.anaconda.com/download) installation as these libraries should be installed already. If not, set up an environment with:  
conda env create -f alphafold_tools.yml  
conda activate alphafold_tools   

Add the following to your ~/.bashrc file: 
```
export PATH=<full_path_to>/alphafold_tools:$PATH  
```
Or if on Windows add <full_path_to>/alphafold_tools to your system PATH environment variable.  

Then it can be run using:
```
structure_sequence_alignment.py  
```
(Note you may have to do "conda activate alphafold_tools" first depending on how you set it up)
