#!/usr/bin/env python

__author__ = "Hannah Holland-Moritz"
__email__ = "hannah.hollandmoritz@colorado.edu"
__version__ = "0.0.1"

"""
Rewrite silva taxonomy fasta file to be complient with USEARCH sinTAX command.
This program depends on biopython; to install use 

pip install biopython

learn more about biopython here: https://biopython.org/
"""

import argparse
from Bio import SeqIO
import re


def main():
    # Create the help messages and input requirements for the script

    parser = argparse.ArgumentParser(description=\
    '''
    Modifies SILVA taxonomy fasta file to be complient with the USEARCH sinTAX command.
    This program depends on biopython; to install it use
    
    pip install biopython
    ''',
    formatter_class=argparse.RawTextHelpFormatter)
    # what follows is some fancy code so that the required arguments are printed
    # above the optional arguments in the --help display.
    parser._action_groups.pop()
    req = parser.add_argument_group('required arguments')
    opt = parser.add_argument_group('optional arguments')

    # now the real arguments
    req.add_argument('-i', '--input_fp', required=True, type=str,
                     help=
     '''
     The input fasta file from SILVA, the expected header format is
     >SILVA_SAMPLE_ACCESSION_ID TAXONOMY;FIELDS;SEPARATED;BY;SEMICOLONS.
     ''')

    req.add_argument('-o', '--output_fp',
        required=True, type=str,
                     help=
     '''
     The output file path.
     ''')

    req.add_argument('-t', '--taxonomy_map_fp',
                     required=True, type=str,
                     help=
     '''
     The silva taxonomy map file path. This is the file named with the format
     tax_slv_(ls)su_VERSION.map
     ''')

    opt.add_argument('-l', '--length', type=int,
                     help=
     '''
     The number of sequences in input file. This argument
     is not required but may speed up processing for long files
     if it is included. If you choose to use this argument,
     ***PLEASE MAKE SURE IT IS CORRECT***, otherwise, the program
     will take the time to read through your entire fasta file first
     before failing which will be a very sad day for you if your 
     fasta file is massive.
     ''')

    # Now the actual code

    # First extract the arguments  from the command line and save the input and
    # output file and the length, if it is given
    args = parser.parse_args()

    input_fp = args.input_fp
    output_fp = args.output_fp
    map_fp = args.taxonomy_map_fp
    length = args.length

    # create a variable to hold modified sequences
    if length is not None:
        mod_records = [None] * length
        i = 0
    else:
        mod_records = []

    # Get the taxonomy dictionary (before the loop so that it only has to be created once)
    tax_dict = create_taxonomy_dictionary(map_fp)
    print("Taxa dictionary created.")

    print("Modifying sequence headers...")

    # read through the fasta file, for every record run the modify_fasta_header function
    for record in SeqIO.parse(input_fp, "fasta"):
        # extract header from record
        header=record.description

        # modify the header to match the sinTAX format, add ";" to sequence id
        seq_id_out, header_out = modify_fasta_header(header, tax_dict)

        # overwrite the record description and id with the new header and new sequence id
        record.description = header_out
        record.id = seq_id_out

        # save the record to the modified records list;
        # if list is not pre-allocated in size, simply append the records
        # otherwise, incrementally replace value of "None" in the records list
        # with the record.
        if length is None:
            mod_records.append(record)
        else:
            mod_records[i] = record
            i += 1
        # print(mod_records) # debugging

    # If the input length doesn't match the actual number of sequences
    # exit with an error; Note, this is a stupid place to exit from the program,
    # since it exists after doing all the work, but if these two don't match,
    # the program will exit anyway, since writing the output file will fail.
    # This is just to make it fail with an error that makes sense.
    if length is not None and len(mod_records) != i:
        raise RuntimeError('It looks like the number of sequences you have specified '
                           'in the --length option does not match the actual number of '
                           'sequences in the input fasta file. Please recount the '
                           'sequences and try again.')


    print("Writing modified headers to file...")

    # Write modified records to file
    with open(output_fp, "w") as handle:
        SeqIO.write(mod_records, handle, "fasta")
        #record.description = header_out


# Helper functions #


def create_taxonomy_dictionary(tax_map_fp):
    """
    create_taxonomy_dictionary(tax_map_fp)
    this function takes the tax.map file from silva
    and creates a dictionary from it where the keys
    are taxonomy levels and the values are the taxonomy
    numbers (i.e. 10000 for domain, 98 for genera, other numbers in between).
    The taxonomy levels that have numbers > 10000 or are labelled
    with a 1 or a 0 are for the most part eukaryotic kingdoms/domains
    for example, Fungi and Metazoa are 1, Opisthokonts and SAR are 10001
    Monilophytes (ferns) and Bilatarians (animals with bilateran symmetry)
    are 0. If the taxonomy level doesn't conform to the 7-level system,
    the value is "non-conforming". The function returns a dictionary
    with the taxonomy name as the key and the sinTax abbreviation as the value.
    """
    tax_map_handle = open(tax_map_fp, "r")

    # initialize dictionary
    tax_dict = {}

    # read the lines, split the lines
    for line in tax_map_handle:
        line = line.strip()
        sep_line = line.split("\t")
        taxonomy_name = sep_line[1]
        taxonomy_number = sep_line[3]
        # Translate the taxonomy numbers into a 6-level system
        if taxonomy_number == "98":
            tax_abr = "g:"
        elif taxonomy_number == "5":
            tax_abr = "f:"
        elif taxonomy_number == "4":
            tax_abr = "o:"
        elif taxonomy_number == "3":
            tax_abr = "c:"
        elif taxonomy_number == "2":
            tax_abr = "p:"
        elif taxonomy_number == "10000":
            tax_abr = "d:"
        # the taxonomy levels that have numbers > 10000 or are labelled
        # with a 1 or a 0 are for the most part eukaryotic kingdoms/domains
        # for example, Fungi and Metazoa are 1, Opisthokonts and SAR are 10001
        # Monilophytes (ferns) and Bilatarians (animals with bilateran symmetry)
        # are 0.
        else:
            tax_abr = "non-conforming:"
        tax_dict[taxonomy_name] = tax_abr

    return tax_dict


def multiple_replace(text, adict):

    """
    multiple_replace(text, adict)
    takes a text string and a dictionary of replacements
    allows multiple replacement of characters so that those characters
    only need to be specified in a dictionary once
    depends on the re module
    for more information on this function, see
    https://www.safaribooksonline.com/library/view/python-cookbook-2nd/0596007973/ch01s19.html
    """
    rx = re.compile('|'.join(map(re.escape, adict)))

    def one_xlat(match):
        return adict[match.group(0)]

    final_replace = rx.sub(one_xlat, text)

    return final_replace


def match_taxa(tax_dict, split_taxonomy):

    """
    match_taxa(tax_dict, split_taxonomy)
    this function will take the taxonomy dictionary and the
    split taxa string as input and assign it the proper
    d:/p:/c:/o:/f:/g:. If a split is the last one in the
    list, it will be given the species name. Non-conforming
    taxonomy levels are removed from the list. It will return a
    list of the taxonomy matches modified to include  d:/p:/c:/o:/f:/g:/s:
    #
    Note:
    this step would probably be better if vectorized, but
    I don't know how to do that easily, plus I want to keep
    the script dependencies small. So it is in a loop instead.
    """

    for taxa_level in list(range(len(split_taxonomy))):

        # save the current taxonomy level as a new variable
        current_tax_level = split_taxonomy[taxa_level]

        # then if the level is above species, look up it's abbreviation
        if taxa_level is not (len(split_taxonomy) - 1):
            # for abbreviations that conform to the sinTax style,
            # overwrite the split_taxonomy list with the new
            # style, if they are non-conforming, mark them by
            # labeling them "non-conforming".

            if tax_dict[current_tax_level] not in ["non-conforming:"]:
                split_taxonomy[taxa_level] = tax_dict[current_tax_level] + current_tax_level
            else:
                split_taxonomy[taxa_level] = "non-conforming"

        # if a taxa is the species level (i.e. the last level in the string)
        #  put an "s:" before it.
        else:
            split_taxonomy[taxa_level] = "s:" + current_tax_level

    # remove non-conforming taxonomy levels from the list
    if "non-conforming" in split_taxonomy:
        while "non-conforming" in split_taxonomy:
            split_taxonomy.remove("non-conforming")
        modified_split_taxonomy = split_taxonomy
    else:
        modified_split_taxonomy = split_taxonomy

    return modified_split_taxonomy


def modify_fasta_header(fasta_header, tax_dict):

    """
    modify_fasta_headers(fasta_header)
    This function takes a fasta header string as input and a taxonomy dictionary and
    1) splits the string into the sequence id and the taxonomy string
    2) splits the string based on semicolon placement
    3) places d:, p:, c:, o:, f:, g:, and s: at the appropriate places before each
    taxonomy level with the use of a helper function, match_taxa()
    4) removes special characters from the taxonomy strings by replacing them with underscores
    5) rebuilds the string and appends it to the sequence ids
    6) returns the sequence id with a ";" after it and the modified taxonomy header
    """

    # First split the header string on the first space and save the
    # sequence ID and taxonomy [second half]
    taxonomy_string = fasta_header.split(" ", 1)[1]
    sequence_id = fasta_header.split(" ", 1)[0]

    # Then split the string based on semicolons
    split_taxonomy = taxonomy_string.split(";")

    # for each split_taxonomy, check which taxa level it is, and rename it with the proper
    # d:/p:/c:/o:/f:/g: heading from the taxa_dict (species are not included in that list)
    # so the last split in the list that does not match any taxa level
    # will be given the species name.
    modified_taxa_list = match_taxa(tax_dict, split_taxonomy)

    # Then replace special characters with underscoes
    replacement_dict = {
        " " : "_",
        "(" : "",
        ")" : "",
        "/" : ""
    }

    modified_taxa_list_no_spc_char = []

    for taxa in range(len(modified_taxa_list)):
        #print(taxa)
        taxa_name = modified_taxa_list[taxa]
        #print(taxa_name)
        taxa_str_no_spc_char = multiple_replace(taxa_name, replacement_dict)
        #print(taxa_str_no_spc_char)
        modified_taxa_list[taxa] = taxa_str_no_spc_char
        #print(modified_taxa_list[taxa])

    # Then rebuild the string from the list by joining taxa list
    # together with a comma as the separator; also add "tax=" to the beginning
    full_taxonomy_string = "tax=" + ",".join(modified_taxa_list)

    # now append the modified string to the sequence id

    modified_header = full_taxonomy_string + ";"
    modified_seq_id = sequence_id + ";"

    return modified_seq_id, modified_header


if __name__ == "__main__":
    main()
