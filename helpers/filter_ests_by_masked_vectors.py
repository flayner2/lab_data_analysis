import os
import pandas as pd
from Bio import SeqIO
from filter_seqs import process_seq_by_class
from masked_seqs_stats import find_x_regions_and_calculate_stats


# Path to vector-masked sequence files
path_to_seqs = "/home/maycon/Documents/LAB/eusociality/local_data/crossmatch_filtered"

# DataFrame columns
xgroups_data = pd.DataFrame(
    columns=[
        "seq_id",
        "taxon",
        "seq_len",
        "seq_xgroup_count",
        "xgroup_len",
        "dist_from_3",
        "dist_from_5",
        "seq_class",
    ]
)

# Creating the dataframe with the x-regions data
for seq_file in os.listdir(path_to_seqs):
    if seq_file.endswith(".screen"):
        taxon_name = seq_file.split(".")[0].split("_ests")[0]
        file_path = os.path.join(path_to_seqs, seq_file)

        for seq in SeqIO.parse(file_path, "fasta"):
            # Only calculate x-region stats if there are x-regions in the sequence
            if "X" in seq.seq:
                current_seqlist = find_x_regions_and_calculate_stats(seq, taxon_name)

                if len(current_seqlist) > 0:
                    for seq_dict in current_seqlist:
                        xgroups_data = xgroups_data.append(seq_dict, ignore_index=True)

        # List of columns that should be numeric
        nums_list = [
            "seq_len",
            "seq_xgroup_count",
            "xgroup_len",
            "dist_from_3",
            "dist_from_5",
        ]

        # Converting those columns to a numeric data type
        for col in nums_list:
            xgroups_data[col] = pd.to_numeric(
                xgroups_data[col], downcast="integer", errors="coerce"
            )

# The sequence classes that sould be filtered by x-groups
classes_to_filter = [1, 3, 6, 7]
# Path to the sequence files without vector-masking (without 'X's)
path_to_clean_seqs = (
    "/home/maycon/Documents/LAB/eusociality/local_data/blast_removal_results"
)

for x_file in os.listdir(path_to_seqs):
    if x_file.endswith(".screen"):
        for clean_file in os.listdir(path_to_clean_seqs):
            # Find the vector-masked and the clean file for the same species
            if clean_file.endswith(".fasta") and (
                clean_file.split("_")[0] == x_file.split("_")[0]
            ):
                x_input_file = os.path.join(path_to_seqs, x_file)
                clean_input_file = os.path.join(path_to_clean_seqs, clean_file)

                org_name = "_".join(x_file.split(".")[0].split("_")[:2])
                out_name = f"{org_name}_ests_len_blast_vector_filtered.fasta"
                out_dir = (
                    "/home/maycon/Documents/LAB/eusociality/local_data/"
                    "len_blast_vec_filtered"
                )
                out_path = os.path.join(out_dir, out_name)

                print(org_name)

                with open(out_path, "a") as out_file:
                    # Load all the clean sequences to a list
                    clean_seqs = list(SeqIO.parse(clean_input_file, "fasta"))

                    for seq in SeqIO.parse(x_input_file, "fasta"):
                        # Match the DataFrame entry with the same id as the
                        # current sequence
                        selector = xgroups_data["seq_id"] == seq.id
                        selected = xgroups_data[selector].reset_index()

                        # Means the sequence is not vector-masked and should
                        # be kept as-is
                        if selected.empty:
                            print(f"{seq.id}, not masked")

                            out_file.write(seq.format("fasta"))
                        # Otherwise it should be filtered according to its
                        # sequence class
                        else:
                            seq_class = selected.iloc[0]["seq_class"]

                            print(f"{seq.id}: {seq_class}")

                            # Sequences of class 4 are removed
                            if seq_class == 4:
                                continue
                            # Of class 1, 3, 6 and 7 are trimmed
                            elif seq_class in classes_to_filter:
                                seq.seq = process_seq_by_class(seq.seq, seq_class)
                                out_file.write(seq.format("fasta"))
                            # The rest is kept as-is, without the masked vectors
                            # That's why we write the sequences from the clean_seq
                            else:
                                for clean_seq in clean_seqs:
                                    if clean_seq.id == seq.id:
                                        out_file.write(clean_seq.format("fasta"))
