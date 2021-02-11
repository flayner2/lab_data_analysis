import os
import pandas as pd
from Bio import SeqIO
from filter_seqs import process_seq_by_class
from masked_seqs_stats import find_x_regions_and_calculate_stats

path_to_seqs = "/home/maycon/Documents/LAB/eusociality/local_data/crossmatch_filtered"

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

for seq_file in os.listdir(path_to_seqs):
    if seq_file.endswith(".screen"):
        taxon_name = seq_file.split(".")[0].split("_ests")[0]
        file_path = os.path.join(path_to_seqs, seq_file)

        for seq in SeqIO.parse(file_path, "fasta"):
            if "X" in seq.seq:
                current_seqlist = find_x_regions_and_calculate_stats(seq, taxon_name)

                if len(current_seqlist) > 0:
                    for seq_dict in current_seqlist:
                        xgroups_data = xgroups_data.append(seq_dict, ignore_index=True)

        nums_list = [
            "seq_len",
            "seq_xgroup_count",
            "xgroup_len",
            "dist_from_3",
            "dist_from_5",
        ]

        for col in nums_list:
            xgroups_data[col] = pd.to_numeric(
                xgroups_data[col], downcast="integer", errors="coerce"
            )

classes_to_filter = [1, 3, 6, 7]
path_to_clean_seqs = (
    "/home/maycon/Documents/LAB/eusociality/local_data/blast_removal_results"
)

for x_file in os.listdir(path_to_seqs):
    if x_file.endswith(".screen"):
        for clean_file in os.listdir(path_to_clean_seqs):
            if clean_file.endswith(".fasta") and (
                clean_file.split("_")[0] == x_file.split("_")[0]
            ):
                x_input_file = os.path.join(path_to_seqs, x_file)
                clean_input_file = os.path.join(path_to_clean_seqs, clean_file)

                org_name = "_".join(x_file.split(".")[0].split("_")[:2])
                out_name = f"{org_name}_ests_len_blast_vector_filtered.fasta"
                out_dir = "/home/maycon/Documents/LAB/eusociality/local_data/len_blast_vec_filtered"
                out_path = os.path.join(out_dir, out_name)

                print(org_name)

                with open(out_path, "a") as out_file:
                    clean_seqs = list(SeqIO.parse(clean_input_file, "fasta"))

                    for seq in SeqIO.parse(x_input_file, "fasta"):
                        selector = xgroups_data["seq_id"] == seq.id
                        selected = xgroups_data[selector].reset_index()

                        if selected.empty:
                            print(f"{seq.id}, not masked")

                            out_file.write(seq.format("fasta"))
                        else:
                            seq_class = selected.iloc[0]["seq_class"]

                            print(f"{seq.id}: {seq_class}")

                            if seq_class == 4:
                                continue
                            elif seq_class in classes_to_filter:
                                seq.seq = process_seq_by_class(seq.seq, seq_class)
                                out_file.write(seq.format("fasta"))
                            else:
                                for clean_seq in clean_seqs:
                                    if clean_seq.id == seq.id:
                                        out_file.write(clean_seq.format("fasta"))
