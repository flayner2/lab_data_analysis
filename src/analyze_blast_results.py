# Standard lib imports
import os
import sys

# Third-party imports
from Bio.Blast import NCBIXML


def parse_blast_xml(taxon: str, files_path: str, result_type: str):
    results = {}

    for file in os.listdir(files_path):
        if taxon in file and file.endswith(".xml") and result_type in file:
            curr_file = os.path.join(files_path, file)

            with open(curr_file, "r") as handle:
                blast_records = NCBIXML.parse(handle)

                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        hit_counter = 0

                        for hsp in alignment.hsps:
                            hit = f"hit_{hit_counter}"
                            query = blast_record.query.split(".")[-1]

                            if query not in results:
                                results[query] = {}
                                results[query]["count"] = 0
                                results[query]["hits"] = {}

                            if hit not in results[query]["hits"]:
                                results[query]["hits"][hit] = {}

                            results[query]["count"] += 1
                            results[query]["hits"][hit]["score"] = hsp.score
                            results[query]["hits"][hit]["expect"] = hsp.expect
                            results[query]["hits"][hit]["start"] = hsp.sbjct_start
                            results[query]["hits"][hit]["end"] = (
                                hsp.sbjct_start + hsp.align_length
                            )
                            hit_counter += 1

    return results


def main() -> None:
    if len(sys.argv) < 2:
        sys.exit(
            (
                "ERROR: Run this script with: `python analyze_blast_results.py)"
                "{path_to_results_directory}`"
            )
        )

    files_path = sys.argv[1]
    taxa = ["Apis_mellifera", "Polistes_canadensis", "Solenopsis_invicta"]
    result_types = ["cap3", "phrap"]

    try:
        for taxon in taxa:
            for result_type in result_types:
                results = parse_blast_xml(taxon, files_path, result_type)

                out_file = f"{taxon}_{result_type}_multialign.txt"

                with open(out_file, "a+") as result_file:
                    header = f"{taxon} | {result_type}\n"
                    result_file.write(header)

                    for result in results.keys():
                        if results[result]["count"] > 1:
                            sub_header = f">{result}:{results[result]['count']}\n"
                            result_file.write(sub_header)

                            for hit in results[result]["hits"]:
                                hit_str = (
                                    f"\t$ {hit}\t"
                                    f"{results[result]['hits'][hit]['score']}\t"
                                    f"{results[result]['hits'][hit]['expect']}\t"
                                    f"{results[result]['hits'][hit]['start']}\t"
                                    f"{results[result]['hits'][hit]['end']}\n"
                                )
                                result_file.write(hit_str)

    except FileNotFoundError:
        sys.exit("ERROR: Invalid path. Please provide a valid path to a directory")


if __name__ == "__main__":
    main()
