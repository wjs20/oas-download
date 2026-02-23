#!/usr/bin/env python

import argparse
import csv
import gzip
import io
import json
import sys
from itertools import islice
from pathlib import Path

description = """
The arrangement of columns is inconsistent between OAS data units.
This script will identify all possible columns, create a canonical sort order and write the data units back out with consistent column ordering.
"""


def parse_metadata(metadata):
    return json.loads(metadata.replace('""', '"').strip('"').strip().rstrip('"'))


def get_metadata_and_header(path):
    with gzip.open(path) as f:
        metadata, header = islice(f, 2)
        metadata = parse_metadata(metadata.decode("utf-8"))
        header = header.decode("utf-8").strip().split(",")
    return metadata, header


def list_gzipped_csvs(path):
    return [p for p in list(Path(path).iterdir()) if p.name.endswith(".csv.gz")]


def get_user_args():
    parser = argparse.ArgumentParser("Normalize", description=description)
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        help="Directory to find raw gzipped csv OAS data units",
    )
    parser.add_argument(
        "-o", "--output", type=Path, help="Directory to write normalized csvs"
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help="Directory to write the metadata.json file",
        default=Path(),
    )
    parser.add_argument(
        "--headers",
        type=Path,
        help="Directory to write the headers.json file",
        default=Path(),
    )
    return parser.parse_args()


def main():
    args = get_user_args()
    if not args.metadata.is_dir() and args.headers.is_dir():
        raise SystemExit("metadata and headers arguments must be directories")
    paths = list_gzipped_csvs(args.input)
    metadatas, headers = [], []
    for path in paths:
        metadata, header = get_metadata_and_header(path)
        metadatas.append(metadata)
        headers.extend(header)
    normalized_headers = sorted(set(headers))
    (args.metadata / "metadata.json").write_text(json.dumps(metadatas))
    (args.headers / "headers.json").write_text(json.dumps(normalized_headers))
    args.output.mkdir(exist_ok=True, parents=True)
    for in_path in paths:
        out_path = args.output / in_path.name

        with out_path.open("wb") as f:
            with gzip.open(f, "wb") as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="") as tf:
                    w = csv.DictWriter(
                        tf, fieldnames=normalized_headers, extrasaction="ignore"
                    )
                    w.writeheader()
                    print(f"Writing {in_path} to {out_path}", file=sys.stderr)
                    with gzip.open(in_path, "rt", newline="", encoding="utf-8") as in_f:
                        _ = next(in_f)  # drop metadata
                        r = csv.DictReader(in_f)
                        for row in r:
                            w.writerow({k: row.get(k, "") for k in normalized_headers})


if __name__ == "__main__":
    main()
