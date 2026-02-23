#!/usr/bin/env python
import argparse
import json
import polars as pl
from pathlib import Path
from typing import List, Dict
import logging

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

TYPE_MAPPING = {
    "string": pl.String,
    "boolean": pl.Boolean,
    "integer": pl.Int64,
    "number": pl.Float64,
    "json string": pl.String,
    "enum [0, 1, 2]": pl.Int8,
    'enum ["imgt", "kabat", "chothia", "martin"]': pl.String,
}

RAW_TYPE_DICT = {
    "ANARCI_numbering": "json string",
    "ANARCI_status": "string",
    "additional_validation_flags_heavy": "json string",
    "c_alignment_end": "integer",
    "c_alignment_start": "integer",
    "c_call": "string",
    "c_cigar": "string",
    "c_germline_alignment": "string",
    "c_germline_end": "integer",
    "c_germline_start": "integer",
    "c_identity": "number",
    "c_score": "number",
    "c_sequence_alignment": "string",
    "c_sequence_end": "integer",
    "c_sequence_start": "integer",
    "c_support": "number",
    "c_region": "string",
    "cdr1": "string",
    "cdr1_aa": "string",
    "cdr1_end": "integer",
    "cdr1_start": "integer",
    "cdr2": "string",
    "cdr2_aa": "string",
    "cdr2_end": "integer",
    "cdr2_start": "integer",
    "cdr3": "string",
    "cdr3_aa": "string",
    "cdr3_end": "integer",
    "cdr3_start": "integer",
    "complete_vdj": "boolean",
    "d_alignment_end": "integer",
    "d_alignment_start": "integer",
    "d_call": "string",
    "d_cigar": "string",
    "d_germline_alignment": "string",
    "d_germline_alignment_aa": "string",
    "d_germline_end": "integer",
    "d_germline_start": "integer",
    "d_identity": "number",
    "d_score": "number",
    "d_sequence_alignment": "string",
    "d_sequence_alignment_aa": "string",
    "d_sequence_end": "integer",
    "d_sequence_start": "integer",
    "d_support": "number",
    "exc": "string",
    "fwr1": "string",
    "fwr1_aa": "string",
    "fwr1_end": "integer",
    "fwr1_start": "integer",
    "fwr2": "string",
    "fwr2_aa": "string",
    "fwr2_end": "integer",
    "fwr2_start": "integer",
    "fwr3": "string",
    "fwr3_aa": "string",
    "fwr3_end": "integer",
    "fwr3_start": "integer",
    "fwr4": "string",
    "fwr4_aa": "string",
    "fwr4_end": "integer",
    "fwr4_start": "integer",
    "germline_alignment": "string",
    "germline_alignment_aa": "string",
    "j_alignment_end": "integer",
    "j_alignment_start": "integer",
    "j_call": "string",
    "j_cigar": "string",
    "j_frame": "enum [0, 1, 2]",
    "j_frameshift": "boolean",
    "j_germline_alignment": "string",
    "j_germline_alignment_aa": "string",
    "j_germline_end": "integer",
    "j_germline_start": "integer",
    "j_identity": "number",
    "j_score": "number",
    "j_sequence_alignment": "string",
    "j_sequence_alignment_aa": "string",
    "j_sequence_end": "integer",
    "j_sequence_start": "integer",
    "j_support": "number",
    "junction": "string",
    "junction_aa": "string",
    "junction_aa_length": "integer",
    "junction_length": "integer",
    "locus": "string",
    "numbering_scheme": 'enum ["imgt", "kabat", "chothia", "martin"]',
    "np2": "string",
    "np2_length": "integer",
    "np1": "string",
    "np1_length": "integer",
    "positional_scheme_mapping": "json string",
    "productive": "boolean",
    "rev_comp": "boolean",
    "scheme_residue_mapping": "json string",
    "sequence": "string",
    "sequence_aa": "string",
    "sequence_aa_scheme_cigar": "string",
    "sequence_alignment": "string",
    "sequence_alignment_aa": "string",
    "sequence_header": "string",
    "stop_codon": "boolean",
    "v_alignment_end": "integer",
    "v_alignment_start": "integer",
    "v_call": "string",
    "v_cigar": "string",
    "v_frame": "enum [0, 1, 2]",
    "v_frameshift": "boolean",
    "v_germline_alignment": "string",
    "v_germline_alignment_aa": "string",
    "v_germline_end": "integer",
    "v_germline_start": "integer",
    "v_identity": "number",
    "v_score": "number",
    "v_sequence_alignment": "string",
    "v_sequence_alignment_aa": "string",
    "v_sequence_end": "integer",
    "v_sequence_start": "integer",
    "v_support": "number",
    "vj_in_frame": "boolean",
    "Redundancy": "string",
}


def parse_user_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("column_file", type=Path)
    return parser.parse_args()


def build_schema(columns: List[str], paired: bool = False) -> Dict[str, pl.DataType]:
    """
    Build the schema from an input columns list. Used to specify the column types in constructed polars dataframes.

    Args:
        columns: list of columns found in csvs to be converted to parquets.
        paired: boolean indicating whether the schema will be used to construct parquets for the paired antibody database,
        or for the unpaired antibody database. If paired, each column is appended with '_heavy' or '_light'.
    """
    if paired:
        type_mapping = {
            k + suffix: v
            for suffix in ("_heavy", "_light")
            for k, v in RAW_TYPE_DICT.items()
        }
    else:
        type_mapping = TYPE_MAPPING
    return {col: type_mapping[RAW_TYPE_DICT[col]] for col in columns}

    column_names = read_column_names(Path(column_file))
    schema = build_schema(column_names)
    for file in input_dir.glob("*.csv.gz"):
        output_path = (output_dir / file.name).with_suffix(".parquet")
        try:
            df = pl.read_csv(file, schema=schema, ignore_errors=True)
            df.write_parquet(output_path)
            logging.info(f"Processed {file.name} successfully.")
        except Exception as e:
            logging.debug(f"Failed to process {file}: {e}")


if __name__ == "__main__":
    args = parse_user_args()
    main(args.input_dir, args.output_dir, args.column_file)
