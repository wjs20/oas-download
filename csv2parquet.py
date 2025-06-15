#!/usr/bin/env python
import gzip
import polars as pl
from pathlib import Path
from typing import List, Dict
import logging
import json

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
    "sequence_header_heavy": "string",
    "sequence_heavy": "string",
    "sequence_aa_heavy": "string",
    "numbering_scheme_heavy": 'enum ["imgt", "kabat", "chothia", "martin"]',
    "locus_heavy": "string",
    "stop_codon_heavy": "boolean",
    "vj_in_frame_heavy": "boolean",
    "v_frameshift_heavy": "boolean",
    "j_frameshift_heavy": "boolean",
    "productive_heavy": "boolean",
    "rev_comp_heavy": "boolean",
    "complete_vdj_heavy": "boolean",
    "v_call_heavy": "string",
    "d_call_heavy": "string",
    "j_call_heavy": "string",
    "c_call_heavy": "string",
    "v_frame_heavy": "enum [0, 1, 2]",
    "j_frame_heavy": "enum [0, 1, 2]",
    "sequence_alignment_heavy": "string",
    "germline_alignment_heavy": "string",
    "sequence_alignment_aa_heavy": "string",
    "germline_alignment_aa_heavy": "string",
    "v_alignment_start_heavy": "integer",
    "v_alignment_end_heavy": "integer",
    "d_alignment_start_heavy": "integer",
    "d_alignment_end_heavy": "integer",
    "j_alignment_start_heavy": "integer",
    "j_alignment_end_heavy": "integer",
    "c_alignment_start_heavy": "integer",
    "c_alignment_end_heavy": "integer",
    "v_sequence_alignment_heavy": "string",
    "v_sequence_alignment_aa_heavy": "string",
    "v_germline_alignment_heavy": "string",
    "v_germline_alignment_aa_heavy": "string",
    "d_sequence_alignment_heavy": "string",
    "d_germline_alignment_heavy": "string",
    "j_sequence_alignment_heavy": "string",
    "j_sequence_alignment_aa_heavy": "string",
    "j_germline_alignment_heavy": "string",
    "j_germline_alignment_aa_heavy": "string",
    "c_sequence_alignment_heavy": "string",
    "c_germline_alignment_heavy": "string",
    "fwr1_heavy": "string",
    "fwr1_aa_heavy": "string",
    "cdr1_heavy": "string",
    "cdr1_aa_heavy": "string",
    "fwr2_heavy": "string",
    "fwr2_aa_heavy": "string",
    "cdr2_heavy": "string",
    "cdr2_aa_heavy": "string",
    "fwr3_heavy": "string",
    "fwr3_aa_heavy": "string",
    "cdr3_heavy": "string",
    "cdr3_aa_heavy": "string",
    "fwr4_heavy": "string",
    "fwr4_aa_heavy": "string",
    "junction_heavy": "string",
    "junction_aa_heavy": "string",
    "junction_length_heavy": "integer",
    "junction_aa_length_heavy": "integer",
    "v_score_heavy": "number",
    "d_score_heavy": "number",
    "j_score_heavy": "number",
    "c_score_heavy": "number",
    "v_cigar_heavy": "string",
    "d_cigar_heavy": "string",
    "j_cigar_heavy": "string",
    "c_cigar_heavy": "string",
    "v_support_heavy": "number",
    "d_support_heavy": "number",
    "j_support_heavy": "number",
    "c_support_heavy": "number",
    "v_identity_heavy": "number",
    "d_identity_heavy": "number",
    "j_identity_heavy": "number",
    "c_identity_heavy": "number",
    "v_sequence_start_heavy": "integer",
    "v_sequence_end_heavy": "integer",
    "d_sequence_start_heavy": "integer",
    "d_sequence_end_heavy": "integer",
    "j_sequence_start_heavy": "integer",
    "j_sequence_end_heavy": "integer",
    "c_sequence_start_heavy": "integer",
    "c_sequence_end_heavy": "integer",
    "v_germline_start_heavy": "integer",
    "v_germline_end_heavy": "integer",
    "d_germline_start_heavy": "integer",
    "d_germline_end_heavy": "integer",
    "j_germline_start_heavy": "integer",
    "j_germline_end_heavy": "integer",
    "c_germline_start_heavy": "integer",
    "c_germline_end_heavy": "integer",
    "fwr1_start_heavy": "integer",
    "fwr1_end_heavy": "integer",
    "cdr1_start_heavy": "integer",
    "cdr1_end_heavy": "integer",
    "fwr2_start_heavy": "integer",
    "fwr2_end_heavy": "integer",
    "cdr2_start_heavy": "integer",
    "cdr2_end_heavy": "integer",
    "fwr3_start_heavy": "integer",
    "fwr3_end_heavy": "integer",
    "cdr3_start_heavy": "integer",
    "cdr3_end_heavy": "integer",
    "fwr4_start_heavy": "integer",
    "fwr4_end_heavy": "integer",
    "sequence_aa_scheme_cigar_heavy": "string",
    "scheme_residue_mapping_heavy": "json string",
    "positional_scheme_mapping_heavy": "json string",
    "exc_heavy": "string",
    "additional_validation_flags_heavy": "json string",
    "sequence_header_light": "string",
    "sequence_light": "string",
    "sequence_aa_light": "string",
    "numbering_scheme_light": 'enum ["imgt", "kabat", "chothia", "martin"]',
    "locus_light": "string",
    "stop_codon_light": "boolean",
    "vj_in_frame_light": "boolean",
    "v_frameshift_light": "boolean",
    "j_frameshift_light": "boolean",
    "productive_light": "boolean",
    "rev_comp_light": "boolean",
    "complete_vdj_light": "boolean",
    "v_call_light": "string",
    "d_call_light": "string",
    "j_call_light": "string",
    "c_call_light": "string",
    "v_frame_light": "enum [0, 1, 2]",
    "j_frame_light": "enum [0, 1, 2]",
    "sequence_alignment_light": "string",
    "germline_alignment_light": "string",
    "sequence_alignment_aa_light": "string",
    "germline_alignment_aa_light": "string",
    "v_alignment_start_light": "integer",
    "v_alignment_end_light": "integer",
    "d_alignment_start_light": "integer",
    "d_alignment_end_light": "integer",
    "j_alignment_start_light": "integer",
    "j_alignment_end_light": "integer",
    "c_alignment_start_light": "integer",
    "c_alignment_end_light": "integer",
    "v_sequence_alignment_light": "string",
    "v_sequence_alignment_aa_light": "string",
    "v_germline_alignment_light": "string",
    "v_germline_alignment_aa_light": "string",
    "d_sequence_alignment_light": "string",
    "d_germline_alignment_light": "string",
    "j_sequence_alignment_light": "string",
    "j_sequence_alignment_aa_light": "string",
    "j_germline_alignment_light": "string",
    "j_germline_alignment_aa_light": "string",
    "c_sequence_alignment_light": "string",
    "c_germline_alignment_light": "string",
    "fwr1_light": "string",
    "fwr1_aa_light": "string",
    "cdr1_light": "string",
    "cdr1_aa_light": "string",
    "fwr2_light": "string",
    "fwr2_aa_light": "string",
    "cdr2_light": "string",
    "cdr2_aa_light": "string",
    "fwr3_light": "string",
    "fwr3_aa_light": "string",
    "cdr3_light": "string",
    "cdr3_aa_light": "string",
    "fwr4_light": "string",
    "fwr4_aa_light": "string",
    "junction_light": "string",
    "junction_aa_light": "string",
    "junction_length_light": "integer",
    "junction_aa_length_light": "integer",
    "v_score_light": "number",
    "d_score_light": "number",
    "j_score_light": "number",
    "c_score_light": "number",
    "v_cigar_light": "string",
    "d_cigar_light": "string",
    "j_cigar_light": "string",
    "c_cigar_light": "string",
    "v_support_light": "number",
    "d_support_light": "number",
    "j_support_light": "number",
    "c_support_light": "number",
    "v_identity_light": "number",
    "d_identity_light": "number",
    "j_identity_light": "number",
    "c_identity_light": "number",
    "v_sequence_start_light": "integer",
    "v_sequence_end_light": "integer",
    "d_sequence_start_light": "integer",
    "d_sequence_end_light": "integer",
    "j_sequence_start_light": "integer",
    "j_sequence_end_light": "integer",
    "c_sequence_start_light": "integer",
    "c_sequence_end_light": "integer",
    "v_germline_start_light": "integer",
    "v_germline_end_light": "integer",
    "d_germline_start_light": "integer",
    "d_germline_end_light": "integer",
    "j_germline_start_light": "integer",
    "j_germline_end_light": "integer",
    "c_germline_start_light": "integer",
    "c_germline_end_light": "integer",
    "fwr1_start_light": "integer",
    "fwr1_end_light": "integer",
    "cdr1_start_light": "integer",
    "cdr1_end_light": "integer",
    "fwr2_start_light": "integer",
    "fwr2_end_light": "integer",
    "cdr2_start_light": "integer",
    "cdr2_end_light": "integer",
    "fwr3_start_light": "integer",
    "fwr3_end_light": "integer",
    "cdr3_start_light": "integer",
    "cdr3_end_light": "integer",
    "fwr4_start_light": "integer",
    "fwr4_end_light": "integer",
    "sequence_aa_scheme_cigar_light": "string",
    "scheme_residue_mapping_light": "json string",
    "positional_scheme_mapping_light": "json string",
    "exc_light": "string",
    "additional_validation_flags_light": "json string",
    "ANARCI_numbering_heavy": "json string",
    "ANARCI_numbering_light": "json string",
}


def read_column_names(path: Path) -> List[str]:
    with open(path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def build_schema(columns: List[str]) -> Dict[str, pl.DataType]:
    return {
        col: TYPE_MAPPING[RAW_TYPE_DICT[col]] for col in columns if col in RAW_TYPE_DICT
    }


def process_file(
    input_path: Path, output_path: Path, schema: Dict[str, pl.DataType], metadata
):
    try:
        with gzip.open(input_path, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
            if len(lines[1].split(",")) != 198:
                raise ValueError("Trucated columns")
            csv_data = "".join(lines[1:])  # skip malformed header

        df = pl.read_csv(
            csv_data.encode("utf-8"), schema_overrides=schema, ignore_errors=True
        )
        df = df.select(list(schema.keys()))
        df = df.with_columns(
            [
                pl.col("ANARCI_numbering_light")
                .str.replace_all("'", '"')
                .str.replace_all(" ", ""),
                pl.col("ANARCI_numbering_heavy")
                .str.replace_all("'", '"')
                .str.replace_all(" ", ""),
            ]
        )
        df = df.with_columns(
            *(pl.lit(value).alias(field) for field, value in metadata.items())
        )
        df.write_parquet(output_path)
        logging.info(f"Processed {input_path.name} successfully.")
    except Exception as e:
        logging.error(f"Failed to process {input_path.name}: {e}")
        return str(e)
    return None


def main(input_dir: Path, output_dir: Path, column_file: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    column_names = read_column_names(Path(column_file))
    schema = build_schema(column_names)
    metadata = json.loads(Path("metadata.json").read_text())
    log = {}
    for file in input_dir.glob("*.csv.gz"):
        out_file = output_dir / (file.with_suffix(".parquet"))
        error = process_file(file, out_file, schema, metadata[file.name])
        if error:
            log[file.name] = error

    if log:
        log_file = output_dir / "failed_files.log"
        with open(log_file, "w") as f:
            for fname, errmsg in log.items():
                f.write(f"{fname}: {errmsg}\n")
        logging.info(f"Some files failed. See {log_file}")
    else:
        logging.info("All files processed successfully.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("column_file", type=Path)
    args = parser.parse_args()
    main(args.input_dir, args.output_dir, args.column_file)
