#!/usr/bin/env python3
"""Run the ACTB sequence retrieval, BLASTN search, and annotation workflow."""

from __future__ import annotations

import csv
import json
import logging
import os
from datetime import datetime, timezone
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML


# ---------------------------------------------------------------------
# PARAMETERS: edit these values in one place when changing the analysis
# ---------------------------------------------------------------------
ACCESSION = "NM_001101.5"
BLAST_PROGRAM = "blastn"
BLAST_DATABASE = "nt"
EXPECT_THRESHOLD = 10.0
MAX_HITS = 10
MEGABLAST = True

# All paths are anchored to this script, so the pipeline works from
# any current working directory after a clean repository checkout.
PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
INTERMEDIATE_DIR = RESULTS_DIR / "intermediate"
FINAL_DIR = RESULTS_DIR / "final"

FASTA_FILE = DATA_DIR / "query_sequence.fasta"
BLAST_XML_FILE = INTERMEDIATE_DIR / "blast_results.xml"
ANNOTATION_CSV_FILE = FINAL_DIR / "annotation_table.csv"
RUN_METADATA_FILE = RESULTS_DIR / "run_metadata.json"
LOG_FILE = RESULTS_DIR / "pipeline.log"


def configure_logging() -> None:
    """Write progress messages to both the terminal and a persistent log."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(LOG_FILE, encoding="utf-8"),
        ],
        force=True,
    )


def fetch_sequence():
    """Step 1: retrieve the query FASTA record from NCBI Entrez."""
    email = os.environ.get("NCBI_EMAIL")
    if not email:
        raise RuntimeError(
            "NCBI_EMAIL is not set. Set it to your email address before running."
        )

    Entrez.email = email
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    logging.info("Step 1/3: fetching %s from NCBI Nucleotide", ACCESSION)

    with Entrez.efetch(
        db="nucleotide",
        id=ACCESSION,
        rettype="fasta",
        retmode="text",
    ) as handle:
        record = SeqIO.read(handle, "fasta")

    SeqIO.write(record, FASTA_FILE, "fasta")
    logging.info("Saved FASTA intermediate: %s (%d nt)", FASTA_FILE, len(record.seq))
    return record


def run_blast(record) -> None:
    """Step 2: submit the retrieved sequence to NCBI BLASTN and retain raw XML."""
    INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)
    logging.info(
        "Step 2/3: submitting %s search against %s",
        BLAST_PROGRAM.upper(),
        BLAST_DATABASE,
    )

    result_handle = NCBIWWW.qblast(
        program=BLAST_PROGRAM,
        database=BLAST_DATABASE,
        sequence=record.format("fasta"),
        expect=EXPECT_THRESHOLD,
        hitlist_size=MAX_HITS,
        megablast=MEGABLAST,
        format_type="XML",
    )
    try:
        xml_text = result_handle.read()
    finally:
        result_handle.close()

    if not xml_text.strip():
        raise RuntimeError("NCBI BLAST returned an empty XML response.")

    BLAST_XML_FILE.write_text(xml_text, encoding="utf-8")
    logging.info("Saved raw BLAST XML intermediate: %s", BLAST_XML_FILE)


def create_annotation_table(query_length: int) -> int:
    """Step 3: parse retained BLAST XML and write a structured CSV result."""
    FINAL_DIR.mkdir(parents=True, exist_ok=True)
    logging.info("Step 3/3: parsing BLAST XML and creating annotation table")

    rows = []
    with BLAST_XML_FILE.open(encoding="utf-8") as handle:
        for result in NCBIXML.parse(handle):
            for alignment in result.alignments:
                for hsp in alignment.hsps:
                    identity = (
                        hsp.identities / hsp.align_length * 100
                        if hsp.align_length
                        else 0.0
                    )
                    coverage = (
                        hsp.align_length / query_length * 100
                        if query_length
                        else 0.0
                    )
                    rows.append(
                        {
                            "subject_accession": alignment.accession,
                            "subject_title": alignment.title,
                            "alignment_length": hsp.align_length,
                            "percent_identity": round(identity, 2),
                            "query_coverage_approx_percent": round(
                                min(coverage, 100.0), 2
                            ),
                            "e_value": hsp.expect,
                            "bit_score": round(hsp.bits, 2),
                        }
                    )

    rows.sort(key=lambda row: (row["e_value"], -row["bit_score"]))
    columns = [
        "subject_accession",
        "subject_title",
        "alignment_length",
        "percent_identity",
        "query_coverage_approx_percent",
        "e_value",
        "bit_score",
    ]

    with ANNOTATION_CSV_FILE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    logging.info("Saved %d annotation row(s): %s", len(rows), ANNOTATION_CSV_FILE)
    return len(rows)


def save_run_metadata(record, row_count: int) -> None:
    """Record key settings and outputs so a run can be inspected later."""
    metadata = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "query_accession": ACCESSION,
        "query_sequence_id": record.id,
        "query_length_nt": len(record.seq),
        "blast_program": BLAST_PROGRAM,
        "blast_database": BLAST_DATABASE,
        "expect_threshold": EXPECT_THRESHOLD,
        "maximum_hits": MAX_HITS,
        "megablast": MEGABLAST,
        "annotation_rows": row_count,
        "intermediate_files": [
            str(FASTA_FILE.relative_to(PROJECT_DIR)),
            str(BLAST_XML_FILE.relative_to(PROJECT_DIR)),
        ],
        "final_output": str(ANNOTATION_CSV_FILE.relative_to(PROJECT_DIR)),
    }
    RUN_METADATA_FILE.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    logging.info("Saved run metadata: %s", RUN_METADATA_FILE)


def main() -> int:
    configure_logging()
    logging.info("Starting end-to-end bioinformatics pipeline")
    try:
        record = fetch_sequence()
        run_blast(record)
        row_count = create_annotation_table(len(record.seq))
        save_run_metadata(record, row_count)
    except Exception:
        logging.exception("Pipeline failed")
        return 1

    logging.info("Pipeline completed successfully")
    logging.info("Final annotation table: %s", ANNOTATION_CSV_FILE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
