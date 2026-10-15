from pathlib import Path

from helpers.fastq import (
    check_gc,
    check_length,
    check_quality,
    check_fastq_record,
    write_fastq_record,
    read_fastq_record,
)

from helpers.bio_files import (
    convert_multiline_fasta_to_oneline,
    parse_blast_output,
    select_genes_from_gbk_to_fasta,
)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple | int | float = (0, 100),
    length_bounds: tuple | int = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filter FASTQ reads by GC content, length and mean Phred33 quality.

    Args:
        input_fastq: Path to the input FASTQ file.
        output_fastq: Name of the output file saved in the filtered directory.
        gc_bounds: Allowed GC percentage range or a single upper bound.
        length_bounds: Allowed read length range or a single upper bound.
        quality_threshold: Minimum acceptable mean Phred33 quality.

    Returns:
        None. Filtered reads are written to a FASTQ file.

    Raises:
        ValueError: If a FASTQ record is incomplete or invalid.
    """

    output_dir = Path("filtered")
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / Path(output_fastq).name

    with open(input_fastq, "r") as file:
        with open(output_path, "w") as output:
            while True:
                record = read_fastq_record(file)
                if record is None:
                    break

                header, sequence, separator, quality = record

                sequence_clean = sequence.strip()
                quality_clean = quality.strip()

                gc_ok = check_gc(sequence_clean, gc_bounds)
                length_ok = check_length(sequence_clean, length_bounds)
                quality_ok = check_quality(quality_clean, quality_threshold)

                if gc_ok and length_ok and quality_ok:
                    write_fastq_record(
                        output,
                        header,
                        sequence,
                        separator,
                        quality,
                    )
