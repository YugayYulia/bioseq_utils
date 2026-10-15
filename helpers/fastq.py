def check_gc(sequence: str, gc_bounds: tuple | int | float) -> bool:
    """
    Check whether the GC content of a sequence is within the given bounds.

    Args:
        sequence: Nucleotide sequence.
        gc_bounds: Allowed GC percentage range or a single upper bound.

    Returns:
        True if GC content is within bounds, otherwise False.
    """
    if not sequence:
        return False

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)

    gc_count = sequence.upper().count("G") + sequence.upper().count("C")
    gc_percent = gc_count / len(sequence) * 100

    return gc_bounds[0] <= gc_percent <= gc_bounds[1]


def check_length(sequence: str, length_bounds: tuple | int) -> bool:
    """
    Check whether the sequence length is within the given bounds.

    Args:
        sequence: Nucleotide sequence.
        length_bounds: Allowed length range or a single upper bound.

    Returns:
        True if sequence length is within bounds, otherwise False.
    """
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    length = len(sequence)

    return length_bounds[0] <= length <= length_bounds[1]


def check_quality(quality: str, quality_threshold: float) -> bool:
    """
    Check whether the mean Phred33 quality meets the threshold.

    Args:
        quality: String containing ASCII-encoded quality scores.
        quality_threshold: Minimum acceptable mean Phred33 quality.

    Returns:
        True if mean quality meets the threshold, otherwise False.
    """
    if not quality:
        return False

    scores = [ord(symbol) - 33 for symbol in quality]
    mean_quality = sum(scores) / len(scores)

    return mean_quality >= quality_threshold


def check_fastq_record(
    header: str,
    sequence: str,
    separator: str,
    quality: str,
) -> bool:
    """
    Validate the basic structure of a FASTQ record.

    Args:
        header: Read identifier line.
        sequence: Nucleotide sequence.
        separator: Separator line.
        quality: Phred33 quality string.

    Returns:
        True if the record is structurally valid, otherwise False.
    """
    if not header.startswith("@"):
        return False

    if not separator.startswith("+"):
        return False

    if len(sequence.strip()) != len(quality.strip()):
        return False

    return True


def write_fastq_record(
    output,
    header: str,
    sequence: str,
    separator: str,
    quality: str,
) -> None:
    """
    Write a FASTQ record to an open output file.

    Args:
        output: Open file object for writing.
        header: Read identifier line.
        sequence: Nucleotide sequence line.
        separator: Separator line.
        quality: Phred33 quality line.

    Returns:
        None.
    """
    output.write(header)
    output.write(sequence)
    output.write(separator)
    output.write(quality)


def read_fastq_record(file) -> tuple | None:
    """
    Read one FASTQ record from an open file.

    Args:
        file: Open FASTQ file for reading.

    Returns:
        Tuple containing four FASTQ lines,
        or None if the end of the file is reached.

    Raises:
        ValueError: If the FASTQ record is incomplete.
    """
    header = file.readline()

    if not header:
        return None

    sequence = file.readline()
    separator = file.readline()
    quality = file.readline()

    if not sequence or not separator or not quality:
        raise ValueError("Incomplete FASTQ record")

    return header, sequence, separator, quality
