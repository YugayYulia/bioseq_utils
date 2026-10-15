from pathlib import Path


def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: str | None = None,
) -> None:
    """
    Convert multiline FASTA sequences into single-line sequences.

    Args:
        input_fasta: Path to the input FASTA file.
        output_fasta: Optional path to the output FASTA file.
            If None, a filename ending with "_oneline.fasta" is generated.

    Returns:
        None. Converted sequences are written to a file.
    """

    if output_fasta is None:
        input_path = Path(input_fasta)
        output_fasta = str(input_path.with_name(input_path.stem + "_oneline.fasta"))

    with open(input_fasta, "r") as file:
        with open(output_fasta, "w") as output:
            sequence_parts = []

            for line in file:
                line = line.strip()

                if line.startswith(">"):
                    if sequence_parts:
                        output.write("".join(sequence_parts) + "\n")
                        sequence_parts = []

                    output.write(line + "\n")
                else:
                    sequence_parts.append(line)

            if sequence_parts:
                output.write("".join(sequence_parts) + "\n")


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Extract the best BLAST hit description for each query.

    Args:
        input_file: Path to the BLAST text report.
        output_file: Path to the output text file.

    Returns:
        None. Best hit descriptions are written to a file.
    """
    descriptions = []
    in_alignments = False
    waiting_for_hit = False

    with open(input_file, "r") as file:
        for line in file:
            if "Sequences producing significant alignments:" in line:
                in_alignments = True
                waiting_for_hit = False
                continue

            if in_alignments and line.strip().startswith("Description"):
                waiting_for_hit = True
                continue

            if waiting_for_hit:
                if not line.strip():
                    continue

                description = line.split("  ")[0].strip()
                descriptions.append(description)

                in_alignments = False
                waiting_for_hit = False

    descriptions.sort(key=str.casefold)

    with open(output_file, "w") as output:
        for description in descriptions:
            output.write(description + "\n")


def read_gbk(input_gbk: str) -> list[dict]:
    """
    Read CDS annotations from a GenBank file.

    Args:
        input_gbk: Path to the input GenBank file.

    Returns:
        List of dictionaries with gene names, locus tags
        and protein sequences.
    """
    cds_list = []
    cds_lines = []
    in_cds = False

    with open(input_gbk, "r") as file:
        for line in file:
            if line.startswith("     CDS "):
                if cds_lines:
                    cds_list.append(parse_cds(cds_lines))

                cds_lines = [line]
                in_cds = True

            elif line.startswith("     ") and not line.startswith(
                "                     "
            ):
                if in_cds:
                    cds_list.append(parse_cds(cds_lines))
                    cds_lines = []
                    in_cds = False

            elif in_cds:
                cds_lines.append(line)

        if cds_lines:
            cds_list.append(parse_cds(cds_lines))

    return cds_list


def extract_qualifier(line: str) -> str:
    """
    Extract a qualifier value from a GenBank annotation line.

    Args:
        line: GenBank qualifier line.

    Returns:
        Extracted qualifier value.
    """
    return line.split('="', 1)[1].rstrip().strip('"')


def extract_translation(lines: list[str]) -> str:
    """
    Extract a protein sequence from GenBank CDS annotation lines.

    Args:
        lines: Lines belonging to one CDS feature.

    Returns:
        Protein sequence as a single string.
    """
    sequence_parts = []
    reading_translation = False

    for line in lines:
        line = line.strip()

        if line.startswith('/translation="'):
            reading_translation = True
            part = line.split('"', 1)[1]

        elif reading_translation:
            part = line

        else:
            continue

        if '"' in part:
            sequence_parts.append(part.split('"', 1)[0])
            break

        sequence_parts.append(part)

    return "".join(sequence_parts)


def parse_cds(lines: list[str]) -> dict:
    """
    Extract gene information from one CDS annotation.

    Args:
        lines: Lines belonging to one CDS feature.

    Returns:
        Dictionary with gene name, locus tag and protein sequence.
    """
    gene = None
    locus_tag = None

    for line in lines:
        line = line.strip()

        if line.startswith('/gene="'):
            gene = extract_qualifier(line)

        elif line.startswith('/locus_tag="'):
            locus_tag = extract_qualifier(line)

    translation = extract_translation(lines)

    return {
        "gene": gene,
        "locus_tag": locus_tag,
        "translation": translation,
    }


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: str | list[str] | tuple[str, ...] | set[str],
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = "selected_genes.fasta",
) -> None:
    """
    Save protein sequences of genes neighboring the target genes.

    Args:
        input_gbk: Path to the input GenBank file.
        genes: Target gene name or collection of gene names.
        n_before: Number of CDS before each target gene.
        n_after: Number of CDS after each target gene.
        output_fasta: Path to the output FASTA file.

    Returns:
        None. Neighboring protein sequences are saved to FASTA.
    """
    if n_before < 1 or n_after < 1:
        raise ValueError("n_before and n_after must be at least 1")

    if isinstance(genes, str):
        genes = [genes]

    target_genes = set(genes)
    cds_list = read_gbk(input_gbk)
    target_indices = set()

    for i, cds in enumerate(cds_list):
        if cds["gene"] in target_genes:
            target_indices.add(i)

    if not target_indices:
        raise ValueError("Target genes not found")

    neighbor_indices = set()

    for i in target_indices:
        start = max(0, i - n_before)
        end = min(len(cds_list), i + n_after + 1)

        for j in range(start, end):
            if j not in target_indices:
                neighbor_indices.add(j)

    with open(output_fasta, "w") as output:
        for i in sorted(neighbor_indices):
            cds = cds_list[i]

            if not cds["translation"]:
                continue

            gene_name = cds["gene"] or "unknown"
            locus_tag = cds["locus_tag"] or "unknown"

            output.write(f">{locus_tag} {gene_name}\n")
            output.write(cds["translation"] + "\n")
