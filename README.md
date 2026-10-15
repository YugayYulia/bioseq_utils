# bioseq_utils

A Python project for working with biological sequence files.

## Functions

The project contains three main functions:

- `filter_fastq()` - filters FASTQ reads by GC content, sequence length, and quality.
- `convert_multiline_fasta_to_oneline()` - converts multiline FASTA sequences into single-line sequences.
- `parse_blast_output()` - extracts the best hits from a BLAST report and sorts them alphabetically.
- `select_genes_from_gbk_to_fasta()` - extracts protein sequences of genes neighboring selected genes from a GenBank file.
## Usage

### 1. FASTQ filtering

```python
from bioseq_utils import filter_fastq

filter_fastq(
    "test_data/test.fastq",
    "result.fastq",
    gc_bounds=(20, 80),
    length_bounds=(5, 10),
    quality_threshold=25,
)
```

The filtered reads are saved in the `filtered` folder.

### 2. FASTA conversion

```python
from bioseq_utils import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline("test_data/test.fasta")
```

The result is saved as `test_oneline.fasta` in the same folder as the input file.

### 3. BLAST results

```python
from bioseq_utils import parse_blast_output

parse_blast_output(
    "test_data/example_blast_results.txt", "test_data/blast_best_hits.txt"
)
```

The function saves the best hit descriptions in alphabetical order, one per line.

### 4. GenBank processing

```python
from bioseq_utils import select_genes_from_gbk_to_fasta

select_genes_from_gbk_to_fasta(
    "test_data/example_gbk.gbk",
    "phrB",
    n_before=1,
    n_after=1,
    output_fasta="test_data/phrB_neighbors.fasta",
)
```

The function saves protein sequences of neighboring genes to a FASTA file, excluding the target gene.

## Author

Yulia Yugay
