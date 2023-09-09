# UniAligner

### Version 0.1

UniAligner (formerly, TandemAligner) is the first parameter-free algorithm for sequence alignment that introduces a _sequence-dependent alignment scoring_ 
that automatically changes for any pair of compared sequences.
Classical alignment approaches, such as the Smith-Waterman algorithm, that work well for most sequences,
fail to construct biologically adequate alignments of extra-long tandem repeats (ETRs), such as human centromeres and immunoglobulin loci.
This limitation was overlooked in the previous studies since the sequences of the centromeres and other ETRs across multiple genomes only became available recently.

UniAligner addresses this limitation.
As an input, it takes two strings in fasta format and outputs the alignment in CIGAR format.

## Building from source

To build UniAligner from source, one needs Linux, GNU make and C++-20 compatible compiler (typically takes less than 5 minutes):

```
  git clone git@github.com:seryrzu/unialigner.git
  cd unialigner/unialigner
  make -j4
```

Then UniAligner will be available at `./build/bin/unialigner`.

To launch a test launch: `make test_launch` (typically instant).

## Synopsis

To align sequences `first.fasta` and `second.fasta` and save the alignment into the directory `output_dir`:

```
  ./build/bin/unialigner --first first.fasta --second second.fasta -o output_dir
```

CIGAR string will be saved to `output_dir/cigar.txt`

## Citation

TandemAligner: a new parameter-free framework for fast sequence alignment. Andrey V. Bzikadze, Pavel A. Pevzner, _bioRxiv_

## Contact
Please report any problems to the [issue tracker](https://github.com/seryrzu/unialigner/issues).
Alternatively, you can write directly to [abzikadze@ucsd.edu](mailto:abzikadze@ucsd.edu).
