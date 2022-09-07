# TandemAligner

### Version 0.1

TandemAligner is the first parameter-free algorithm for sequence alignment that introduces a _sequence-dependent alignment scoring_ 
that automatically changes for any pair of compared sequences.
Classical alignment approaches, such as the Smith-Waterman algorithm, that work well for most sequences,
fail to construct biologically adequate alignments of extra-long tandem repeats (ETRs), such as human centromeres and immunoglobulin loci.
This limitation was overlooked in the previous studies since the sequences of the centromeres and other ETRs across multiple genomes only became available recently.

TandemAligner addresses this limitation.
As an input, it takes two strings in fasta format and outputs the alignment in CIGAR format.

## Building from source

To build TandemAligner from source, one needs GNU make and C++-20 compatible compiler:

```
  git clone git@github.com:seryrzu/tandem_aligner.git
  cd tandem_aligner/tandem_aligner
  make -j4
```

Then TandemAligner will be available at `./build/bin/tandem_aligner`.

To launch a test launch: `make test_launch`.

## Synopsis

To align sequences `first.fasta` and `second.fasta` and save the alignment into the directory `output_dir`:

```
  ./build/bin/tandem_aligner --first first.fasta --second second.fasta -o output_dir
```

CIGAR string will be saved to `output_dir/cigar.txt`

## Citation

TandemAligner: a new parameter-free framework for fast sequence alignment. Andrey V. Bzikadze, Pavel A. Pevzner, _bioRxiv_

## Contact
Please report any problems to the [issue tracker](https://github.com/seryrzu/tandem_aligner/issues).
Alternatively, you can write directly to [abzikadze@ucsd.edu](mailto:abzikadze@ucsd.edu).
