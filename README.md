# kaldi-lattice-utils

This repository builds a searchable probabilistic index of mathematical expressions
recognized in images. A math expression recognition (MER) model can produce
several plausible readings of the same image. A Kaldi recognition lattice
represents those alternatives and their scores. The indexer extracts likely
full expressions and subexpressions from the lattice so a search can find
an expression on its own or as part of a larger formula.

Each index entry identifies the source lattice, the recognized tokens, their
estimated probability, and the frame range where the expression appears.
The indexer uses LaTeX-style delimiters to identify nested groups and can
split expressions around operators. It extends Joan Puigcerver's
[Kaldi lattice indexing tools](kwsbin2/README.md).

## Pipeline

1. Normalize LaTeX expressions and generate structural annotations.
   [latex-math-normalization](https://github.com/d-anitei/latex-math-normalization)
2. Train a math expression recognition (MER) model on the annotated data.
   [PyLaia-math](https://github.com/d-anitei/PyLaia-math)
3. Run the model on expression images and convert its posteriorgrams into Kaldi
   recognition lattices. [PyLaia-math](https://github.com/d-anitei/PyLaia-math)
4. Build the index from a gzipped recognition lattice archive and its matching
   Kaldi symbol table. [build-math-index.sh](util/build-math-index.sh)

## Build and run

Set `KALDI_ROOT` to your built Kaldi checkout, then run these commands from the
repository root. The wrapper also requires Bash, GNU awk, gzip, and Python 3.

The `examples/` directory contains a sample recognition lattice archive and
its matching symbol table. Use them to try the indexing stage.

```bash
make -C kwsbin2 lattice-char-index-position-math
PATH="$PWD/kwsbin2:$PATH" ./util/build-math-index.sh \
  examples/ts-lm-2.lat.gz examples/ts-lm-2.words.txt
```

The symbol table maps the numeric labels in the lattice to expression tokens.
The wrapper uses this mapping to decode the index, converts log scores to
probabilities, and filters out expressions with mismatched or unpaired
delimiters. Its `--min-posterior=1e-10` cutoff prunes low-probability region
starts and hypotheses. For direct use and other options, run
`./kwsbin2/lattice-char-index-position-math --help`.

## Index format

`TMP-IDX/idx.gz` contains one hypothesis per line with six space-separated columns.

| Column | Meaning |
| --- | --- |
| 1 | Lattice key (utterance ID) |
| 2 | Expression tokens separated by `¬` |
| 3 | Posterior probability |
| 4 | Start frame |
| 5 | End frame |
| 6 | Hypothesis ordinal position within the lattice. Multiple hypotheses can share a position. |

The sample lattice produces this row.

```text
0005174_page003.0005174_page003-1-emb Z¬Z¬_¬{¬4¬} 5.50884e-05 0 14 2
```

This hypothesis belongs to region 2 of the named lattice and spans frames 0
through 14. The `¬` characters mark token boundaries within the expression.

## License

[Apache License 2.0](LICENSE)
