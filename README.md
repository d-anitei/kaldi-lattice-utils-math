# kaldi-lattice-utils

Utilities for building probabilistic indexes from Kaldi lattices, with an
additional math-aware character lattice indexer for mathematical expressions.

This repository is based on Joan Puigcerver's Kaldi lattice indexing tools and
extends them with `lattice-char-index-position-math`, a command-line tool that
extracts n-best expression/sub-expression hypotheses from character lattices.
The math indexer understands LaTeX-style delimiter scopes such as `{ ... }` and
`{@ ... }@`, and can add extra segmentation points around operator symbols such
as `@bin...`, `@rel...`, and `@punct...`.

## Main Tool

### `lattice-char-index-position-math`

Builds a position-level index from character lattices. Each output hypothesis
contains:

- an expression encoded as underscore-separated character labels,
- an expression position,
- initial and final frame indexes,
- a log score.

Required delimiter options:

- `--symbol-table=words.txt`
- `--open-delimiters="{ {@"`
- `--close-delimiters="} }@"`
- `--operator-segment-prefixes="@bin @rel @punct"`

The default delimiter values match the expected math vocabulary, so in normal
use only `--symbol-table` is required explicitly.

## Usage

```bash
lattice-char-index-position-math --symbol-table=words.txt [options] \
  <lattice-rspecifier> <index-wspecifier>
```

Example:

```bash
lattice-char-index-position-math --print-args=true --verbose=2 \
  --nbest=25 \
  --min-start-state-posterior=1e-10 \
  --symbol-table=words.txt \
  --open-delimiters="{ {@" \
  --close-delimiters="} }@" \
  --operator-segment-prefixes="@bin @rel @punct" \
  "ark:gunzip -c lat_0001.gz |" ark,t:-
```

## Symbol Table Format

The symbol table is an OpenFst/Kaldi text symbol table:

```text
@bin+ 26
@bin- 27
@punct, 66
@rel= 71
{ 494
{@ 495
} 496
}@ 497
```

Open and close delimiters are paired by position:

```bash
--open-delimiters="{ {@"
--close-delimiters="} }@"
```

This means:

- `{` pairs with `}`
- `{@` pairs with `}@`

Operator segmentation labels are all labels whose symbol names start with one
of the prefixes passed through `--operator-segment-prefixes`.

## Output

The output is written with Kaldi table writers. With `ark,t:-`, each line starts
with the lattice key followed by tuples:

```text
utt-id expr position start-frame end-frame log-score ; ...
```

Expression labels are printed as integer character labels joined with `_`.
A downstream vocabulary table can be used to map these labels back to symbols.

## Build

This project is built inside a Kaldi source tree or against an existing Kaldi
checkout. The `kaldi.mk` file must point to your Kaldi installation.

Build all tools:

```bash
make
```

Build only the math indexer:

```bash
make -C kwsbin2 lattice-char-index-position-math
```

Install binaries:

```bash
make install
```

## Other Tools

The repository also includes the original lattice indexing utilities:

- `lattice-word-index-utterance`
- `lattice-word-index-segment`
- `lattice-word-index-position`
- `lattice-char-index-segment`
- `lattice-char-index-position`

See [kwsbin2/README.md](kwsbin2/README.md) for examples of the original word and
character lattice indexers.

## License

The original code is distributed under the Apache License, Version 2.0. See the
license headers in the source files for details.
