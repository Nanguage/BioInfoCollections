# Pair-end reads align

Align R1 and R2, print the align results.

## Build

```bash
$ cargo build --release
```

## Usage

For example:

```bash
$ target/release/pair_end_align ./data/test_R1.fq ./data/test_R2.fq | less -S
```