# pilercr-parser

A parser for the output of the [PILER-CR](https://www.drive5.com/pilercr/) CRISPR annotation tool.

PILER-CR v1.06 (at least) reports incorrect coordinates if any of the repeat sequences contains gaps. This parser will correct those errors, and also provides the repeat sequence of each repeat-spacer (which is given only as a difference pattern to the consensus in the PILER-CR output).

### Installation

Add the following to Cargo.toml:

`pilercr-parser = 1.0.2`

### Usage

```rust
use std::fs::File;
use std::io::{BufReader, Read};

fn main() {
    let file = File::open("examples/example.txt").unwrap();
    let mut reader = BufReader::new(file);
    let mut input = String::new();
    reader.read_to_string(&mut input).unwrap();
    let arrays = pilercr_parser::parse(&input).unwrap();
    for array in arrays {
        println!(
            "{} has {} arrays",
            array.accession,
            array.repeat_spacers.len()
        );
    }
}
```
