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
