#![deny(warnings, missing_docs)]
//! Parses the output produced by PILER-CR (<https://www.drive5.com/pilercr/>), a CRISPR array
//! annotation tool.
//!
//! PILER-CR v1.06 (at least) reports incorrect coordinates if any of the repeat sequences contains gaps.
//! This parser will correct those errors, and also provides the repeat sequence of each repeat-spacer
//! (which is given only as a difference pattern to the consensus in the PILER-CR output).
//!
//! ## Example
//!
//! ```rust
//! use std::fs::File;
//! use std::io::{BufReader, Read};
//!
//! let file = File::open("examples/example.txt").unwrap();
//! let mut reader = BufReader::new(file);
//! let mut input = String::new();
//! reader.read_to_string(&mut input).unwrap();
//! let arrays = pilercr_parser::parse(&input).unwrap();
//! for array in arrays {
//!     println!(
//!         "{} has {} arrays",
//!         array.accession,
//!         array.repeat_spacers.len()
//!     );
//! }
//! ```

use nom::{
    bytes::complete::tag,
    character::complete::{
        alpha0, char, digit0, digit1, line_ending, multispace0, multispace1, not_line_ending,
    },
    error::Error,
    multi::{many0, many1},
    number::complete::float,
    sequence::{pair, tuple},
    Err, IResult, InputTakeAtPosition,
};

#[derive(Debug, PartialEq)]
/// Represents the information of a repeat-spacer as reflected in the PILER-CR output.
/// The coordinates here can be incorrect (we will correct them later) and the
/// repeat sequence has not yet been constructed.
struct RawRepeatSpacer<'a> {
    /// Zero-indexed, inclusive start coordinate.
    start: usize,
    /// Zero-indexed, exclusive end coordinate.
    end: usize,
    /// A pattern representing the difference between this repeat and the consensus repeat.
    repeat_diff: &'a str,
    /// Sequence of the spacer.
    spacer: &'a str,
}

#[derive(Debug, PartialEq)]
/// A single repeat-spacer.
pub struct RepeatSpacer<'a> {
    /// Zero-indexed, inclusive start coordinate.
    pub start: usize,
    /// Zero-indexed, exclusive end coordinate.
    pub end: usize,
    /// Zero-indexed, inclusive start coordinate of the spacer.
    pub spacer_start: usize,
    /// Zero-indexed, exclusive end coordinate of the spacer.
    pub spacer_end: usize,
    /// Zero-indexed, inclusive start coordinate of the repeat.
    pub repeat_start: usize,
    /// Zero-indexed, exclusive end coordinate of the repeat.
    pub repeat_end: usize,
    /// Sequence of the repeat.
    pub repeat: String,
    /// Sequence of the spacer.
    pub spacer: &'a str,
}

#[derive(Debug, PartialEq)]
/// A single CRISPR array.
pub struct Array<'a> {
    /// Accession of the contig/genome.
    pub accession: &'a str,
    /// The Nth CRISPR array in the PILER-CR output.
    pub order: usize,
    /// Zero-indexed, inclusive start coordinate.
    pub start: usize,
    /// Zero-indexed, exclusive end coordinate.
    pub end: usize,
    /// The consensus repeat sequence for this array. May include gaps.
    pub consensus_repeat_sequence: &'a str,
    /// The repeat-spacers in this array.
    pub repeat_spacers: Vec<RepeatSpacer<'a>>,
}

/// Parses the output of PILER-CR for a single contig/genome.
pub fn parse(input: &str) -> Result<Vec<Array>, Err<Error<&str>>> {
    let result = tuple((skip_header, many0(parse_array)))(input);
    match result {
        Ok((_, (_, arrays))) => Ok(arrays),
        Err(e) => Err(e),
    }
}

/// Gets space-delimited text.
fn not_space(input: &str) -> IResult<&str, &str> {
    input.split_at_position_complete(char::is_whitespace)
}

/// Skips the lines at the beginning of the PILER-CR output.
fn skip_header(input: &str) -> IResult<&str, ()> {
    let result = tuple((
        skip_one_line,
        skip_one_line,
        skip_empty_line,
        skip_one_line,
        skip_empty_line,
        skip_empty_line,
        skip_empty_line,
        skip_one_line,
        skip_empty_line,
        skip_empty_line,
        skip_empty_line,
    ))(input);
    match result {
        Ok((remainder, _)) => Ok((remainder, ())),
        Err(e) => Err(e),
    }
}

/// Parses a single repeat spacer. These may have incorrect coordinates, and the repeat sequence
/// has not yet been determined.
fn parse_raw_repeat_spacer(input: &str) -> IResult<&str, RawRepeatSpacer> {
    let result = tuple((
        multispace0,
        digit1,
        multispace1,
        digit1,
        multispace1,
        float,
        multispace1,
        digit0,
        multispace0,
        alpha0,
        multispace1,
        not_space,
        multispace1,
        alpha0,
    ))(input);
    match result {
        Ok((remainder, data)) => {
            let repeat_diff = data.11;
            let spacer = data.13;
            let start = data.1.parse::<usize>().unwrap() - 1;
            let end = start + repeat_diff.len() + spacer.len();
            let raw_repeat_spacer = RawRepeatSpacer {
                start,
                end,
                repeat_diff,
                spacer,
            };
            Ok((remainder, raw_repeat_spacer))
        }
        Err(e) => Err(e),
    }
}

/// Skips a line with text.
fn skip_one_line(input: &str) -> IResult<&str, ()> {
    let result = pair(not_line_ending, line_ending)(input);
    match result {
        Ok((remaining, _)) => Ok((remaining, ())),
        Err(e) => Err(e),
    }
}

/// Skips an empty line.
fn skip_empty_line(input: &str) -> IResult<&str, ()> {
    let result = line_ending(input);
    match result {
        Ok((remaining, _)) => Ok((remaining, ())),
        Err(e) => Err(e),
    }
}

/// Gets the consensus sequence from the last line of an array and discards everything else.
fn parse_array_summary_line(input: &str) -> IResult<&str, &str> {
    let result = tuple((
        multispace0,
        digit1,
        multispace1,
        digit1,
        multispace1,
        digit1,
        multispace1,
        not_space,
        line_ending,
    ))(input);
    match result {
        Ok((remainder, data)) => Ok((remainder, data.7)),
        Err(e) => Err(e),
    }
}

/// Parses a single CRISPR array.
fn parse_array(input: &str) -> IResult<&str, Array> {
    let result = tuple((
        tag("Array "),
        digit1,
        line_ending,
        char('>'),
        not_space,
        line_ending,
        skip_empty_line,
        skip_one_line,
        skip_one_line,
        many1(parse_raw_repeat_spacer),
        skip_one_line,
        skip_one_line,
        parse_array_summary_line,
        skip_empty_line,
        skip_empty_line,
    ))(input);
    match result {
        Err(e) => Err(e),
        Ok((remainder, data)) => {
            let order = data.1.parse::<usize>().unwrap() - 1;
            let accession = data.4;
            let raw_repeat_spacers = data.9;
            let consensus_repeat_sequence = data.12;
            let repeat_spacers =
                convert_raw_rs_to_final_rs(consensus_repeat_sequence, &raw_repeat_spacers);
            let start = repeat_spacers.first().unwrap().start;
            let end = repeat_spacers.last().unwrap().end;
            Ok((
                remainder,
                Array {
                    start,
                    end,
                    order,
                    accession,
                    consensus_repeat_sequence,
                    repeat_spacers,
                },
            ))
        }
    }
}

/// Due to a bug in PILER-CR, coordinates don't take gaps in repeat sequences into account.
/// We correct those coordinates here. Additionally, each repeat has a difference pattern instead
/// of an actual sequence, so we determine what the true repeat sequence is.
fn convert_raw_rs_to_final_rs<'a>(
    consensus_repeat: &'a str,
    raw_repeat_spacers: &[RawRepeatSpacer<'a>],
) -> Vec<RepeatSpacer<'a>> {
    let mut output = vec![];
    let mut total_gap_count = 0usize;
    for raw in raw_repeat_spacers {
        assert_eq!(raw.repeat_diff.len(), consensus_repeat.len());
        let repeat = raw
            .repeat_diff
            .chars()
            .zip(consensus_repeat.chars())
            .filter(|(r, _)| *r != '-')
            .map(|(r, c)| if r == '.' { c } else { r })
            .collect::<String>();
        let gap_count = raw.repeat_diff.matches('-').count();
        let rs = RepeatSpacer {
            start: raw.start - total_gap_count,
            end: raw.end - total_gap_count - gap_count,
            repeat_start: raw.start - total_gap_count,
            repeat_end: raw.start - total_gap_count + repeat.len(),
            spacer_start: raw.start - total_gap_count + repeat.len(),
            spacer_end: raw.end - total_gap_count - gap_count,
            repeat,
            spacer: raw.spacer,
        };
        output.push(rs);
        total_gap_count += gap_count;
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_raw_repeat_spacer() {
        let input = "       462      36   100.0      29  CTTTCTGAAG    ....................................    CGTGCTCGCTTTGAATTTGTAGAACCCGA";
        let expected = RawRepeatSpacer {
            start: 461,
            end: 461 + 36 + 29,
            repeat_diff: "....................................",
            spacer: "CGTGCTCGCTTTGAATTTGTAGAACCCGA",
        };
        let (_, actual) = parse_raw_repeat_spacer(input).unwrap();
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_parse_array_summary_line() {
        let input = "        22      36              29                GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC\n";
        let expected = "GTTGTGGTTTGATGTAGGAATCAAAAGATATACAAC";
        let (_, actual) = parse_array_summary_line(input).unwrap();
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_convert_raw_rs_to_final_rs() {
        let consensus = "AAGTTTCCGTCCCCTTTCGGGGAATCATTTAGAAAAT--A";
        let raws = vec![
            RawRepeatSpacer {
                start: 3831,
                end: 3906,
                repeat_diff: "..A..................................CC.",
                spacer: "GAATTACATCGTATGCCAATACGCAGTTGCTTTT",
            },
            RawRepeatSpacer {
                start: 3831,
                end: 3906,
                repeat_diff: "GG............-......................--.",
                spacer: "ATCACATTCA",
            },
        ];
        let rs = convert_raw_rs_to_final_rs(consensus, &raws);
        assert_eq!(rs.len(), 2);
        assert_eq!(rs[0].repeat, "AAATTTCCGTCCCCTTTCGGGGAATCATTTAGAAAATCCA");
        assert_eq!(rs[1].repeat, "GGGTTTCCGTCCCCTTCGGGGAATCATTTAGAAAATA");
    }

    #[test]
    fn test_parse_array() {
        let input = "Array 5
>MGYG000273829_14

       Pos  Repeat     %id  Spacer  Left flank    Repeat                                  Spacer
==========  ======  ======  ======  ==========    ====================================    ======
     16576      36   100.0      30  AAACAGTTCT    ....................................    ACGAACTTAGTACCCTTTTCTGGGCGGCAT
     16642      36   100.0      30  TGGGCGGCAT    ....................................    CCGCAGGTGCTACCGCTGTTATACTCTGTT
     16708      36   100.0      30  ATACTCTGTT    ....................................    CGTAAATCGTTGGCGAAACGCTACCAACTG
     16774      36   100.0      30  CTACCAACTG    ....................................    CCTCGGTCTGCTCTAACAGATCCCCCAAGT
     16840      36   100.0      30  TCCCCCAAGT    ....................................    ACAGAGAAAGAAAGAGAGATTAACGACTAC
     16906      36   100.0      30  TAACGACTAC    ....................................    TGAAACGGAGTGGACAGGTAAAGGAATGGG
     16972      36   100.0      30  AAGGAATGGG    ....................................    TGCGGTCCCTTGGTTCCGTCAACAACATCA
     17038      36   100.0      30  AACAACATCA    ....................................    TGTCCTATTCCCTTTTATGCTGCGTGTATA
     17104      36   100.0      30  TGCGTGTATA    ....................................    AATACAAGCATAAAGAACGAACCGCAACGG
     17170      36   100.0          ACCGCAACGG    ....................................    AGGGAA
==========  ======  ======  ======  ==========    ====================================
        10      36              30                GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT


";
        let (_, actual) = parse_array(input).unwrap();
        assert_eq!(actual.repeat_spacers.len(), 10);
        assert_eq!(actual.accession, "MGYG000273829_14");
        assert_eq!(
            actual.consensus_repeat_sequence,
            "GCTGTAGTTCCCGGTTATTACTTGGTATGTTATAAT"
        );
        assert_eq!(actual.repeat_spacers[0].start, 16575);
        assert_eq!(actual.repeat_spacers[9].start, 17169);
    }

    #[test]
    fn test_parse_array_with_gaps() {
        let input = "Array 18
>MGYG000232241_150

       Pos  Repeat     %id  Spacer  Left flank    Repeat                                      Spacer
==========  ======  ======  ======  ==========    ========================================    ======
      3832      40    92.5      34  CATATAGCAA    ..A..................................CC.    GAATTACATCGTATGCCAATACGCAGTTGCTTTT
      3906      40    97.5      41  AGTTGCTTTT    .....................................---    TGTACTACTATGCGGTATTCCATCTGAAGGATGGCGGCTAC
      3987      40    92.5          TGGCGGCTAC    GG............-......................--.    ATCACATTCA
==========  ======  ======  ======  ==========    ========================================
         3      40              37                AAGTTTCCGTCCCCTTTCGGGGAATCATTTAGAAAAT--A


";
        let expected = Array {
            accession: "MGYG000232241_150",
            consensus_repeat_sequence: "AAGTTTCCGTCCCCTTTCGGGGAATCATTTAGAAAAT--A",
            start: 3831,
            end: 4030,
            order: 17,
            repeat_spacers: vec![
                RepeatSpacer {
                    start: 3831,
                    end: 3905,
                    repeat_start: 3831,
                    repeat_end: 3871,
                    spacer_start: 3871,
                    spacer_end: 3905,
                    spacer: "GAATTACATCGTATGCCAATACGCAGTTGCTTTT",
                    repeat: "AAATTTCCGTCCCCTTTCGGGGAATCATTTAGAAAATCCA".to_string(),
                },
                RepeatSpacer {
                    start: 3905,
                    end: 3983,
                    repeat_start: 3905,
                    repeat_end: 3942,
                    spacer_start: 3942,
                    spacer_end: 3983,
                    spacer: "TGTACTACTATGCGGTATTCCATCTGAAGGATGGCGGCTAC",
                    repeat: "AAGTTTCCGTCCCCTTTCGGGGAATCATTTAGAAAAT".to_string(),
                },
                RepeatSpacer {
                    start: 3983,
                    end: 4030,
                    repeat_start: 3983,
                    repeat_end: 4020,
                    spacer_start: 4020,
                    spacer_end: 4030,
                    spacer: "ATCACATTCA",
                    repeat: "GGGTTTCCGTCCCCTTCGGGGAATCATTTAGAAAATA".to_string(),
                },
            ],
        };
        let (_, actual) = parse_array(input).unwrap();
        assert_eq!(expected, actual);
    }
}
