use nom::{
    bytes::complete::tag,
    character::complete::{
        alpha0, char, digit0, digit1, line_ending, multispace0, multispace1, not_line_ending,
    },
    multi::many1,
    number::complete::float,
    sequence::{pair, tuple},
    IResult, InputTakeAtPosition,
};

#[derive(Debug, PartialEq)]
struct RawRepeatSpacer<'a> {
    start: usize,
    end: usize,
    repeat_diff: &'a str,
    spacer: &'a str,
}

#[derive(Debug, PartialEq)]
struct RepeatSpacer<'a> {
    start: usize,
    end: usize,
    repeat: String,
    spacer: &'a str,
}

#[derive(Debug, PartialEq)]
struct Array<'a> {
    accession: &'a str,
    order: usize,
    start: usize,
    end: usize,
    consensus_repeat_sequence: &'a str,
    repeat_spacers: Vec<RepeatSpacer<'a>>,
}

fn not_space(input: &str) -> IResult<&str, &str> {
    input.split_at_position_complete(char::is_whitespace)
}

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
                    spacer: "GAATTACATCGTATGCCAATACGCAGTTGCTTTT",
                    repeat: "AAATTTCCGTCCCCTTTCGGGGAATCATTTAGAAAATCCA".to_string(),
                },
                RepeatSpacer {
                    start: 3905,
                    end: 3983,
                    spacer: "TGTACTACTATGCGGTATTCCATCTGAAGGATGGCGGCTAC",
                    repeat: "AAGTTTCCGTCCCCTTTCGGGGAATCATTTAGAAAAT".to_string(),
                },
                RepeatSpacer {
                    start: 3983,
                    end: 4030,
                    spacer: "ATCACATTCA",
                    repeat: "GGGTTTCCGTCCCCTTCGGGGAATCATTTAGAAAATA".to_string(),
                },
            ],
        };
        let (_, actual) = parse_array(input).unwrap();
        assert_eq!(expected, actual);
    }
}
