use clap::Parser;
use genome_graph::bigraph::traitgraph::index::GraphIndex;
use genome_graph::bigraph::traitgraph::interface::{
    ImmutableGraphContainer, NavigableGraph, Neighbor,
};
use genome_graph::bigraph::traitgraph::traitsequence::interface::Sequence;
use genome_graph::compact_genome::implementation::DefaultSequenceStore;
use genome_graph::compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use genome_graph::compact_genome::interface::sequence::GenomeSequence;
use genome_graph::compact_genome::interface::sequence_store::SequenceStore;
use genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_edge_centric;
use genome_graph::types::PetBCalm2EdgeGraph;
use log::{info, warn, LevelFilter};
use simplelog::{ColorChoice, CombinedLogger, Config, TermLogger, TerminalMode};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::PathBuf;

#[derive(Parser, Debug)]
struct Cli {
    /// The input file containing a node-centric de Bruijn graph.
    /// The file should be in bcalm2 format.
    #[clap(long)]
    input: PathBuf,

    /// The k-mer size used to generate the de Bruijn graph.
    #[clap(short)]
    k: usize,

    /// The output file where the arc-centric de Bruijn graph should be written to.
    #[clap(long)]
    output: PathBuf,

    /// The desired log level.
    #[clap(long, default_value = "Info")]
    log_level: LevelFilter,
}

pub fn initialise_logging(log_level: LevelFilter) {
    CombinedLogger::init(vec![TermLogger::new(
        log_level,
        Config::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )])
    .unwrap();

    info!("Logging initialised successfully");
}

fn node_to_arc_centric_dbg(k: usize, input: &mut impl BufRead, output: &mut impl Write) {
    info!("Reading graph");
    let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::new();
    let graph: PetBCalm2EdgeGraph<_> =
        read_bigraph_from_bcalm2_as_edge_centric(input, &mut sequence_store, k).unwrap();

    info!("Writing graph...");
    output_arc_centric_dbg(&graph, &sequence_store, k, output);
}

fn output_arc_centric_dbg(
    graph: &PetBCalm2EdgeGraph<
        <DefaultSequenceStore<DnaAlphabet> as SequenceStore<DnaAlphabet>>::Handle,
    >,
    sequence_store: &DefaultSequenceStore<DnaAlphabet>,
    k: usize,
    output: &mut impl Write,
) {
    writeln!(output, "{}", graph.node_count()).unwrap();
    for n1 in graph.node_indices() {
        let mut neighbors: Vec<_> = graph.out_neighbors(n1).collect();
        neighbors.sort_unstable_by_key(|neighbor| neighbor.node_id);

        for Neighbor {
            node_id: n2,
            edge_id,
        } in neighbors
        {
            let edge_data = graph.edge_data(edge_id);
            let kmer_count = edge_data.length - (k - 1);
            if edge_data.total_abundance % kmer_count != 0 {
                let sequence = sequence_store.get(&edge_data.sequence_handle);
                let sequence = &sequence[..(k + 10).min(sequence.len())];
                warn!(
                    "Found edge with non-integer average abundance: {}",
                    sequence.as_string()
                );
            }

            let n1 = n1.as_usize();
            let n2 = n2.as_usize();
            let weight = edge_data.total_abundance / kmer_count;
            write!(output, "{n1} {n2} {weight} ").unwrap();

            let sequence = sequence_store.get(&edge_data.sequence_handle);
            if edge_data.forwards {
                for character in sequence.iter() {
                    write!(output, "{}", character).unwrap();
                }
            } else {
                for character in sequence.reverse_complement_iter() {
                    write!(output, "{}", character).unwrap();
                }
            }
            writeln!(output).unwrap();
        }
    }
}

fn main() {
    let cli = Cli::parse();
    initialise_logging(cli.log_level);

    info!(
        "Loading graph from {:?} with k = {} and writing to {:?}",
        cli.input, cli.k, cli.output
    );
    let mut input = BufReader::new(File::open(&cli.input).unwrap());
    let mut output = BufWriter::new(File::create(&cli.output).unwrap());
    node_to_arc_centric_dbg(cli.k, &mut input, &mut output);

    info!("Success!");
}

#[cfg(test)]
mod tests {
    use crate::node_to_arc_centric_dbg;
    use std::io::BufReader;

    #[test]
    fn test_complex_file() {
        let mut file = BufReader::new(
            ">0 LN:i:14 KC:i:21 km:f:21.0   L:-:2:+  L:+:2:+
ATCGATCGATCGAT
>1 LN:i:14 KC:i:20 km:f:20.0   L:-:2:-  L:+:2:-
CGATCGATCGATCG
>2 LN:i:14 KC:i:43 km:f:43.0   L:+:1:+ L:+:1:- L:+:3:+  L:-:0:+ L:-:0:-
TCGATCGATCGATC
>3 LN:i:16 KC:i:3 km:f:1.0   L:-:2:-
CGATCGATCGATCAGT"
                .as_bytes(),
        );

        let expected = "6
0 1 21 ATCGATCGATCGAT
0 1 21 ATCGATCGATCGAT
1 2 43 TCGATCGATCGATC
2 3 20 CGATCGATCGATCG
2 3 20 CGATCGATCGATCG
2 4 1 CGATCGATCGATCAGT
3 0 43 GATCGATCGATCGA
5 3 1 ACTGATCGATCGATCG
";

        let mut output = Vec::new();
        node_to_arc_centric_dbg(14, &mut file, &mut output);
        println!("{}", String::from_utf8(output.clone()).unwrap());

        assert_eq!(expected.as_bytes(), output);
    }

    #[test]
    fn test_complex_circularised_file() {
        let mut file = BufReader::new(
            ">0 LN:i:14 KC:i:20 km:f:20.0   L:-:1:-  L:+:1:-
CGATCGATCGATCG
>1 LN:i:14 KC:i:43 km:f:43.0   L:+:0:+ L:+:0:- L:+:4:+ L:+:5:+  L:-:2:+ L:-:2:- L:-:5:-
TCGATCGATCGATC
>2 LN:i:14 KC:i:21 km:f:21.0   L:-:1:+  L:+:1:+
ATCGATCGATCGAT
>3 LN:i:27 KC:i:14 km:f:1.0   L:-:4:-  L:+:4:-
GATCGATCGATCAGTGATCGATCGATC
>4 LN:i:14 KC:i:2 km:f:2.0   L:-:1:-  L:+:3:+ L:+:3:-
CGATCGATCGATCA
>5 LN:i:26 KC:i:13 km:f:1.0   L:-:1:-  L:+:1:+
CGATCGATCGATCTCGATCGATCGAT"
                .as_bytes(),
        );

        let expected = "6
0 1 20 CGATCGATCGATCG
0 1 20 CGATCGATCGATCG
0 2 1 CGATCGATCGATCTCGATCGATCGAT
0 4 2 CGATCGATCGATCA
1 3 43 GATCGATCGATCGA
2 0 43 TCGATCGATCGATC
3 1 1 ATCGATCGATCGAGATCGATCGATCG
3 2 21 ATCGATCGATCGAT
3 2 21 ATCGATCGATCGAT
4 5 1 GATCGATCGATCACTGATCGATCGATC
4 5 1 GATCGATCGATCAGTGATCGATCGATC
5 1 2 TGATCGATCGATCG
";

        let mut output = Vec::new();
        node_to_arc_centric_dbg(14, &mut file, &mut output);
        println!("{}", String::from_utf8(output.clone()).unwrap());

        assert_eq!(expected.as_bytes(), output);
    }

    #[test]
    fn test_bcalm2_multiplicities() {
        let strings = vec![
            "ACTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC",
            "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA",
        ];
        let strings: Vec<_> = strings
            .iter()
            .flat_map(|&string| {
                [
                    string.to_owned(),
                    string
                        .chars()
                        .rev()
                        .map(|c| match c {
                            'A' => 'T',
                            'C' => 'G',
                            'G' => 'C',
                            'T' => 'A',
                            c => panic!("{c}"),
                        })
                        .collect(),
                ]
            })
            .collect();

        let requirements = [
            ("ATCGATCGATCGAT", 42),
            ("TCGATCGATCGATC", 43),
            ("CGATCGATCGATCG", 40),
            ("CGATCGATCGATCAGT", 1),
        ];

        for (kmer, count) in requirements {
            for offset in 0..kmer.len() - 13 {
                let kmer = &kmer[offset..offset + 14];
                let mut actual_count = 0;
                for string in &strings {
                    let mut string = string.as_str();
                    while let Some(index) = string.find(kmer) {
                        actual_count += 1;
                        string = &string[(index + 1).min(string.len())..];
                    }
                }
                assert_eq!(
                    count, actual_count,
                    "kmer {kmer} should match {count} times but matches {actual_count} times"
                );
            }
        }
    }

    #[test]
    fn test_bcalm2_circularised_multiplicities() {
        let strings = vec![
            "ACTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCACTGATCGATCGA",
            "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAGATCGATCGATCG",
        ];
        let strings: Vec<_> = strings
            .iter()
            .flat_map(|&string| {
                [
                    string.to_owned(),
                    string
                        .chars()
                        .rev()
                        .map(|c| match c {
                            'A' => 'T',
                            'C' => 'G',
                            'G' => 'C',
                            'T' => 'A',
                            c => panic!("{c}"),
                        })
                        .collect(),
                ]
            })
            .collect();
        for string in &strings {
            println!("{string}");
        }

        let requirements = [
            ("ATCGATCGATCGAT", 42),
            ("CGATCGATCGATCG", 40),
            ("TCGATCGATCGATC", 43),
            ("CGATCGATCGATCA", 2),
            ("GATCGATCGATCAGTGATCGATCGATC", 1),
            ("CGATCGATCGATCTCGATCGATCGAT", 1),
        ];

        for (kmer, count) in requirements {
            for offset in 0..kmer.len() - 13 {
                let kmer = &kmer[offset..offset + 14];
                let mut actual_count = 0;
                for string in &strings {
                    let mut string = string.as_str();
                    while let Some(index) = string.find(kmer) {
                        actual_count += 1;
                        string = &string[(index + 1).min(string.len())..];
                    }
                }
                assert_eq!(
                    count, actual_count,
                    "kmer {kmer} should match {count} times but matches {actual_count} times"
                );
            }
        }
    }
}
