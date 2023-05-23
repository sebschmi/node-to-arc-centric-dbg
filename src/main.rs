use clap::Parser;
use genome_graph::bigraph::traitgraph::index::GraphIndex;
use genome_graph::bigraph::traitgraph::interface::{
    ImmutableGraphContainer, NavigableGraph, Neighbor,
};
use genome_graph::bigraph::traitgraph::traitsequence::interface::Sequence;
use genome_graph::compact_genome::implementation::DefaultSequenceStore;
use genome_graph::compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use genome_graph::compact_genome::interface::sequence_store::SequenceStore;
use genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_edge_centric_from_file;
use genome_graph::types::PetBCalm2EdgeGraph;
use log::{info, warn, LevelFilter};
use simplelog::{ColorChoice, CombinedLogger, Config, TermLogger, TerminalMode};
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
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

fn main() {
    let cli = Cli::parse();
    initialise_logging(cli.log_level);

    let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::new();
    let graph: PetBCalm2EdgeGraph<_> =
        read_bigraph_from_bcalm2_as_edge_centric_from_file(&cli.input, &mut sequence_store, cli.k)
            .unwrap();

    let mut output = BufWriter::new(File::create(&cli.output).unwrap());
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
            let kmer_count = edge_data.length - (cli.k - 1);
            if edge_data.total_abundance % kmer_count != 0 {
                warn!("Found edge with non-integer average abundance");
            }

            let n1 = n1.as_usize();
            let n2 = n2.as_usize();
            let weight = edge_data.total_abundance / kmer_count;
            write!(output, "{n1} {n2} {weight} ").unwrap();

            let sequence = sequence_store.get(&edge_data.sequence_handle);
            for character in sequence.iter() {
                write!(output, "{}", character).unwrap();
            }
            writeln!(output).unwrap();
        }
    }

    println!("Hello, world!");
}
