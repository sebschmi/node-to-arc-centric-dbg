use clap::Parser;
use genome_graph::compact_genome::implementation::DefaultSequenceStore;
use genome_graph::compact_genome::interface::alphabet::dna_alphabet::DnaAlphabet;
use genome_graph::io::bcalm2::read_bigraph_from_bcalm2_as_edge_centric_from_file;
use genome_graph::types::PetBCalm2EdgeGraph;
use log::{info, LevelFilter};
use simplelog::{ColorChoice, CombinedLogger, Config, TermLogger, TerminalMode};
use std::path::PathBuf;

#[derive(Parser, Debug)]
struct Cli {
    /// The input file containing a node-centric de Bruijn graph.
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
    let _graph: PetBCalm2EdgeGraph<_> =
        read_bigraph_from_bcalm2_as_edge_centric_from_file(&cli.input, &mut sequence_store, cli.k)
            .unwrap();

    println!("Hello, world!");
}
