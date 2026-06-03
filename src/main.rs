use needletail::{parse_fastx_file, FastxReader};
use clap::Parser;
use rayon::prelude::*;

#[derive(Parser)]
#[command(about = "Compute length-weighted mean and median coverage from a FASTA/FASTQ file")]
struct Args {
    #[arg(long = "in", help = "Input FASTA/FASTQ file")]
    input: String,
    #[arg(long, short = 't', help = "Number of threads (default: all available)")]
    threads: Option<usize>,
}

fn find_ka_index(header: &[u8]) -> Option<usize> {
    std::str::from_utf8(header)
        .ok()?
        .split_whitespace()
        .position(|token| token.starts_with("ka:f:"))
}

fn extract_coverage_at(header: &[u8], index: usize) -> Option<f64> {
    std::str::from_utf8(header)
        .ok()?
        .split_whitespace()
        .nth(index)?
        .strip_prefix("ka:f:")?
        .parse()
        .ok()
}

struct RecordIter {
    reader: Box<dyn FastxReader>,
    ka_index: Option<usize>,
}

impl Iterator for RecordIter {
    type Item = (f64, usize); // (coverage, sequence_length)

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let seqrec = self.reader.next()?.expect("Invalid record");
            let cov = if let Some(idx) = self.ka_index {
                extract_coverage_at(seqrec.id(), idx)
            } else {
                let idx = find_ka_index(seqrec.id())?;
                self.ka_index = Some(idx);
                extract_coverage_at(seqrec.id(), idx)
            };
            if let Some(c) = cov {
                return Some((c, seqrec.num_bases()));
            }
        }
    }
}

fn sequence_median(sorted_pairs: &[(f64, usize)]) -> f64 {
    let n = sorted_pairs.len();
    if n % 2 == 1 {
        sorted_pairs[n / 2].0
    } else {
        (sorted_pairs[n / 2 - 1].0 + sorted_pairs[n / 2].0) / 2.0
    }
}

// Length-weighted median: coverage at which cumulative base count crosses 50% of total.
fn weighted_median(sorted_pairs: &[(f64, usize)], total_len: usize) -> f64 {
    let half = (total_len + 1) / 2;
    let mut cumulative = 0usize;
    for &(cov, len) in sorted_pairs {
        cumulative += len;
        if cumulative >= half {
            return cov;
        }
    }
    sorted_pairs.last().map(|&(c, _)| c).unwrap_or(0.0)
}

fn main() {
    let args = Args::parse();

    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to build thread pool");
    }

    let reader = parse_fastx_file(&args.input).expect("Cannot open input file");

    let mut pairs: Vec<(f64, usize)> = RecordIter { reader, ka_index: None }
        .par_bridge()
        .collect();

    if pairs.is_empty() {
        eprintln!("No sequences with ka:f: coverage tag found.");
        std::process::exit(1);
    }

    let total_len: usize = pairs.par_iter().map(|&(_, len)| len).sum();
    let weighted_sum: f64 = pairs.par_iter().map(|&(cov, len)| cov * len as f64).sum();

    pairs.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    println!("Number of sequences:          {}", pairs.len());
    println!("Total length:                 {} bp", total_len);
    println!("Mean coverage:                {:.3}", weighted_sum / total_len as f64);
    println!("Median coverage (sequences):  {:.3}", sequence_median(&pairs));
    println!("Median coverage (bases):      {:.3}", weighted_median(&pairs, total_len));
}
