use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use kseq::parse_path;
use serde::Deserialize;
use structopt::StructOpt;
use indicatif::{ProgressBar, ProgressStyle};



#[derive(Debug, StructOpt)]
struct Opt {
    #[structopt(short = "fq", long = "fastq", help = "FASTQ file")]
    fastq: Option<String>,

    #[structopt(short = "u", long = "unblocks", help = "Unblocks file")]
    unblocks: Option<String>,

    #[structopt(short = "c", long = "toml", help = "Channels TOML file")]
    toml: String,

    #[structopt(short = "n", long = "num_reads", help = "number of reads in the fastq file")]
    num_reads: u64,
}

#[derive(Debug, Deserialize)]
struct TomlData {
    conditions: HashMap<String, Condition>,
}

#[derive(Debug, Deserialize)]
struct Condition {
    name: String,
    channels: Vec<i32>,
}

fn channels_toml(toml_file: &str) -> io::Result<HashMap<i32, String>> {
    let toml_data: TomlData = toml::from_slice(&std::fs::read(toml_file)?)?;
    let mut channels = HashMap::new();

    for (_condition, cond_data) in toml_data.conditions {
        for &channel in &cond_data.channels {
            channels.insert(channel, cond_data.name.clone());
        }
    }

    Ok(channels)
}

fn run(opts: Opt) -> io::Result<()> {

    let total = opts.num_reads;
    // Parse the TOML file
    let channels = channels_toml(&opts.toml)?;
    let unblocked = if let Some(unblock_file) = opts.unblocks {
        std::fs::read_to_string(unblock_file)?.lines().map(String::from).collect::<Vec<_>>()
    } else {
        vec![]
    };

    let mut as_unblocked = vec![];
    let mut as_sequenced = vec![];
    let mut control_sequenced = vec![];
    let bar = ProgressBar::new(total);
    bar.set_style(ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}[{eta_precise}]")
        .unwrap()
        .progress_chars("##-"));

    if let Some(fastq_file) = opts.fastq {
        let mut total_reads = 0;
        let mut processed_reads = 0;
        
        let mut records = parse_path(fastq_file).unwrap();

        // Create an iterator that processes each record
        while let Some(record) = records.iter_record().unwrap()  {
            bar.inc(1);
            // Split the description by whitespace
            let desc: Vec<&str> = record.des().split_whitespace().collect();

        // Find the element that starts with "ch="
            if let Some(ch_element) = desc.iter().find(|&&s| s.starts_with("ch=")) {
                let ch = &ch_element["ch=".len()..];
                // convert str to i32
                let ch = ch.parse::<i32>().unwrap_or(0);

            let id = record.head().to_string();
        //for (id, ch) in get_read_info_kseq(&fastq_file)? {
            total_reads += 1;
            if unblocked.contains(&id) {
                if channels.get(&ch).map_or(false, |name| name != "control condtion 1") {
                    as_unblocked.push(id);
                } else {
                    control_sequenced.push(id);
                }
            } else {
                if channels.get(&ch).map_or(false, |name| name != "control condtion 1") {
                    as_sequenced.push(id);
                } else {
                    control_sequenced.push(id);
                }
            }
            processed_reads += 1;
        }

        
    }
    
    println!("Total reads processed: {}", total_reads);
    println!("Reads after filtering: {}", processed_reads);
    
    write_to_file("AS_unblocked.txt", &as_unblocked)?;
    write_to_file("AS_sequenced.txt", &as_sequenced)?;
    write_to_file("control_sequenced.txt", &control_sequenced)?;

    Ok(())
    } else {
        Ok(())
    }
}


fn write_to_file<P: AsRef<Path>>(path: P, data: &[String]) -> io::Result<()> {
    let mut file = File::create(path)?;
    for line in data {
        writeln!(file, "{}", line)?;
    }
    Ok(())
}

fn main() -> io::Result<()> {
    let opts = Opt::from_args();
    run(opts)
}
