use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use bio::io::fastq;
use flate2::read::GzDecoder;
use serde::Deserialize;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
struct Opt {
    #[structopt(short = "fq", long = "fastq", help = "FASTQ file")]
    fastq: Option<String>,

    #[structopt(short = "u", long = "unblocks", help = "Unblocks file")]
    unblocks: Option<String>,

    #[structopt(short = "c", long = "toml", help = "Channels TOML file")]
    toml: String,
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
    let file = File::open(toml_file)?;
    let toml_data: TomlData = toml::from_slice(&std::fs::read(file)?)?;
    let mut channels = HashMap::new();

    for (condition, cond_data) in toml_data.conditions {
        for &channel in &cond_data.channels {
            channels.insert(channel, cond_data.name.clone());
        }
    }

    Ok(channels)
}

fn get_read_info(fq: &str) -> io::Result<impl Iterator<Item = (String, i32, String)>> {
    let file = File::open(fq)?;
    let gz = GzDecoder::new(file);
    let reader = BufReader::new(gz);
    let fastq_reader = fastq::Reader::new(reader);

    let iter = fastq_reader.records().filter_map(|rec| {
        if let Ok(record) = rec {
            let desc = record.desc().unwrap_or("").split_whitespace();
            let mut desc_dict = HashMap::new();

            for d in desc {
                if let Some((key, value)) = d.split_once('=') {
                    desc_dict.insert(key, value);
                }
            }

            if let (Some(ch), Some(start_time)) = (desc_dict.get("ch"), desc_dict.get("start_time")) {
                return Some((record.id().to_string(), ch.parse::<i32>().unwrap_or(0), start_time.to_string()));
            }
        }
        None
    });

    Ok(iter)
}

fn run(opts: Opt) -> io::Result<()> {
    let channels = channels_toml(&opts.toml)?;
    let unblocked = if let Some(unblock_file) = opts.unblocks {
        std::fs::read_to_string(unblock_file)?.lines().map(String::from).collect::<Vec<_>>()
    } else {
        vec![]
    };

    let mut as_unblocked = vec![];
    let mut as_sequenced = vec![];
    let mut control_sequenced = vec![];

    if let Some(fastq_file) = opts.fastq {
        for (id, ch, _start_time) in get_read_info(&fastq_file)? {
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
        }
    }

    write_to_file("AS_unblocked.txt", &as_unblocked)?;
    write_to_file("AS_sequenced.txt", &as_sequenced)?;
    write_to_file("control_sequenced.txt", &control_sequenced)?;

    Ok(())
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
