use onecode::OneFile;
use std::collections::HashMap;
use std::io::{self, Write};
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "oneview-rs")]
#[command(about = "View alignments from ONE format files", long_about = None)]
struct Args {
    /// Input .1aln file path
    #[arg(value_name = "FILE")]
    input: String,

    /// Alignment number to read (0-indexed)
    #[arg(short, long, value_name = "NUM")]
    alignment: Option<usize>,
    
    /// Print only file metadata (sequences and trace spacing)
    #[arg(short, long)]
    metadata: bool,

    /// Emit alignments in PAF format
    #[arg(long)]
    paf: bool,
}

#[derive(Debug, Default)]
struct AlignmentData {
    query_name: String,
    query_length: i64,
    query_start: i64,
    query_end: i64,
    target_name: String,
    target_length: i64,
    target_start: i64,
    target_end: i64,
    strand: char,
    differences: i64,
    tracepoints: Vec<i64>,
    trace_diffs: Vec<i64>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    if args.metadata && args.paf {
        return Err("Cannot combine --metadata with --paf output".into());
    }

    let output_format = if args.paf {
        OutputFormat::Paf
    } else {
        OutputFormat::Human
    };
    
    let (metadata, trace_spacing) = get_file_metadata(&args.input)?;
    
    match (args.metadata, args.alignment) {
        (true, _) => {
            // Only metadata
            print_metadata(&metadata, trace_spacing, &args.input)?;
        }
        (false, Some(idx)) => {
            // Only specific alignment
            read_single_alignment(
                &args.input,
                idx,
                &metadata,
                trace_spacing,
                output_format,
            )?;
        }
        (false, None) => {
            // Default: metadata + all alignments
            if output_format == OutputFormat::Human {
                print_metadata(&metadata, trace_spacing, &args.input)?;
                writeln!(io::stdout(), "\n=== ALIGNMENTS ===\n")?;
            }
            read_all_alignments(
                &args.input,
                &metadata,
                trace_spacing,
                output_format,
            )?;
        }
    }
    Ok(())
}

struct FileMetadata {
    seq_names: HashMap<i64, String>,
    seq_lengths: HashMap<i64, i64>,
}

fn get_file_metadata(path: &str) -> Result<(FileMetadata, i64), Box<dyn std::error::Error>> {
    let mut file = OneFile::open_read(path, None, None, 1)?;
    let metadata = FileMetadata {
        seq_names: file.get_all_sequence_names(),
        seq_lengths: file.get_all_sequence_lengths(),
    };
    
    // Get trace spacing
    let mut trace_spacing = 100; // default
    loop {
        match file.read_line() {
            't' => { trace_spacing = file.int(0); break; }
            'A' | '\0' => break,
            _ => {}
        }
    }
    
    Ok((metadata, trace_spacing))
}

fn print_metadata(
    metadata: &FileMetadata,
    trace_spacing: i64,
    path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    
    writeln!(handle, "=== METADATA ===\n")?;
    writeln!(handle, "File: {}", path)?;
    writeln!(handle, "Trace spacing: {}", trace_spacing)?;
    writeln!(handle, "Number of sequences: {}", metadata.seq_names.len())?;
    
    // Sort sequences by ID for consistent output
    let mut sequences: Vec<_> = metadata.seq_names.iter().collect();
    sequences.sort_by_key(|&(&id, _)| id);
    
    writeln!(handle, "\nSequences:")?;
    for (&id, name) in sequences {
        let length = metadata.seq_lengths.get(&id).copied().unwrap_or(0);
        writeln!(handle, "  {}) {} (length: {})", id, name, length)?;
    }
    
    // Count alignments
    write!(handle, "\nCounting alignments...")?;
    handle.flush()?;
    
    let mut file = OneFile::open_read(path, None, None, 1)?;
    let mut count = 0;
    loop {
        match file.read_line() {
            '\0' => break,
            'A' => count += 1,
            _ => {}
        }
    }
    writeln!(handle, "\rTotal alignments: {}    ", count)?;
    
    Ok(())
}

#[derive(Copy, Clone, Eq, PartialEq)]
enum OutputFormat {
    Human,
    Paf,
}

fn read_single_alignment(
    path: &str,
    idx: usize,
    metadata: &FileMetadata,
    trace_spacing: i64,
    format: OutputFormat,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = OneFile::open_read(path, None, None, 1)?;

    // Require O(1) access via binary index
    if file.goto('A', (idx + 1) as i64).is_err() {
        return Err(format!(
            "Cannot access alignment {} directly. Binary index not available for this file.\n\
             Please ensure the file has an associated .1idx index file.",
            idx
        ).into());
    }

    eprintln!("Using O(1) binary index to jump to alignment {}", idx);
    file.read_line(); // Read the 'A' line we jumped to
    let (aln, _) = parse_alignment(&mut file, metadata)?;

    print_alignment(&aln, trace_spacing, format)?;
    Ok(())
}

fn read_all_alignments(
    path: &str,
    metadata: &FileMetadata,
    trace_spacing: i64,
    format: OutputFormat,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = OneFile::open_read(path, None, None, 1)?;

    let mut current_line = file.read_line();
    loop {
        match current_line {
            '\0' => break,
            'A' => {
                let (aln, next_line) = parse_alignment(&mut file, metadata)?;
                print_alignment(&aln, trace_spacing, format)?;
                current_line = next_line;
            }
            _ => {
                current_line = file.read_line();
            }
        }
    }
    Ok(())
}

fn parse_alignment(
    file: &mut OneFile,
    metadata: &FileMetadata,
) -> Result<(AlignmentData, char), Box<dyn std::error::Error>> {
    // Read alignment coordinates from current 'A' line
    let query_id = file.int(0);
    let target_id = file.int(3);

    let query_name = metadata
        .seq_names
        .get(&query_id)
        .cloned()
        .ok_or_else(|| format!("Query sequence with ID {} not found in metadata", query_id))?;
    let target_name = metadata
        .seq_names
        .get(&target_id)
        .cloned()
        .ok_or_else(|| format!("Target sequence with ID {} not found in metadata", target_id))?;
    let query_length = metadata
        .seq_lengths
        .get(&query_id)
        .copied()
        .ok_or_else(|| {
            format!(
                "Query sequence length for ID {} not found in metadata",
                query_id
            )
        })?;
    let target_length = metadata
        .seq_lengths
        .get(&target_id)
        .copied()
        .ok_or_else(|| {
            format!(
                "Target sequence length for ID {} not found in metadata",
                target_id
            )
        })?;

    let mut aln = AlignmentData {
        query_name,
        query_length,
        query_start: file.int(1),
        query_end: file.int(2),
        target_name,
        target_length,
        target_start: file.int(4),
        target_end: file.int(5),
        strand: '+',
        ..Default::default()
    };

    // Read associated lines
    let next_line = loop {
        let line_type = file.read_line();
        match line_type {
            'R' => aln.strand = '-',
            'D' => aln.differences = file.int(0),
            'T' => aln.tracepoints = file.int_list().map(|v| v.to_vec()).unwrap_or_default(),
            'X' => aln.trace_diffs = file.int_list().map(|v| v.to_vec()).unwrap_or_default(),
            'A' | 'a' | 'g' | 'S' | '^' | '\0' => break line_type,
            _ => {}
        }
    };

    if matches!(aln.strand, '-' | '\'') {
        let orig_start = aln.target_start;
        let orig_end = aln.target_end;
        // Reverse-complement target coordinates so start/end reflect forward strand
        aln.target_start = aln.target_length - orig_end;
        aln.target_end = aln.target_length - orig_start;
    }

    Ok((aln, next_line))
}

fn print_alignment(aln: &AlignmentData, trace_spacing: i64, format: OutputFormat) -> io::Result<()> {
    match format {
        OutputFormat::Human => print_alignment_human(aln, trace_spacing),
        OutputFormat::Paf => print_alignment_paf(aln),
    }
}

fn print_alignment_human(aln: &AlignmentData, trace_spacing: i64) -> io::Result<()> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    
    writeln!(handle, "Query: {}:{}-{}, query total length: {}", 
        aln.query_name, aln.query_start, aln.query_end, aln.query_length)?;
    writeln!(handle, "Target: {}:{}-{}, target total length: {}", 
        aln.target_name, aln.target_start, aln.target_end, aln.target_length)?;
    writeln!(handle, "Strand: {}", aln.strand)?;
    writeln!(handle, "Differences: {}", aln.differences)?;
    writeln!(handle, "Trace spacing: {}", trace_spacing)?;
    
    print_trace_data(&mut handle, "Tracepoints", &aln.tracepoints)?;
    print_trace_data(&mut handle, "Trace diffs", &aln.trace_diffs)?;

    writeln!(handle)?;
    Ok(())
}

fn print_alignment_paf(aln: &AlignmentData) -> io::Result<()> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    let query_span = (aln.query_end - aln.query_start).max(0);
    let target_span = (aln.target_end - aln.target_start).max(0);
    let block_length = query_span.max(target_span);
    let matches = (block_length - aln.differences).max(0);
    let mapq = 255;

    write!(
        handle,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        aln.query_name,
        aln.query_length,
        aln.query_start,
        aln.query_end,
        aln.strand,
        aln.target_name,
        aln.target_length,
        aln.target_start,
        aln.target_end,
        matches,
        block_length,
        mapq
    )?;

    write!(handle, "\tdf:i:{}", aln.differences)?;

    let tp_pairs = if !aln.trace_diffs.is_empty() && !aln.tracepoints.is_empty() {
        let pair_count = aln.trace_diffs.len().min(aln.tracepoints.len());
        aln.trace_diffs
            .iter()
            .zip(aln.tracepoints.iter())
            .take(pair_count)
            .map(|(&diff, &tp)| format!("{},{}", diff, tp))
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };
    if !tp_pairs.is_empty() {
        write!(handle, "\ttp:Z:{}", tp_pairs.join(";"))?;
    }

    writeln!(handle)?;
    Ok(())
}

fn print_trace_data(handle: &mut io::StdoutLock, label: &str, data: &[i64]) -> io::Result<()> {
    if !data.is_empty() {
        writeln!(handle, "{}: {} values", label, data.len())?;
        write!(handle, "  ")?;
        for (i, &val) in data.iter().enumerate() {
            if i > 0 { write!(handle, " ")?; }
            write!(handle, "{}", val)?;
        }
        writeln!(handle)?;
    } else {
        writeln!(handle, "{}: 0 values", label)?;
    }
    Ok(())
}
