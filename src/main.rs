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
    
    let (metadata, trace_spacing) = get_file_metadata(&args.input)?;
    
    match (args.metadata, args.alignment) {
        (true, _) => {
            // Only metadata
            print_metadata(&metadata, trace_spacing, &args.input)?;
        }
        (false, Some(idx)) => {
            // Only specific alignment
            read_single_alignment(&args.input, idx, &metadata, trace_spacing)?;
        }
        (false, None) => {
            // Default: metadata + all alignments
            print_metadata(&metadata, trace_spacing, &args.input)?;
            writeln!(io::stdout(), "\n=== ALIGNMENTS ===\n")?;
            read_all_alignments(&args.input, &metadata, trace_spacing)?;
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

fn read_single_alignment(
    path: &str,
    idx: usize,
    metadata: &FileMetadata,
    trace_spacing: i64,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = OneFile::open_read(path, None, None, 1)?;
    
    // Try O(1) access first, fall back to sequential scan
    let aln = if file.goto('A', (idx + 1) as i64).is_ok() {
        eprintln!("Using O(1) binary index to jump to alignment {}", idx);
        file.read_line(); // Read the 'A' line we jumped to
        parse_alignment(&mut file, metadata)?
    } else {
        eprintln!("Binary index unavailable, scanning sequentially to alignment {}...", idx);
        find_nth_alignment(&mut file, idx, metadata)?
    };
    
    print_alignment(&aln, trace_spacing)?;
    Ok(())
}

fn find_nth_alignment(
    file: &mut OneFile,
    target_idx: usize,
    metadata: &FileMetadata,
) -> Result<AlignmentData, Box<dyn std::error::Error>> {
    let mut current_idx = 0;
    
    loop {
        match file.read_line() {
            '\0' => return Err(format!("Alignment {} not found", target_idx).into()),
            'A' => {
                if current_idx == target_idx {
                    return parse_alignment(file, metadata);
                }
                current_idx += 1;
            }
            _ => {}
        }
    }
}

fn read_all_alignments(
    path: &str,
    metadata: &FileMetadata,
    trace_spacing: i64,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = OneFile::open_read(path, None, None, 1)?;
    
    loop {
        match file.read_line() {
            '\0' => break,
            'A' => {
                let aln = parse_alignment(&mut file, metadata)?;
                print_alignment(&aln, trace_spacing)?;
            }
            _ => {}
        }
    }
    Ok(())
}

fn parse_alignment(
    file: &mut OneFile,
    metadata: &FileMetadata,
) -> Result<AlignmentData, Box<dyn std::error::Error>> {
    // Read alignment coordinates from current 'A' line
    let query_id = file.int(0);
    let target_id = file.int(3);
    
    let mut aln = AlignmentData {
        query_name: metadata.seq_names.get(&query_id)
            .cloned().unwrap_or_else(|| "unknown".to_string()),
        query_length: metadata.seq_lengths.get(&query_id).copied().unwrap_or(0),
        query_start: file.int(1),
        query_end: file.int(2),
        target_name: metadata.seq_names.get(&target_id)
            .cloned().unwrap_or_else(|| "unknown".to_string()),
        target_length: metadata.seq_lengths.get(&target_id).copied().unwrap_or(0),
        target_start: file.int(4),
        target_end: file.int(5),
        strand: '+',
        ..Default::default()
    };
    
    // Read associated lines
    loop {
        match file.read_line() {
            'R' => aln.strand = '-',
            'D' => aln.differences = file.int(0),
            'T' => aln.tracepoints = file.int_list().map(|v| v.to_vec()).unwrap_or_default(),
            'X' => aln.trace_diffs = file.int_list().map(|v| v.to_vec()).unwrap_or_default(),
            'A' | 'a' | 'g' | '\0' => break,
            _ => {}
        }
    }
    
    Ok(aln)
}

fn print_alignment(aln: &AlignmentData, trace_spacing: i64) -> io::Result<()> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    
    writeln!(handle, "Query: {} (len={}) {}:{}-{}", 
        aln.query_name, aln.query_length, aln.query_name, 
        aln.query_start, aln.query_end)?;
    writeln!(handle, "Target: {} (len={}) {}:{}-{}", 
        aln.target_name, aln.target_length, aln.target_name, 
        aln.target_start, aln.target_end)?;
    writeln!(handle, "Strand: {}", aln.strand)?;
    writeln!(handle, "Differences: {}", aln.differences)?;
    writeln!(handle, "Trace spacing: {}", trace_spacing)?;
    
    print_trace_data(&mut handle, "Tracepoints", &aln.tracepoints)?;
    print_trace_data(&mut handle, "Trace diffs", &aln.trace_diffs)?;
    
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