use onecode::OneFile;
use std::collections::HashMap;
use std::fs;
use std::io::{self, Write};
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "oneview-rs")]
#[command(about = "View alignments from ONE format files", long_about = None)]
struct Args {
    /// Input .1aln file path
    #[arg(value_name = "FILE")]
    input: String,

    /// Alignment number to read (0-indexed). If not specified, prints all alignments.
    #[arg(short, long, value_name = "NUM")]
    alignment: Option<usize>,
}

struct AlignmentIndex {
    /// Line numbers where each 'A' record starts
    line_numbers: Vec<i64>,
    trace_spacing: i64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    
    let file_path = &args.input;
    let index_path = format!("{}.idx", file_path);
    
    // Check if index exists, if not create it
    let index = if std::path::Path::new(&index_path).exists() {
        eprintln!("Loading existing index from {}", index_path);
        load_index(&index_path)?
    } else {
        eprintln!("Building alignment index...");
        let index = build_alignment_index(file_path)?;
        save_index(&index, &index_path)?;
        eprintln!("Saved index with {} alignments to {}", 
            index.line_numbers.len(), index_path);
        index
    };
    
    if let Some(aln_idx) = args.alignment {
        // Read specific alignment by reopening file and seeking
        if aln_idx >= index.line_numbers.len() {
            eprintln!("Error: Alignment {} out of range (max {})", 
                aln_idx, index.line_numbers.len() - 1);
            std::process::exit(1);
        }
        
        // Get sequence metadata from a separate file handle
        let mut metadata_file = OneFile::open_read(file_path, None, None, 1)?;
        let seq_names = metadata_file.get_all_sequence_names();
        let seq_lengths = metadata_file.get_all_sequence_lengths();
        drop(metadata_file); // Close this file handle
        
        // Open a fresh file handle for reading the alignment
        let mut file = OneFile::open_read(file_path, None, None, 1)?;
        let aln = read_alignment_at_line(&mut file, index.line_numbers[aln_idx], 
                                         &seq_names, &seq_lengths)?;
        print_alignment_stdout(&aln, index.trace_spacing)?;
    } else {
        // Read all alignments sequentially in one pass
        let mut file = OneFile::open_read(file_path, None, None, 1)?;
        let seq_names = file.get_all_sequence_names();
        let seq_lengths = file.get_all_sequence_lengths();
        
        // Reopen file for sequential reading to start from beginning
        drop(file);
        let mut file = OneFile::open_read(file_path, None, None, 1)?;
        
        // Read alignments sequentially
        read_all_alignments_sequential(&mut file, &seq_names, &seq_lengths, index.trace_spacing)?;
    }

    Ok(())
}

/// Read all alignments sequentially (efficient for full file scan)
fn read_all_alignments_sequential(
    file: &mut OneFile,
    seq_names: &HashMap<i64, String>,
    seq_lengths: &HashMap<i64, i64>,
    trace_spacing: i64,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut current_alignment: Option<AlignmentData> = None;
    
    loop {
        let line_type = file.read_line();
        
        if line_type == '\0' {
            // EOF - print last alignment if any
            if let Some(aln) = current_alignment.take() {
                print_alignment_stdout(&aln, trace_spacing)?;
            }
            break;
        }
        
        match line_type {
            'A' => {
                // Print previous alignment if we have one
                if let Some(aln) = current_alignment.take() {
                    print_alignment_stdout(&aln, trace_spacing)?;
                }
                
                // Start new alignment
                let query_id = file.int(0);
                let query_start = file.int(1);
                let query_end = file.int(2);
                let target_id = file.int(3);
                let target_start = file.int(4);
                let target_end = file.int(5);
                
                let query_name = seq_names.get(&query_id)
                    .map(|s| s.as_str())
                    .unwrap_or("unknown")
                    .to_string();
                let target_name = seq_names.get(&target_id)
                    .map(|s| s.as_str())
                    .unwrap_or("unknown")
                    .to_string();
                
                let query_length = seq_lengths.get(&query_id).copied().unwrap_or(0);
                let target_length = seq_lengths.get(&target_id).copied().unwrap_or(0);
                
                current_alignment = Some(AlignmentData {
                    query_name,
                    query_length,
                    query_start,
                    query_end,
                    target_name,
                    target_length,
                    target_start,
                    target_end,
                    strand: '+',
                    differences: 0,
                    tracepoints: Vec::new(),
                    trace_diffs: Vec::new(),
                });
            }
            'R' => {
                if let Some(ref mut aln) = current_alignment {
                    aln.strand = '-';
                }
            }
            'D' => {
                if let Some(ref mut aln) = current_alignment {
                    aln.differences = file.int(0);
                }
            }
            'T' => {
                if let Some(ref mut aln) = current_alignment {
                    if let Some(tracepoints) = file.int_list() {
                        aln.tracepoints = tracepoints.to_vec();
                    }
                }
            }
            'X' => {
                if let Some(ref mut aln) = current_alignment {
                    if let Some(trace_diffs) = file.int_list() {
                        aln.trace_diffs = trace_diffs.to_vec();
                    }
                }
            }
            _ => {
                // Skip other line types
            }
        }
    }
    
    Ok(())
}

/// Save index to a simple text file
fn save_index(index: &AlignmentIndex, path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut content = String::new();
    content.push_str(&format!("{}\n", index.trace_spacing));
    content.push_str(&format!("{}\n", index.line_numbers.len()));
    
    for &line_num in &index.line_numbers {
        content.push_str(&format!("{}\n", line_num));
    }
    
    fs::write(path, content)?;
    Ok(())
}

/// Load index from text file
fn load_index(path: &str) -> Result<AlignmentIndex, Box<dyn std::error::Error>> {
    let content = fs::read_to_string(path)?;
    let mut lines = content.lines();
    
    let trace_spacing = lines.next()
        .ok_or("Missing trace_spacing in index")?
        .parse()?;
    
    let num_alignments: usize = lines.next()
        .ok_or("Missing num_alignments in index")?
        .parse()?;
    
    let mut line_numbers = Vec::with_capacity(num_alignments);
    
    for line in lines {
        let line_num: i64 = line.trim().parse()?;
        line_numbers.push(line_num);
    }
    
    if line_numbers.len() != num_alignments {
        return Err(format!("Expected {} alignments, found {}", 
            num_alignments, line_numbers.len()).into());
    }
    
    Ok(AlignmentIndex {
        line_numbers,
        trace_spacing,
    })
}

/// Build an index by scanning the entire file and recording line numbers of 'A' records
fn build_alignment_index(file_path: &str) -> Result<AlignmentIndex, Box<dyn std::error::Error>> {
    let mut file = OneFile::open_read(file_path, None, None, 1)?;
    
    // Read trace spacing from header
    let mut trace_spacing = 100;
    let mut line_numbers = Vec::new();
    
    loop {
        let line_type = file.read_line();
        
        if line_type == '\0' {
            break; // EOF
        }
        
        if line_type == 't' {
            trace_spacing = file.int(0);
        } else if line_type == 'A' {
            // Record this alignment's line number (AFTER reading)
            line_numbers.push(file.line_number());
        }
    }
    
    Ok(AlignmentIndex {
        line_numbers,
        trace_spacing,
    })
}

/// Read alignment at a specific line number
fn read_alignment_at_line(
    file: &mut OneFile,
    target_line: i64,
    seq_names: &HashMap<i64, String>,
    seq_lengths: &HashMap<i64, i64>,
) -> Result<AlignmentData, Box<dyn std::error::Error>> {
    // Scan forward until we've read the target line
    loop {
        let line_type = file.read_line();
        
        if line_type == '\0' {
            return Err("Reached EOF before target line".into());
        }
        
        // Check if we just read the target line
        if file.line_number() == target_line {
            if line_type != 'A' {
                return Err(format!("Expected 'A' at line {}, got '{}'", target_line, line_type).into());
            }
            
            // We've read the 'A' line, now parse it
            let query_id = file.int(0);
            let query_start = file.int(1);
            let query_end = file.int(2);
            let target_id = file.int(3);
            let target_start = file.int(4);
            let target_end = file.int(5);
            
            let query_name = seq_names.get(&query_id)
                .map(|s| s.as_str())
                .unwrap_or("unknown")
                .to_string();
            let target_name = seq_names.get(&target_id)
                .map(|s| s.as_str())
                .unwrap_or("unknown")
                .to_string();
            
            let query_length = seq_lengths.get(&query_id).copied().unwrap_or(0);
            let target_length = seq_lengths.get(&target_id).copied().unwrap_or(0);
            
            let mut aln = AlignmentData {
                query_name,
                query_length,
                query_start,
                query_end,
                target_name,
                target_length,
                target_start,
                target_end,
                strand: '+',
                differences: 0,
                tracepoints: Vec::new(),
                trace_diffs: Vec::new(),
            };
            
            // Read associated lines (R, D, T, X)
            loop {
                let line_type = file.read_line();
                
                if line_type == '\0' || line_type == 'A' || line_type == 'a' || line_type == 'g' {
                    break;
                }
                
                match line_type {
                    'R' => aln.strand = '-',
                    'D' => aln.differences = file.int(0),
                    'T' => {
                        if let Some(tracepoints) = file.int_list() {
                            aln.tracepoints = tracepoints.to_vec();
                        }
                    }
                    'X' => {
                        if let Some(trace_diffs) = file.int_list() {
                            aln.trace_diffs = trace_diffs.to_vec();
                        }
                    }
                    _ => {}
                }
            }
            
            return Ok(aln);
        }
    }
}

#[derive(Debug)]
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

fn print_alignment_stdout(aln: &AlignmentData, trace_spacing: i64) -> io::Result<()> {
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
    writeln!(handle, "Tracepoints: {} values", aln.tracepoints.len())?;
    
    if !aln.tracepoints.is_empty() {
        write!(handle, "  ")?;
        for (i, &tp) in aln.tracepoints.iter().enumerate() {
            if i > 0 {
                write!(handle, " ")?;
            }
            write!(handle, "{}", tp)?;
        }
        writeln!(handle)?;
    }
    
    if !aln.trace_diffs.is_empty() {
        writeln!(handle, "Trace diffs: {} values", aln.trace_diffs.len())?;
        write!(handle, "  ")?;
        for (i, &diff) in aln.trace_diffs.iter().enumerate() {
            if i > 0 {
                write!(handle, " ")?;
            }
            write!(handle, "{}", diff)?;
        }
        writeln!(handle)?;
    }
    
    writeln!(handle)?; // Empty line between alignments
    
    Ok(())
}