use onecode::OneFile;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut file = OneFile::open_read("/home/guarracino/Desktop/Garrison/impg/x.n5.1aln", None, None, 1)?;

    // Get all sequence names and lengths once (efficient for multiple lookups)
    let seq_names = file.get_all_sequence_names();
    let seq_lengths = file.get_all_sequence_lengths();
    println!("Found {} sequences", seq_names.len());

    // Read the global trace point spacing
    let mut trace_spacing = 100; // Default
    loop {
        let line_type = file.read_line();
        if line_type == '\0' { break; }
        
        if line_type == 't' {
            trace_spacing = file.int(0);
            println!("Trace point spacing: {}", trace_spacing);
        }
        
        if line_type == 'g' || line_type == 'a' {
            break;
        }
    }
    
    // Now read alignments - restart from beginning
    file = OneFile::open_read("/home/guarracino/Desktop/Garrison/impg/x.n5.1aln", None, None, 1)?;
    
    // Track current alignment data
    let mut current_alignment: Option<AlignmentData> = None;
    
    loop {
        let line_type = file.read_line();
        if line_type == '\0' { break; }

        match line_type {
            'A' => {
                // Print previous alignment if any
                if let Some(aln) = current_alignment.take() {
                    print_alignment(&aln, trace_spacing);
                }
                
                // Read new alignment
                let query_id = file.int(0);
                let query_start = file.int(1);
                let query_end = file.int(2);
                let target_id = file.int(3);
                let target_start = file.int(4);
                let target_end = file.int(5);
                
                let query_name = seq_names.get(&query_id)
                    .map(|s| s.as_str())
                    .unwrap_or("unknown");
                let target_name = seq_names.get(&target_id)
                    .map(|s| s.as_str())
                    .unwrap_or("unknown");
                
                // Get FULL sequence lengths from GDB
                let query_length = seq_lengths.get(&query_id).copied().unwrap_or(0);
                let target_length = seq_lengths.get(&target_id).copied().unwrap_or(0);
                
                current_alignment = Some(AlignmentData {
                    query_name: query_name.to_string(),
                    query_length,
                    query_start,
                    query_end,
                    target_name: target_name.to_string(),
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
            _ => {}
        }
    }
    
    // Print last alignment
    if let Some(aln) = current_alignment {
        print_alignment(&aln, trace_spacing);
    }

    Ok(())
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

fn print_alignment(aln: &AlignmentData, trace_spacing: i64) {
    println!("\nAlignment:");
    println!("  Query:  {} (len={}) {}:{}-{}", 
        aln.query_name, aln.query_length, aln.query_name, 
        aln.query_start, aln.query_end);
    println!("  Target: {} (len={}) {}:{}-{}", 
        aln.target_name, aln.target_length, aln.target_name, 
        aln.target_start, aln.target_end);
    println!("  Strand: {}", aln.strand);
    println!("  Differences: {}", aln.differences);
    println!("  Trace spacing: {}", trace_spacing);
    println!("  Tracepoints: {} values", aln.tracepoints.len());
    
    if !aln.tracepoints.is_empty() {
        print!("    ");
        for (i, &tp) in aln.tracepoints.iter().enumerate() {
            if i < 10 {
                print!("{} ", tp);
            }
        }
        if aln.tracepoints.len() > 10 {
            println!("...");
        } else {
            println!();
        }
    }
    
    if !aln.trace_diffs.is_empty() {
        println!("  Trace diffs: {} values", aln.trace_diffs.len());
        print!("    ");
        for (i, &diff) in aln.trace_diffs.iter().enumerate() {
            if i < 10 {
                print!("{} ", diff);
            }
        }
        if aln.trace_diffs.len() > 10 {
            println!("...");
        } else {
            println!();
        }
    }
}