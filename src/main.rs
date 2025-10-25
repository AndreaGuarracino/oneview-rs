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
    // Query genome (gdb1 - first reference, or embedded if self-alignment)
    query_seq_names: HashMap<i64, String>,
    query_seq_lengths: HashMap<i64, i64>,
    query_contig_offsets: HashMap<i64, (i64, i64)>,

    // Target genome (gdb2 - second reference, or embedded skeleton)
    target_seq_names: HashMap<i64, String>,
    target_seq_lengths: HashMap<i64, i64>,
    target_contig_offsets: HashMap<i64, (i64, i64)>,
}

fn get_file_metadata(path: &str) -> Result<(FileMetadata, i64), Box<dyn std::error::Error>> {
    let mut file = OneFile::open_read(path, None, None, 1)?;

    // Check if there are reference paths to external GDB files
    let references = file.get_references();

    // Initialize empty metadata structures
    let mut query_seq_names = HashMap::new();
    let mut query_seq_lengths = HashMap::new();
    let mut query_contig_offsets = HashMap::new();
    let mut target_seq_names = HashMap::new();
    let mut target_seq_lengths = HashMap::new();
    let mut target_contig_offsets = HashMap::new();

    // The embedded GDB skeleton (if present) is the target genome (gdb2)
    let embedded_names = file.get_all_sequence_names();
    let embedded_lengths = file.get_all_sequence_lengths();
    let embedded_offsets = file.get_all_contig_offsets();

    // Process references:
    // - First reference (count=1) is the QUERY genome (gdb1/A-read)
    // - Second reference (count=2) is the TARGET genome (gdb2/B-read)
    // - If there's an embedded skeleton, it's for the TARGET (gdb2/B-read)

    let mut has_external_query = false;
    let mut has_external_target = false;

    for (ref_idx, (ref_path, ref_count)) in references.iter().enumerate() {
        if ref_path.is_empty() {
            continue;
        }

        // Skip if this is not a genome reference (count > 2 might be other metadata)
        if *ref_count > 2 {
            continue;
        }

        let is_query = *ref_count == 1;    // First reference is query (A-read)
        let is_target = *ref_count == 2;   // Second reference is target (B-read)

        eprintln!("Processing reference {}: {} (count: {}, type: {})",
            ref_idx + 1, ref_path, ref_count,
            if is_query { "query" } else if is_target { "target" } else { "unknown" });

        // Try to load genome metadata
        let query_path = ref_path;

        // Helper function to strip fasta extensions
        let strip_fasta_ext = |p: &str| -> String {
            let path_str = p.to_string();
            // Try to remove common fasta extensions
            for ext in &[".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna"] {
                if path_str.ends_with(ext) {
                    return path_str[..path_str.len() - ext.len()].to_string();
                }
            }
            path_str
        };

        // Get alignment file directory for relative path resolution
        let aln_dir = std::path::Path::new(path)
            .parent()
            .unwrap_or_else(|| std::path::Path::new("."));

        // Try multiple strategies to find the GDB file
        let mut gdb_path: Option<String> = None;

        // Strategy 1: Try as absolute path (as-is)
        if std::path::Path::new(query_path).exists() {
            gdb_path = Some(query_path.clone());
        }

        // Strategy 2: Try adding .1gdb or .gdb extension to the original path
        if gdb_path.is_none() {
            for ext in &[".1gdb", ".gdb"] {
                let with_ext = format!("{}{}", query_path, ext);
                if std::path::Path::new(&with_ext).exists() {
                    gdb_path = Some(with_ext);
                    break;
                }
            }
        }

        // Strategy 3: Strip fasta extension and try with GDB extensions (absolute path)
        if gdb_path.is_none() {
            let base_path = strip_fasta_ext(query_path);
            if base_path != *query_path {
                for ext in &[".1gdb", ".gdb"] {
                    let with_ext = format!("{}{}", base_path, ext);
                    if std::path::Path::new(&with_ext).exists() {
                        gdb_path = Some(with_ext);
                        break;
                    }
                }
            }
        }

        // Strategy 4: Try relative path with fasta extension stripped and GDB extension added
        // (DO THIS BEFORE trying as-is, to avoid finding the fasta file itself)
        if gdb_path.is_none() {
            let base_path = strip_fasta_ext(query_path);
            if base_path != *query_path {
                for ext in &[".1gdb", ".gdb"] {
                    let relative_with_ext = aln_dir.join(format!("{}{}", base_path, ext));
                    if relative_with_ext.exists() {
                        gdb_path = Some(relative_with_ext.to_string_lossy().to_string());
                        break;
                    }
                }
            }
        }

        // Strategy 5: Try relative path with .1gdb or .gdb extension
        if gdb_path.is_none() {
            for ext in &[".1gdb", ".gdb"] {
                let relative_with_ext = aln_dir.join(format!("{}{}", query_path, ext));
                if relative_with_ext.exists() {
                    gdb_path = Some(relative_with_ext.to_string_lossy().to_string());
                    break;
                }
            }
        }

        // Strategy 6: Try relative to alignment file directory (as-is)
        // (Only as last resort, to avoid finding non-GDB files)
        if gdb_path.is_none() {
            let relative_path = aln_dir.join(query_path);
            if relative_path.exists() {
                gdb_path = Some(relative_path.to_string_lossy().to_string());
            }
        }

        let gdb_path = if let Some(found_path) = gdb_path {
            found_path
        } else {
            eprintln!("Warning: Could not find query GDB file for reference: {}", query_path);
            eprintln!("Tried:");
            eprintln!("  - {}", query_path);
            eprintln!("  - {}.1gdb / {}.gdb", query_path, query_path);
            let base = strip_fasta_ext(query_path);
            if base != *query_path {
                eprintln!("  - {}.1gdb / {}.gdb", base, base);
            }
            eprintln!("  - Relative paths in {}", aln_dir.display());
            eprintln!("Query contig-to-scaffold mappings will not be available");
            query_path.clone()
        };

        // Try to load the GDB metadata
        if let Ok((ref_names, ref_lengths, ref_offsets)) = OneFile::read_gdb_metadata(&gdb_path) {
            if is_query {
                query_seq_names = ref_names;
                query_seq_lengths = ref_lengths;
                query_contig_offsets = ref_offsets;
                has_external_query = true;
                eprintln!("Loaded query genome metadata from: {} ({} sequences)", gdb_path, query_seq_names.len());
            } else if is_target {
                target_seq_names = ref_names;
                target_seq_lengths = ref_lengths;
                target_contig_offsets = ref_offsets;
                has_external_target = true;
                eprintln!("Loaded target genome metadata from: {} ({} sequences)", gdb_path, target_seq_names.len());
            }
        } else {
            eprintln!("Warning: Failed to load GDB metadata from: {}", gdb_path);
        }
    }

    // If we didn't load external target, use embedded skeleton
    if !has_external_target && !embedded_names.is_empty() {
        target_seq_names = embedded_names;
        target_seq_lengths = embedded_lengths;
        target_contig_offsets = embedded_offsets;
        eprintln!("Using embedded skeleton for target genome ({} sequences)", target_seq_names.len());
    }

    // If this is a self-alignment (no external query), use target for query too
    if !has_external_query && !target_seq_names.is_empty() {
        query_seq_names = target_seq_names.clone();
        query_seq_lengths = target_seq_lengths.clone();
        query_contig_offsets = target_contig_offsets.clone();
        eprintln!("Self-alignment detected: using target genome for query");
    }

    if query_seq_names.is_empty() && target_seq_names.is_empty() {
        eprintln!("Warning: No sequence metadata found in file or external references");
    }

    let metadata = FileMetadata {
        query_seq_names,
        query_seq_lengths,
        query_contig_offsets,
        target_seq_names,
        target_seq_lengths,
        target_contig_offsets,
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
    writeln!(handle, "Query sequences: {}", metadata.query_seq_names.len())?;
    writeln!(handle, "Target sequences: {}", metadata.target_seq_names.len())?;

    // Sort query sequences by ID for consistent output
    if !metadata.query_seq_names.is_empty() {
        let mut query_sequences: Vec<_> = metadata.query_seq_names.iter().collect();
        query_sequences.sort_by_key(|&(&id, _)| id);

        writeln!(handle, "\nQuery Sequences:")?;
        for &(&id, ref name) in query_sequences.iter().take(10) {
            let length = metadata.query_seq_lengths.get(&id).copied().unwrap_or(0);
            writeln!(handle, "  {}) {} (length: {})", id, name, length)?;
        }
        if query_sequences.len() > 10 {
            writeln!(handle, "  ... and {} more", query_sequences.len() - 10)?;
        }
    }

    // Sort target sequences by ID for consistent output
    if !metadata.target_seq_names.is_empty() {
        let mut target_sequences: Vec<_> = metadata.target_seq_names.iter().collect();
        target_sequences.sort_by_key(|&(&id, _)| id);

        writeln!(handle, "\nTarget Sequences:")?;
        for &(&id, ref name) in target_sequences.iter().take(10) {
            let length = metadata.target_seq_lengths.get(&id).copied().unwrap_or(0);
            writeln!(handle, "  {}) {} (length: {})", id, name, length)?;
        }
        if target_sequences.len() > 10 {
            writeln!(handle, "  ... and {} more", target_sequences.len() - 10)?;
        }
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
        .query_seq_names
        .get(&query_id)
        .cloned()
        .ok_or_else(|| format!("Query sequence with ID {} not found in metadata", query_id))?;
    let target_name = metadata
        .target_seq_names
        .get(&target_id)
        .cloned()
        .ok_or_else(|| format!("Target sequence with ID {} not found in metadata", target_id))?;
    let query_length = metadata
        .query_seq_lengths
        .get(&query_id)
        .copied()
        .ok_or_else(|| {
            format!(
                "Query sequence length for ID {} not found in metadata",
                query_id
            )
        })?;
    let target_length = metadata
        .target_seq_lengths
        .get(&target_id)
        .copied()
        .ok_or_else(|| {
            format!(
                "Target sequence length for ID {} not found in metadata",
                target_id
            )
        })?;

    let (query_offset, _) = metadata
        .query_contig_offsets
        .get(&query_id)
        .copied()
        .ok_or_else(|| {
            format!(
                "Contig offset for query sequence ID {} not found in metadata",
                query_id
            )
        })?;

    let (target_offset, target_contig_len) = metadata
        .target_contig_offsets
        .get(&target_id)
        .copied()
        .ok_or_else(|| {
            format!(
                "Contig offset for target sequence ID {} not found in metadata",
                target_id
            )
        })?;

    let query_contig_start = file.int(1);
    let query_contig_end = file.int(2);
    let mut target_contig_start = file.int(4);
    let mut target_contig_end = file.int(5);

    let mut aln = AlignmentData {
        query_name,
        query_length,
        query_start: 0,
        query_end: 0,
        target_name,
        target_length,
        target_start: 0,
        target_end: 0,
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
        let orig_start = target_contig_start;
        let orig_end = target_contig_end;
        // Reverse-complement target coordinates so start/end reflect forward strand
        target_contig_start = target_contig_len - orig_end;
        target_contig_end = target_contig_len - orig_start;
    }

    aln.query_start = add_offset(query_offset, query_contig_start)?;
    aln.query_end = add_offset(query_offset, query_contig_end)?;
    aln.target_start = add_offset(target_offset, target_contig_start)?;
    aln.target_end = add_offset(target_offset, target_contig_end)?;

    Ok((aln, next_line))
}

fn add_offset(offset: i64, position: i64) -> Result<i64, Box<dyn std::error::Error>> {
    offset
        .checked_add(position)
        .ok_or_else(|| "Coordinate overflow when applying contig offset".into())
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

    // Match ALNtoPAF calculation (when not computing CIGAR):
    let block_length = query_span + target_span;
    let matches = ((block_length - aln.differences) / 2).max(0);
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
