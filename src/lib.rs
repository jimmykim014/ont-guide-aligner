//! Rust counting core for ont_guide_aligner (`ont_guide_aligner._core`).
//!
//! `count_file` is the read loop: it finds guides next to the learned flanks, assigns each read
//! to at most one guide, and optionally rescues reads with a 1 bp indel in the guide. Python
//! keeps everything else (argument parsing, library loading, flank learning, output).
//!
//! The 1-mismatch lookup table follows the approach of guide-counter's `build_lookup`
//! (Fulcrum Genomics, MIT license, https://github.com/fulcrumgenomics/guide-counter), as does
//! decompressing on a background thread; no code was copied from it.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{sync_channel, Receiver};
use std::thread;
use std::time::Duration;

use ahash::{AHashMap, AHashSet};
use flate2::read::MultiGzDecoder;
use memchr::memmem::Finder;
use pyo3::prelude::*;

type Lookup = AHashMap<Vec<u8>, (u32, bool)>;

/// Lookup value for a variant shared by more than one guide: a read with it fits several
/// guides, so if nothing else matches it is counted as ambiguous rather than unmatched.
const AMBIGUOUS: u32 = u32::MAX;

/// sequence -> (guide index, is_exact). One-mismatch variants shared by >1 guide map to
/// AMBIGUOUS; exact sequences are inserted last so they always win.
fn build_lookup(seqs: &[Vec<u8>], allow_mismatch: bool) -> Lookup {
    let mut lookup: Lookup = AHashMap::default();
    if allow_mismatch {
        let mut clash: AHashSet<Vec<u8>> = AHashSet::default();
        for (i, s) in seqs.iter().enumerate() {
            let i = i as u32;
            for p in 0..s.len() {
                for &b in b"ACGT" {
                    if b != s[p] {
                        let mut v = s.clone();
                        v[p] = b;
                        if let Some(&(prev, _)) = lookup.get(&v) {
                            if prev != i {
                                clash.insert(v.clone());
                            }
                        }
                        lookup.insert(v, (i, false));
                    }
                }
            }
        }
        for v in clash {
            lookup.insert(v, (AMBIGUOUS, false));
        }
    }
    for (i, s) in seqs.iter().enumerate() {
        lookup.insert(s.clone(), (i as u32, true));
    }
    lookup
}

/// (guide with one base deleted) -> guide index, for spans 1 bp shorter than a guide. Variants
/// shared by >1 guide map to AMBIGUOUS.
fn build_deletion_lookup(seqs: &[Vec<u8>]) -> AHashMap<Vec<u8>, u32> {
    let mut lookup: AHashMap<Vec<u8>, u32> = AHashMap::default();
    let mut clash: AHashSet<Vec<u8>> = AHashSet::default();
    let mut own: AHashSet<Vec<u8>> = AHashSet::default();
    for (i, s) in seqs.iter().enumerate() {
        let i = i as u32;
        own.clear();
        for p in 0..s.len() {
            let mut v = Vec::with_capacity(s.len() - 1);
            v.extend_from_slice(&s[..p]);
            v.extend_from_slice(&s[p + 1..]);
            own.insert(v); // a homopolymer gives the same variant more than once
        }
        for v in own.drain() {
            if let Some(&prev) = lookup.get(&v) {
                if prev != i {
                    clash.insert(v.clone());
                }
            }
            lookup.insert(v, i);
        }
    }
    for v in clash {
        lookup.insert(v, AMBIGUOUS);
    }
    lookup
}

/// Error for input that isn't a well-formed 4-line FASTQ. Reads are taken by position (line 2
/// of each record), so a damaged record would otherwise shift every read after it.
fn not_fastq(why: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "not a valid FASTQ file: {} (it may be truncated or corrupted)",
            why
        ),
    )
}

/// Prefix an I/O error with the file it came from.
fn with_path(e: io::Error, path: &str) -> io::Error {
    io::Error::new(e.kind(), format!("{}: {}", path, e))
}

/// Reader fed by a background decompression thread, so gzip inflation and matching run on
/// different cores.
struct ChannelReader {
    rx: Receiver<io::Result<Vec<u8>>>,
    buf: Vec<u8>,
    pos: usize,
}

impl Read for ChannelReader {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        while self.pos >= self.buf.len() {
            match self.rx.recv() {
                Ok(Ok(chunk)) => {
                    self.buf = chunk;
                    self.pos = 0;
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => return Ok(0), // producer finished
            }
        }
        let n = out.len().min(self.buf.len() - self.pos);
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

/// Open a FASTQ or FASTQ.gz; returns the reader and whether the file is gzip.
fn open_reader(path: &str) -> io::Result<(Box<dyn BufRead + Send>, bool)> {
    // Opened once and peeked, so a pipe works too. Gzip is recognised by its magic bytes, so
    // names like .FASTQ.GZ or .bgz (or a misnamed .fastq) work.
    let mut file = BufReader::with_capacity(1 << 20, File::open(path)?);
    if !file.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        return Ok((Box::new(file), false));
    }
    let (tx, rx) = sync_channel::<io::Result<Vec<u8>>>(8);
    thread::spawn(move || {
        let mut dec = MultiGzDecoder::new(file);
        loop {
            let mut chunk = vec![0u8; 1 << 20];
            let mut filled = 0;
            // fill the chunk (or hit EOF/error) before sending, to keep channel traffic low
            let res = loop {
                match dec.read(&mut chunk[filled..]) {
                    Ok(0) => break Ok(true),
                    Ok(n) => {
                        filled += n;
                        if filled == chunk.len() {
                            break Ok(false);
                        }
                    }
                    Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                    Err(e) => break Err(e),
                }
            };
            match res {
                Ok(eof) => {
                    chunk.truncate(filled);
                    if !chunk.is_empty() && tx.send(Ok(chunk)).is_err() {
                        return; // consumer gone
                    }
                    if eof {
                        return;
                    }
                }
                Err(e) => {
                    // the caller adds the path
                    let msg = format!("failed to decompress: {}", e);
                    let _ = tx.send(Err(io::Error::new(e.kind(), msg)));
                    return;
                }
            }
        }
    });
    let reader = BufReader::with_capacity(
        1 << 20,
        ChannelReader {
            rx,
            buf: Vec::new(),
            pos: 0,
        },
    );
    Ok((Box::new(reader), true))
}

/// Characters Python's str.strip() removes, restricted to ASCII (FASTQ is ASCII).
#[inline]
fn is_py_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c | 0x1c..=0x1f)
}

#[inline]
fn strip(s: &[u8]) -> &[u8] {
    let start = s.iter().position(|&b| !is_py_space(b)).unwrap_or(s.len());
    let end = s
        .iter()
        .rposition(|&b| !is_py_space(b))
        .map_or(start, |e| e + 1);
    &s[start..end]
}

#[inline]
fn rc_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other, // N and anything else unchanged, as Python's str.translate
    }
}

/// All start positions of `needle` in `hay`, including overlapping ones (like repeated
/// str.find(sub, i + 1) in Python; memmem's own iterator skips overlaps).
#[inline]
fn find_all(finder: &Finder, hay: &[u8], out: &mut Vec<usize>) {
    let mut from = 0;
    while from <= hay.len() {
        match finder.find(&hay[from..]) {
            Some(i) => {
                out.push(from + i);
                from += i + 1;
            }
            None => break,
        }
    }
}

/// Per-read matching state shared by all chunks of one sample.
struct Matcher<'a> {
    lookup: &'a Lookup,
    n_guides: usize,
    glen: usize,
    up: &'a [u8],
    down: &'a [u8],
    min_len: usize,
    max_len: usize,
    /// Indel rescue: Some(deletion lookup) when enabled.
    del_lookup: Option<&'a AHashMap<Vec<u8>, u32>>,
}

#[derive(Default, Clone, Copy)]
struct Tallies {
    total: u64,
    filtered: u64,
    assigned: u64,
    both_flank: u64,
    mismatch: u64,
    ambiguous: u64,
    unmatched: u64,
    indel: u64,
}

impl Tallies {
    fn add(&mut self, o: &Tallies) {
        self.total += o.total;
        self.filtered += o.filtered;
        self.assigned += o.assigned;
        self.both_flank += o.both_flank;
        self.mismatch += o.mismatch;
        self.ambiguous += o.ambiguous;
        self.unmatched += o.unmatched;
        self.indel += o.indel;
    }
}

impl Matcher<'_> {
    /// Guides exactly one insertion or deletion away from `span` (the bases between the two
    /// flanks), added to `cand` without duplicates.
    fn indel_guides(
        &self,
        span: &[u8],
        del_lookup: &AHashMap<Vec<u8>, u32>,
        buf: &mut Vec<u8>,
        cand: &mut Vec<u32>,
    ) {
        let glen = self.glen;
        let add = |j: u32, cand: &mut Vec<u32>| {
            if !cand.contains(&j) {
                cand.push(j);
            }
        };
        if span.len() + 1 == glen {
            if let Some(&j) = del_lookup.get(span) {
                add(j, cand);
            }
        } else if span.len() == glen + 1 {
            for k in 0..=glen {
                buf.clear();
                buf.extend_from_slice(&span[..k]);
                buf.extend_from_slice(&span[k + 1..]);
                // exact library sequences are the lookup entries with is_exact = true
                if let Some(&(j, true)) = self.lookup.get(&buf[..]) {
                    add(j, cand);
                }
            }
        }
    }

    /// The read loop. `reader` must start at the first line of a FASTQ record. `stop` is set
    /// when the user presses Ctrl-C; the loop then returns early.
    fn run(&self, reader: &mut dyn BufRead, stop: &AtomicBool) -> io::Result<(Vec<u64>, Tallies)> {
        let (glen, up, down) = (self.glen, self.up, self.down);
        let flen = up.len();
        let up_f = Finder::new(up);
        let down_f = Finder::new(down);
        let mut counts = vec![0u64; self.n_guides];
        let mut t = Tallies::default();

        let mut line: Vec<u8> = Vec::with_capacity(1024);
        let mut read: Vec<u8> = Vec::with_capacity(1024);
        let mut rc: Vec<u8> = Vec::with_capacity(1024);
        let mut starts: Vec<usize> = Vec::with_capacity(16);
        // per orientation: guide starts after `up`, and positions of `down` (for indel rescue)
        let mut ups: [Vec<usize>; 2] = [Vec::with_capacity(8), Vec::with_capacity(8)];
        let mut downs: [Vec<usize>; 2] = [Vec::with_capacity(8), Vec::with_capacity(8)];
        let mut hits: Vec<usize> = Vec::with_capacity(16);
        let mut cand: Vec<u32> = Vec::with_capacity(4);
        let mut buf: Vec<u8> = Vec::with_capacity(64);
        // guide index -> best (number of intact flanks, is_exact); tiny, so a Vec beats a map
        let mut support: Vec<(u32, (u8, bool))> = Vec::with_capacity(8);
        let mut line_no: u64 = 0;
        let mut seq_len = 0;

        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            let pos = line_no % 4;
            if pos == 0 {
                if strip(&line).is_empty() {
                    continue; // blank line between records
                }
                if line[0] != b'@' {
                    return Err(not_fastq("a read's first line doesn't start with '@'"));
                }
            } else if pos == 2 && line.first() != Some(&b'+') {
                return Err(not_fastq(
                    "the line after a read's sequence doesn't start with '+'",
                ));
            }
            line_no += 1;
            if pos == 3 && strip(&line).len() != seq_len {
                return Err(not_fastq(
                    "a read's quality line is not the same length as its sequence",
                ));
            }
            if pos != 1 {
                continue;
            }
            seq_len = strip(&line).len();
            read.clear();
            read.extend(strip(&line).iter().map(|b| b.to_ascii_uppercase()));
            t.total += 1;
            if t.total.is_multiple_of(4096) && stop.load(Ordering::Relaxed) {
                return Err(io::Error::new(io::ErrorKind::Interrupted, "interrupted"));
            }
            let n = read.len();
            if n < self.min_len || n > self.max_len {
                t.filtered += 1;
                continue;
            }
            support.clear();
            // a window 1 mismatch from two guides; only decides ambiguous vs unmatched
            let mut near_two = false;
            rc.clear();
            rc.extend(read.iter().rev().map(|&b| rc_base(b)));
            for (o, seq) in [&read[..], &rc[..]].into_iter().enumerate() {
                starts.clear();
                hits.clear();
                find_all(&up_f, seq, &mut hits);
                ups[o].clear();
                ups[o].extend(hits.iter().map(|&i| i + flen));
                starts.extend_from_slice(&ups[o]);
                downs[o].clear();
                find_all(&down_f, seq, &mut downs[o]);
                // skip `down` hits too close to the read start: the guide would start before it
                starts.extend(downs[o].iter().filter(|&&i| i >= glen).map(|&i| i - glen));
                for &s in &starts {
                    if s + glen > seq.len() {
                        continue; // window runs past the end of the read
                    }
                    let Some(&(j, is_exact)) = self.lookup.get(&seq[s..s + glen]) else {
                        continue;
                    };
                    if j == AMBIGUOUS {
                        near_two = true;
                        continue;
                    }
                    let up_ok = s >= flen && &seq[s - flen..s] == up;
                    let down_ok =
                        s + glen + flen <= seq.len() && &seq[s + glen..s + glen + flen] == down;
                    let score = (up_ok as u8 + down_ok as u8, is_exact);
                    match support.iter_mut().find(|(g, _)| *g == j) {
                        Some((_, best)) => {
                            if score > *best {
                                *best = score;
                            }
                        }
                        None => {
                            if score > (0, false) {
                                support.push((j, score));
                            }
                        }
                    }
                }
            }

            if support.is_empty() {
                if let Some(del_lookup) = self.del_lookup {
                    cand.clear();
                    for (o, seq) in [&read[..], &rc[..]].into_iter().enumerate() {
                        for &s in &ups[o] {
                            for &d in &downs[o] {
                                // both flanks, glen +/- 1 bp apart (d > s always holds then)
                                if d == s + glen + 1 || d + 1 == s + glen {
                                    self.indel_guides(&seq[s..d], del_lookup, &mut buf, &mut cand);
                                }
                            }
                        }
                    }
                    if cand.len() == 1 && cand[0] != AMBIGUOUS {
                        counts[cand[0] as usize] += 1;
                        t.assigned += 1;
                        t.both_flank += 1;
                        t.indel += 1;
                        continue;
                    }
                    if !cand.is_empty() {
                        t.ambiguous += 1;
                        continue;
                    }
                }
                if near_two {
                    t.ambiguous += 1;
                } else {
                    t.unmatched += 1;
                }
                continue;
            }
            let top = support.iter().map(|&(_, sc)| sc).max().unwrap();
            let mut winners = support.iter().filter(|&&(_, sc)| sc == top);
            let first = winners.next().unwrap().0;
            if winners.next().is_none() {
                counts[first as usize] += 1;
                t.assigned += 1;
                t.both_flank += (top.0 == 2) as u64;
                t.mismatch += (!top.1) as u64;
            } else {
                t.ambiguous += 1;
            }
        }
        if !line_no.is_multiple_of(4) {
            return Err(not_fastq("the file ends in the middle of a read"));
        }
        Ok((counts, t))
    }
}

/// Does a complete gzip member start at `buf[0]`, and does it decode to the start of a FASTQ
/// record ('@' header line, '+' third line)? Decoding the whole member also checks its CRC,
/// so a chance 1f 8b 08 inside compressed data is rejected.
fn member_starts_record(buf: &[u8]) -> bool {
    let mut dec = flate2::read::GzDecoder::new(buf);
    let mut out = Vec::new();
    if dec.read_to_end(&mut out).is_err() {
        return false;
    }
    let mut lines = out.split(|&b| b == b'\n');
    let l1 = lines.next().unwrap_or(&[]);
    let _l2 = lines.next();
    let l3 = lines.next().unwrap_or(&[]);
    l1.first() == Some(&b'@') && l3.first() == Some(&b'+')
}

/// Byte offsets where the file can be split for parallel decoding: the first gzip member at or
/// after each of `parts - 1` evenly spaced points that starts a FASTQ record. Files written
/// as one gzip stream (most non-ONT data) have no such points, and are decoded in one piece.
fn split_points(path: &str, parts: usize) -> io::Result<Vec<u64>> {
    const WINDOW: usize = 4 << 20;
    let mut f = File::open(path)?;
    let size = f.metadata()?.len();
    let mut points = vec![0u64];
    // Smaller files are decoded in one piece. The env var is for tests only, so CI can exercise
    // the parallel path without a 64 MB file.
    let min_size: u64 = std::env::var("ONT_GUIDE_ALIGNER_MIN_SPLIT_BYTES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(64 << 20);
    if parts < 2 || size < min_size {
        return Ok(points);
    }
    let mut buf = vec![0u8; WINDOW];
    for k in 1..parts {
        let target = size * k as u64 / parts as u64;
        if target <= *points.last().unwrap() {
            continue;
        }
        f.seek(SeekFrom::Start(target))?;
        let mut filled = 0;
        while filled < WINDOW {
            match f.read(&mut buf[filled..])? {
                0 => break,
                n => filled += n,
            }
        }
        let win = &buf[..filled];
        let mut from = 0;
        while let Some(i) = memchr::memmem::find(&win[from..], &[0x1f, 0x8b, 0x08]) {
            let at = from + i;
            if member_starts_record(&win[at..]) {
                points.push(target + at as u64);
                break;
            }
            from = at + 1;
        }
    }
    Ok(points)
}

/// The byte offsets a .gz would be split at for `parts` threads; for tests.
#[pyfunction]
fn gz_split_points(path: &str, parts: usize) -> PyResult<Vec<u64>> {
    Ok(split_points(path, parts)?)
}

/// Count one FASTQ. Returns (counts per guide, tallies: total, length_filtered, assigned,
/// both_flank, mismatch, ambiguous, unmatched, indel). `threads` > 1 decodes and counts a
/// multi-member .gz in parallel pieces. `allow_indel` turns on the 1 bp indel rescue (ignored
/// with `exact_only`).
#[pyfunction]
#[pyo3(signature = (path, seqs, glen, up, down, min_len, max_len, exact_only, threads=1, allow_indel=false))]
#[allow(clippy::too_many_arguments)]
fn count_file(
    py: Python<'_>,
    path: &str,
    seqs: Vec<String>,
    glen: usize,
    up: &str,
    down: &str,
    min_len: usize,
    max_len: usize,
    exact_only: bool,
    threads: usize,
    allow_indel: bool,
) -> PyResult<(Vec<u64>, HashMap<String, u64>)> {
    let seqs: Vec<Vec<u8>> = seqs.into_iter().map(String::into_bytes).collect();
    let lookup = build_lookup(&seqs, !exact_only);
    let del_lookup = if allow_indel && !exact_only {
        Some(build_deletion_lookup(&seqs))
    } else {
        None
    };
    let m = Matcher {
        lookup: &lookup,
        n_guides: seqs.len(),
        glen,
        up: up.as_bytes(),
        down: down.as_bytes(),
        min_len,
        max_len,
        del_lookup: del_lookup.as_ref(),
    };

    let n_guides = seqs.len();

    // All file I/O and counting run with the GIL released. Counting runs in worker threads
    // while this thread checks for Ctrl-C: Python only handles signals on the main thread.
    let stop = AtomicBool::new(false);
    let (counts, t) = py.detach(|| -> PyResult<(Vec<u64>, Tallies)> {
        let (reader, gz) = open_reader(path).map_err(|e| with_path(e, path))?;
        // a large multi-member .gz is split into pieces decoded and counted in parallel
        let split = if gz {
            split_points(path, threads).map_err(|e| with_path(e, path))?
        } else {
            vec![0]
        };
        thread::scope(|sc| {
            let (m, stop) = (&m, &stop);
            let waiter = thread::current(); // woken by each worker when it finishes
            let handles: Vec<_> = if split.len() < 2 {
                let (waiter, mut reader) = (waiter.clone(), reader);
                vec![sc.spawn(move || -> io::Result<(Vec<u64>, Tallies)> {
                    let r = m.run(&mut *reader, stop);
                    waiter.unpark();
                    r
                })]
            } else {
                drop(reader);
                let size = std::fs::metadata(path)
                    .map_err(|e| with_path(e, path))?
                    .len();
                let mut bounds = split.clone();
                bounds.push(size);
                bounds
                    .windows(2)
                    .map(|w| {
                        let (lo, hi) = (w[0], w[1]);
                        let waiter = waiter.clone();
                        sc.spawn(move || -> io::Result<(Vec<u64>, Tallies)> {
                            let count_piece = || -> io::Result<(Vec<u64>, Tallies)> {
                                let mut f = File::open(path)?;
                                f.seek(SeekFrom::Start(lo))?;
                                let piece = BufReader::with_capacity(1 << 20, f.take(hi - lo));
                                let mut reader =
                                    BufReader::with_capacity(1 << 20, MultiGzDecoder::new(piece));
                                m.run(&mut reader, stop)
                            };
                            let r = count_piece();
                            waiter.unpark();
                            r
                        })
                    })
                    .collect()
            };
            while !handles.iter().all(|h| h.is_finished()) {
                thread::park_timeout(Duration::from_millis(50));
                if let Err(e) = Python::attach(|py| py.check_signals()) {
                    stop.store(true, Ordering::Relaxed);
                    for h in handles {
                        let _ = h.join();
                    }
                    return Err(e); // usually KeyboardInterrupt
                }
            }
            let mut counts = vec![0u64; n_guides];
            let mut t = Tallies::default();
            for h in handles {
                let (c, tt) = h
                    .join()
                    .expect("counting thread panicked")
                    .map_err(|e| with_path(e, path))?;
                for (a, b) in counts.iter_mut().zip(c) {
                    *a += b;
                }
                t.add(&tt);
            }
            Ok((counts, t))
        })
    })?;

    let tallies: HashMap<String, u64> = [
        ("total", t.total),
        ("length_filtered", t.filtered),
        ("assigned", t.assigned),
        ("both_flank", t.both_flank),
        ("mismatch", t.mismatch),
        ("ambiguous", t.ambiguous),
        ("unmatched", t.unmatched),
        ("indel", t.indel),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_string(), v))
    .collect();
    Ok((counts, tallies))
}

/// Number of guides that have at least one other guide exactly 1 substitution away (Hamming
/// distance 1). Guides are compared byte for byte as given (Python passes upper-case A/C/G/T of
/// equal length).
#[pyfunction]
fn close_guide_count(py: Python<'_>, seqs: Vec<String>) -> usize {
    py.detach(|| {
        let index: HashMap<&[u8], usize> = seqs
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_bytes(), i))
            .collect();
        let mut close = vec![false; seqs.len()];
        let mut buf = Vec::new();
        for (i, s) in seqs.iter().enumerate() {
            if close[i] {
                continue;
            }
            buf.clear();
            buf.extend_from_slice(s.as_bytes());
            'search: for p in 0..buf.len() {
                let orig = buf[p];
                for &b in b"ACGT" {
                    if b == orig {
                        continue;
                    }
                    buf[p] = b;
                    if let Some(&j) = index.get(buf.as_slice()) {
                        close[i] = true;
                        close[j] = true;
                        buf[p] = orig;
                        break 'search;
                    }
                }
                buf[p] = orig;
            }
        }
        close.iter().filter(|&&c| c).count()
    })
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(count_file, m)?)?;
    m.add_function(wrap_pyfunction!(close_guide_count, m)?)?;
    m.add_function(wrap_pyfunction!(gz_split_points, m)?)?;
    Ok(())
}
