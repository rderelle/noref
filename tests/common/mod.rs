use flate2::read::MultiGzDecoder;
use std::fs;
use std::io::Read;
use std::path::Path;

const OUTPUT_FILES: &[&str] = &[
    "noref_raw_alignment.fas",
    "noref_filtered_alignment.fas",
    "noref_raw_missing.tsv",
    "noref_filtered_missing.tsv",
    "noref_raw_NJ.tree",
    "noref_filtered_NJ.tree",
];

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("Failed to read {:?}: {}", path, e))
}

/// Compare all output files with the expected ones.
pub fn compare_test_outputs(expected_dir: &Path, work_dir: &Path) {
    for name in OUTPUT_FILES {
        assert_eq!(
            read(&expected_dir.join(name)),
            read(&work_dir.join(name)),
            "{name} differs"
        );
    }
}

pub fn read_outputs(dir: &Path) -> Vec<String> {
    OUTPUT_FILES
        .iter()
        .map(|name| read(&dir.join(name)))
        .collect()
}

/// Read a FASTA file, decompressing `.gz` files.
pub fn read_fasta_text(path: &Path) -> String {
    let bytes = fs::read(path).unwrap();
    if path.extension().is_some_and(|e| e == "gz") {
        let mut text = String::new();
        MultiGzDecoder::new(bytes.as_slice())
            .read_to_string(&mut text)
            .unwrap();
        text
    } else {
        String::from_utf8(bytes).unwrap()
    }
}
