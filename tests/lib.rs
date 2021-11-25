use itertools::izip;
use itertools::Itertools;
use std::fs;
use std::process::Command;

#[test]
fn test_report() {
    assert!(Command::new("bash")
        .arg("-c")
        .arg("target/debug/pf -q tests/resources/example.fastq")
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
}
