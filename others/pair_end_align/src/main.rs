extern crate bio;
extern crate clap;

use std::fs::File;
use std::str;
use std::cmp;

use clap::{Arg, App};
use bio::alignment::pairwise::*;
use bio::io::fastq;

fn ann_line(start: usize, end: usize, c: u8) -> Vec<u8> {
    let mut vec = vec![b' '; start];
    vec.push(c);
    let vec2 = vec![b' '; end - start];
    vec.extend(vec2);
    vec.push(c);
    vec
}


fn main() {
    let matches = App::new("PET Extract")
        .arg(Arg::with_name("fq1")
             .takes_value(true)
             .help("Fastq file of reads 1."))
        .arg(Arg::with_name("fq2")
             .takes_value(true)
             .help("Fastq file of reads 2."))
        .get_matches();

    let fq1_path = matches.value_of("fq1").expect("Please input fq1.");
    let fq2_path = matches.value_of("fq2").expect("Please input fq2.");
    let fq1_file = File::open(fq1_path).unwrap();
    let fq1 = fastq::Reader::new(fq1_file);
    let fq2_file = File::open(fq2_path).unwrap();
    let fq2 = fastq::Reader::new(fq2_file);

    let mut records1 = fq1.records();
    let mut records2 = fq2.records();

    let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};

    loop {
        let r1 = match records1.next() {
            Some(r) => match r {
                Ok(r_) => r_,
                Err(e) => panic!("{:?}", e),
            },
            None => break
        };
        let r2 = match records1.next() {
            Some(r) => match r {
                Ok(r_) => r_,
                Err(e) => panic!("{:?}", e),
            },
            None => break
        };

        let vec1 = r1.seq();
        let vec2 = r2.seq();
        let seq1 = &vec1[..(vec1.len()-1)];
        let seq2 = &vec2[..(vec2.len()-1)];
        let seq2_revc = bio::alphabets::dna::revcomp(&vec2[..(vec2.len()-1)]);

        let mut aligner = Aligner::with_capacity(seq1.len(), seq2_revc.len(), -5, -1, &score);
        let alignment = aligner.local(&seq1, &seq2_revc);
        let common = &seq1[alignment.xstart..alignment.xend];
        let common_seq = str::from_utf8(common).unwrap();

        let m = (alignment.ystart as i32) - (alignment.xstart as i32);
        let prefix = vec![b' '; m.abs() as usize];
        let rep = str::from_utf8(&prefix).unwrap();

        println!("");
        if m > 0 {
            let mut pre_ann = ann_line(alignment.ystart, alignment.yend, b'v');
            let line_ann = str::from_utf8(&pre_ann).unwrap();
            println!("{}", line_ann);
            println!("{}{}", rep, str::from_utf8(seq1).unwrap());
            println!("{}", str::from_utf8(&seq2_revc).unwrap());
        } else {
            println!("{}", str::from_utf8(seq1).unwrap());
            println!("{}{}", rep, str::from_utf8(&seq2_revc).unwrap());
            let mut pre_ann = ann_line(alignment.xstart, alignment.xend, b'^');
            let line_ann = str::from_utf8(&pre_ann).unwrap();
            println!("{}", line_ann);
        }

    }

    
    let x =     b"ACCCTGGATGGGGG";
    let y = b"AAAAACCGTTGAT";
    let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
    let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
    let alignment = aligner.local(x, y);
    println!("{:?}", alignment);
}
