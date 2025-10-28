use crate::{DIM, Matrix};
use crate::flag::Flag;

use std::io::{Read, Write};

pub fn load_flag_file(mut file: std::fs::File) -> Vec<Flag> {
    // parse numbers
    let mut string: String = String::new();
    file.read_to_string(&mut string).unwrap();
    let mut nums = string
        .split(|c: char| c.is_whitespace())     // isolate each number
        .flat_map(|n| n.parse::<f32>());  // parse into numbers

    // convert this data into the vertices of the simplices
    let mut fundamental_simplices: Vec<Flag> = Vec::new();
    for c in std::iter::from_fn(move || {
        let chunk: Vec<f32> = nums.by_ref().take(DIM*DIM).collect(); // take DIM^2 at a time
        if chunk.is_empty() {None} else {Some(chunk)}
    }) {
        fundamental_simplices.push(Flag { vertices: Matrix::from_vec(c) })
    }
    fundamental_simplices
}

pub fn save_flags_to_file(filename: String, flags: &Vec<Flag>) -> Result<(), std::io::Error> {
    let mut f = std::fs::File::create_new(filename)?;
    for flag in flags.iter() {
        let vertices = flag.vertices;
        for v in vertices.iter() {
            write!(f, "{} ", *v).unwrap();
        }
        writeln!(f).unwrap();
    }
    Ok(())
}

pub fn coxmat_to_name(m: [[u8; DIM]; DIM]) -> String {
    let mut name = String::new();
    for r in 0..DIM {
        for c in 0..r+1 {
            let n = m[r][c];
            let c = char::from_digit(n as u32, 36).expect("invalid CD (too big)");
            name.push(c);
        }
    }
    name.push_str(".flag");
    name
}