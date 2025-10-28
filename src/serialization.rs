use crate::{DIM, Matrix};
use crate::flag::{Flag, FlagGraph};

use std::io::{Read, Write};

pub fn load_flag_file(mut file: std::fs::File) -> FlagGraph {
    // parse numbers
    let mut string: String = String::new();
    file.read_to_string(&mut string).unwrap();
    let mut nums = string
        .split(|c: char| c.is_whitespace())     // isolate each number
        .flat_map(|n| n.parse::<f32>());  // parse into numbers

    // convert this data into the vertices of the flags/simplices
    let mut flag_graph: FlagGraph = FlagGraph {flags: Vec::new(), edges: Vec::new()};
    for c in std::iter::from_fn(move || {
        let chunk: Vec<f32> = nums.by_ref().take(DIM*DIM).collect(); // take DIM^2 at a time
        if chunk.is_empty() {None} else {Some(chunk)}
    }) {
        flag_graph.flags.push(Flag { vertices: Matrix::from_vec(c) })
    }
    flag_graph
}

pub fn save_flags_to_file(filename: String, flag_graph: &FlagGraph) -> Result<(), std::io::Error> {
    let mut f = std::fs::File::create_new(filename)?;
    for flag in flag_graph.flags.iter() {
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