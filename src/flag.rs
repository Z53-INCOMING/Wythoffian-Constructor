use crate::{DIM, Matrix, Vector};
use std::io::{Read, Write};

pub struct FlagGraph {
    pub flags: Vec<Flag>,
    pub edges: Vec<(usize, usize, u32)> // indices into flag list, + number of shared vertices
}

#[derive(Clone)]
pub struct Flag {
    pub vertices: Matrix, // points of the fundamental simplex
}

impl FlagGraph {
    /// takes one flag and a set of mirrors, and returns all flags and their connections
    pub fn generate(start_flag: Flag, mirrors: Matrix) -> Self {
        // start finding all the other fundamental simplices
        let mut flags: Vec<Flag> = vec![start_flag.clone()];

        // algorithm state (starting simplex, current reflection)
        let mut stack: Vec<(Flag, usize)> = vec![(start_flag.clone(), 0)];
        // depth-first search n-tree of potential symmetries
        'wythoffian_outer: while stack.len() > 0 {
            // get flag and mirror id from top of stack
            let (base_flag, mirror) = stack.last_mut().unwrap();
            let new_flag = base_flag.reflect(mirrors.column(*mirror).into());

            // if this flag has been visited before, move on
            for f in flags.iter() {
                if new_flag.compare(&f) == DIM as u32 {
                    *mirror += 1;
                    if *mirror == DIM {stack.pop();}
                    continue 'wythoffian_outer;
                }
            }

            // otherwise, proceed by moving one step deeper
            flags.push(new_flag.clone());
            *mirror += 1;
            if *mirror == DIM {stack.pop();}
            stack.push((new_flag, 0));
        }

        Self {
            flags,
            edges: Vec::new(),
        }
    }

    pub fn serialize(&self, filename: String) -> Result<(), std::io::Error> {
        let mut f = std::fs::File::create_new(filename)?;
        for flag in self.flags.iter() {
            let vertices = flag.vertices;
            for v in vertices.iter() {
                write!(f, "{} ", *v)?;
            }
            writeln!(f)?;
        }
        Ok(())
    }

    pub fn deserialize(filename: String) -> Result<Self, std::io::Error> {
        // open input file
        let mut f = std::fs::File::open(filename)?;
        // get text from file
        let mut string: String = String::new();
        f.read_to_string(&mut string)?;
        // input parsing iterator
        let mut nums = string
            .split(|c: char| c.is_whitespace())     // isolate each number
            .flat_map(|n| n.parse::<f32>());  // parse into numbers

        // convert parsed input numbers to vertices of flags/simplices
        let mut flag_graph: FlagGraph = FlagGraph {flags: Vec::new(), edges: Vec::new()};
        for c in std::iter::from_fn(move || {
            let chunk: Vec<f32> = nums.by_ref().take(DIM*DIM).collect(); // take DIM^2 at a time
            if chunk.is_empty() {None} else {Some(chunk)}
        }) {
            flag_graph.flags.push(Flag { vertices: Matrix::from_vec(c) })
        }
        Ok(flag_graph)
    }
}

impl Flag {
    /// reflects the entire flag along some vector.
    /// the reflection fixes the hyperplane orthogonal to vector v.
    pub fn reflect(&self, v: Vector) -> Self {
        let mut vertices = self.vertices.clone();
        for mut c in vertices.column_iter_mut() {
            c -= 2. * (c.dot(&v)) / v.magnitude_squared() * v;
        }
        Self {vertices}
    }

    /// generates a vertex based on which nodes of the CD are ringed.
    pub fn rings_to_point(&self, rings: [bool; DIM]) -> Vector {
        self.vertices
            .column_iter()
            .zip(rings.iter())
            .map(|(m, r)| if *r {m.into()} else {Vector::zeros()})
            .sum::<Vector>()
            .normalize()
    }

    /// returns the number of vertices two simplices share
    pub fn compare(&self, other: &Self) -> u32 {
        self.vertices.column_iter().map(|self_col| 
            if other.vertices.column_iter().find(|other_col| 
                (self_col - other_col).magnitude_squared() < 0.00001
            ).is_some() {1} else {0}
        ).sum()
    }

    /// generates one fundamental simplex/flag from a matrix whose columns are reflection vectors
    pub fn from_mirrors(m: Matrix) -> Self {
        Self { 
            // basically the cofactor matrix lol
            vertices: Matrix::from_fn(|r, c| {
                (-1f32).powi(r as i32) * m.remove_row(r).remove_column(c).determinant()
            })
        }
    }
}