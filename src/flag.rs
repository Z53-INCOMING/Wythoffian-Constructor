use crate::{DIM, Matrix, Vector};

#[derive(Clone)]
pub struct Flag {
    pub vertices: Matrix, // points of the fundamental simplex
}

impl Flag {
    /// takes one flag and a set of mirrors, and returns all flags
    pub fn find_all(&self, mirrors: Matrix) -> Vec<Self> {
        // start finding all the other fundamental simplices
        let mut flags: Vec<Flag> = vec![self.clone()];

        // algorithm state (starting simplex, current reflection)
        let mut stack: Vec<(Flag, usize)> = vec![(self.clone(), 0)];
        // depth-first search n-tree of potential symmetries
        'wythoffian_outer: while stack.len() > 0 {
            /*for (i, (_, m)) in stack.iter().enumerate() {
                if i == stack.len() - 1 {print!("{m}");}
                else {print!("{}", m-1);}
            }println!();*/

            // get flag and mirror id from top of stack
            let (base_flag, mirror) = stack.last_mut().unwrap();
            let new_flag = base_flag.reflect(mirrors.column(*mirror).into());

            // if this flag has been visited before, move on
            for f in flags.iter() {
                if new_flag.compare(&f) == DIM {
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

        flags
    }

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
    pub fn compare(&self, other: &Self) -> usize {
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