use std::f32::consts::PI;

use nalgebra::{ArrayStorage, Const};

// number of dimensions for polytope
const DIM: usize = 3;
// Coxeter diagram's matrix and ringed nodes
const COXMAT: [[i8; DIM]; DIM] = [
    [1, 5, 2],
    [5, 1, 3],
    [2, 3, 1],
];
const RINGS: [bool; DIM] = [false, false, true];

type Vector = nalgebra::Matrix<f32, Const<DIM>, Const<1>, ArrayStorage<f32, DIM, 1>>;
type Matrix = nalgebra::Matrix<f32, Const<DIM>, Const<DIM>, ArrayStorage<f32, DIM, DIM>>;

#[derive(Clone)]
struct FundamentalSimplex {
    vertices: Matrix, // points of the actual simplex
}

fn main() {
    // matrix for which a_ij is the dot product of vector i and vector j
    let dot_matrix: Matrix = Matrix::from_fn(
        |r, c| {
            if r == c {1.} // diagonals should be 1 always
            else {f32::cos(PI / COXMAT[r][c] as f32)} // other entries are angles based on coxeter matrix
        }
    );

    // columns of Cholesky decomposition upper-triangular martix are normal to mirrors
    let mirrors: Matrix = dot_matrix.cholesky().expect("invalid coxeter diagram").l().transpose();
    // Find a starting vertex of the resulting polytope
    let starting_simplex = FundamentalSimplex::from_mirrors(mirrors);
    
    // start finding all the flags
    let mut fundamental_simplices: Vec<FundamentalSimplex> = vec![starting_simplex.clone()];

    // algorithm state (starting simplex, current reflection)
    let mut stack: Vec<(FundamentalSimplex, usize)> = vec![(starting_simplex, 0)];
    // depth-first search n-tree of potential symmetries
    'wythoffian_outer: while stack.len() > 0 {
        /*for (i, (_, m)) in stack.iter().enumerate() {
            if i == stack.len() - 1 {print!("{m}");}
            else {print!("{}", m-1);}
        }println!();*/

        // get simplex and mirror id from top of stack
        let (base_simplex, mirror) = stack.last_mut().unwrap();
        let new_simplex = base_simplex.reflect(mirrors.column(*mirror).into());

        // if this simplex has been visited before, move on
        for f in fundamental_simplices.iter() {
            if new_simplex.compare(&f) == DIM {
                *mirror += 1;
                if *mirror == DIM {stack.pop();}
                continue 'wythoffian_outer;
            }
        }

        // otherwise, proceed by moving one step deeper
        fundamental_simplices.push(new_simplex.clone());
        *mirror += 1;
        if *mirror == DIM {stack.pop();}
        stack.push((new_simplex, 0));
    }

    // find polytope elements
    let mut vertices: Vec<Vector> = Vec::new();
    for s in fundamental_simplices {
        let v = s.rings_to_point(RINGS);
        if vertices.iter().find(|f| (*f - v).magnitude_squared() < 0.0001).is_none() {
            vertices.push(v);
        }
    }
    
    print!("V = [");
    for v in vertices.iter().take(vertices.len()-1) {
        print!("(");
        for r in v.iter().take(DIM-1) {
            print!("{}, ", r);
        }
        print!("{}), ", v[DIM-1]);
    }
    let v = vertices.last().unwrap();
    print!("(");
    for r in v.iter().take(DIM-1) {
        print!("{}, ", r);
    }
    println!("{})] ", v[DIM-1]);

}

impl FundamentalSimplex {
    fn from_mirrors(m: Matrix) -> Self {
        Self { 
            // basically the cofactor matrix lol
            vertices: Matrix::from_fn(|r, c| {
                (-1f32).powi(r as i32) * m.remove_row(r).remove_column(c).determinant()
            }) 
        }
    }
    fn reflect(&self, v: Vector) -> Self {
        let mut vertices = self.vertices.clone();
        for mut c in vertices.column_iter_mut() {
            c -= 2. * (c.dot(&v)) / v.magnitude_squared() * v;
        }
        Self {vertices}
    }
    fn rings_to_point(&self, rings: [bool; DIM]) -> Vector {
        self.vertices
            .column_iter()
            .zip(rings.iter())
            .map(|(m, r)| if *r {m.into()} else {Vector::zeros()})
            .sum::<Vector>()
            .normalize()
    }
    fn compare(&self, other: &Self) -> usize {
        self.vertices.column_iter().map(|self_col| 
            if other.vertices.column_iter().find(|other_col| 
                (self_col - other_col).magnitude_squared() < 0.00001
            ).is_some() {1} else {0}
        ).sum()
    }
}