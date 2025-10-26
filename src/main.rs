use std::{f32::consts::PI, io::{Read, Write}};
use nalgebra::{ArrayStorage, Const};
use macroquad::prelude::*;

// number of dimensions for polytope
const DIM: usize = 4;
// Coxeter diagram's matrix and ringed nodes
const COXMAT: [[u8; DIM]; DIM] = [
    [1, 3, 2, 2],
    [3, 1, 4, 2],
    [2, 4, 1, 3],
    [2, 2, 3, 1],
];
const RINGS: [bool; DIM] = [true, false, false, false];

const SCALE: f32 = 800.;

type Vector = nalgebra::Matrix<f32, Const<DIM>, Const<1>, ArrayStorage<f32, DIM, 1>>;
type Matrix = nalgebra::Matrix<f32, Const<DIM>, Const<DIM>, ArrayStorage<f32, DIM, DIM>>;

#[derive(Clone)]
struct FundamentalSimplex {
    vertices: Matrix, // points of the actual simplex
}

#[macroquad::main("wireframe")]
async fn main() {
    // compute fundamental simplices / flags
    let fundamental_simplices = if let Ok(mut file) = std::fs::File::open(coxmat_to_name(COXMAT)) {
        println!("found cached flag file");
        
        // parse numbers
        let mut string: String = String::new();
        file.read_to_string(&mut string).unwrap();
        let mut nums = string
            .split(|c: char| c.is_whitespace())     // isolate each number
            .flat_map(|n| n.parse::<f32>());  // parse into numbers

        // convert this data into the vertices of the simplices
        let mut fundamental_simplices: Vec<FundamentalSimplex> = Vec::new();
        for c in std::iter::from_fn(move || {
            let chunk: Vec<f32> = nums.by_ref().take(DIM*DIM).collect(); // take DIM^2 at a time
            if chunk.is_empty() {None} else {Some(chunk)}
        }) {
            fundamental_simplices.push(FundamentalSimplex { vertices: Matrix::from_vec(c) })
        }
        fundamental_simplices
    } else {
        println!("no cached flag file found, generating...");

        // matrix for which a_ij is the dot product of vector i and vector j
        let dot_matrix: Matrix = Matrix::from_fn(
            |r, c| {
                if r == c {1.} // diagonals should be 1 always
                else {f32::cos(PI / COXMAT[r][c] as f32)} // other entries are angles based on coxeter matrix
            }
        );

        // find mirrors, then one fundamental simplex from those. use these to generate all fundamental simplices / flags
        let mirrors: Matrix = dot_matrix.cholesky().expect("invalid coxeter diagram").l().transpose();
        let starting_simplex = FundamentalSimplex::from_mirrors(mirrors);
        let fundamental_simplices = starting_simplex.find_all(mirrors);

        // cache the result
        let mut f = std::fs::File::create_new(coxmat_to_name(COXMAT)).expect("couldn't cache output (???)");
        for v in fundamental_simplices.iter() {
            let v = v.vertices;
            for n in v.iter() {
                write!(f, "{} ", *n).unwrap();
            }
            writeln!(f).unwrap();
        }

        fundamental_simplices
    };
    println!("{} flags computed", fundamental_simplices.len());
    

    // find polytope elements
    // vertex list without duplicates
    let mut vertices: Vec<Vector> = Vec::new();
    let mut flag_vertex: Vec<usize> = Vec::new();
    for s in fundamental_simplices.iter() {
        let v = s.rings_to_point(RINGS);
        if let Some((i, _)) = vertices.iter().enumerate().find(|(_, f)| (*f - v).magnitude_squared() < 0.0001) {
            flag_vertex.push(i)
        }
        else {
            vertices.push(v);
            flag_vertex.push(vertices.len() - 1);
        }
    }
    println!("{} vertices computed", vertices.len());

    //edge list without self-connections
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for (s, i) in fundamental_simplices.iter().zip(flag_vertex.iter()) {
        fundamental_simplices.iter().zip(flag_vertex.iter()).filter(|(f, _)| s.compare(f) == DIM - 1).for_each(|(_, j)| {
            let (i, j) = {
                if i < j {(*i, *j)}
                else {(*j, *i)}
            };
            if i != j && !edges.contains(&(i, j)) {
                edges.push((i, j));
            }
        })
    }
    println!("{} edges computed", edges.len());

    // Draw polytope to screen
    println!("rendering");
    let mut width: f32 = 2.;
    let mut scale = SCALE;
    loop {
        clear_background(BLACK);

        for (a, b) in edges.iter() {
            let a_perspective = (vertices[*a].z + 2.)*(vertices[*a].w + 2.);
            let b_perspective = (vertices[*b].z + 2.)*(vertices[*b].w + 2.);
            draw_line(scale*vertices[*a].x/a_perspective + screen_width()/2., scale*vertices[*a].y/a_perspective + screen_height()/2., scale*vertices[*b].x/b_perspective + screen_width()/2., scale*vertices[*b].y/b_perspective + screen_height()/2., width, WHITE);
        }

        if is_key_down(KeyCode::Up) {scale += screen_width()/60.}
        if is_key_down(KeyCode::Down) {scale -= screen_width()/60.}
        if is_key_down(KeyCode::Left) {width += 0.1}
        if is_key_down(KeyCode::Right) {width -= 0.1}

        for v in vertices.iter_mut() {
            *v = Matrix::new(
                f32::cos(1./60.), 0., 0., -f32::sin(1./60.), 
                0., f32::cos(1./60.), -f32::sin(1./60.), 0.,
                0., f32::sin(1./60.), f32::cos(1./60.), 0.,
                f32::sin(1./60.), 0., 0., f32::cos(1./60.),
            ) * *v;
        }

        next_frame().await
    }
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
    fn find_all(&self, mirrors: Matrix) -> Vec<Self> {
        // start finding all the other fundamental simplices
        let mut fundamental_simplices: Vec<FundamentalSimplex> = vec![self.clone()];

        // algorithm state (starting simplex, current reflection)
        let mut stack: Vec<(FundamentalSimplex, usize)> = vec![(self.clone(), 0)];
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

        fundamental_simplices
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

fn coxmat_to_name(m: [[u8; DIM]; DIM]) -> String {
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