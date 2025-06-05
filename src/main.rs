// filepath: /rust-1d-fem/rust-1d-fem/src/main.rs

fn main() {
  let x_min = -1.0; // Calculation domain minimum
  let x_max = 1.0; // Calculation domain maximum
  let node_total = 4; // Number of nodes (>=2)
  let ele_total = node_total - 1; // Number of elements
  let func_f = 1.0; // Constant function f
  let alpha = 1.0; // Dirichlet boundary condition at the left end
  let beta = -1.0; // Neumann boundary condition at the right end

  let node_x_glo = (0..node_total).map(|i| x_min + i as f64 * (x_max - x_min) / (node_total - 1) as f64).collect::<Vec<f64>>();

  let mut mat_a_glo = vec![vec![0.0; node_total]; node_total];
  let mut vec_b_glo = vec![0.0; node_total];

  // Element matrix and vector initialization
  let mut mat_a_ele = vec![vec![vec![0.0; 3]; 3]; ele_total];
  let mut vec_b_ele = vec![vec![0.0; 3]; ele_total];

  // Calculate element lengths
  let length: Vec<f64> = (0..ele_total).map(|e| (node_x_glo[e + 1] - node_x_glo[e]).abs()).collect();

  // Calculate local matrices
  for e in 0..ele_total {
    for i in 0..2 {
      for j in 0..2 {
        mat_a_ele[e][i][j] = ((-1i32).pow((i + 1) as u32) * (-1i32).pow((j + 1) as u32)) as f64 / length[e];
      }
      vec_b_ele[e][i] = -func_f * length[e] / 2.0;
    }
  }

  // Assemble global matrix
  for e in 0..ele_total {
    for i in 0..2 {
      for j in 0..2 {
        mat_a_glo[e + i][e + j] += mat_a_ele[e][i][j];
      }
      vec_b_glo[e + i] += vec_b_ele[e][i];
    }
  }

  // Apply boundary conditions
  boundary(0, alpha, 0.0, &mut mat_a_glo, &mut vec_b_glo);
  boundary(node_total - 1, f64::INFINITY, beta, &mut mat_a_glo, &mut vec_b_glo);

  // Solve the linear system
  let unknown_vec_u = solve(&mat_a_glo, &vec_b_glo);

  // Output results
  println!("Unknown vector u: {:?}", unknown_vec_u);
  println!(
    "Max u: {}, Min u: {}",
    unknown_vec_u.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    unknown_vec_u.iter().cloned().fold(f64::INFINITY, f64::min)
  );
}

// Function to apply boundary conditions
fn boundary(node_num_glo: usize, dirichlet: f64, neumann: f64, mat_a_glo: &mut Vec<Vec<f64>>, vec_b_glo: &mut Vec<f64>) {
  if dirichlet != f64::INFINITY {
    vec_b_glo.iter_mut().for_each(|b| *b -= dirichlet * mat_a_glo[node_num_glo][node_num_glo]);
    vec_b_glo[node_num_glo] = dirichlet;
    for j in 0..mat_a_glo.len() {
      mat_a_glo[node_num_glo][j] = 0.0;
      mat_a_glo[j][node_num_glo] = 0.0;
    }
    mat_a_glo[node_num_glo][node_num_glo] = 1.0;
  }

  if neumann != f64::INFINITY {
    vec_b_glo[node_num_glo] += neumann;
  }
}

// Function to solve the linear system
fn solve(_mat_a_glo: &Vec<Vec<f64>>, vec_b_glo: &Vec<f64>) -> Vec<f64> {
  // Implement a solver for the linear system here (e.g., Gaussian elimination)
  // Placeholder for the actual implementation
  vec![0.0; vec_b_glo.len()] // Replace with actual solution
}
