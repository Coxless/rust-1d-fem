// This file contains the core logic of the program, including the implementation of the finite element method, boundary condition handling, matrix assembly, and solving the linear system.

use ndarray::Zip;
use ndarray::{Array1, Array2};
use ndarray_linalg::solve::Inverse;
use std::f64;

pub struct FiniteElement {
  pub node_total: usize,
  pub ele_total: usize,
  pub alpha: f64,
  pub beta: f64,
  pub func_f: f64,
  pub mat_a_glo: Array2<f64>,
  pub vec_b_glo: Array1<f64>,
}

impl FiniteElement {
  pub fn new(node_total: usize, alpha: f64, beta: f64, func_f: f64) -> Self {
    let ele_total = node_total - 1;
    let mat_a_glo = Array2::<f64>::zeros((node_total, node_total));
    let vec_b_glo = Array1::<f64>::zeros(node_total);

    FiniteElement {
      node_total,
      ele_total,
      alpha,
      beta,
      func_f,
      mat_a_glo,
      vec_b_glo,
    }
  }

  pub fn boundary(&mut self, node_num_glo: usize, dirichlet: Option<f64>, neumann: Option<f64>) {
    if let Some(value) = dirichlet {
      let row = self.mat_a_glo.row(node_num_glo).to_owned();
      Zip::from(&mut self.vec_b_glo).and(&row).for_each(|b, &a| *b -= value * a);
      self.vec_b_glo[node_num_glo] = value;
      self.mat_a_glo.row_mut(node_num_glo).fill(0.0);
      self.mat_a_glo.column_mut(node_num_glo).fill(0.0);
      self.mat_a_glo[(node_num_glo, node_num_glo)] = 1.0;
    }

    if let Some(value) = neumann {
      self.vec_b_glo[node_num_glo] += value;
    }
  }

  pub fn assemble_global_matrix(&mut self) {
    // Implementation for assembling the global matrix from element matrices
  }

  pub fn solve(&self) -> Array1<f64> {
    let inv = self.mat_a_glo.inv().unwrap();
    inv.dot(&self.vec_b_glo)
  }

  pub fn print_results(&self, unknown_vec_u: &Array1<f64>) {
    println!("Unknown vector u = {:?}", unknown_vec_u);
    println!(
      "Max u = {}, Min u = {}",
      unknown_vec_u.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
      unknown_vec_u.iter().cloned().fold(f64::INFINITY, f64::min)
    );
  }
}
