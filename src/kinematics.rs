use crate::math_utils::{ad, exponential_map};
use na::{DMatrix, Matrix3, Matrix4, Matrix6, Vector3, Vector6};
use std::f64::consts::PI;

pub fn ur5_params() -> (Matrix4<f64>, DMatrix<f64>) {
    let l1 = 0.5;
    let l2 = 0.45;
    let h1 = 0.18;
    let h2 = 0.1;
    let w1 = 0.15;
    let w2 = 0.09;

    let initial_position = Matrix4::new(
        -1.0,
        0.0,
        0.0,
        l1 + l2,
        0.0,
        0.0,
        1.0,
        w1 + w2,
        0.0,
        1.0,
        0.0,
        h1 - h2,
        0.0,
        0.0,
        0.0,
        1.0,
    );

    let screw = DMatrix::from_row_slice(
        6,
        6,
        &[
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            -h1,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            -h1,
            0.0,
            l1,
            0.0,
            1.0,
            0.0,
            -h1,
            0.0,
            l1 + l2,
            0.0,
            0.0,
            -1.0,
            -w1,
            l1 + l2,
            0.0,
            0.0,
            1.0,
            0.0,
            h2 - h1,
            0.0,
            l1 + l2,
        ],
    );

    (initial_position, screw)
}

pub fn compute_t(
    screw: &DMatrix<f64>,
    theta: &Vector6<f64>,
    dof: usize,
    m: &Matrix4<f64>,
) -> Matrix4<f64> {
    let mut t = Matrix4::identity();
    for i in 0..dof {
        t = t * exponential_map(
            &screw.column(i).fixed_slice::<6, 1>(0, 0).into_owned(),
            theta[i],
        );
    }
    t * m
}

pub fn root_finding(
    theta_0: Vector6<f64>,
    theta_d: Vector6<f64>,
    try_max: usize,
    dof: usize,
    initial_position: Matrix4<f64>,
    screw: &DMatrix<f64>,
) -> (Vector6<f64>, Matrix4<f64>, usize) {
    let mut n_try = 1;
    let tol = 0.0001;
    let mut theta = theta_0;
    let mut e = compute_e(&theta_d, &theta, dof, &initial_position, screw);

    while n_try < try_max && e.norm() > tol {
        let ja = compute_jacobian(&theta, screw, dof);
        let inverse_jacobian = ja.try_inverse().unwrap();
        theta = theta + inverse_jacobian * (theta_d - theta);
        e = compute_e(&theta_d, &theta, dof, &initial_position, screw);
        n_try += 1;
    }
    (theta, e, n_try)
}

fn compute_e(
    theta_d: &Vector6<f64>,
    theta_0: &Vector6<f64>,
    dof: usize,
    initial_position: &Matrix4<f64>,
    screw: &DMatrix<f64>,
) -> Matrix4<f64> {
    let t_0 = compute_t(screw, theta_0, dof, initial_position);
    let t_d = compute_t(screw, theta_d, dof, initial_position);
    t_d - t_0
}

fn compute_jacobian(theta: &Vector6<f64>, screw: &DMatrix<f64>, dof: usize) -> Matrix6<f64> {
    let mut g_local = Vec::new();
    let mut exp_map_local = Vec::new();

    for i in 0..dof {
        g_local.push(denavit_transformation(theta[i], i));
        exp_map_local.push(
            g_local[i]
                * exponential_map(
                    &screw.column(i).fixed_slice::<6, 1>(0, 0).into_owned(),
                    theta[i],
                ),
        );
    }

    let mut g = Vec::new();
    let mut g_accum = Matrix4::identity();
    for i in 0..dof {
        g_accum = g_accum * exp_map_local[i];
        g.push(g_accum);
    }

    let mut j_s = Vec::new();
    for i in 0..dof {
        j_s.push(ad(&g[i]) * screw.column(i).fixed_slice::<6, 1>(0, 0).into_owned());
    }

    let p_k = Vector3::zeros();
    let p_0_extended = g[5] * p_k.push(1.0);
    let p_0 = p_0_extended.fixed_slice::<3, 1>(0, 0);

    let p_0_cpo = Matrix3::new(0.0, -p_0.z, p_0.y, p_0.z, 0.0, -p_0.x, -p_0.y, p_0.x, 0.0);

    let j_g = Matrix6::new(
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        -p_0_cpo[(0, 1)],
        p_0_cpo[(0, 2)],
        1.0,
        0.0,
        0.0,
        p_0_cpo[(0, 1)],
        0.0,
        -p_0_cpo[(0, 0)],
        0.0,
        1.0,
        0.0,
        -p_0_cpo[(0, 2)],
        p_0_cpo[(0, 0)],
        0.0,
        0.0,
        0.0,
        1.0,
    );

    let r = g[5].fixed_slice::<3, 3>(0, 0);
    let r_pitch =
        (-r[(2, 0)] / (r[(2, 2)] * r[(2, 2)] + r[(2, 1)] * r[(2, 1)]).sqrt()).atan2(r[(2, 2)]);
    let r_yaw = (r[(1, 0)] / r[(0, 0)]).atan2(r[(0, 0)]);

    let b = Matrix3::new(
        (r_pitch.cos() * r_yaw.cos()).cos(),
        -r_yaw.sin(),
        0.0,
        (r_pitch.cos() * r_yaw.cos()).sin(),
        r_yaw.cos(),
        0.0,
        -r_pitch.sin(),
        0.0,
        1.0,
    );

    let j_a = Matrix6::new(
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        b[(0, 0)],
        b[(0, 1)],
        b[(0, 2)],
        0.0,
        0.0,
        0.0,
        b[(1, 0)],
        b[(1, 1)],
        b[(1, 2)],
        0.0,
        0.0,
        0.0,
        b[(2, 0)],
        b[(2, 1)],
        b[(2, 2)],
    );

    j_a * j_g
}

fn denavit_transformation(theta: f64, index: usize) -> Matrix4<f64> {
    let d = [0.1, 0.0, 0.0, 0.12, 0.1, 0.06];
    let alpha = [PI / 2.0, 0.0, 0.0, PI / 2.0, -PI / 2.0, 0.0];
    let r = [0.0, -0.5, -0.45, 0.0, 0.0, 0.0];

    let c_theta = theta.cos();
    let s_theta = theta.sin();
    let c_alpha = alpha[index].cos();
    let s_alpha = alpha[index].sin();

    let r_z = Matrix4::new(
        c_theta, -s_theta, 0.0, 0.0, s_theta, c_theta, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
        1.0,
    );

    let t_z = Matrix4::new(
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, d[index], 0.0, 0.0, 0.0, 1.0,
    );

    let t_x = Matrix4::new(
        1.0, 0.0, 0.0, r[index], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    );

    let r_x = Matrix4::new(
        1.0, 0.0, 0.0, 0.0, 0.0, c_alpha, -s_alpha, 0.0, 0.0, s_alpha, c_alpha, 0.0, 0.0, 0.0, 0.0,
        1.0,
    );

    r_z * t_z * t_x * r_x
}
