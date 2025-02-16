mod kinematics;
mod math_utils;

extern crate nalgebra as na;

use kinematics::{compute_t, root_finding, ur5_params};
use na::{Quaternion, UnitQuaternion, Vector3, Vector6};
use std::f64::consts::PI;

fn main() {
    let dof = 6;
    let (initial_position, screw) = ur5_params();

    // Initial Position
    let initial_position_vec = Vector3::new(0.2, 0.4, 0.2);
    let initial_orientation = UnitQuaternion::new_normalize(Quaternion::from_parts(
        0.0,
        Vector3::new(-0.7071068, 0.0, 0.7071068),
    ));

    // Goal Position
    let goal_position = Vector3::new(0.38, 0.125, 1.108);
    let goal_orientation = UnitQuaternion::new_normalize(Quaternion::from_parts(
        0.0,
        Vector3::new(0.7071068, 0.7071068, 0.0),
    ));

    // Convert positions and orientations to joint angles (theta_0 and theta_d)
    let theta_0 = Vector6::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0); // Placeholder
    let theta_d = Vector6::new(0.0, -PI / 2.0, 0.0, 0.0, PI / 2.0, 0.0); // Placeholder

    let t_0 = compute_t(&screw, &theta_0, dof, &initial_position);
    let t_d = compute_t(&screw, &theta_d, dof, &initial_position);

    println!("Initial position: \n{}", t_0);
    println!("Desired position: \n{}", t_d);

    let (theta, _, _) = root_finding(theta_0, theta_d, 20, dof, initial_position, &screw);
    let t_theta = compute_t(&screw, &theta, dof, &initial_position);

    println!("Solution: \n{}", theta);
    println!("Transformation: \n{}", t_theta);
}
