use na::{Matrix3, Matrix4, Matrix6, Vector3, Vector6};

pub fn cross_product_operator(vector: &Vector3<f64>) -> Matrix3<f64> {
    let w1 = vector.x;
    let w2 = vector.y;
    let w3 = vector.z;
    Matrix3::new(0.0, -w3, w2, w3, 0.0, -w1, -w2, w1, 0.0)
}

pub fn exponential_map(action_axis: &Vector6<f64>, theta: f64) -> Matrix4<f64> {
    let linear = action_axis.fixed_slice::<3, 1>(3, 0).into_owned();
    let angular = action_axis.fixed_slice::<3, 1>(0, 0).into_owned();

    let exp_rot = exponential_form_rotation(angular, theta);
    let exp_tras = exponential_form_translation(linear, angular, theta);

    Matrix4::new(
        exp_rot[(0, 0)],
        exp_rot[(0, 1)],
        exp_rot[(0, 2)],
        exp_tras.x,
        exp_rot[(1, 0)],
        exp_rot[(1, 1)],
        exp_rot[(1, 2)],
        exp_tras.y,
        exp_rot[(2, 0)],
        exp_rot[(2, 1)],
        exp_rot[(2, 2)],
        exp_tras.z,
        0.0,
        0.0,
        0.0,
        1.0,
    )
}

pub fn exponential_form_rotation(angular: Vector3<f64>, theta: f64) -> Matrix3<f64> {
    let c = cross_product_operator(&angular) * cross_product_operator(&angular);
    Matrix3::identity() + theta.sin() * cross_product_operator(&angular) + (1.0 - theta.cos()) * c
}

pub fn exponential_form_translation(
    linear: Vector3<f64>,
    angular: Vector3<f64>,
    theta: f64,
) -> Vector3<f64> {
    let angular_mat = cross_product_operator(&angular);
    let c = angular_mat * angular_mat;
    (theta * Matrix3::identity() + (1.0 - theta.cos()) * angular_mat + (theta - theta.sin()) * c)
        * linear
}

pub fn ad(mat4: &Matrix4<f64>) -> Matrix6<f64> {
    let rot = mat4.fixed_slice::<3, 3>(0, 0);
    let tras = mat4.fixed_slice::<3, 1>(0, 3);
    let ad = Matrix6::new(
        rot[(0, 0)],
        rot[(0, 1)],
        rot[(0, 2)],
        0.0,
        0.0,
        0.0,
        rot[(1, 0)],
        rot[(1, 1)],
        rot[(1, 2)],
        0.0,
        0.0,
        0.0,
        rot[(2, 0)],
        rot[(2, 1)],
        rot[(2, 2)],
        0.0,
        0.0,
        0.0,
        0.0,
        -tras.z,
        tras.y,
        rot[(0, 0)],
        rot[(0, 1)],
        rot[(0, 2)],
        tras.z,
        0.0,
        -tras.x,
        rot[(1, 0)],
        rot[(1, 1)],
        rot[(1, 2)],
        -tras.y,
        tras.x,
        0.0,
        rot[(2, 0)],
        rot[(2, 1)],
        rot[(2, 2)],
    );
    ad
}
