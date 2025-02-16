# Inverse Kinematics and Motion Planning in Robotics

## Overview
This project implements inverse kinematics (IK) for a 6-degree-of-freedom (DOF) robotic arm using numerical methods. The code uses the **Newton-Raphson root-finding algorithm** to iteratively compute joint angles that achieve a desired end-effector position and orientation.

## Degrees of Freedom (DOF) in Robotics
In robotics, **DOF** refers to the number of independent movements a robotic system can perform. A robotic arm with **6 DOF** can move in three translational directions (X, Y, Z) and rotate about three axes (roll, pitch, yaw). This makes it capable of reaching any pose in 3D space.

A **6 DOF robotic manipulator**, like the UR5 used in this project, can:
- Move forward and backward
- Move left and right
- Move up and down
- Rotate around its base
- Tilt its wrist up and down
- Rotate its wrist

## What is Inverse Kinematics (IK)?
**Inverse Kinematics** is the process of calculating the joint angles required for a robotic arm to reach a desired end-effector position and orientation.

### Why Do We Calculate Inverse Kinematics?
Unlike **forward kinematics**, which computes the position of the end-effector given the joint angles, inverse kinematics solves the opposite problem: finding the necessary joint angles for a desired end-effector pose.

## Implementation in This Project
This project computes inverse kinematics for a **UR5 robotic arm** by:
1. Defining the **robot parameters** (link lengths, offsets, and screw axes) in `kinematics.rs`.
2. Using **exponential mapping** to calculate transformation matrices.
3. Implementing **numerical root-finding** to iteratively solve for joint angles (`root_finding` function).

### Motion Planning and Practical Examples
Motion planners such as **CHOMP, RRT, and Trajectory Optimization** rely on inverse kinematics to generate feasible paths. For example:
- **MoveIt! framework** in ROS uses IK solvers to plan smooth motions.
- **Trajectory optimization** algorithms adjust robot paths while satisfying constraints.

In this project, a motion trajectory is planned by:
1. Setting an **initial pose** and **goal pose**.
2. Computing the **desired transformation matrix**.
3. Iteratively solving for joint angles using **root finding**.
4. Evaluating the **final transformation matrix** to verify correctness.

## How to Run the Project
### Dependencies
- Rust programming language
- `nalgebra` library for matrix operations

### Running the Code
```sh
cargo build
```

```sh
cargo run
```
