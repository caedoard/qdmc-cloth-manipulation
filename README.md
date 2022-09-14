## Matlab Code

This repository contains the Matlab files for the simulation of QDMC and linear MPC applied to cloth manipulation. THe simulations can be run by simply executing the files `qdmc_simulation.m` and `lmpc_simulation.m` in Matlab.

## QDMC ROS Nodes

Furthermore, this repository contains the C++ implementation of the ROS nodes used in our real-world experiments. The code was adapted from the implementation available at [https://github.com/Alados5/mpc_node.git], by Adrià Luque Acera.

The nodes are:
- `qdmc_rt_node.cpp` and `qdmc_opti_node.cpp` for running the standard QDMC algorithm without disturbances or chance constraints;
- `qdmc_rt_node_wind_disturbance.cpp` and `qdmc_opti_node_wind_disturbance.cpp` for running the standard QDMC algorithm with a known wind profile incorporated in the prediction model;
- `qdmc_opti_node_stochastic.cpp` for replacing the deterministic constraints of standard QDMC with our chance constrained formulation (used for obstacle avoidance).

**Contact author:** Edoardo Caldarelli, ecaldarelli@iri.upc.edu

