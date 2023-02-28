## Matlab Code

This repository contains the Matlab files for the simulation of QDMC and linear MPC applied to cloth manipulation. The simulations can be run by simply executing the files `qdmc_simulation.m`, `lmpc_simulation.m`, `qdmc_obstacle_avoidance_chance.m`, `lmpc_obstacle_avoidance_chance.m` in Matlab.

## QDMC ROS Nodes

Furthermore, this repository contains the C++ implementation of the ROS nodes used in our real-world experiment. The code was adapted from the implementation available at [https://github.com/Alados5/mpc_node.git], by Adrià Luque Acera.

The nodes are `qdmc_rt_node.cpp` and `qdmc_opti_node.cpp` for running the standard QDMC algorithm.

**Contact author:** Edoardo Caldarelli, ecaldarelli@iri.upc.edu

