# Percolation Thresholds and Connectivity in Quantum Networks

This repository contains the code and data associated with the paper:

**Percolation Thresholds and Connectivity in Quantum Networks**  
*Author(s): Andrea De Girolamo, Giuseppe Magnifico, Cosmo Lupo*  
*Preprint: [arXiv:2502.17191](https://arxiv.org/abs/2502.17191)*

## 📄 Abstract

We study entanglement percolation in qubit-based planar quantum network models of arbitrary topology, where neighboring nodes are initially connected by pure states with quenched disorder in their entanglement. To address this, we develop a physics-informed heuristic algorithm designed to find a sequence of entanglement swapping and distillation operations to connect any pair of distant nodes.

The algorithm combines locally optimal percolation strategies between nodes at a maximum distance of one swapping operation. If this fails to produce a maximally entangled state, it looks for alternative paths surrounding intermediate states within the process. We analytically find and numerically verify thresholds in quantum percolation, which depend on the initial network configuration and entanglement, and are associated with specific percolation strategies.

We classify these strategies based on the *connectivity*, a quantity that relates the entanglement in the final state and the level of integrity of the network at the end of the process. We find distinct regimes of quantum percolation, clearly separated by percolation thresholds that vary significantly with the underlying topology.

## 🚀 Project Overview

This codebase implements:

- A simulation environment for quantum network topologies
- Modeling of initial entangled link states with **quenched disorder** in entanglement
- A **heuristic algorithm** for entanglement percolation using:
  - Local entanglement swapping strategies
  - Exploration of alternative paths + heuristics
