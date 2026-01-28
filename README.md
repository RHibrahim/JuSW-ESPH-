# JuSW-ESPH: A High-Performance Eulerian SPH Solver for Shallow Water Equations

JuSW-ESPH is a high-performance Eulerian Smoothed Particle Hydrodynamics (SPH) solver written in Julia, designed for accurate and efficient simulation of shallow water flows. It is built for speed, robustness, and clarity, making it suitable for research, benchmarking, and free-surface flow simulations.

✨ **Key Features**

-	Eulerian SPH formulation for shallow water equations

-	HLL Riemann solver for stable and accurate flux computation

-	Well-balanced scheme preserving lake-at-rest conditions

-	Leap-frog time stepping: explicit, stable, and efficient

-	CPU multithreading (Julia Threads)

-	GPU acceleration via CUDA.jl

-	Modular code structure for easy extension and experimentation

-	Benchmark-driven simulations with configurable parameters

⚙️ **Requirements**

- Julia ≥ 1.9

- For GPU version:

  - NVIDIA GPU
  - CUDA.jl

▶️ **Usage**

1️⃣ Configure the Simulation

Each benchmark comes with a text input file describing the simulation setup.

Before running a simulation:

1. Open src/parameters.jl
2. Adjust the parameters according to the selected benchmark input file.
   
2️⃣ Run the Solver

Launch Julia and execute:   include("SPH_RUN.jl"). 

That’s it ✅

🚀 **CPU vs GPU Versions**

This repository provides two versions of the solver:

- CPU version : 
  Optimized with Julia multithreading for shared-memory systems

- GPU version : 
  Fully accelerated using CUDA for large-scale, high-resolution simulations

Both versions share the same physical model and numerical scheme, ensuring consistent results across architectures.
