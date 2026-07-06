# Coding Exercise: Structural Topology Optimization from Scratch

Welcome to this coding exercise on structural topology optimization! In this project, you will study the mathematics of material design under physical constraints, implement a 2D finite element method (FEM) solver in Python, and optimize a load-bearing structure to minimize structural compliance.

Topology optimization is a powerful engineering method that determines the optimal layout of material within a given design domain to support specified loads and boundary conditions. The resulting structures often have highly organic, bio-mimetic shapes (like bones or tree roots).

## Learning Objectives

By completing this exercise, you will learn to:
1. **Formulate** compliance minimization problems mathematically with density design variables and volume constraints.
2. **Translate** physical equilibrium equations (linear elasticity PDEs) into a discrete finite element model (FEM).
3. **Analyze** the mathematical properties (specifically convexity and local minima) of the Solid Isotropic Material with Penalization (SIMP) interpolation model.
4. **Implement** highly efficient numerical routines in Python using sparse matrices, vectorized element matrices, and adjoint sensitivity analysis.
5. **Optimize** structures using gradient-based mathematical programming (specifically the Method of Moving Asymptotes, MMA).

---

## Repository Structure

- [exercise.md](file:///Users/j.morlier/.gemini/antigravity/scratch/structural_optimization_exercise/exercise.md): The assignment description, covering optimization formulation, physical assumptions, convexity analysis, numerical efficiency questions, and programming instructions.
- [solutions.md](file:///Users/j.morlier/.gemini/antigravity/scratch/structural_optimization_exercise/solutions.md): The analytical answers to the theoretical questions and reference implementation code.
- [template.py](file:///Users/j.morlier/.gemini/antigravity/scratch/structural_optimization_exercise/template.py): A skeleton script with marked `TODO` blocks where you will write your implementation.
- [solution.py](file:///Users/j.morlier/.gemini/antigravity/scratch/structural_optimization_exercise/solution.py): A fully functioning, optimized reference solver featuring caching wrappers for numpy arrays.
- [test_cases.py](file:///Users/j.morlier/.gemini/antigravity/scratch/structural_optimization_exercise/test_cases.py): A unit-testing suite to verify your completed implementation.

---

## Getting Started

### 1. Prerequisites

This exercise requires Python 3.7+ and several standard scientific computing libraries. Install the dependencies using `pip`:

```bash
pip install numpy scipy matplotlib autograd nlopt
```

*Note: On some systems, compiling `nlopt` requires a C compiler. In Google Colab or conda environments, it can be installed via `pip install nlopt` or `conda install -c conda-forge nlopt`.*

### 2. Running the Reference Solution

To run the reference solution and see structural optimization in action:

```bash
python solution.py
```

This will run an 80-step optimization of an MBB beam (a classic structural benchmark) and display the resulting organic-looking support structure.

### 3. Completing the Assignment

1. Read the theoretical questions in [exercise.md](file:///Users/j.morlier/.gemini/antigravity/scratch/structural_optimization_exercise/exercise.md) and write down your answers.
2. Open [template.py](file:///Users/j.morlier/.gemini/antigravity/scratch/structural_optimization_exercise/template.py) and complete the functions marked with `# TODO`.
3. Test your code using the test suite:
   ```bash
   python test_cases.py
   ```
4. Check your analytical answers and code against [solutions.md](file:///Users/j.morlier/.gemini/antigravity/scratch/structural_optimization_exercise/solutions.md).
