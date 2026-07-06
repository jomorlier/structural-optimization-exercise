# Homework Assignment: Structural Topology Optimization from Scratch

In this assignment, you will explore the mathematical foundations of structural topology optimization, analyze its convexity properties, and implement a highly efficient 2D solver in Python based on the Solid Isotropic Material with Penalization (SIMP) method.

---

## 1. Problem Formulation & Physics

Consider a 2D design domain discretized into a grid of $N_x \times N_y$ rectangular finite elements. The goal is to distribute a fixed amount of material within this domain to design a structure that supports an external load $\mathbf{f}$ while minimizing its deformability (compliance).

Let $x_i \in [x_{\min}, 1]$ be the design variable representing the physical density of finite element $i$, where $x_i = 1$ denotes solid material, $x_i = x_{\min}$ (a small positive value, e.g., $10^{-9}$) represents a void (to prevent the global stiffness matrix from becoming singular), and $\mathbf{x} = [x_1, x_2, \dots, x_N]^T$ is the vector of all element densities.

Under the **SIMP** model, the Young's modulus $E_i$ of element $i$ is interpolated as:
$$E_i(x_i) = E_{\min} + x_i^p (E_0 - E_{\min})$$
where $E_0$ is the modulus of the solid material, $E_{\min}$ is the modulus of the void, and $p \ge 1$ is the penalization power. This elasticity interpolation scales the element stiffness matrix:
$$\mathbf{K}_i(x_i) = E_i(x_i) \mathbf{K}_0$$
where $\mathbf{K}_0$ is the baseline stiffness matrix of a solid 4-node quad element. The global stiffness matrix $\mathbf{K}(\mathbf{x})$ is assembled by summing the contributions of all element stiffness matrices:
$$\mathbf{K}(\mathbf{x}) = \sum_{i=1}^N \mathbf{K}_i(x_i)$$

### Question 1: Physical Interpretation of Compliance
The compliance of the structure is defined as:
$$C(\mathbf{x}) = \mathbf{f}^T \mathbf{u}$$
where $\mathbf{u}$ is the nodal displacement vector obtained from the static equilibrium equation:
$$\mathbf{K}(\mathbf{x}) \mathbf{u} = \mathbf{f}$$
1. Show that the compliance can also be written as $C(\mathbf{x}) = \mathbf{u}^T \mathbf{K}(\mathbf{x}) \mathbf{u}$.
2. Explain the physical meaning of compliance. What does minimizing compliance correspond to? Why is it twice the total strain energy of the structure?

### Question 2: Optimization Setup
Write down the complete mathematical optimization problem for compliance minimization. Clearly state:
1. The objective function.
2. The design variables and their bounds.
3. The equality constraint representing the physics of the system.
4. The inequality constraint limiting the total volume of material to a fraction $f$ of the total design domain volume $V_0$.

---

## 2. Convexity Analysis

Convexity is a fundamental property of optimization problems. An optimization problem is convex if the objective function is convex, the inequality constraint functions are convex, and the equality constraint functions are affine. Convex problems have the nice property that any local minimum is a global minimum.

### Question 3: The Case $p = 1$ (Linear Elasticity Interpolation)
Suppose $p = 1$ and $E_{\min} = 0$. The Young's modulus is linear in $x_i$: $E_i(x_i) = x_i E_0$.
1. Show that the global stiffness matrix $\mathbf{K}(\mathbf{x})$ is a linear function of $\mathbf{x}$.
2. The compliance objective can be written as $C(\mathbf{x}) = \mathbf{f}^T \mathbf{K}(\mathbf{x})^{-1} \mathbf{f}$. Using the fact that the function $\mathbf{A} \mapsto \mathbf{c}^T \mathbf{A}^{-1} \mathbf{c}$ is convex on the set of symmetric positive-definite matrices $\mathbf{A}$, prove that the compliance objective $C(\mathbf{x})$ is a convex function of $\mathbf{x}$ for $\mathbf{x} > \mathbf{0}$.
3. Is the overall optimization problem convex when $p = 1$? Explain.

### Question 4: The Case $p = 3$ (SIMP Penalization)
In practice, we use $p \ge 3$ (typically $p = 3$).
1. Explain why we use $p = 3$ instead of $p = 1$. (Hint: What kind of designs do we get with $p=1$, and what does the penalization power do to intermediate densities like $x_i = 0.5$?)
2. Is the optimization problem convex when $p = 3$? Discuss the implications of this for gradient-based optimization algorithms (e.g., convergence, dependency on initial designs, local minima).
3. What is a **continuation method** and how does it help alleviate issues caused by non-convexity in SIMP?

---

## 3. Computational Efficiency

Implementing structural optimization requires solving physical simulations inside an optimization loop. Writing efficient code is critical.

### Question 5: Dense vs. Sparse Solvers
At each iteration of the optimizer, we must solve $\mathbf{K}(\mathbf{x}) \mathbf{u} = \mathbf{f}$ for the displacements $\mathbf{u}$.
1. If the design domain is $100 \times 100$ elements, how many degrees of freedom (DOFs) does the system have?
2. Compare the computational complexity (in terms of Big-O notation) of solving this system using a dense matrix solver (e.g., `numpy.linalg.solve`) versus a sparse direct solver (e.g., nested dissection or Cholesky/LU factorization on sparse banded matrices). Why is the sparse solver essential?

### Question 6: Adjoint Sensitivity Analysis
To update the design variables using gradient-based algorithms (like MMA), we need the sensitivities (derivatives) of the compliance w.r.t the densities, i.e., $\frac{\partial C}{\partial x_i}$.
1. Suppose we have $N$ elements. If we compute these derivatives using finite differences (e.g., perturbing each $x_i$ by $\epsilon$ and re-solving the equilibrium equations), how many times must we solve the sparse linear system of equations per optimization step?
2. Show that by using the **adjoint method**, we can obtain all $N$ derivatives by solving the linear system **exactly once**. Prove that the compliance gradient is:
   $$\frac{\partial C}{\partial x_i} = - \mathbf{u}_i^T \frac{\partial \mathbf{K}_i}{\partial x_i} \mathbf{u}_i = - p x_i^{p-1} (E_0 - E_{\min}) \mathbf{u}_i^T \mathbf{K}_0 \mathbf{u}_i$$
   where $\mathbf{u}_i$ is the vector of displacements of the 8 degrees of freedom associated with element $i$.
   *(Hint: Start from $C = \mathbf{f}^T \mathbf{u}$, differentiate with respect to $x_i$, and substitute the derivative of $\mathbf{K}\mathbf{u} = \mathbf{f}$.)*

### Question 7: Vectorization and NumPy's `einsum`
Evaluating compliance requires computing $\mathbf{u}_i^T \mathbf{K}_0 \mathbf{u}_i$ for every element $i = 1, \dots, N$.
In Python, writing a loop over all elements to compute these quadratic forms is extremely slow.
Explain how the code in `solution.py` uses `numpy.einsum` to perform this calculation in a single vectorized step:
```python
ke_u = anp.einsum('ij,jkl->ikl', ke, u_selected)
ce = anp.einsum('ijk,ijk->jk', u_selected, ke_u)
```
Specifically:
1. What is the shape of `u_selected` and what does each axis represent?
2. Explain the contractions performed by the two `einsum` calls.

### Question 8: Hashable LRU Caching
In `solution.py`, the author introduces a custom NumPy array LRU cache decorator:
```python
@ndarray_safe_lru_cache(maxsize=1)
def _get_solver(a_entries, a_indices, size, sym_pos):
    ...
```
1. Why does standard Python `@functools.lru_cache` fail when applied to functions with NumPy array arguments?
2. What are the major performance bottlenecks in the optimization loop that do not change their structure (only their values), and how does caching intermediate structural arrays (like boundary condition indexes or solver matrices) speed up execution?

---

## 4. Programming Assignment

Your practical task is to complete the template code in [template.py](file:///Users/j.morlier/.gemini/antigravity/scratch/structural_optimization_exercise/template.py).

### Files Provided:
- `template.py`: Code skeleton.
- `solution.py`: Reference solution with optimization caching.
- `test_cases.py`: Script to test your progress.

### Tasks to Complete in `template.py`:
1. **`young_modulus`**: Implement the SIMP Young's modulus interpolation: $E(x) = E_{\min} + x^p (E_0 - E_{\min})$.
2. **`compliance`**: Compute the global compliance of the structure using vectorized operations.
3. **`get_stiffness_matrix`**: Formulate the stiffness matrix $\mathbf{K}_0$ for a 2D square element.
4. **`_get_dof_indices`**: Implement the index-mapping array for node degrees of freedom to assemble the global sparse stiffness matrix.
5. **`solve_coo`**: Define a custom autograd-primitive solver for the sparse system of equations, and implement its vector-Jacobian product (VJP) for automatic differentiation.

Run `python test_cases.py` to verify each component as you code. Good luck!
