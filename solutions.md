# Solution Manual: Structural Topology Optimization from Scratch

This document provides detailed analytical answers to the theoretical questions in [exercise.md](file:///Users/j.morlier/.gemini/antigravity/scratch/structural_optimization_exercise/exercise.md) and presents the complete reference code implementation.

---

## 1. Problem Formulation & Physics

### Question 1: Physical Interpretation of Compliance
1. **Derivation:**
   We are given the compliance definition $C(\mathbf{x}) = \mathbf{f}^T \mathbf{u}$ and the static equilibrium equation:
   ```math
\mathbf{K}(\mathbf{x}) \mathbf{u} = \mathbf{f}
```
   Since the external load vector $\mathbf{f}$ is equal to $\mathbf{K}(\mathbf{x}) \mathbf{u}$, we can substitute this into the compliance equation:
   ```math
C(\mathbf{x}) = (\mathbf{K}(\mathbf{x}) \mathbf{u})^T \mathbf{u} = \mathbf{u}^T \mathbf{K}(\mathbf{x})^T \mathbf{u}
```
   In linear elasticity, the global stiffness matrix $\mathbf{K}(\mathbf{x})$ is symmetric because it is assembled from symmetric element stiffness matrices ($\mathbf{K}_i^T = \mathbf{K}_i$). Thus, $\mathbf{K}(\mathbf{x})^T = \mathbf{K}(\mathbf{x})$.
   Substituting this symmetry property yields:
   ```math
C(\mathbf{x}) = \mathbf{u}^T \mathbf{K}(\mathbf{x}) \mathbf{u}
```

2. **Physical Meaning:**
   - **Compliance** is a measure of the overall structural deformability under the applied load. It represents the external work done by the applied forces: $W = \mathbf{f}^T \mathbf{u}$.
   - **Stiffness** is the inverse of compliance. Minimizing compliance is mathematically identical to maximizing the structural stiffness under the given loading conditions.
   - The total strain energy stored in an elastic structure is given by:
     ```math
U_e = \frac{1}{2} \mathbf{u}^T \mathbf{K}(\mathbf{x}) \mathbf{u}
```
     Comparing this to compliance, we see that $C(\mathbf{x}) = 2 U_e$. Compliance is exactly twice the total elastic strain energy of the structure.

### Question 2: Optimization Setup
The nested formulation of the topology optimization problem is written as:
```math
\begin{aligned}
\min_{\mathbf{x}} \quad & C(\mathbf{x}) = \mathbf{f}^T \mathbf{u}(\mathbf{x}) \\
\text{subject to} \quad & \mathbf{K}(\mathbf{x}) \mathbf{u}(\mathbf{x}) = \mathbf{f} \\
& \frac{\sum_{i=1}^N x_i v_i}{V_0} \le f \\
& x_{\min} \le x_i \le 1, \quad i = 1, \dots, N
\end{aligned}
```

Where:
- **Objective function:** Structural compliance $C(\mathbf{x}) = \mathbf{f}^T \mathbf{u}(\mathbf{x})$.
- **Design variables:** Element densities $\mathbf{x} = [x_1, \dots, x_N]^T$.
- **Lower/Upper bounds:** $x_i \in [x_{\min}, 1]$.
- **Physical constraint (equality):** State equation $\mathbf{K}(\mathbf{x}) \mathbf{u} = \mathbf{f}$ (implicitly defines $\mathbf{u}(\mathbf{x})$).
- **Resource constraint (inequality):** Total volume fraction limit $\frac{V(\mathbf{x})}{V_0} \le f$.

---

## 2. Convexity Analysis

### Question 3: The Case $p = 1$ (Linear Interpolation)
1. **Linearity of $\mathbf{K}(\mathbf{x})$:**
   With $p=1$ and $E_{\min}=0$, the element Young's modulus is $E_i(x_i) = x_i E_0$.
   The global stiffness matrix is the sum of element matrices:
   ```math
\mathbf{K}(\mathbf{x}) = \sum_{i=1}^N \mathbf{K}_i(x_i) = \sum_{i=1}^N x_i E_0 \mathbf{K}_0^{(i)}
```
   where $\mathbf{K}_0^{(i)}$ is the baseline solid element stiffness matrix expanded to global dimensions. Since $\mathbf{K}(\mathbf{x})$ is a linear combination of the scalars $x_i$, it is a linear (affine) function of $\mathbf{x}$.

2. **Convexity of Compliance:**
   Let $g(\mathbf{A}) = \mathbf{f}^T \mathbf{A}^{-1} \mathbf{f}$ be a function mapping a symmetric positive-definite (SPD) matrix $\mathbf{A}$ to a scalar.
   The function $g(\mathbf{A})$ is convex on the cone of SPD matrices. To see this, we look at the Hessian or show that for any $\mathbf{A}, \mathbf{B} \succ 0$ and $\theta \in [0, 1]$:
   ```math
\theta \mathbf{f}^T \mathbf{A}^{-1} \mathbf{f} + (1-\theta) \mathbf{f}^T \mathbf{B}^{-1} \mathbf{f} \ge \mathbf{f}^T (\theta \mathbf{A} + (1-\theta) \mathbf{B})^{-1} \mathbf{f}
```
   Since the global stiffness matrix $\mathbf{K}(\mathbf{x})$ is a linear function of $\mathbf{x}$, and for $\mathbf{x} > \mathbf{0}$ the matrix $\mathbf{K}(\mathbf{x})$ is SPD, the composition $C(\mathbf{x}) = g(\mathbf{K}(\mathbf{x}))$ of the convex function $g$ with the linear function $\mathbf{K}$ is convex in $\mathbf{x}$.

3. **Overall Convexity:**
   In the nested formulation, the variables $\mathbf{u}$ are eliminated by writing them as $\mathbf{u}(\mathbf{x}) = \mathbf{K}(\mathbf{x})^{-1} \mathbf{f}$.
   - The objective $C(\mathbf{x})$ is convex.
   - The volume constraint $\frac{1}{V_0} \sum_i x_i v_i - f \le 0$ is linear, which is convex.
   - The bound constraints $x_{\min} \le x_i \le 1$ define a convex set.
   
   Therefore, in the **nested formulation**, the optimization problem is **convex**.
   *Note: If we do not eliminate $\mathbf{u}$ and optimize over both $\mathbf{x}$ and $\mathbf{u}$ simultaneously, the equilibrium equation $\mathbf{K}(\mathbf{x})\mathbf{u} = \mathbf{f}$ acts as an equality constraint. Since $\mathbf{K}(\mathbf{x})\mathbf{u}$ is bilinear in $\mathbf{x}$ and $\mathbf{u}$, this equality constraint is non-affine, making the simultaneous formulation non-convex.*

### Question 4: The Case $p = 3$ (SIMP Penalization)
1. **Why we use $p = 3$:**
   If we solve the convex case ($p=1$), the optimizer allocates intermediate densities (porous material, e.g., $x_i = 0.5$) over large regions. Such structures are extremely expensive or impossible to manufacture.
   By setting $p = 3$, we penalize intermediate densities. For example, at $x_i = 0.5$:
   - The element consumes $50\%$ of its solid volume.
   - Its stiffness is only $0.5^3 = 12.5\%$ of its solid stiffness.
   This interpolation makes intermediate densities highly inefficient. The optimizer is forced to drive densities towards the bounds ($x_i \approx 0$ or $x_i = 1$) to get maximum stiffness per unit mass, resulting in clean solid-void designs.

2. **Non-Convexity & Implications:**
   For $p > 1$, the Young's modulus interpolation is non-convex, which makes the objective function $C(\mathbf{x})$ non-convex.
   - The optimization landscape contains many local minima.
   - Gradient-based optimization algorithms (like MMA) will converge to a local minimum.
   - The final structural design is sensitive to the initial guess, mesh size, filter parameters, and the optimization path.

3. **Continuation Methods:**
   A continuation method solves a sequence of optimization problems where $p$ is gradually increased (e.g., $p = 1.0 \to 1.5 \to 2.0 \to 3.0$).
   By starting with $p=1$, the algorithm solves a convex problem first to find the optimal global distribution of stiffness. As $p$ is incrementally raised, the design is gradually forced to become discrete, preventing the optimizer from getting trapped in poor local minima early in the process.

---

## 3. Computational Efficiency

### Question 5: Dense vs. Sparse Solvers
1. **Degrees of Freedom:**
   For $100 \times 100$ elements, there are $101 \times 101 = 10,201$ nodes.
   Since each node has 2 degrees of freedom (horizontal and vertical displacements), the system size is:
   ```math
\text{DOFs} = 2 \times 10,201 = 20,402 \text{ equations}
```

2. **Computational Complexity:**
   - **Dense Solvers ($O(M^3)$):** Solving a system of $M = 20,402$ equations using a dense LU/Cholesky solver requires $M^3 \approx 8.5 \times 10^{12}$ floating-point operations. Doing this on every optimization iteration is computationally prohibitive.
   - **Sparse Solvers ($O(M^{1.5})$ to $O(M^2)$):** In finite element grids, nodes only interact with neighboring nodes. The stiffness matrix is highly sparse and banded (only $0.05\%$ of entries are non-zero). Direct sparse solvers (like CHOLMOD or SuperLU) use graph-reordering techniques to minimize fill-in, solving the system in a fraction of a second.

### Question 6: Adjoint Sensitivity Analysis
1. **Finite Differences:**
   To compute the gradient of compliance w.r.t $N$ density variables using finite differences:
   ```math
\frac{\partial C}{\partial x_i} \approx \frac{C(\mathbf{x} + \epsilon \mathbf{e}_i) - C(\mathbf{x})}{\epsilon}
```
   We must perturb each element $i$ and solve the equilibrium equations. This requires **$N$ sparse matrix solves** per optimization step. For $N = 10,000$, this is impossible to run.

2. **Adjoint Derivation:**
   Differentiating the compliance objective $C = \mathbf{f}^T \mathbf{u}$ w.r.t $x_i$:
   ```math
\frac{d C}{d x_i} = \mathbf{f}^T \frac{\partial \mathbf{u}}{\partial x_i}
```
   Differentiating the state equation $\mathbf{K}\mathbf{u} = \mathbf{f}$ w.r.t $x_i$ (where $\mathbf{f}$ is constant):
   ```math
\frac{\partial \mathbf{K}}{\partial x_i} \mathbf{u} + \mathbf{K} \frac{\partial \mathbf{u}}{\partial x_i} = \mathbf{0} \implies \frac{\partial \mathbf{u}}{\partial x_i} = - \mathbf{K}^{-1} \frac{\partial \mathbf{K}}{\partial x_i} \mathbf{u}
```
   Substitute this back into the compliance derivative:
   ```math
\frac{d C}{d x_i} = - \mathbf{f}^T \mathbf{K}^{-1} \frac{\partial \mathbf{K}}{\partial x_i} \mathbf{u}
```
   Define the adjoint variable $\boldsymbol{\lambda}$ as the solution to:
   ```math
\mathbf{K}^T \boldsymbol{\lambda} = \mathbf{f}
```
   Since $\mathbf{K}$ is symmetric, $\mathbf{K}^T = \mathbf{K}$, which implies $\boldsymbol{\lambda} = \mathbf{K}^{-1} \mathbf{f} = \mathbf{u}$ (the problem is self-adjoint).
   Substituting $\boldsymbol{\lambda} = \mathbf{u}$ yields:
   ```math
\frac{d C}{d x_i} = - \mathbf{u}^T \frac{\partial \mathbf{K}}{\partial x_i} \mathbf{u}
```
   Since $\mathbf{K}(\mathbf{x}) = \sum_j \mathbf{K}_j(x_j)$, the derivative $\frac{\partial \mathbf{K}}{\partial x_i}$ is only non-zero for element $i$:
   ```math
\frac{d C}{d x_i} = - \mathbf{u}_i^T \frac{\partial \mathbf{K}_i}{\partial x_i} \mathbf{u}_i
```
   Using SIMP interpolation $\mathbf{K}_i(x_i) = \left( E_{\min} + x_i^p(E_0 - E_{\min}) \right) \mathbf{K}_0$:
   ```math
\frac{\partial \mathbf{K}_i}{\partial x_i} = p x_i^{p-1} (E_0 - E_{\min}) \mathbf{K}_0
```
   Substituting this:
   ```math
\frac{\partial C}{\partial x_i} = - p x_i^{p-1} (E_0 - E_{\min}) \mathbf{u}_i^T \mathbf{K}_0 \mathbf{u}_i
```

   **Efficiency:** The adjoint method only requires **1 sparse solve** per iteration (to find $\mathbf{u}$). Once $\mathbf{u}$ is known, all $N$ gradients are calculated via cheap local dot products, reducing the gradient computation time to almost zero.

### Question 7: Vectorization and NumPy's `einsum`
1. **Shape of `u_selected`:** `(8, nely, nelx)`.
   - Axis 0 (size 8): The 8 degrees of freedom associated with the 4 nodes of a quadrilateral element (horizontal and vertical displacements).
   - Axis 1 (`nely`): The vertical element grid coordinate.
   - Axis 2 (`nelx`): The horizontal element grid coordinate.

2. **einsum Contractions:**
   - `ke_u = anp.einsum('ij,jkl->ikl', ke, u_selected)`:
     Multiplies the $8 \times 8$ element stiffness matrix `ke` (indices `ij`) by the displacements `u_selected` (indices `jkl`) for all elements. The index `j` is contracted (summed over), yielding `ke_u` of shape `(8, nely, nelx)` (indices `ikl`). This represents $\mathbf{K}_0 \mathbf{u}_e$ for every element.
   - `ce = anp.einsum('ijk,ijk->jk', u_selected, ke_u)`:
     Performs an element-wise dot product between the displacements `u_selected` and the resulting forces `ke_u` along the degree-of-freedom axis (index `i`), yielding `ce` of shape `(nely, nelx)` (indices `jk`). This calculates the quadratic form $\mathbf{u}_e^T \mathbf{K}_0 \mathbf{u}_e$ for all elements simultaneously without any Python loops.

### Question 8: Hashable LRU Caching
1. **Why standard cache fails:**
   NumPy arrays are mutable and do not implement `__hash__`. When passed as arguments to standard `@functools.lru_cache`, Python raises a `TypeError: unhashable type: 'numpy.ndarray'`. The `ndarray_safe_lru_cache` solves this by wrapping arrays in a hashable class that computes a hash based on the flattened string representation of the array.

2. **Sparsity Patterns and Adjoint reuse:**
   - The degree-of-freedom boundary mapping and sparsity locations (`x_list`, `y_list`) do not change across iterations. Caching these index setups saves significant index assembly time.
   - **Crucially**, during each optimization step, we solve the system $\mathbf{K}\mathbf{u} = \mathbf{f}$ (forward pass) and then evaluate the automatic gradient which calls `solve_coo` again (backward pass) on the same stiffness matrix $\mathbf{K}$. By caching `_get_solver` with `maxsize=1`, the backward pass reuses the already factorized LU/Cholesky solver from the forward pass, making the gradient calculation practically free.

---

## 4. Completed Solution Code (`solution.py`)

Below is the complete reference implementation containing the vectorized solver and array caching:

```python
# A Tutorial on Structural Optimization | Sam Greydanus | 2022
import time, functools, nlopt
import numpy as np                                                # for dense matrix ops
import matplotlib.pyplot as plt                                   # for plotting
import autograd, autograd.core, autograd.extend, autograd.tracer  # for adjoints
import autograd.numpy as anp      
import scipy, scipy.ndimage, scipy.sparse, scipy.sparse.linalg    # sparse matrices

##### Problem setup #####
class ObjectView(object):
    def __init__(self, d): self.__dict__ = d
    
def get_args(normals, forces, density=0.4):  # Manage the problem setup parameters
  width = normals.shape[0] - 1
  height = normals.shape[1] - 1
  fixdofs = np.flatnonzero(normals.ravel())
  alldofs = np.arange(2 * (width + 1) * (height + 1))
  freedofs = np.sort(list(set(alldofs) - set(fixdofs)))
  params = {
      # material properties
      'young': 1, 'young_min': 1e-9, 'poisson': 0.3, 'g': 0,
      # constraints
      'density': density, 'xmin': 0.001, 'xmax': 1.0,
      # input parameters
      'nelx': width, 'nely': height, 'mask': 1, 'penal': 3.0, 'filter_width': 1,
      'freedofs': freedofs, 'fixdofs': fixdofs, 'forces': forces.ravel(),
      # optimization parameters
      'opt_steps': 80, 'print_every': 10}
  return ObjectView(params)

def mbb_beam(width=80, height=25, density=0.4, y=1, x=0):  # textbook beam example
  normals = np.zeros((width + 1, height + 1, 2))
  normals[-1, -1, y] = 1
  normals[0, :, x] = 1
  forces = np.zeros((width + 1, height + 1, 2))
  forces[0, 0, y] = -1
  return normals, forces, density

def young_modulus(x, e_0, e_min, p=3):
  return e_min + x ** p * (e_0 - e_min)

def physical_density(x, args, volume_contraint=False, use_filter=True):
  x = args.mask * x.reshape(args.nely, args.nelx)  # reshape from 1D to 2D
  return gaussian_filter(x, args.filter_width) if use_filter else x  # maybe filter

def mean_density(x, args, volume_contraint=False, use_filter=True):
  return anp.mean(physical_density(x, args, volume_contraint, use_filter)) / anp.mean(args.mask)

##### Caching for NumPy arrays (gives a 2x speedup) #####
class _WrappedArray:
  """Hashable wrapper for NumPy arrays."""
  def __init__(self, value):
    self.value = value

  def __eq__(self, other):
    return np.array_equal(self.value, other.value)

  def __hash__(self):
    return hash(repr(self.value.ravel()))

def ndarray_safe_lru_cache(maxsize=128):
  """An ndarray compatible version of functools.lru_cache."""
  def decorator(func):
    @functools.lru_cache(maxsize)
    def cached_func(*args, **kwargs):
      args = tuple(a.value if isinstance(a, _WrappedArray) else a for a in args)
      kwargs = {k: v.value if isinstance(v, _WrappedArray) else v for k, v in kwargs.items()}
      return func(*args, **kwargs)

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
      args = tuple(_WrappedArray(a) if isinstance(a, np.ndarray) else a for a in args)
      kwargs = {k: _WrappedArray(v) if isinstance(v, np.ndarray) else v for k, v in kwargs.items()}
      return cached_func(*args, **kwargs)
    return wrapper
  return decorator

##### Optimization objective + physics of elastic materials #####
def objective(x, args, volume_contraint=False, use_filter=True):
  kwargs = dict(penal=args.penal, e_min=args.young_min, e_0=args.young)
  x_phys = physical_density(x, args, volume_contraint=volume_contraint, use_filter=use_filter)
  ke     = get_stiffness_matrix(args.young, args.poisson)  # stiffness matrix
  u      = displace(x_phys, ke, args.forces, args.freedofs, args.fixdofs, **kwargs)
  c      = compliance(x_phys, u, ke, **kwargs)
  return c

def compliance(x_phys, u, ke, *, penal=3, e_min=1e-9, e_0=1):
  nely, nelx = x_phys.shape
  ely, elx = anp.meshgrid(range(nely), range(nelx))  # x, y coords for the index map

  n1 = (nely+1)*(elx+0) + (ely+0)  # nodes
  n2 = (nely+1)*(elx+1) + (ely+0)
  n3 = (nely+1)*(elx+1) + (ely+1)
  n4 = (nely+1)*(elx+0) + (ely+1)
  all_ixs = anp.array([2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n3, 2*n3+1, 2*n4, 2*n4+1])
  u_selected = u[all_ixs]  # select from u matrix

  ke_u = anp.einsum('ij,jkl->ikl', ke, u_selected)  # compute x^penal * U.T @ ke @ U
  ce = anp.einsum('ijk,ijk->jk', u_selected, ke_u)
  C = young_modulus(x_phys, e_0, e_min, p=penal) * ce.T
  return anp.sum(C)

def get_stiffness_matrix(e, nu):  # e=young's modulus, nu=poisson coefficient
  k = anp.array([1/2-nu/6, 1/8+nu/8, -1/4-nu/12, -1/8+3*nu/8,
                -1/4+nu/12, -1/8-nu/8, nu/6, 1/8-3*nu/8])
  return e/(1-nu**2)*anp.array([[k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
                               [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
                               [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
                               [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
                               [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
                               [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
                               [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
                               [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]])
  
def get_k(stiffness, ke):
  # Constructs sparse stiffness matrix k (used in the displace fn)
  nely, nelx = stiffness.shape
  ely, elx = anp.meshgrid(range(nely), range(nelx))  # x, y coords
  ely, elx = ely.reshape(-1, 1), elx.reshape(-1, 1)

  n1 = (nely+1)*(elx+0) + (ely+0)
  n2 = (nely+1)*(elx+1) + (ely+0)
  n3 = (nely+1)*(elx+1) + (ely+1)
  n4 = (nely+1)*(elx+0) + (ely+1)
  edof = anp.array([2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n3, 2*n3+1, 2*n4, 2*n4+1])
  edof = edof.T[0]
  x_list = anp.repeat(edof, 8)  # flat list pointer of each node in an element
  y_list = anp.tile(edof, 8).flatten()  # flat list pointer of each node in elem

  # make the global stiffness matrix K
  kd = stiffness.T.reshape(nelx*nely, 1, 1)
  value_list = (kd * anp.tile(ke, kd.shape)).flatten()
  return value_list, y_list, x_list

def displace(x_phys, ke, forces, freedofs, fixdofs, *, penal=3, e_min=1e-9, e_0=1):
  # Displaces the load x using finite element techniques (solve_coo=most of runtime)
  stiffness = young_modulus(x_phys, e_0, e_min, p=penal)
  k_entries, k_ylist, k_xlist = get_k(stiffness, ke)
  index_map, keep, indices = _get_dof_indices(freedofs, fixdofs, k_ylist, k_xlist)
  
  u_nonzero = solve_coo(k_entries[keep], indices, forces[freedofs], sym_pos=True)
  u_values = anp.concatenate([u_nonzero, anp.zeros(len(fixdofs))])
  return u_values[index_map]
  
##### Sparse matrix (COO) helper functions #####
@ndarray_safe_lru_cache(1)
def _get_dof_indices(freedofs, fixdofs, k_xlist, k_ylist):  # sparse matrix helper function
  index_map = inverse_permutation(anp.concatenate([freedofs, fixdofs]))
  keep = anp.isin(k_xlist, freedofs) & anp.isin(k_ylist, freedofs)
  i = index_map[k_ylist][keep]
  j = index_map[k_xlist][keep]
  return index_map, keep, anp.stack([i, j])

def inverse_permutation(indices):  # reverses an index operation
  inverse_perm = np.zeros(len(indices), dtype=anp.int64)
  inverse_perm[indices] = np.arange(len(indices), dtype=anp.int64)
  return inverse_perm

@ndarray_safe_lru_cache(1)
def _get_solver(a_entries, a_indices, size, sym_pos):
  a = scipy.sparse.coo_matrix((a_entries, a_indices), shape=(size,)*2).tocsc()
  return scipy.sparse.linalg.splu(a).solve

##### Autograd custom gradients #####
@autograd.primitive
def solve_coo(a_entries, a_indices, b, sym_pos=False):
  solver = _get_solver(a_entries, a_indices, b.size, sym_pos)
  return solver(b)

def grad_solve_coo_entries(ans, a_entries, a_indices, b, sym_pos=False):
  def jvp(grad_ans):
    lambda_ = solve_coo(a_entries, a_indices if sym_pos else a_indices[::-1],
                        grad_ans, sym_pos)
    i, j = a_indices
    return -lambda_[i] * ans[j]
  return jvp
autograd.extend.defvjp(solve_coo, grad_solve_coo_entries,
                       lambda: print('err: gradient undefined'),
                       lambda: print('err: gradient not implemented'))

@autograd.extend.primitive
def gaussian_filter(x, width): # 2D gaussian blur/filter
  return scipy.ndimage.gaussian_filter(x, width, mode='reflect')

def _gaussian_filter_vjp(ans, x, width): # gives the gradient of orig. function w.r.t. x
  del ans, x  # unused
  return lambda g: gaussian_filter(g, width)
autograd.extend.defvjp(gaussian_filter, _gaussian_filter_vjp)

##### Main optimization function #####
def fast_stopt(args, x=None, verbose=True):
  if x is None:
    x = anp.ones((args.nely, args.nelx)) * args.density  # init mass
  reshape = lambda x: x.reshape(args.nely, args.nelx)
  objective_fn = lambda x: objective(reshape(x), args)
  constraint = lambda params: mean_density(reshape(params), args) - args.density

  def wrap_autograd_func(func, losses=None, frames=None):
    def wrapper(x, grad):
      if grad.size > 0:
        value, grad[:] = autograd.value_and_grad(func)(x)
      else:
        value = func(x)
      if losses is not None:
        losses.append(value)
      if frames is not None:
        frames.append(reshape(x).copy())
        if verbose and len(frames) % args.print_every == 0:
          print('step {}, loss {:.2e}, t={:.2f}s'.format(len(frames), value, time.time()-dt))
      return value
    return wrapper

  losses, frames = [], [] ; dt = time.time()
  print('Optimizing a problem with {} nodes'.format(len(args.forces)))
  opt = nlopt.opt(nlopt.LD_MMA, x.size)
  opt.set_lower_bounds(0.0) ; opt.set_upper_bounds(1.0)
  opt.set_min_objective(wrap_autograd_func(objective_fn, losses, frames))
  opt.add_inequality_constraint(wrap_autograd_func(constraint), 1e-8)
  opt.set_maxeval(args.opt_steps + 1)
  opt.optimize(x.flatten())
  return np.array(losses), reshape(frames[-1]), np.array(frames)

if __name__ == '__main__':
  args = get_args(*mbb_beam())
  losses, x, mbb_frames = fast_stopt(args=args, verbose=True)
  plt.figure(dpi=100) ; print('\nFinal MBB beam design:')
  plt.imshow(np.concatenate([x[:,::-1],x], axis=1)) ; plt.show()
```
