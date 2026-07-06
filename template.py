# Structural Topology Optimization Code Template | Student Version
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
  """
  TODO: Task 1 - Implement the SIMP Young's modulus interpolation equation:
  E(x) = E_min + (x^p) * (E_0 - E_min)
  """
  # Please replace the line below with your implementation:
  return None

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
  """
  TODO: Task 2 - Calculate the compliance C of the design structure.
  Specifically, you must:
  1. Construct the grid of nodes (n1, n2, n3, n4) corresponding to each element.
  2. Select the nodal displacements 'u_selected' for each element.
  3. Vectorize the computation of elements compliance 'ce' using numpy's einsum:
     ce_e = u_e^T * ke * u_e
  4. Scale by physical stiffness (SIMP interpolation of young modulus) and sum to return the total compliance.
  """
  nely, nelx = x_phys.shape
  # Create a meshgrid ely, elx for coordinates:
  ely, elx = anp.meshgrid(range(nely), range(nelx))

  # 1. Compute the node indices n1, n2, n3, n4 of the 4 nodes of each element
  # Hint:
  # - Node n1 corresponds to (elx+0, ely+0)
  # - Node n2 corresponds to (elx+1, ely+0)
  # - Node n3 corresponds to (elx+1, ely+1)
  # - Node n4 corresponds to (elx+0, ely+1)
  # Look at the node-numbering scheme in the assignment to map them to 1D node indices.
  # Replace the following dummy equations:
  n1 = None
  n2 = None
  n3 = None
  n4 = None

  # 2. Gather all degrees of freedom (all_ixs) of the elements:
  # Since each node n has DOFs (2*n, 2*n+1), construct the array:
  # all_ixs = anp.array([2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n3, 2*n3+1, 2*n4, 2*n4+1])
  # Replace the dummy equation:
  all_ixs = None

  # 3. Slice displacement vector 'u' to find element displacements 'u_selected':
  # Replace the dummy equation:
  u_selected = None

  # 4. Perform einsum contractions:
  # - Contract 'ke' and 'u_selected' to compute K_0 * u_e for all elements (shape: 8, nely, nelx).
  # - Contract 'u_selected' and the result to compute u_e^T * K_0 * u_e (shape: nely, nelx).
  # Replace the dummy equations:
  ke_u = None
  ce = None

  # 5. Compute Young's modulus at x_phys, multiply by ce, and sum the result:
  # Replace the dummy equation:
  return None

def get_stiffness_matrix(e, nu):
  """
  TODO: Task 3 - Formulate the element stiffness matrix ke of a 2D square quad element.
  Use Young's modulus 'e' and Poisson's ratio 'nu'.
  The stiffness matrix is represented by a symmetric 8x8 matrix based on coefficients:
  k = [1/2-nu/6, 1/8+nu/8, -1/4-nu/12, -1/8+3*nu/8, -1/4+nu/12, -1/8-nu/8, nu/6, 1/8-3*nu/8]
  """
  # Please replace the dummy array below:
  return None

def get_k(stiffness, ke):
  # Constructs sparse stiffness matrix entries and positions
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
  # Displaces the load x using finite element techniques
  stiffness = young_modulus(x_phys, e_0, e_min, p=penal)
  k_entries, k_ylist, k_xlist = get_k(stiffness, ke)
  index_map, keep, indices = _get_dof_indices(freedofs, fixdofs, k_ylist, k_xlist)
  
  u_nonzero = solve_coo(k_entries[keep], indices, forces[freedofs], sym_pos=True)
  u_values = anp.concatenate([u_nonzero, anp.zeros(len(fixdofs))])
  return u_values[index_map]
  
##### Sparse matrix (COO) helper functions #####
def _get_dof_indices(freedofs, fixdofs, k_xlist, k_ylist):
  """
  TODO: Task 4 - Implement the indexing setup to filter the sparse stiffness matrix
  for free degrees of freedom (ignoring the fixed boundary DOFs).
  1. Build an index permutation map that goes from all degrees of freedom to [freedofs, fixdofs].
  2. Create a boolean mask 'keep' selecting entries of 'k_xlist' and 'k_ylist' that belong to 'freedofs'.
  3. Look up mapped coordinate indexes 'i' and 'j' using 'index_map'.
  """
  # Please replace the dummy equations below:
  index_map = None
  keep = None
  i = None
  j = None
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
"""
TODO: Task 5 - Define a custom solver primitive with Autograd.
Since we solve a sparse matrix equation in each step, standard autograd cannot propagate
gradients through scipy.sparse solvers automatically. You must implement the forward pass
and a custom vector-Jacobian product (VJP) for solve_coo.
"""

@autograd.primitive
def solve_coo(a_entries, a_indices, b, sym_pos=False):
  # Forward solve:
  # Use _get_solver (sparse LU solver) to solve Ax = b.
  # Replace the dummy code below:
  return None

def grad_solve_coo_entries(ans, a_entries, a_indices, b, sym_pos=False):
  # Backward pass (VJP w.r.t matrix a_entries):
  # Complete the vector-Jacobian product. The output is a function that maps
  # incoming gradients 'grad_ans' to sensitivities w.r.t 'a_entries'.
  # Under adjoint analysis, the sensitivity w.r.t an entry A_ij is:
  # dC/dA_ij = - lambda_i * u_j, where lambda satisfies A^T * lambda = grad_ans, and u = ans.
  def jvp(grad_ans):
    # 1. Compute adjoint variable lambda_ by solving A^T * lambda = grad_ans.
    # Note: If sym_pos is False, solve with a_indices[::-1] (transposed indices).
    lambda_ = None
    
    # 2. Gather indices:
    i, j = a_indices
    
    # 3. Return the element-wise multiplication: -lambda_[i] * ans[j]
    return None
  return jvp

# Register the VJP w.r.t the first argument (a_entries):
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
  try:
    opt.optimize(x.flatten())
  except Exception as e:
    print(f"Optimization halted. Did you implement the missing sections correctly? Error: {e}")
    return np.array([]), x, np.array([])
  return np.array(losses), reshape(frames[-1]), np.array(frames)

if __name__ == '__main__':
  args = get_args(*mbb_beam())
  losses, x, mbb_frames = fast_stopt(args=args, verbose=True)
  if len(losses) > 0:
    plt.figure(dpi=100) ; print('\nFinal MBB beam design:')
    plt.imshow(np.concatenate([x[:,::-1],x], axis=1)) ; plt.show()
