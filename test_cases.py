# Unit Tests for Structural Topology Optimization
# Run with: python test_cases.py
import sys
import unittest
import numpy as np
import scipy.sparse
import autograd

# Import from template.py:
try:
  import template
except ImportError:
  print("Error: Could not import template.py. Make sure you are in the correct directory.")
  sys.exit(1)

class TestStructuralOptimization(unittest.TestCase):
    def setUp(self):
        # Setup parameters for a tiny 2x2 design space
        self.width = 2
        self.height = 2
        self.density = 0.5
        
        # Normals (fixed nodes)
        self.normals = np.zeros((self.width + 1, self.height + 1, 2))
        self.normals[0, :, 0] = 1  # fix left wall horizontally
        self.normals[0, 0, 1] = 1  # fix bottom-left node vertically
        
        # Forces (loads)
        self.forces = np.zeros((self.width + 1, self.height + 1, 2))
        self.forces[-1, 0, 1] = -1.0  # apply vertical downward load at bottom-right node
        
        # Args view
        self.args = template.get_args(self.normals, self.forces, self.density)
        self.args.penal = 3.0
        self.args.young_min = 1e-9
        self.args.young = 1.0
        self.args.poisson = 0.3
        
        # Uniform initial layout
        self.x = np.ones((self.args.nely, self.args.nelx)) * self.args.density

    def test_task1_young_modulus(self):
        print("Testing Task 1: young_modulus...")
        e_0 = 1.0
        e_min = 1e-9
        p = 3.0
        
        # Try different values of x
        val1 = template.young_modulus(1.0, e_0, e_min, p)
        self.assertAlmostEqual(val1, 1.0, places=6)
        
        val2 = template.young_modulus(0.0, e_0, e_min, p)
        self.assertAlmostEqual(val2, 1e-9, places=9)
        
        val3 = template.young_modulus(0.5, e_0, e_min, p)
        expected = e_min + (0.5**3.0) * (e_0 - e_min)
        self.assertAlmostEqual(val3, expected, places=7)
        print("Task 1 passed!")

    def test_task2_compliance(self):
        print("Testing Task 2: compliance...")
        # Since compliance requires displace, stiffness matrix and index mapping,
        # we will test if it computes a scalar value correctly when using the template.
        ke = template.get_stiffness_matrix(self.args.young, self.args.poisson)
        if ke is None:
            self.skipTest("Skipping compliance test: Task 3 (get_stiffness_matrix) is not implemented.")
            
        u = template.displace(self.x, ke, self.args.forces, self.args.freedofs, self.args.fixdofs, 
                              penal=self.args.penal, e_min=self.args.young_min, e_0=self.args.young)
        if u is None or np.all(u == 0):
            self.skipTest("Skipping compliance test: displace/solve_coo is not implemented.")
            
        c = template.compliance(self.x, u, ke, penal=self.args.penal, e_min=self.args.young_min, e_0=self.args.young)
        self.assertIsNotNone(c, "Compliance should not be None")
        self.assertTrue(c > 0, f"Compliance should be positive, got {c}")
        print("Task 2 compliance returned a valid positive scalar:", c)

    def test_task3_stiffness_matrix(self):
        print("Testing Task 3: get_stiffness_matrix...")
        e = 1.0
        nu = 0.3
        ke = template.get_stiffness_matrix(e, nu)
        
        self.assertIsNotNone(ke, "Stiffness matrix is None.")
        self.assertEqual(ke.shape, (8, 8), "Stiffness matrix must be 8x8.")
        
        # Symmetry check
        np.testing.assert_allclose(ke, ke.T, rtol=1e-10, atol=1e-10, err_msg="Stiffness matrix must be symmetric.")
        
        # Row sums check (rigid body modes have zero net internal force)
        row_sums = np.sum(ke, axis=1)
        # Note: some rows should sum to zero for translation invariance, but let's check eigenvectors instead
        eigenvalues = np.linalg.eigvalsh(ke)
        # 2D quadrilateral elements have 3 rigid body modes (translation x, y, and rotation).
        # This means the stiffness matrix must have exactly 3 eigenvalues that are zero.
        zero_eigenvalues = np.sum(np.abs(eigenvalues) < 1e-10)
        self.assertEqual(zero_eigenvalues, 3, f"Element stiffness matrix must have exactly 3 rigid body modes (zero eigenvalues), found {zero_eigenvalues}.")
        self.assertTrue(np.all(eigenvalues >= -1e-10), "Stiffness matrix must be positive semi-definite (eigenvalues >= 0).")
        print("Task 3 passed!")

    def test_task4_dof_indices(self):
        print("Testing Task 4: _get_dof_indices...")
        fixdofs = self.args.fixdofs
        freedofs = self.args.freedofs
        
        # Generate dummy k lists
        k_xlist = np.array([0, 1, 2, 3, 4, 5, 0, 1, 2, 3])
        k_ylist = np.array([1, 2, 3, 4, 5, 0, 2, 3, 4, 5])
        
        index_map, keep, indices = template._get_dof_indices(freedofs, fixdofs, k_xlist, k_ylist)
        
        self.assertIsNotNone(index_map, "index_map is None.")
        self.assertEqual(len(index_map), len(freedofs) + len(fixdofs), "index_map must cover all degrees of freedom.")
        
        # Test permutation inverse
        combined = np.concatenate([freedofs, fixdofs])
        perm_check = index_map[combined]
        np.testing.assert_array_equal(perm_check, np.arange(len(combined)), "index_map does not form a correct permutation mapping.")
        
        self.assertIsNotNone(keep, "keep boolean array is None.")
        # entries in keep should represent whether both coordinates are free dofs
        expected_keep = np.isin(k_xlist, freedofs) & np.isin(k_ylist, freedofs)
        np.testing.assert_array_equal(keep, expected_keep, "keep mask does not correctly filter elements belonging to freedofs.")
        
        self.assertIsNotNone(indices, "indices array is None.")
        self.assertEqual(indices.ndim, 2)
        self.assertEqual(indices.shape[0], 2)
        print("Task 4 passed!")

    def test_task5_solve_coo_and_autograd(self):
        print("Testing Task 5: solve_coo and VJP gradient...")
        # Create a simple 3x3 positive definite matrix A in COO format:
        # A = [[4, 1, 0],
        #      [1, 5, 2],
        #      [0, 2, 6]]
        row = np.array([0, 0, 1, 1, 1, 2, 2])
        col = np.array([0, 1, 0, 1, 2, 1, 2])
        data = np.array([4.0, 1.0, 1.0, 5.0, 2.0, 2.0, 6.0])
        b = np.array([1.0, 2.0, 3.0])
        
        # Test forward pass
        x_solve = template.solve_coo(data, np.stack([row, col]), b, sym_pos=True)
        self.assertIsNotNone(x_solve, "solve_coo returned None.")
        
        A_dense = np.array([[4.0, 1.0, 0.0],
                            [1.0, 5.0, 2.0],
                            [0.0, 2.0, 6.0]])
        x_expected = np.linalg.solve(A_dense, b)
        np.testing.assert_allclose(x_solve, x_expected, rtol=1e-7, atol=1e-7, err_msg="solve_coo forward pass is incorrect.")
        
        # Test backward pass (VJP)
        # Define a function: f(data) = sum(solve_coo(data, indices, b))
        indices = np.stack([row, col])
        def wrapper(d):
            return autograd.numpy.sum(template.solve_coo(d, indices, b, sym_pos=True))
            
        try:
            grad_autograd = autograd.grad(wrapper)(data)
        except Exception as e:
            self.fail(f"Autograd raised an error during backward pass: {e}")
            
        # Compute gradient using finite differences:
        eps = 1e-6
        grad_fd = np.zeros_like(data)
        for idx in range(len(data)):
            data_plus = data.copy()
            data_plus[idx] += eps
            val_plus = wrapper(data_plus)
            
            data_minus = data.copy()
            data_minus[idx] -= eps
            val_minus = wrapper(data_minus)
            
            grad_fd[idx] = (val_plus - val_minus) / (2.0 * eps)
            
        np.testing.assert_allclose(grad_autograd, grad_fd, rtol=1e-4, atol=1e-4, 
                                   err_msg="solve_coo backward pass (VJP) is incorrect compared to finite differences.")
        print("Task 5 passed!")

if __name__ == '__main__':
    unittest.main()
