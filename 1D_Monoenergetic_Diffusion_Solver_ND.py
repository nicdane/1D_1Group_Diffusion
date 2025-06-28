import numpy as np
import matplotlib.pyplot as plt

# --- Physical and Numerical Constants ---
D = 1.0                 # Diffusion coefficient
SIGMA_A = 0.6           # Absorption cross section (left material)
SIGMA_B = 1.5           # Absorption cross section (right material)
SIGMA_F = 0.1           # Fission cross section (assumed same everywhere)
NU = 2.5                # Average number of neutrons per fission
L = 100                 # Total length of system
N = 50                  # Number of discretization points (higher resolution)
h = L / N               # Mesh spacing (cell width)
material_boundary = int(0.5 * N)  # Middle of domain (where material changes)

# --- Matrix Assembly: Diffusion + Absorption System ---
def create_matrix_A(N, D, SIGMA_A, SIGMA_B, material_boundary, h):
    """Creates the A matrix for the diffusion equation with spatially varying absorption."""
    A = np.zeros((N, N))

    for i in range(N):
        # Use correct absorption value based on spatial location
        sigma = SIGMA_A if i < material_boundary else SIGMA_B

        if i == 0:
            # Reflecting boundary condition on the left
            A[i, i] = 2 * D / h**2 + sigma
            A[i, i + 1] = -2 * D / h**2
        elif i == N - 1:
            # Vacuum boundary condition on the right
            A[i, i] = D / h**2 + sigma
            A[i, i - 1] = -2 * D / h**2
        else:
            # Standard internal cell discretization
            A[i, i - 1] = -D / h**2
            A[i, i]     = 2 * D / h**2 + sigma
            A[i, i + 1] = -D / h**2
    return A

# --- Matrix Assembly: Fission Source Matrix ---
def create_matrix_F(N, SIGMA_F, NU):
    """Creates the diagonal fission matrix F = NU * SIGMA_F * I"""
    return np.diag(np.full(N, NU * SIGMA_F))

# --- Eigenvalue Solver Using Inverse Power Method (Stable) ---
def gauss_seidel_eigenvalue(A, F, tol=1e-5, max_iter=100):
    """Iteratively solves A*phi = (1/k)*F*phi using power iteration."""
    N = A.shape[0]
    phi = np.ones(N)  # Initial guess for flux
    k = 1.0           # Initial guess for eigenvalue

    for _ in range(max_iter):
        # Solve for new phi using linear solve
        phi_new = np.linalg.solve(A, (1 / k) * F @ phi)

        # Normalize flux vector
        phi_new /= np.linalg.norm(phi_new)

        # Update eigenvalue estimate
        k_new = np.dot(phi_new, F @ phi_new) / np.dot(phi_new, A @ phi_new)

        # Check for convergence
        if np.abs(k_new - k) < tol:
            return 1 / k_new, phi_new

        phi = phi_new
        k = k_new

    return 1 / k, phi

# --- Eigenvalue Test Using Direct Matrix Inversion ---
def eigenvalue_test(A, F):
    """Computes eigenvalues/eigenvectors from the generalized eigenvalue problem."""
    F_inv = np.linalg.inv(F)
    M = F_inv @ A
    eigenvalues, eigenvectors = np.linalg.eig(M)
    return eigenvalues, eigenvectors

# --- Visualization: Flux Distribution ---
def plot_flux_distributions(eigenvectors, N):
    """Plots the flux distribution using the first eigenmode."""
    x = np.linspace(0, N - 1, N)
    flux = np.abs(eigenvectors[:, 0])
    flux /= np.max(flux)

    plt.plot(x, flux, label='Normalized Flux')
    plt.xlabel('Cell Index')
    plt.ylabel('Neutron Flux (Relative)')
    plt.title('Neutron Flux Distribution (First Eigenvector)')
    plt.grid(True)
    plt.legend()
    plt.show()

# --- Visualization: Matrix Heatmap ---
def plot_matrix_heatmap(matrix, title='Matrix Heatmap'):
    """Visualizes any matrix as a heatmap (useful for debugging)."""
    plt.figure(figsize=(8, 6))
    plt.imshow(matrix, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.title(title)
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.show()

def plot_symmetric_flux(flux, N):
    """Plots the flux mirrored across the y-axis to visualize symmetry."""
    flux = np.abs(flux)
    flux /= np.max(flux)

    x_positive = np.linspace(0, N - 1, N)
    x_negative = -np.flip(x_positive)
    x_full = np.concatenate((x_negative, x_positive))
    
    flux_mirrored = np.flip(flux)
    flux_full = np.concatenate((flux_mirrored, flux))

    plt.plot(x_full, flux_full, label='Symmetric Flux')
    plt.xlabel('Spatial Coordinate')
    plt.ylabel('Neutron Flux (Normalized)')
    plt.title('Symmetric Neutron Flux Profile')
    plt.grid(True)
    plt.legend()
    plt.show()
    

# === MAIN EXECUTION ===

# Assemble system matrices
A = create_matrix_A(N, D, SIGMA_A, SIGMA_B, material_boundary, h)
F = create_matrix_F(N, SIGMA_F, NU)

# Solve using direct eigenvalue method (matrix inversion)
eigenvalues, eigenvectors = eigenvalue_test(A, F)

# Solve using iterative Gauss-Seidel-like method
k_inv, flux = gauss_seidel_eigenvalue(A, F)

# Output results
print(f"\nConverged eigenvalue (1/k): {k_inv:.4f}")
print(f"Flux distribution:\n{flux}")

plot_symmetric_flux(flux, N)

# Show matrix structure of A
plot_matrix_heatmap(A, title='Heatmap of Matrix A')
