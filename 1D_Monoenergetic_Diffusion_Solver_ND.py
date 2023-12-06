import numpy as np
import matplotlib.pyplot as plt

# Constants
D = 1.0        # Diffusion coefficient
SIGMA_A = 0.6  # Absorption cross section for left material
SIGMA_F = 0.1  # Fission cross section
SIGMA_B = 1.5  # Absorption cross section for the right material
material_boundary = 5  # Index where the material changes
NU = 2.50      # Average number of neutrons produced per fission
L = 100        # total length of cells
N = 9         # Number of discretization points (basically resolution)
h = L / N      # Mesh spacing
S = 1.0        #Source strength

# Function to create A matrix (physical system)
def create_matrix_A(N, D, SIGMA_A, SIGMA_B, material_boundary, h):
    A = np.zeros((N, N))

    for i in range(N):
        # Determine the absorption cross section for the current point
        sigma = SIGMA_A * S if i < material_boundary else SIGMA_B

        if i == 0:
            # Reflecting boundary at the left
            A[i, i] = 2 * D / h**2 + sigma
            A[i, i + 1] = -2 * D / h**2
        elif i == N - 1:
            # Vacuum boundary at the right
            A[i, i] = -D / h** + sigma #not sure if this sigma is right but it helped fix some vacuum boundary matrix issues
            A[i, i - 1] = 2 * D / h**2 
        else:
            # Internal cells
            A[i, i - 1] = -D / h**2
            A[i, i] = 2 * D / h**2 + sigma
            A[i, i + 1] = -D / h**2

    return A


# Function to create the F matrix using fission cross section vals
def create_matrix_F(N, SIGMA_F, NU):
    F = np.diag(np.full(N, NU * SIGMA_F))
    return F

#iterative solver to converge flux + k values -- modified from HW2 to fit within context of problem -- still a bit broken, will take some time to update
#returns 1/k to be converted
def gauss_seidel_eigenvalue(A, F, tol=1e-4, max_iter=100):
    # Initial guess for the flux vector phi and eigenvalue lambda (1/k)
    phi = np.ones(N)
    lambda_old = 0.5

    for _ in range(max_iter):
        phi_old = phi.copy()

        # Gauss-Seidel iteration for A*phi = lambda*F*phi
        for i in range(N):
            sigma = np.dot(A[i, :i], phi[:i]) + np.dot(A[i, i+1:], phi_old[i+1:])
            phi[i] = (lambda_old * np.dot(F[i], phi_old)) / A[i, i]

        # Normalize phi
        phi /= np.linalg.norm(phi)

        # Update lambda (1/k)
        lambda_new = np.dot(phi, np.dot(F, phi)) / np.dot(phi, np.dot(A, phi))

        # Convergence check
        if np.linalg.norm(phi - phi_old, ord=np.inf) < tol:
            return 1/lambda_new, phi

        lambda_old = lambda_new

    return 1/lambda_new, phi

# Using inverted matrix method to quickly see flux vals This DOES NOT use Gauss-Seidel method
#NB: Some plots resulting from this subroutine were used during presentation. Kept subroutine still intact as it's still 'most' accurate
def eigenvalue_test(A, F):

    phi = np.ones(N)
    F_inv = np.linalg.inv(F)
    eigenvalue_matrix = F_inv.dot(A)

    # calculate eigenvalues
    eigenvalues, eigenvectors = np.linalg.eig(eigenvalue_matrix)

    return eigenvalues, eigenvectors

# helpful for quickly visualizing 1D system without having to closely inspect numbers
def plot_matrix_heatmap(matrix, title='Matrix Heatmap'):

    plt.figure(figsize=(8, 6))
    plt.imshow(matrix, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.title(title)
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.show()

#again this is an older subroutine but my gauss seidel method appears to be performing worse as opposed to inverse method
def plot_flux_distributions(eigenvectors, N):

    x = range(N)  # Discretization points

    x_positive = np.linspace(0, N-1, N)  # Original discretization points
    x_negative = -np.flip(x_positive) - 0.5  # Mirrored discretization points
    x_full = np.concatenate((x_negative, x_positive))

    #loop is redundant now...modified so it wouldnt overplot several eigenvalues at once, index 0 is only one that matters
    for i in range(eigenvectors.shape[1]):
        flux_full = np.concatenate((np.flip(eigenvectors[:, 0]), eigenvectors[:, 0]))

        plt.plot(x_full, flux_full) #plots what you first see when running code -- supposed to emulate flux distribution w/ reflecting boundary

    plt.xlabel('Discretization Points of System')
    plt.ylabel('Neutron Flux - Relative')
    plt.title('Neutron Flux Distribution')
    plt.legend()
    plt.show()

#MAIN ROUTINE

# Creating matrices A and F
A = create_matrix_A(N, D, SIGMA_A, SIGMA_B, material_boundary, h)
F = create_matrix_F(N, SIGMA_F, NU)

# Solving for eigenvalue
eigenvalues, eigenvectors = eigenvalue_test(A, F)

solution = gauss_seidel_eigenvalue(A, F)
print("\n ",solution, "\n") #first output here -- these are directly correlating to A matrix -- something went wrong

# Printing the results for matrices A and F
print(A)
print(F)

# After solving for eigenvalues and eigenvectors
plot_flux_distributions(eigenvectors, N)
print('\n \n')

# Main execution
A = create_matrix_A(N, D, SIGMA_A, SIGMA_B, material_boundary, h)
F = create_matrix_F(N, SIGMA_F, NU)

k, flux = gauss_seidel_eigenvalue(A, F)

print(f"Converged eigenvalue (1/k): {k}")
print(f"Flux distribution: {flux}")
print('\n \n \n')

#This will pop up after you close the flux plot. Bottom cell second from right vacuum appears to still be faulty in A matrix creation
plot_matrix_heatmap(A, title='Heatmap of Matrix A')

#EXTRA STUFF
# Warning >> Power Method is very broken. You will get a 0 value K...NOT GOOD! Terrible neutron economy >__<
"""def power_method(A, F, tol=1e-4, max_iter=100):
    N = len(A)  
    phi = np.ones(N)  # Initial guess for the flux
    k_old = 1.0  # Initial guess for the eigenvalue -- shouldn't matter but 1.0 is standard

    for _ in range(max_iter):
        # Calculate F * phi
        F_phi = np.dot(F, phi)

        # Solve A * new_phi = F * phi
        new_phi = gauss_seidel(A, F_phi)

        # Normalize new_phi and prevent overflow
        norm = np.linalg.norm(new_phi)
        if norm == 0 or np.isinf(norm) or np.isnan(norm):
            break  # Break if norm is zero, infinity, or nan
        new_phi /= norm

        # some checks for weird cases (infinity or nan values to break method -- for some situations the code just won't exit
        k_new = k_old * norm
        if np.isinf(k_new) or np.isnan(k_new):
            break  # Break if k_new is infinity or NaN

        # Check for convergence
        if np.abs(k_new - k_old) / np.abs(k_old) < tol:
            return 1 / k_new, new_phi

        phi = new_phi
        k_old = k_new

    return 1 / k_new, phi

# Main execution using Power Method
A = create_matrix_A(N, D, SIGMA_A, SIGMA_B, material_boundary, h)
F = create_matrix_F(N, SIGMA_F, NU)

k, flux = power_method(A, F)

print(f"Converged eigenvalue (1/k): {k}")
print(f"Flux distribution: {flux}")
plt.plot(flux)"""

k, flux = power_method(A, F)

print(f"Converged eigenvalue (1/k): {k}")
print(f"Flux distribution: {flux}")
plt.plot(flux)"""
