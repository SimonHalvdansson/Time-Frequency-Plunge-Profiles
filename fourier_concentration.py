import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.optimize import curve_fit, OptimizeWarning
import warnings
from gplearn.genetic import SymbolicRegressor
from gplearn.functions import make_function
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# Define all your existing functions here

def create_indicator_mask(N, start, end):
    """
    Create an indicator mask of length N with ones from start to end (inclusive)
    and zeros elsewhere. Assumes that end >= start (no wrap-around).
    """
    mask = np.zeros(N)
    mask[start:end+1] = 1
    return mask

def construct_fourier_matrix(N):
    """
    Construct the NxN discrete Fourier transform matrix.
    """
    n = np.arange(N)
    k = n.reshape((N, 1))
    omega = np.exp(-2j * np.pi * k * n / N)
    F = omega / np.sqrt(N)
    return F

def fitting_function(k, a, b):
    """
    The function to fit: 0.5 * erfc((k - b) / a)
    """
    return 0.5 * erfc((k - b) / a)

def compute_eigenvalues(N, E_start, E_end, F_start, F_end):
    """
    Compute and return the sorted eigenvalues of the Fourier concentration operator.
    """
    # Create indicator masks
    mask_E = create_indicator_mask(N, E_start, E_end)
    mask_F = create_indicator_mask(N, F_start, F_end)

    # Compute lengths |E| and |F|
    E_length = int(np.sum(mask_E))  # Ensure it's an integer
    F_length = int(np.sum(mask_F))

    # Construct Fourier transform matrix and its conjugate transpose
    F_matrix = construct_fourier_matrix(N)
    F_conj = F_matrix.conj().T

    # Construct the concentration operator A = M_F * F^* * M_E * F * M_F
    # Efficient Construction of A without explicit diagonal matrices

    # Step 1: Compute F * M_F by multiplying each column of F by mask_F
    FMF = F_matrix * mask_F  # Broadcasting mask_F over columns

    # Step 2: Compute M_E * (F * M_F) by multiplying each row by mask_E
    MEMFMF = FMF * mask_E[:, np.newaxis]  # Broadcasting mask_E over rows

    # Step 3: Compute F^* * (M_E * F * M_F)
    F_conj_MEMFMF = F_conj @ MEMFMF

    # Step 4: Multiply by M_F on the left by element-wise multiplication
    A = F_conj_MEMFMF * mask_F[:, np.newaxis]  # Broadcasting mask_F over rows

    # Ensure the operator is Hermitian
    if not np.allclose(A, A.conj().T, atol=1e-10):
        print("Warning: Operator A is not Hermitian.")

    # Compute eigenvalues only using numpy.linalg.eigvalsh
    eigenvalues = np.linalg.eigvalsh(A)

    # Sort eigenvalues in decreasing order
    eigenvalues_sorted = np.sort(eigenvalues)[::-1]

    # Define k indices
    k_indices = np.arange(1, N+1)  # k from 1 to N

    return k_indices, eigenvalues_sorted, E_length, F_length

def fit_eigenvalues(k_indices, eigenvalues_sorted, E_length, F_length):
    """
    Fit the sorted eigenvalues to the function 0.5 * erfc((k - b) / a).
    """
    N = len(k_indices)

    # Initial guesses for a and b
    # a controls the width of the transition; b controls the shift
    # These guesses can be adjusted based on data characteristics
    initial_guess = [N / 10, N / 2]  # Example initial guess

    # Suppress warnings temporarily
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)
        try:
            popt, pcov = curve_fit(
                fitting_function, 
                k_indices, 
                eigenvalues_sorted, 
                p0=initial_guess,
                bounds=([0, 0], [np.inf, N])  # a > 0, 0 <= b <= N
            )
            a_fit, b_fit = popt
            print(f"Fitted parameters:\na = {a_fit}\nb = {b_fit}")

            # Compute additional quantities
            ratio_b = b_fit / (E_length * F_length)*N
            ratio_a = a_fit / (E_length + F_length)*N
            print(f"b / (|E| * |F|) = {ratio_b}")
            print(f"a / (|E| + |F|) = {ratio_a}")

        except RuntimeError as e:
            print("Error in curve fitting:", e)
            return None, None, None

    # Compute the fitted function using the optimized parameters
    fitted_f = fitting_function(k_indices, a_fit, b_fit)
    return a_fit, b_fit, fitted_f

def get_optimal_a_b(N, E_start, E_end, F_start, F_end):
    """
    Compute eigenvalues and fit them to obtain optimal a and b.
    """
    k_indices, eigenvalues_sorted, E_length, F_length = compute_eigenvalues(
        N, E_start, E_end, F_start, F_end
    )
    a_fit, b_fit, fitted_f = fit_eigenvalues(k_indices, eigenvalues_sorted, E_length, F_length)
    return a_fit, b_fit, k_indices, eigenvalues_sorted, fitted_f

def get_random_limits(N):
    """
    Generate random start and end indices for E and F masks ensuring that end > start.
    """
    # Ensure that start < end for E
    E_start = np.random.randint(0, N - 10)  # Avoid very small intervals
    E_end = np.random.randint(E_start + 10, N)  # Ensure at least 10 elements

    # Ensure that start < end for F
    F_start = np.random.randint(0, N - 10)
    F_end = np.random.randint(F_start + 10, N)

    return E_start, E_end, F_start, F_end

def experiment_a_dependency(N, E_sizes, F_sizes):
    """
    Conduct an experiment to determine how the fitted parameter 'a' depends on |E| and |F|.
    """
    a_matrix = np.zeros((len(E_sizes), len(F_sizes)))
    a_matrix[:] = np.nan  # Initialize with NaNs for cases where fitting fails

    for i, E_size in enumerate(E_sizes):
        for j, F_size in enumerate(F_sizes):
            E_start = 0
            E_end = E_start + E_size - 1
            F_start = 0
            F_end = F_start + F_size - 1

            # Validate bounds
            if E_end >= N or F_end >= N:
                print(f"Skipping |E|={E_size} and |F|={F_size} due to index out of bounds.")
                continue

            try:
                # Compute eigenvalues
                k_indices, eigenvalues_sorted, E_length, F_length = compute_eigenvalues(
                    N, E_start, E_end, F_start, F_end
                )
                
                # Fit to obtain 'a'
                a_fit, b_fit, fitted_f = fit_eigenvalues(k_indices, eigenvalues_sorted, E_length, F_length)
                
                if a_fit is not None:
                    a_matrix[i, j] = a_fit
                else:
                    print(f"Fitting failed for |E|={E_size}, |F|={F_size}.")
            except Exception as e:
                print(f"Error for |E|={E_size}, |F|={F_size}: {e}")
                continue

    return a_matrix

def plot_a_dependency(E_sizes, F_sizes, a_matrix):
    """
    Plot the dependency of 'a' on |E| and |F|.
    """
    # Mask NaN values for plotting
    a_matrix_plot = np.copy(a_matrix)
    a_matrix_plot[np.isnan(a_matrix_plot)] = 0  # Replace NaNs with zeros or use masked arrays

    # Plotting a_heatmap
    plt.figure(figsize=(10, 8))
    X, Y = np.meshgrid(F_sizes, E_sizes)
    cp = plt.contourf(X, Y, a_matrix_plot, levels=50, cmap='viridis')
    plt.colorbar(cp, label='Fitted a')
    plt.xlabel('|F|')
    plt.ylabel('|E|')
    plt.title('Dependency of Fitted Parameter a on |E| and |F|')
    plt.show()

    # Plotting a vs |E| for selected |F|
    plt.figure(figsize=(12, 6))
    num_lines = 10  # Number of lines to plot
    F_indices = np.linspace(0, len(F_sizes)-1, num_lines, dtype=int)
    for idx in F_indices:
        plt.plot(E_sizes, a_matrix[:, idx], label=f'|F|={F_sizes[idx]}')
    plt.xlabel('|E|')
    plt.ylabel('Fitted a')
    plt.title('Fitted a vs |E| for Various |F|')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plotting a vs |F| for selected |E|
    plt.figure(figsize=(12, 6))
    E_indices = np.linspace(0, len(E_sizes)-1, num_lines, dtype=int)
    for idx in E_indices:
        plt.plot(F_sizes, a_matrix[idx, :], label=f'|E|={E_sizes[idx]}')
    plt.xlabel('|F|')
    plt.ylabel('Fitted a')
    plt.title('Fitted a vs |F| for Various |E|')
    plt.legend()
    plt.grid(True)
    plt.show()

def perform_symbolic_regression(E_sizes, F_sizes, a_matrix):
    """
    Perform symbolic regression to find a functional relationship between |E|, |F|, and a.
    """
    # Prepare the dataset
    E_grid, F_grid = np.meshgrid(E_sizes, F_sizes, indexing='ij')
    E_flat = E_grid.flatten()
    F_flat = F_grid.flatten()
    a_flat = a_matrix.flatten()

    # Remove entries where a is NaN or zero (if any)
    valid_indices = ~np.isnan(a_flat) & (a_flat > 0)
    E_valid = E_flat[valid_indices]
    F_valid = F_flat[valid_indices]
    a_valid = a_flat[valid_indices]

    # Create feature matrix
    X = np.vstack((E_valid, F_valid)).T
    y = a_valid

    # Split into training and testing datasets
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    # Initialize Symbolic Regressor
    est_gp = SymbolicRegressor(
        population_size=5000,
        generations=20,
        tournament_size=20,
        stopping_criteria=0.01,
        const_range=(0, 1000),
        init_depth=(2, 6),
        function_set=['add', 'sub', 'mul', 'div', 'sqrt', 'log', 'abs', 'neg', 'inv'],
        metric='mse',
        parsimony_coefficient=0.001,
        max_samples=0.9,
        verbose=1,
        n_jobs=-1,
        random_state=42
    )

    # Fit the model
    est_gp.fit(X_train, y_train)

    # Predict on test set
    y_pred = est_gp.predict(X_test)

    # Evaluate the model
    mse = mean_squared_error(y_test, y_pred)
    print(f"Symbolic Regression Mean Squared Error on Test Set: {mse}")

    # Display the best program
    print("Best Program:", est_gp._program)

    # Plotting actual vs predicted
    plt.figure(figsize=(8, 6))
    plt.scatter(y_test, y_pred, alpha=0.5)
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--')
    plt.xlabel('Actual a')
    plt.ylabel('Predicted a')
    plt.title('Actual vs Predicted a (Symbolic Regression)')
    plt.grid(True)
    plt.show()

    return est_gp

def main():
    # Set parameters
    N = 1024*2

    # Define ranges for |E| and |F|
    # Choose step sizes to balance computational load and resolution
    E_sizes = np.arange(10, N//2 + 1, 50)  # From 10 to N/2 in steps of 50
    F_sizes = np.arange(10, N//2 + 1, 50)  # From 10 to N/2 in steps of 50

    print(f"Running experiment with |E| from {E_sizes[0]} to {E_sizes[-1]} and |F| from {F_sizes[0]} to {F_sizes[-1]}.")

    # Conduct the experiment
    a_matrix = experiment_a_dependency(N, E_sizes, F_sizes)

    # Plot the results
    plot_a_dependency(E_sizes, F_sizes, a_matrix)

    # Perform Symbolic Regression
    est_gp = perform_symbolic_regression(E_sizes, F_sizes, a_matrix)

    # Optionally, use the symbolic regression model to predict 'a' for new |E| and |F|
    # Example:
    # new_E = 200
    # new_F = 300
    # predicted_a = est_gp.predict([[new_E, new_F]])
    # print(f"Predicted a for |E|={new_E} and |F|={new_F}: {predicted_a[0]}")

if __name__ == "__main__":
    main()
