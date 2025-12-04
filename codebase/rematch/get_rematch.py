import numpy as np
from dscribe.kernels import REMatchKernel
from sklearn.preprocessing import normalize

_kernel_cache = {}

def _get_kernel(gamma=1.0, alpha=1.0, threshold=1e-6):
    """
    Reuse REMatch kernels for identical hyperparameters since constructing them
    repeatedly is relatively expensive.
    """
    key = (gamma, alpha, threshold)
    kernel = _kernel_cache.get(key)
    if kernel is None:
        kernel = REMatchKernel(metric="rbf", gamma=gamma, alpha=alpha, threshold=threshold)
        _kernel_cache[key] = kernel
    return kernel

def get_cached_rematch_kernel(gamma=1.0, alpha=1.0, threshold=1e-6):
    """
    Public helper to fetch a cached REMatch kernel for the given hyperparameters.
    This avoids constructing a new REMatchKernel every time a comparison is made.
    """
    return _get_kernel(gamma=gamma, alpha=alpha, threshold=threshold)

def load_soap_from_csv(csv_path):
    """
    read SOAP (.csv), return np
    """
    return np.loadtxt(csv_path, delimiter=',')

def soap_rematch_similarity(soap1, soap2, gamma=1.0, alpha=1.0, threshold=1e-6, kernel=None):
    """
    input two SOAP, return REMatch
    """
    soap1_norm = normalize(soap1)
    soap2_norm = normalize(soap2)

    if kernel is None:
        kernel = _get_kernel(gamma=gamma, alpha=alpha, threshold=threshold)
    kernel_matrix = kernel.create([soap1_norm, soap2_norm])

    return kernel_matrix[0, 1]

def soap_rematch_kernel_matrix(soap_descriptors, gamma=1.0, alpha=1.0, threshold=1e-6, kernel=None):
    """
    Compute the REMatch Gram matrix for a list of SOAP descriptors in one shot.
    """
    normalized = [normalize(soap) for soap in soap_descriptors]
    if kernel is None:
        kernel = _get_kernel(gamma=gamma, alpha=alpha, threshold=threshold)
    return kernel.create(normalized)

def soap_rematch_similarity_from_files(csv_path1, csv_path2, gamma=1.0, alpha=1.0, threshold=1e-6):
    soap1 = load_soap_from_csv(csv_path1)
    soap2 = load_soap_from_csv(csv_path2)
    return soap_rematch_similarity(soap1, soap2, gamma=gamma, alpha=alpha, threshold=threshold)


if __name__ == "__main__":

    file1 = "path_to_soap1.csv"
    file2 = "path_to_soap2.csv"

    sim = soap_rematch_similarity_from_files(file1, file2, gamma=1.0, alpha=1.0)
    print("SOAP REMatch similarity:", sim)
