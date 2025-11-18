import numpy as np
from dscribe.kernels import REMatchKernel
from sklearn.preprocessing import normalize

def load_soap_from_csv(csv_path):
    """
    read SOAP (.csv), return np
    """
    return np.loadtxt(csv_path, delimiter=',')

def soap_rematch_similarity(soap1, soap2, gamma=1.0, alpha=1.0, threshold=1e-6):
    """
    input two SOAP, return REMatch
    """
    soap1_norm = normalize(soap1)
    soap2_norm = normalize(soap2)

    re = REMatchKernel(metric="rbf",gamma=gamma,alpha=alpha,threshold=threshold)

    kernel_matrix = re.create([soap1_norm, soap2_norm])

    return kernel_matrix[0, 1]

def soap_rematch_similarity_from_files(csv_path1, csv_path2, gamma=1.0, alpha=1.0, threshold=1e-6):
    soap1 = load_soap_from_csv(csv_path1)
    soap2 = load_soap_from_csv(csv_path2)
    return soap_rematch_similarity(soap1, soap2, gamma=gamma, alpha=alpha, threshold=threshold)


if __name__ == "__main__":

    file1 = "path_to_soap1.csv"
    file2 = "path_to_soap2.csv"

    sim = soap_rematch_similarity_from_files(file1, file2, gamma=1.0, alpha=1.0)
    print("SOAP REMatch similarity:", sim)
