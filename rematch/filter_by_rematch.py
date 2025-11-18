import os
from rematch.get_rematch import soap_rematch_similarity_from_files

YOUR_PATH = "your_path"
SOAP_DIR = f"{YOUR_PATH}"
LOG_PATH = f"{YOUR_PATH}/rematch_log.txt"
SELECTED_LIST_PATH = f"{YOUR_PATH}/selected_structures.txt"

# Parametrisation
THRESHOLD = 0.9
GAMMA = 1.0
ALPHA = 1.0
THRESHOLD_RE = 1e-6

def main():
    all_files = sorted(
        f for f in os.listdir(SOAP_DIR)
        if f.endswith(".csv")
    )

    n_files = len(all_files)
    if n_files == 0:
        print(f"No .csv files found in {SOAP_DIR}")
        return

    print(f"Found {n_files} SOAP files in {SOAP_DIR}")

    active = [True] * n_files

    with open(LOG_PATH, "w") as log_f:
        log_f.write("fileA\tfileB\tREMatch\n")

        for i in range(n_files):
            if not active[i]:
                continue

            fileA = all_files[i]
            pathA = os.path.join(SOAP_DIR, fileA)

            print(f"Reference: {fileA} (index {i})")

            for j in range(i + 1, n_files):
                if not active[j]:
                    continue

                fileB = all_files[j]
                pathB = os.path.join(SOAP_DIR, fileB)

                sim = soap_rematch_similarity_from_files(
                    pathA,
                    pathB,
                    gamma=GAMMA,
                    alpha=ALPHA,
                    threshold=THRESHOLD_RE,
                )

                log_f.write(f"{fileA}\t{fileB}\t{sim:.8f}\n")
                log_f.flush()

                if sim > THRESHOLD:
                    active[j] = False
                    print(f"  Remove {fileB} (index {j}) due to high similarity: {sim:.4f}")

    selected = [fname for flag, fname in zip(active, all_files) if flag]

    with open(SELECTED_LIST_PATH, "w") as f_sel:
        for fname in selected:
            f_sel.write(fname + "\n")

    print("Filtering done.")
    print(f"Original number of structures: {n_files}")
    print(f"Number of selected structures: {len(selected)}")
    print(f"Pairwise similarity log saved to: {LOG_PATH}")
    print(f"Selected structures list saved to: {SELECTED_LIST_PATH}")


if __name__ == "__main__":
    main()
