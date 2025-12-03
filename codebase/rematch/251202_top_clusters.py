from ase.io import read
import sys
import os
sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules

os.chdir(r"codebase\rematch")

files = [f for f in os.listdir('.') if f.startswith('POSCAR')]
print(files)
clusters = [read(f) for f in files]


#from ase.visualize import view
#view(clusters)

from memory_rematch import get_soap
from get_rematch import soap_rematch_similarity
from codebase.run.evaluate import shave_slab
import timeit

soaps = [get_soap(cluster) for cluster in clusters]
#exec_time = timeit.timeit('soaps = [get_soap(cluster) for cluster in clusters]', globals=globals(), number=10)
#print("Clusters only SOAP Generation:", round(exec_time, 3), "seconds")
#exec_time = timeit.timeit('similarities = [soap_rematch_similarity(soaps[0], soap) for soap in soaps]', globals=globals(), number=10)
#print("Clusters only Similarity:", round(exec_time, 3), "seconds")
similarities = [soap_rematch_similarity(soaps[0], soap) for soap in soaps]

print(similarities)

# shaved slab similarities
clusters_shaved = [shave_slab(cluster, threshold=3)[1] for cluster in clusters]

print(timeit.timeit('get_soap(clusters_shaved)', number=10))
print(timeit.timeit('get_soap(clusters_shaved)', number=10))
soaps = [get_soap(cluster) for cluster in clusters_shaved]
#print([len(cluster) for cluster in clusters_shaved])
soaps = [get_soap(cluster) for cluster in clusters_shaved]
#exec_time = timeit.timeit('soaps = [get_soap(cluster) for cluster in clusters_shaved]', globals=globals(), number=10)
#print("Clusters only SOAP Generation:", round(exec_time, 3), "seconds")
#exec_time = timeit.timeit('similarities = [soap_rematch_similarity(soaps[0], soap) for soap in soaps]', globals=globals(), number=10)
#print("Clusters only Similarity:", round(exec_time, 3), "seconds")
similarities = [soap_rematch_similarity(soaps[0], soap) for soap in soaps]
similarities = [round(sim, 6) for sim in similarities]

[print(i) for i in zip(files, similarities)]

# TODO: Cluster and cluster + surface similarity
# 0.9999+ 4 decimal places  for the shaved slab