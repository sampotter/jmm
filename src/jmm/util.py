import numpy as np
import time

t0 = np.inf

def tic():
    global t0
    t0 = time.time()

def toc():
    global t0
    return time.time() - t0

def otsu(values, bins=100):
    '''Compute a threshold for `values` using Otsu's method, where the
    threshold is computed using the number of bins specified by `bins`.

    '''
    bin_counts, edges = np.histogram(values, bins=bins)
    bin_centers = (edges[:-1] + edges[1:])/2

    def ssd(counts, centers):
        n = counts.sum()
        if n == 0:
            mu = 0
        else:
            mu = np.sum(centers*counts)/n
        return np.sum(counts*((centers - mu)**2))

    # Use Otsu's method to compute the optimum threshold
    ssds = []
    for k in range(1, bin_counts.size):
        left_ssd = ssd(bin_counts[:k], bin_centers[:k])
        right_ssd = ssd(bin_counts[k:], bin_centers[k:])
        ssds.append(left_ssd + right_ssd)

    return bin_centers[np.argmin(ssds)]
