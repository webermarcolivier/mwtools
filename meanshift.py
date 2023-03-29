from sklearn.cluster import MeanShift, estimate_bandwidth
import numpy as np



def clustering_1D_meanshift(x, bandwidth, verbose=0):
    
    if x.ndim != 1:
        raise ValueError("Argument x should be a one-dimensional array.")

    X = np.array(list(zip(x, np.zeros(len(x)))), dtype=np.int)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(X)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)

    centroidList = []
    for k in range(n_clusters_):
        my_members = labels == k
        members = X[my_members, 0]
        if len(members) > 0:
            centroid = members[len(members)//2]
            if verbose >= 2: print("cluster {0}: {1}".format(k, X[my_members, 0]), "centroid:", centroid)
            centroidList.append(centroid)

    centroidList = np.sort(centroidList)
    if verbose >= 2:
        print("len(x):", len(x), "len(centroidList)", len(centroidList))
    return centroidList, labels
