import math


def query_radius(coordinates, radius, sort_results=True):
    try:
        from sklearn.neighbors import KDTree

        tree = KDTree(coordinates)
        return tree.query_radius(coordinates, r=radius, return_distance=True, sort_results=sort_results)
    except ModuleNotFoundError:
        return _query_radius_bruteforce(coordinates, radius, sort_results=sort_results)


def _query_radius_bruteforce(coordinates, radius, sort_results=True):
    all_indices = []
    all_distances = []
    for i, source in enumerate(coordinates):
        neighbors = []
        for j, target in enumerate(coordinates):
            distance = math.dist(source, target)
            if distance <= radius:
                neighbors.append((j, distance))
        if sort_results:
            neighbors.sort(key=lambda item: item[1])
        all_indices.append([idx for idx, _distance in neighbors])
        all_distances.append([distance for _idx, distance in neighbors])
    return all_indices, all_distances
