from kmers.neighbors import query_radius
from kmers.pdb_data import PDBData


def calculate_kmers(pdb_data: PDBData, generate_graph: bool = False):
    """
    Takes in a file of experimental data in the format <AA> <x> <y> <z>
    and constructs proximity based k-mers for each, within 15 angstroms.

    Uses a KDTree-based approach
    """
    residues = pdb_data.residue_list
    coordinates = pdb_data.coordinates

    indices, distances = _nearest_neighbours(residues, coordinates)

    # if e.g. a graph is needed, this is the place to do it
    # something like, maybe this works?
    graph = None
    if generate_graph:
        import networkx as nx

        graph = nx.Graph()
        for residue_idx, (ind, dist) in enumerate(zip(indices, distances)):
            for neighbor_idx, distance in zip(ind, dist):
                if residue_idx != neighbor_idx:  # avoid self-connections
                    graph.add_edge(residue_idx, neighbor_idx, weight=distance)

        # nx.write_graphml(graph, "residues.graphml")

    closest_letters = [''.join([residues[i] for i in ind]) for ind in indices]

    if generate_graph:
        return closest_letters, graph

    return closest_letters


def _nearest_neighbours(_residues, coordinates, search_radius=15):
    """
    Returns 2 lists, each containing lists of indices of residues and distances.
    For every residue, the list will be a sorted list of all residues within the search radius.
    default=15 angstroms
    """
    return query_radius(coordinates, search_radius, sort_results=True)
