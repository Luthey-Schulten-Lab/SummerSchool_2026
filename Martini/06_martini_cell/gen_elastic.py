import json
import sys

import networkx as nx
import numpy as np
import scipy
from networkx.readwrite import json_graph
from tqdm import tqdm
from vermouth import Molecule
from vermouth.forcefield import ForceField
from vermouth.gmx.gro import read_gro
from vermouth.processors.apply_rubber_band import compute_force_constants


def selector(attributes):
    """
    Select backbone SC1 of DNA.
    """
    if attributes["atomname"] in ["BB1", "BB2", "BB3", "SC1"]:
        return True
    return False


def compute_decay(distance, shift, rate, power):
    r"""
    Compute the decay function of the force constant as function to the distance.

    The decay function for the force constant is defined as:

    .. math::

        \exp^{-r(d - s)^p}

    where :math:`r` is the decay rate given by the 'rate' argument,
    :math:`p` is the decay power given by 'power', :math:`s` is a shift
    given by 'shift', and :math:`d` is the distance between the two atoms given
    in 'distance'. If the rate or the power are set to 0, then the decay
    function does not modify the force constant.

    The 'distance' argument can be a scalar or a numpy array. If it is an
    array, then the returned value is an array of decay factors with the same
    shape as the input.
    """
    return np.exp(-rate * np.power((distance - shift), power))


def compute_force_constants(
    distances,
    lower_bound,
    upper_bound,
    decay_factor,
    decay_power,
    base_constant,
    minimum_force,
):
    """
    Compute the force constant of an elastic network bond.

    The force constant can be modified with a decay function, and it can be
    bounded with a minimum threshold, or a distance upper and lower bonds.
    """
    constants = compute_decay(distances, lower_bound, decay_factor, decay_power)
    constants *= base_constant
    constants[constants < minimum_force] = 0
    constants[distances > upper_bound] = 0
    return constants


def elastic_network(coordinates, **kwargs):
    max_dist = 0
    coord_tree = scipy.spatial.ckdtree.cKDTree(coordinates)
    for idx in tqdm(range(0, len(coordinates))):
        distances, pairs = coord_tree.query(
            coordinates[idx],
            k=50,
            distance_upper_bound=kwargs["upper_bound"],
        )
        mask = distances < np.inf
        distances = distances[mask]
        pairs = pairs[mask]
        forces = compute_force_constants(distances, **kwargs)

        for dist, pair, force in zip(distances, pairs, forces):
            if dist > max_dist:
                max_dist = dist

            if force > 0:
                yield (idx, pair, dist, force)


def vermouth_mol_from_json(filename):
    with open(filename) as file_:
        data = json.load(file_)
    molecule_graph = nx.Graph(json_graph.node_link_graph(data))
    ff = ForceField("test")
    mol = Molecule(force_field=ff, nrexcl=1)
    mol.add_nodes_from(molecule_graph.nodes(data=True))
    mol.add_edges_from(molecule_graph.edges)
    return mol


def __main__():

    kwargs = {
        "lower_bound": 0,
        "upper_bound": 1.2,
        "decay_factor": 0,
        "decay_power": 1,
        "base_constant": 500,
        "minimum_force": 13,
    }

    gro = read_gro(sys.argv[1])
    idx_to_node = {}
    missing = []
    coordinates = []
    selection = []
    node_idx = 0
    for node_key, attributes in gro.nodes.items():
        if selector(attributes):
            idx_to_node[node_idx] = node_key
            selection.append(node_idx)
            coordinates.append(attributes.get("position"))
            if coordinates[-1] is None:
                missing.append(node_key)
            node_idx += 1
    if missing:
        raise ValueError(
            "All atoms from the selection must have coordinates. "
            "The following atoms do not have some: {}.".format(" ".join(missing))
        )

    coordinates = np.stack(coordinates)

    template = "{0:7d} {1:7d} 1 {2:3.5F} {3:3.5F}\n"
    with open("elastic_network.itp", "w") as filehandle:
        filehandle.write("[ bonds ]\n")
        filehandle.write("; rubber bands for DNA genome\n")
        filehandle.write("; equivalent to the soft-elastic option in martinize2\n")
        had_pairs = {}
        for idx, jdx, dist, forcek in elastic_network(coordinates, **kwargs):
            if idx != jdx and had_pairs.get(frozenset([idx, jdx]), True):
                filehandle.write(
                    template.format(
                        idx_to_node[idx] + 1, idx_to_node[jdx] + 1, dist, forcek
                    )
                )
                had_pairs[frozenset([idx, jdx])] = False


__main__()
