# -*- coding: utf-8 -*-

"""Utils to be used in notebooks"""

import logging
import os
from collections import Counter, defaultdict
from itertools import product
from typing import Any, Mapping, List, Dict

import matplotlib.pyplot as plt
import pandas as pd
from drug2ways.alternative_graph_traversal import enumerate_paths
from drug2ways.cli_helper import _setup_logging
from drug2ways.pathway import _prepare_json
from drug2ways.rcr import (
    rcr_all_paths,
    validate_paths_with_disease_data,
    pairwise,
)
from drug2ways.wrapper import generate_reduced_graph
from networkx import DiGraph, connected_components
from tqdm import tqdm
from venn import venn

# XREF paths
DATA_DIR = '../data'
KG_DATA_PATH = os.path.join(DATA_DIR, 'kg')

logger = logging.getLogger(__name__)

"""Notebook 1"""


def create_venn_diagram(
    data_dict: Dict,
    plot_title: str,
) -> plt:
    venn(data_dict)
    plt.title(plot_title)
    return plt


def get_stats(
    network_1: Dict,
    network_2: Dict,
    disease_set: Dict = None,
    chemical_set: Dict = None
) -> None:
    """Get KG overlap information."""

    diseases_1 = network_1['disease']
    diseases_2 = network_2['disease']

    chemicals_1 = network_1['chemical']
    chemicals_2 = network_2['chemical']

    if disease_set:
        logger.info(f'No.of elements are: {len(disease_set)}')

        disease_1_intersection = diseases_1.intersection(disease_set)
        logger.info(
            f'No.of elements overlaping with OpenBio network are: '
            f'{len(disease_1_intersection)}'
        )

        disease_2_intersection = diseases_2.intersection(disease_set)
        logger.info(
            f'No.of elements overlaping with Custom network are: '
            f'{len(disease_2_intersection)}'
        )

    if chemical_set:
        chem_1_intersection = chemicals_1.intersection(chemical_set)
        logger.info(
            f'No.of elements overlaping with OpenBio network are: '
            f'{len(chem_1_intersection)}'
        )

        chem_2_intersection = chemicals_2.intersection(chemical_set)
        logger.info(
            f'No.of elements overlaping with Custom network are: '
            f'{len(chem_2_intersection)}'
        )


"""Notebook 2"""


def create_graph_from_df(
    graph_df
) -> DiGraph:
    """Create fully connected graph from dataframe."""
    graph = DiGraph()

    for sub_name, obj_name, relation in graph_df.values:
        # Store edge in the graph
        graph.add_edge(
            sub_name,
            obj_name,
            polarity=relation,
        )

    logger.warning(
        f"Report on the number of relations: "
        f"{dict(Counter(graph_df.polarity))}"
    )

    connected_components_subgraph = [
        component
        for component in sorted(
            connected_components(
                graph.to_undirected()
            ),
            key=len,
            reverse=True
        )
    ]

    final_subgraph = graph.subgraph(connected_components_subgraph[0])

    return final_subgraph


"""Notebook 3"""


def create_df_from_graph(
    graph: DiGraph
) -> pd.DataFrame:
    """Create dataframe from graph."""

    df = pd.DataFrame(columns=['source', 'target', 'polarity'])

    for source, target, relation in tqdm(graph.edges(data=True), desc='Processing data'):
        tmp_df = pd.DataFrame(
            {
                'source': source,
                'target': target,
                'polarity': relation['polarity']
            },
            index=[1, ]
        )
        df = pd.concat([df, tmp_df], ignore_index=True)

    return df


"""Notebook 4"""


def get_path_count(
    directed_graph: DiGraph,
    filtered_paths: list
) -> (int, int):
    """Count number of activatory and inhibitory paths. """
    activatory_paths = 0
    for path in filtered_paths:
        relation_sign = 1
        for node_id in range(len(path) - 1):
            s = path[node_id]
            t = path[node_id + 1]
            relation_sign *= directed_graph[s][t]['polarity']

        if relation_sign == 1:
            activatory_paths += 1

    inhibitory_paths = len(filtered_paths) - activatory_paths

    return activatory_paths, inhibitory_paths


def _validation_paths(
    reduced_graph: DiGraph,
    paths: List[List[int]],
    id2node: Mapping[int, str],
) -> [Dict, set]:
    """Get paths between nodes."""
    results = defaultdict(Counter)

    # Iter over all paths
    for path in paths:
        # Iterate through each node tracking its position in the path
        for index, node in enumerate(path):
            results[index][id2node[node]] += 1

    polarity_dict = {1: '->', -1: '-|'}

    final_paths = set()
    for path in paths:

        reconstructed_path = []

        for index, node in enumerate(path):

            # Avoid crashing
            if index + 1 == len(path):
                continue

            polarity = polarity_dict[
                reduced_graph[node][path[index + 1]]['polarity']
            ]

            if index == 0:
                reconstructed_path.append(node)
                reconstructed_path.append(polarity)
                reconstructed_path.append(path[index + 1])
            else:
                reconstructed_path.append(polarity)
                reconstructed_path.append(path[index + 1])

        final_paths.add(
            tuple(
                id2node[cosa] if cosa in id2node else cosa
                for cosa in reconstructed_path
            )
        )

    final_paths = _prepare_json(final_paths)

    return final_paths


def get_paths(
    graph: DiGraph,
    disease_dict: Dict,
    chemical_dict: Dict,
) -> dict:
    """Get paths in graph."""
    detailed_info = {}

    chemical_els = set(
        node
        for node in graph.nodes()
        if 'pubchem' in node
    )

    disease_els = set(
        node
        for node in graph.nodes()
        if 'mondo' in node
    )

    # Get the reduced version of the graph and the node2id mapping
    target_nodes = list(disease_dict.keys())
    reduced_graph, node2id = generate_reduced_graph(graph, target_nodes)

    id2node = {
        v: k
        for k, v in node2id.items()
    }

    # Iterating different chemical-disease pair
    for chemical, disease in tqdm(
        product(chemical_dict, disease_dict),
        total=len(chemical_dict) * len(disease_dict),
        desc='Getting paths'
    ):
        if chemical in chemical_els and disease in disease_els:
            for lmax in range(3, 8):
                dict_key = f'lmax_{lmax}'

                if dict_key not in detailed_info:
                    detailed_info[dict_key] = []

                paths = enumerate_paths(
                    graph=reduced_graph,
                    source=node2id[chemical],
                    targets=[node2id[disease]],
                    lmax=lmax,
                    cycle_free=True,
                )

                # Skip if there are no paths
                if not paths:
                    detailed_info[dict_key].append(
                        {
                            'source': chemical,
                            'target': disease,
                            'paths': []
                        }
                    )

                else:
                    # Get summarized results for export
                    paths_summary = _validation_paths(
                        reduced_graph=reduced_graph,
                        paths=paths,
                        id2node=id2node,
                    )

                    detailed_info[dict_key].append(
                        {
                            'source': chemical,
                            'target': disease,
                            'paths': paths_summary
                        }
                    )
    return detailed_info


def filter_dataset(
    dataset: Dict,
    graph_df: pd.DataFrame,
) -> dict:
    """Filter dataset based on data in KG"""
    chemical_nodes = set(
        graph_df[
            graph_df['source'].str.contains('pubchem.compound')
        ]['source']
    )
    disease_nodes = set(
        graph_df[
            graph_df['target'].str.contains('mondo')
        ]['target']
    )
    gene_1 = set(
        graph_df[
            graph_df['source'].str.contains('ncbigene')
        ]['source']
    )
    gene_2 = set(
        graph_df[
            graph_df['target'].str.contains('ncbigene')
        ]['target']
    )
    gene_nodes = gene_1.union(gene_2)

    new_dict = {}

    node_list = set()

    for i in dataset:
        if i in chemical_nodes or i in disease_nodes:
            gene_dict = {}
            node_list.add(i)

            for gene in dataset[i]:
                if gene in gene_nodes:
                    gene_dict[gene] = dataset[i][gene]

            new_dict[i] = gene_dict

    # Assert
    for node in node_list:
        assert node in chemical_nodes or node in disease_nodes

    logger.info(f'Reduced dataset from {len(dataset)} to {len(new_dict)}')

    return new_dict


ERRORS_ALLOWED = 1


def get_validated_paths(
    directed_graph: DiGraph,
    source: str,
    target: str,
    all_paths: List[list],
    drug_dict: Dict,
    disease_dict: Dict,
    clinical_pair_dict: Dict,
) -> Dict:
    """Validate paths in KG"""
    _setup_logging(False)

    nodes = set()
    for path in all_paths:
        path_set = set(path)
        nodes.update(path_set)

    # Getting concordant paths with the drug_dict
    filtered_paths = rcr_all_paths(
        directed_graph,
        all_paths,
        drug_dict,
        errors_allowed=ERRORS_ALLOWED
    )

    # Add genes if missing in disease data
    for path in filtered_paths:
        for gene in path[1:]:
            if gene not in disease_dict:
                disease_dict[gene] = 0

    # Check the values for each node in both dictionaries are the opposite
    validate_paths = validate_paths_with_disease_data(
        paths=filtered_paths,
        drug_dict=drug_dict,
        disease_dict=disease_dict,
        errors_allowed=ERRORS_ALLOWED
    )

    # Check if pair in clinical trials dicts
    pair = source + '_' + target

    if pair in clinical_pair_dict:
        trial_pair = True
    else:
        trial_pair = False

    # Count number of activatory and inhibitory paths
    activatory_paths = 0
    for path in validate_paths:
        relation_sign = 1
        for node_id in range(len(path) - 1):
            s = path[node_id]
            t = path[node_id + 1]
            relation_sign *= directed_graph[s][t]['polarity']

        if relation_sign == 1:
            activatory_paths += 1

    inhibitory_paths = len(validate_paths) - activatory_paths

    results = {
        "source": source,
        "target": target,
        "number_of_paths": len(all_paths),
        "number_of_concordant_paths": len(validate_paths),
        "in_clinical_trial": trial_pair,
        "number_of_concordant_activatory_paths": activatory_paths,
        "number_of_concordant_inhibitory_paths": inhibitory_paths,
        "subgraph_size": directed_graph.number_of_nodes(),
        "number_of_unique_nodes": len(nodes),
    }
    return results


"""Protein prioritization"""


def _is_not_concordant(
    graph: DiGraph,
    path: List[str],
    disease_dict: Dict[str, int],
    errors_allowed: int = 0
) -> bool:
    """Calculate if the path is concordant.

    :param graph: original directed graph
    :param path: path to evaluate
    :param disease_dict: dictionary with the fold changes from the disease experiment
    :param errors_allowed: errors allowed in the path
    :return: boolean with the result
    """
    # Calculate the current score
    current_polarity = disease_dict[path[0]]
    # number of errors during evaluation
    current_errors = 0

    for source, target in pairwise(path):

        # Update polarity
        current_polarity = current_polarity * graph.edges[source, target]['polarity']

        target_score = disease_dict[target]

        if current_polarity == target_score:
            # max errors allowed reached
            if current_errors == errors_allowed:
                return False
            # allow for one more error
            current_errors += 1

    return True


def get_protein_paths(
    graph: DiGraph,
    protein_list: List[str],
    target_nodes: List[str],
    lmax: int
) -> Dict:
    """Get paths in graph from protein to diseases."""
    detailed_info = {}

    # Get the reduced version of the graph and the node2id mapping
    reduced_graph, node2id = generate_reduced_graph(graph, target_nodes)

    id2node = {
        v: k
        for k, v in node2id.items()
    }

    # Iterating different chemical-disease pair
    for protein in tqdm(
        protein_list,
        desc='Getting paths'
    ):
        for disease in target_nodes:
            dict_key = f'lmax_{lmax}'

            if dict_key not in detailed_info:
                detailed_info[dict_key] = []

            paths = enumerate_paths(
                graph=reduced_graph,
                source=node2id[protein],
                targets=[node2id[disease]],
                lmax=lmax,
                cycle_free=True,
            )

            # Skip if there are no paths
            if not paths:
                detailed_info[dict_key].append(
                    {
                        'source': protein,
                        'target': disease,
                        'paths': []
                    }
                )

            else:
                # Get summarized results for export
                paths_summary = _validation_paths(
                    reduced_graph=reduced_graph,
                    paths=paths,
                    id2node=id2node,
                )

                detailed_info[dict_key].append(
                    {
                        'source': protein,
                        'target': disease,
                        'paths': paths_summary
                    }
                )
    return detailed_info


def is_a_concordant_protein_target(
    graph: DiGraph,
    all_paths: List[List],
    disease_dict: Dict[str, int],
    errors_allowed: int = 0
) -> List[Any]:
    """Conduct causal reasoning on the paths between drug and disease based on disease experimental data.

    :param graph: original directed graph
    :param all_paths: all paths to evaluate
    :param disease_dict: dictionary with the fold changes from the disease experiment
    :param errors_allowed: errors allowed in the path
    """
    valid_paths = []

    for path in all_paths:

        # remove the last node since it is the disease and it doesnt have any experimental value
        path = path[:-1]

        # Skip path if not all nodes are present in experimental data
        # First node is not considered since it corresponds to the drug which doesnt have experimental value
        if not all(node in disease_dict for node in path):
            continue

        if not _is_not_concordant(graph, path, disease_dict, errors_allowed):
            continue

        valid_paths.append(path)

    return valid_paths


def discover_target(
    directed_graph: DiGraph,
    source: str,
    target: str,
    all_paths: List[List[str]],
    disease_dict: Dict,
) -> Dict:
    """Check if targets have concordant paths with a given disease. """
    _setup_logging(False)

    # Getting concordant paths with the transcriptomic data for the disease
    disease_filtered_paths = is_a_concordant_protein_target(
        directed_graph,
        all_paths,
        disease_dict,
    )

    nodes = set()
    for path in disease_filtered_paths:
        path_set = set(path)
        nodes.update(path_set)

    return {
        "source": source,
        "target": target,
        "number_of_paths": len(all_paths),
        "concordant_paths": disease_filtered_paths,
        "nodes_in_concordant_paths": len(nodes),
    }
