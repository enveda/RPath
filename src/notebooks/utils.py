# -*- coding: utf-8 -*-

"""
    Util functions used in notebook
"""

import os
import json
import logging
import ast
import pandas as pd
from collections import Counter, defaultdict
from tqdm import tqdm
from itertools import product
from typing import Mapping, List

from venn import venn
from networkx import DiGraph, connected_components, has_path
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from drug2ways.alternative_graph_traversal import enumerate_paths
from drug2ways.pathway import _prepare_json
from drug2ways.wrapper import generate_reduced_graph
from drug2ways.rcr import (rcr_all_paths, validate_paths_with_disease_data,
                           disease_rcr_all_paths)
from drug2ways.cli_helper import _setup_logging

# XREF paths
DATA_DIR = '../../data'
KG_DATA_PATH = os.path.join(DATA_DIR, 'kg')

logger = logging.getLogger(__name__)


"""Notebook 1"""


def create_venn_diagram(
    data_dict: dict,
    plot_title: str,
) -> plt:
    venn(data_dict)
    plt.title(plot_title)
    return plt


def get_stats(
    network_1: dict,
    network_2: dict,
    disease_set: dict = None,
    chemical_set: dict = None
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
            index=[1,]
        )
        df = pd.concat([df, tmp_df], ignore_index=True)

    return df


"""Notebook 4"""


def _validation_paths(
    reduced_graph: DiGraph,
    paths: List[List[int]],
    id2node: Mapping[int, str],
) -> [dict, set]:
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
    disease_dict: dict,
    chemical_dict: dict,
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
    dataset: dict,
    graph_df: pd.DataFrame,
) -> dict:
    """Filter dataset based on data in KG"""
    # TODO: Accept graph too.
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
    drug_dict: dict,
    disease_dict: dict,
    clinical_pair_dict: dict,
) -> dict:
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


def get_transcriptomic_paths(
    directed_graph: DiGraph,
    source: str,
    target: str,
    all_paths: List[list],
    drug_dict: dict,
    disease_dict: dict,
    clinical_pair_dict=dict,
) -> [dict, dict]:
    """Validate paths in KG based on the transcriptomic data. """
    _setup_logging(False)

    nodes = set()
    for path in all_paths:
        path_set = set(path)
        nodes.update(path_set)

    # Getting concordant paths with the transcriptomic data
    drug_filtered_paths = rcr_all_paths(
        directed_graph,
        all_paths,
        drug_dict,
        errors_allowed=ERRORS_ALLOWED
    )
    disease_filtered_paths = disease_rcr_all_paths(
        directed_graph,
        all_paths,
        disease_dict,
        errors_allowed=ERRORS_ALLOWED
    )

    # Check if pair in clinical trials
    pair = source + '_' + target
    if pair in clinical_pair_dict:
        trial_pair = True
    else:
        trial_pair = False

    results = {
        "source": source,
        "target": target,
        "number_of_paths": len(all_paths),
        "in_clinical_trial": trial_pair,
        "drug_paths": drug_filtered_paths,
        "disease_paths": disease_filtered_paths,
        "subgraph_size": directed_graph.number_of_nodes(),
        "number_of_unique_nodes": len(nodes),
    }

    return results


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


"""Notebook 6"""


def get_sankey_diagram(
    path_tsv_dir: str,
    source_node: str,
    target_node: str,
    graph_type: str
):
    """ Get sankey graph for source and target. graph_type is either 'custom'
    or 'openbio' """
    paths_df = pd.read_csv(path_tsv_dir, sep='\t')
    query = '(source == "{0}") & (target == "{1}")'.format(source_node, target_node)
    paths_query = paths_df.query(query)
    index = paths_query.index[paths_query['subgraph_name'] == graph_type][0]
    path, polarity = [], []

    ini_list = paths_df.loc[index]['paths']
    ini_pol = paths_df.loc[index]['signs']

    # transform string list into actual list
    if type(ini_list) == str:
        res = ast.literal_eval(ini_list)
        res_pol = ast.literal_eval(ini_pol)
        for n in res:
            path.append(n)
        for n in res_pol:
            polarity.append(n)

    # get set of nodes to use for indexing
    nodes = set()
    for p in path:
        for node in p:
            nodes.add(node)
    nodes = list(nodes)
    # get source target value index corresponding to path connections
    source, target, value = [], [], []
    for paths in path:
        for n in range(len(paths)):
            if n != (len(paths) - 1):
                source.append(nodes.index(paths[n]))
                target.append(nodes.index(paths[n + 1]))
                value.append(1)

    # get colors of nodes
    colors = []
    for n in nodes:
        if n.startswith('mondo'):
            colors.append('blue')
        elif n.startswith('ncbigene'):
            colors.append('green')
        else:
            colors.append('orange')

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=10,
            line=dict(color="black", width=0.3),
            label=nodes,
            color=colors
        ),
        link=dict(
            source=source,
            target=target,
            value=value
        ))])

    fig.update_layout(title_text="Drug Disease Paths between {0} and {1}".format(source_node, target_node),
                      font_size=10)
    fig.show()
