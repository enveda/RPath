# -*- coding: utf-8 -*-

"""
    Util functions used in notebook
"""

import os
import json
import logging
import pandas as pd
from collections import Counter, defaultdict
from tqdm import tqdm
from itertools import product
from typing import Mapping, List

from venn import venn
from networkx import DiGraph, connected_components
import matplotlib.pyplot as plt

from drug2ways.alternative_graph_traversal import enumerate_paths
from drug2ways.pathway import _prepare_json
from drug2ways.wrapper import generate_reduced_graph
from drug2ways.graph_reader import load_graph
from drug2ways.rcr import (rcr_all_paths, validate_paths_with_disease_data,
                           disease_rcr_all_paths)
from drug2ways.cli_helper import _setup_logging

# XREF paths
DATA_DIR = '../../data'
DISEASE_XREF_PATH = os.path.join(DATA_DIR, 'xref', 'pathologies.tsv')
PROTEIN_XREF_PATH = os.path.join(DATA_DIR, 'xref', 'proteins.tsv')
DRUGBANK_XREF_PATH = os.path.join(DATA_DIR, 'xref', 'drugbank_to_pubchem.json')
KG_DATA_PATH = os.path.join(DATA_DIR, 'kg')

logger = logging.getLogger(__name__)

"""Functions related to dataset and network harmonization"""


def get_disease_map():
    """Get disease map"""
    df = pd.read_csv(
        filepath_or_buffer=DISEASE_XREF_PATH,
        sep='\t',
        header=None
    )
    df.columns = ['source', 'target']
    return dict(zip(df['source'], df['target']))


def get_gene_map():
    """Get gene map"""
    df = pd.read_csv(
        filepath_or_buffer=PROTEIN_XREF_PATH,
        sep='\t',
        header=None
    )
    df.columns = ['source', 'target']
    return dict(zip(df['source'], df['target']))


def get_chem_map():
    """Get chemical map"""
    with open(DRUGBANK_XREF_PATH) as f:
        return json.load(f)


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


def harmonize_dataset(
        data_dict: dict,
        threshold: int = 0,
) -> dict:
    """Converts expression data values to +1/-1 based on threshold."""

    for source in tqdm(data_dict, desc='Harmonizing dictionary'):
        for gene, exp_val in data_dict[source].items():
            if threshold == 0:
                if exp_val > threshold:
                    new_exp_val = 1
                else:
                    new_exp_val = -1
            else:
                if exp_val > threshold:
                    new_exp_val = 1
                elif exp_val < -threshold:
                    new_exp_val = -1
                else:
                    new_exp_val = 0

            data_dict[source][gene] = new_exp_val

    return data_dict


"""Notebook 3"""


def normalize_nodes(
        network_df: pd.DataFrame,
        file_name: str
) -> pd.DataFrame:
    """ Harmonizing different namespaces into selected ones. """

    normalized_file_path = os.path.join(
        KG_DATA_PATH,
        'normalized',
        file_name
    )

    if os.path.exists(normalized_file_path):
        os.remove(normalized_file_path)

    # Load mapping dictionary
    mapping_dict = get_disease_map()
    gene_map = get_gene_map()
    chem_map = get_chem_map()

    final_df = pd.DataFrame(columns=['source', 'target', 'polarity'])

    for row in tqdm(
            network_df.values,
            desc='Normalizing graph'
    ):
        tmp = {
            'source': [],
            'target': [],
            'polarity': []
        }
        source_node, target_node, polarity = row

        if source_node.startswith('hgnc:'):
            source_idx = gene_map.get(source_node)

        elif source_node.startswith('drugbank'):
            namespace, idx = source_node.split(':')
            pubchem_idx = chem_map.get(idx)
            source_idx = f'pubchem.compound:{pubchem_idx}'

        else:
            source_idx = source_node

        tmp['source'].append(source_idx)

        # Remove hp nodes
        if target_node.startswith('hp:'):
            continue

        if target_node.startswith('hgnc:'):
            target_idx = gene_map.get(target_node)

        elif target_node.startswith('doid:') or target_node.startswith('umls'):
            target_idx = mapping_dict.get(target_node)

        else:
            target_idx = target_node

        tmp['target'] = target_idx
        tmp['polarity'] = polarity

        tmp_df = pd.DataFrame.from_dict(tmp)
        final_df = pd.concat([final_df, tmp_df], ignore_index=True)

    final_df.dropna(inplace=True)

    return final_df


def filter_graph(
        network_df: pd.DataFrame,
        file_name: str,
        data_dict: dict
):
    """Filter network to have disease and chemicals specific to the dataset."""
    normalized_file_path = os.path.join(KG_DATA_PATH, 'normalized', file_name)

    # Map data to pubchem.compound and mondo
    final_df = normalize_nodes(network_df=network_df, file_name=file_name)

    # Keep only chemical and disease data present in our datasets
    final_nodes = set(data_dict['diseases']).union(set(data_dict['chemicals']))

    chemical_df = final_df[final_df['source'].str.contains('pubchem.compound')]
    chemical_els = set(
        idx
        for idx, x in chemical_df.source.items()
        if x in final_nodes
    )

    disease_df = final_df[final_df['target'].str.contains('mondo')]
    disease_els = set(
        idx
        for idx, x in disease_df.target.items()
        if x in final_nodes
    )

    # Keep all PPI data
    ppi_els = set(final_df[
                      final_df['source'].str.contains('ncbigene:') &
                      final_df['target'].str.contains('ncbigene:')
                      ].index.values.tolist()
                  )

    keep_elements = chemical_els.union(disease_els)
    keep_elements = keep_elements.union(ppi_els)

    # Remove nodes not present in dataset
    final_df = final_df[final_df.index.isin(list(keep_elements))]
    final_df.reset_index(inplace=True)
    final_df.drop('index', axis=1, inplace=True)

    final_df.to_csv(normalized_file_path, sep='\t', index=False)
    return final_df


"""Notebook 5"""


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


def create_subgraph(
        graph_df: pd.DataFrame,
        disease: str,
        chemical: str,
        disease_dict: dict,
        chemical_dict: dict,
) -> pd.DataFrame:
    """Create subgraph based on disease and chemical gene overlap."""
    subgraph_df = pd.DataFrame(columns=['source', 'target', 'polarity'])

    disease_genes = set()
    for i in disease_dict[disease]:
        disease_genes.add(i)

    chemical_genes = set()
    for i in chemical_dict[chemical]:
        chemical_genes.add(i)

    gene_list = disease_genes.union(chemical_genes)

    # No gene overlap found
    if len(gene_list) < 1:
        return subgraph_df

    for source, target, relation in graph_df.values:
        # chemical-gene node
        if (
                source.startswith('pubchem.compound') and
                target.startswith('ncbigene')
        ):
            if target not in gene_list:
                continue

            assert target in gene_list

        # gene-gene node
        elif source.startswith('ncbigene') and target.startswith('ncbigene'):
            if source in gene_list and target not in gene_list:
                continue

            if source not in gene_list and target in gene_list:
                continue

            if source not in gene_list and target not in gene_list:
                continue

            assert source in gene_list and target in gene_list

        # gene-disease node
        elif source.startswith('ncbigene') and target.startswith('mondo'):
            if source not in gene_list:
                continue

            assert source in gene_list

        tmp_dict = {
            'source': [source],
            'target': [target],
            'polarity': [relation]
        }
        tmp_df = pd.DataFrame.from_dict(tmp_dict)
        subgraph_df = pd.concat([subgraph_df, tmp_df], ignore_index=True)

    return subgraph_df


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
        graph_df: pd.DataFrame,
        disease_dict: dict,
        chemical_dict: dict,
        graph_name: str,
        openbio: bool
) -> None:
    """Get paths in graph."""
    detailed_info = {}

    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)

    CACHE_DIR = os.path.join(DATA_DIR, 'lmax-pairs')

    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    chemical_els = set(
        graph_df[
            graph_df['source'].str.contains('pubchem.compound:')
        ]['source']
    )
    disease_els = set(
        graph_df[
            graph_df['target'].str.contains('mondo:')
        ]['target']
    )

    graph = load_graph(graph_df=graph_df)

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
            total=len(chemical_dict) * len(disease_dict)
    ):
        if chemical in chemical_els and disease in disease_els:
            for lmax in range(5, 9):
                dict_key = f'lmax_{lmax}'

                if openbio:
                    file_name = graph_name + '-' + dict_key + '-openbio.json'
                else:
                    file_name = graph_name + '-' + dict_key + '-custom.json'

                if os.path.exists(os.path.join(CACHE_DIR, file_name)):
                    continue

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

    for el in detailed_info:

        if openbio:
            file_name = graph_name + '-' + el + '-openbio.json'
        else:
            file_name = graph_name + '-' + el + '-custom.json'

        with open(os.path.join(CACHE_DIR, file_name), 'w') as f:
            json.dump(detailed_info[el], f, ensure_ascii=False, indent=2)


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


"""Notebook 6"""

ERRORS_ALLOWED = 1


def get_validated_paths(
        directed_graph: DiGraph,
        source: str,
        target: str,
        all_paths: List[list],
        drug_dict: dict,
        disease_dict: dict,
        clinical_pair_dict=dict,
        fda_pairs=set,
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

    # Check if pair in clinical trials
    pair = source + '_' + target
    if pair in clinical_pair_dict:
        trial_pair = True
    else:
        trial_pair = False

    if pair in fda_pairs:
        fda_val = True
    else:
        fda_val = False

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
        "in_fda": fda_val,
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
        chemical_dict: dict,
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
    chemical_filtered_paths = rcr_all_paths(
        directed_graph,
        all_paths,
        chemical_dict,
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
        "chemical_paths": chemical_filtered_paths,
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
