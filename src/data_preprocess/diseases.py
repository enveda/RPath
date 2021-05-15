# -*- coding: utf-8 -*-

"""
    Get disease mapping from various gene expression files
"""

import os
import pickle
import pandas as pd
import logging
import json

from tqdm import tqdm

from resolver import get_local_xrefs
from pyenveda.resolver import canonicalize_name

logger = logging.getLogger(__name__)

disease_map = get_local_xrefs(class_name='pathology')
protein_map = get_local_xrefs(class_name='protein')

DATA_DIR = '../../data'


def get_creed_disease():
    df = pd.read_json(
        os.path.join(DATA_DIR, 'creeds', 'raw', 'data', 'disease_data.json'),
        orient='records',
    )

    df.drop(
        [
            'cell_type',
            'pert_ids',
            'geo_id',
            'platform',
            'curator',
            'geo_id',
            'version',
            'ctrl_ids',
            'id',
        ],
        axis=1,
        inplace=True,
    )

    df = df[df['organism'] == 'human']
    df.drop('organism', axis=1, inplace=True)

    # export cache to pickle each 15 iterations
    CACHE_PATH = os.path.join(DATA_DIR, 'creeds', 'disease_cache.pickle')

    if os.path.exists(CACHE_PATH):
        with open(CACHE_PATH, 'rb') as file:
            cache_dict = pickle.load(file)
    else:
        cache_dict = {}

    counter_cache = 0

    disease_expression_data = cache_dict

    unmapped_proteins = set()
    total_proteins = set()

    for line in tqdm(df.values, desc='Extracting disease information from CREED'):
        (
            doid,
            uml_id,
            down_regulated_genes,
            up_regulated_gene,
            disease_name,
        ) = line

        counter_cache += 1

        if pd.isna(doid) or doid == 'None':
            can_curie = disease_map.get(f'umls:{uml_id}')
        else:
            can_curie = disease_map.get(f'doid:{doid.split(":")[1]}')

        if can_curie in disease_expression_data:
            continue

        if not can_curie:
            continue

        down_regulated_genes.extend(up_regulated_gene)

        logger.info(f'Populating {len(down_regulated_genes)} genes')

        for gene_fold_change_pair in down_regulated_genes:
            hgnc_sym = gene_fold_change_pair[0]
            val = gene_fold_change_pair[1]
            total_proteins.add(hgnc_sym)

            ncbigene_id = protein_map.get(f'hgnc:{hgnc_sym}')

            if ncbigene_id:
                if can_curie not in disease_expression_data:
                    disease_expression_data[can_curie] = {}

                disease_expression_data[can_curie][ncbigene_id] = val

                cache_dict[can_curie][ncbigene_id] = val
            else:
                unmapped_proteins.add(hgnc_sym)

        # export cache to pickle and reset
        if counter_cache == 15:
            with open(CACHE_PATH, 'wb') as handle:
                pickle.dump(cache_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

            counter_cache = 0

    logger.info(f'Extracted {len(disease_expression_data)} diseases from CREED')
    logger.info(f'Unmapped {len(unmapped_proteins)} proteins out of {len(total_proteins)} proteins')

    with open(os.path.join(DATA_DIR, 'creeds', 'normalized', 'disease_expression.json'), 'w') as f:
        json.dump(disease_expression_data, f, ensure_ascii=False, indent=2)


def get_geo_disease():
    DATA_FILE = os.path.join(DATA_DIR, 'geo', 'raw', 'DEgenes')

    expression = {}

    total_proteins = set()
    unmapped_proteins = set()

    for file in tqdm(os.listdir(DATA_FILE)):
        can_id = disease_map.get(f'doid:{file}')

        if not can_id:
            continue

        if can_id not in expression:
            expression[can_id] = {}

        df = pd.read_csv(
            os.path.join(DATA_FILE, file, 'DEgenes.tsv'),
            sep='\t',
            usecols=[
                'HGNC',
                'logFC'
            ]
        )

        ncbigene_ids = []
        for val, hgnc_id in df.values:
            total_proteins.add(hgnc_id)
            ncbigene_id = protein_map.get(f'hgnc:{hgnc_id}')
            if ncbigene_id:
                ncbigene_ids.append(ncbigene_id)
            else:
                unmapped_proteins.add(hgnc_id)
                ncbigene_ids.append(None)
        df['ncbigene_id'] = ncbigene_ids
        df.drop('HGNC', axis=1, inplace=True)
        df.dropna(inplace=True)
        df.set_index('ncbigene_id', inplace=True)

        expression[can_id] = df.to_dict()['logFC']

    logger.info(f'Unmapped {len(unmapped_proteins)} proteins out of {len(total_proteins)} proteins')

    with open(os.path.join(DATA_DIR, 'geo', 'normalized', 'disease_expression.json'), 'w') as f:
        json.dump(expression, f, ensure_ascii=False, indent=2)


def get_open_targets_disease():
    df = pd.read_csv(
        os.path.join(DATA_DIR, 'open_targets', 'raw', 'data', 'open_targets.csv'),
        usecols=[
            'efo.term',
            'gene.symbol',
            'lfc',
        ]
    )

    diseases = df['efo.term'].unique()

    # export cache to pickle each 15 iterations
    CACHE_PATH = os.path.join(DATA_DIR, 'open_targets', 'disease_cache.pickle')

    if os.path.exists(CACHE_PATH):
        with open(CACHE_PATH, 'rb') as file:
            efo_to_mondo = pickle.load(file)
    else:
        efo_to_mondo = {}

    counter_cache = 0

    for i in tqdm(diseases, desc='Normalizing diseases'):
        counter_cache += 1

        if i in efo_to_mondo:
            continue
        efo_to_mondo[i] = canonicalize_name(i, entity_class='pathology').curie

        # export cache to pickle and reset
        if counter_cache == 15:
            with open(CACHE_PATH, 'wb') as handle:
                pickle.dump(efo_to_mondo, handle, protocol=pickle.HIGHEST_PROTOCOL)

            counter_cache = 0

    expression = {}

    total_proteins = set()
    unmapped_protein = set()

    for disease_name, gene, val in tqdm(df.values, desc='Exporting data'):
        if disease_name in efo_to_mondo:
            idx = efo_to_mondo[disease_name]
            if idx not in expression:
                expression[idx] = {}

            ncbigene_id = protein_map.get(f'hgnc:{gene}')
            total_proteins.add(gene)
            if ncbigene_id:
                expression[idx][ncbigene_id] = val
            else:
                unmapped_protein.add(gene)

    logger.info(f'Extracted {len(expression)} diseases from Open Targets dataset')
    logger.info(f'Unmapped {len(unmapped_protein)} protein out of {len(total_proteins)} proteins')

    with open(os.path.join(DATA_DIR, 'open_targets', 'normalized', 'disease_expression.json'), 'w') as f:
        json.dump(expression, f, ensure_ascii=False, indent=2)


def main(
    creed: bool = False,
    geo: bool = False,
    open_target: bool = False,
):

    if creed:
        get_creed_disease()

    if geo:
        get_geo_disease()

    if open_target:
        get_open_targets_disease()


if __name__ == '__main__':
    main()
