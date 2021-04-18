# -*- coding: utf-8 -*-

"""
    Get chemical mapping from various gene expression files
"""

import os
import pandas as pd
import json
import pickle
import logging
from tqdm import tqdm

from pyenveda.resolver import canonicalize_smiles, canonicalize_curie, canonicalize_name
from resolver import get_local_xrefs

logger = logging.getLogger(__name__)

protein_map = get_local_xrefs(class_name='protein')

DATA_DIR = '../../data'


def get_creed_chemicals():
    """
    Get chemical mapping from CREED
    """
    chemical_df = pd.read_json(
        os.path.join(DATA_DIR, 'creeds', 'raw', 'data', 'drug_data.json'),
        orient='records',
    )

    chemical_df.drop(
        [
            'cell_type',
            'pert_ids',
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

    chemical_df = chemical_df[chemical_df['organism'] == 'human']
    chemical_df.drop('organism', axis=1, inplace=True)

    CACHE_PATH = os.path.join(DATA_DIR, 'creeds', 'chemical_cache.pickle')

    if os.path.exists(CACHE_PATH):
        with open(CACHE_PATH, 'rb') as file:
            cache_dict = pickle.load(file)
    else:
        cache_dict = {}

    counter_cache = 0

    expression_data = cache_dict

    unmapped_protein = set()
    total_protein = set()

    for line in tqdm(chemical_df.values, desc='Extracting information from CREED'):
        (
            smiles,
            drugbank_id,
            pubchem_id,
            chem_name,
            down_regulated_genes,
            up_regulated_gene,
        ) = line

        counter_cache += 1

        if pd.isna(pubchem_id) or pubchem_id == 'nan':
            can_curie = canonicalize_curie(f'drugbank:{drugbank_id}')
            if can_curie is None:
                if pd.notna(smiles):
                    can_curie = canonicalize_smiles(smiles)
                if can_curie is None:
                    pubchem_data = canonicalize_name(chem_name, entity_class='chemical')
                    can_curie = pubchem_data.curie
        else:
            pubchem_id = str(pubchem_id)
            can_curie = f"pubchem.compound:{pubchem_id.split('.')[0]}"

        if can_curie is None:
            continue

        if can_curie in expression_data:
            continue

        down_regulated_genes.extend(up_regulated_gene)

        logger.info(f'Populating {len(down_regulated_genes)} genes')

        for gene_fold_change_pair in down_regulated_genes:
            hgnc_sym = gene_fold_change_pair[0]
            val = gene_fold_change_pair[1]

            total_protein.add(hgnc_sym)

            ncbigene_id = protein_map.get(f'hgnc:{hgnc_sym}')

            if ncbigene_id:
                if can_curie not in expression_data:
                    expression_data[can_curie] = {}

                expression_data[can_curie][ncbigene_id] = val

                cache_dict[can_curie][ncbigene_id] = val
            else:
                unmapped_protein.add(hgnc_sym)

        # export cache to pickle and reset
        if counter_cache == 15:
            with open(CACHE_PATH, 'wb') as handle:
                pickle.dump(cache_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

            counter_cache = 0

    logger.info(f'Extracted {len(expression_data)} chemicals from CREED')
    logger.info(f'Unmapped {len(unmapped_protein)} proteins out of {len(total_protein)} protein')

    with open(os.path.join(DATA_DIR, 'creeds', 'normalized', 'chemical_expression.json'), 'w') as f:
        json.dump(expression_data, f, ensure_ascii=False, indent=2)


def get_l1000_chemicals():
    df = pd.read_csv(
        os.path.join(DATA_DIR, 'lc1000', 'raw', 'data', 'l1000.csv'),
        usecols=[
            'pubchem.id',
            'gene.symbol',
            'direction',
        ]
    )

    df['pubchem.id'] = df['pubchem.id'].apply(lambda x: f'pubchem.compound:{x}')

    total_protein = set()
    unmapped_protein = set()

    ncbigene_ids = []
    for i in df['gene.symbol']:
        total_protein.add(i)
        idx = protein_map.get(f'hgnc:{i}')
        if not idx:
            unmapped_protein.add(i)
        ncbigene_ids.append(idx)

    df['ncbigene'] = ncbigene_ids
    df.drop('gene.symbol', axis=1, inplace=True)

    logger.info(f'Unmapped {len(unmapped_protein)} protein out of {len(total_protein)} proteins')

    expression_data = {}

    for pubchem_id, val, gene in tqdm(df.values, desc='Extracting values'):
        if pubchem_id:
            if pubchem_id not in expression_data:
                expression_data[pubchem_id] = {}

            if gene:
                expression_data[pubchem_id][gene] = val

    logger.info(f'Extracted {len(expression_data)} chemicals from L1000 dataset.')

    with open(os.path.join(DATA_DIR, 'lc1000', 'normalized', 'chemical_expression.json'), 'w') as f:
        json.dump(expression_data, f, ensure_ascii=False, indent=2)


def main(
        creed: bool = False,
        l1000: bool = False
):
    if creed:
        get_creed_chemicals()

    if l1000:
        get_l1000_chemicals()


if __name__ == '__main__':
    main()
