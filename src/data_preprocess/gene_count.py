# -*- coding: utf-8 -*-

"""
    Gene level analysis
"""

import os
import pandas as pd
import json
from tqdm import tqdm
import logging
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

DATA_DIR = '../../data'


def get_kg_data():
    """Combine the data of both networks in json. """

    openbio_kg = pd.read_csv(
        os.path.join(DATA_DIR, 'kg', 'normalized', 'openbiolink_kg_normalized.tsv'),
        sep='\t',
        dtype=str,
        usecols=['source', 'target']
    )

    custom_kg = pd.read_csv(
        os.path.join(DATA_DIR, 'kg', 'normalized', 'custom_kg_normalized.tsv'),
        sep='\t',
        dtype=str,
        usecols=['source', 'target']
    )

    data = {
        'openbio': {
            'chemicals': [],
            'diseases': [],
            'genes': []
        },
        'custom': {
            'chemicals': [],
            'diseases': [],
            'genes': []
        },
    }

    for source, target in tqdm(openbio_kg.values, desc='Exporting data from OpenBioLink KG'):
        if ('pubchem.compound:' in source) & ('ncbigene:' in target):
            data['openbio']['chemicals'].append(source)
            data['openbio']['genes'].append(target)
        elif ('mondo:' in target) & ('ncbigene:' in source):
            data['openbio']['diseases'].append(target)
            data['openbio']['genes'].append(source)
        else:
            data['openbio']['genes'].append(source)
            data['openbio']['genes'].append(target)

    for source, target in tqdm(custom_kg.values, desc='Exporting data from Custom KG'):
        if ('pubchem.compound:' in source) & ('ncbigene:' in target):
            data['custom']['chemicals'].append(source)
            data['custom']['genes'].append(target)
        elif ('mondo:' in target) & ('ncbigene:' in source):
            data['custom']['diseases'].append(target)
            data['custom']['genes'].append(source)
        else:
            data['custom']['genes'].append(source)
            data['custom']['genes'].append(target)

    with open(os.path.join(DATA_DIR, 'kg', 'normalized', 'data.json'), 'w') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def gene_count(
        data_file: str,
        data_type: str,
        data_info: dict
):
    """Get statistics on gene data for different datasets. """
    with open(data_file) as f:
        data_dict = json.load(f)

    count_dict = {}

    for i in data_dict:
        if i in data_info[data_type]:
            count_dict[i] = len(data_dict[i])

    folder = data_file.split('/')[-2]

    counter = count_dict.values()

    logger.warning(f'Successfully found {len(counter)} from {len(data_dict)}')
    logger.warning(f'Maximum gene count - {max(counter)}')
    logger.warning(f'Minimum gene count - {min(counter)}')
    logger.warning(f'Average gene count - {sum(counter) / len(counter)}')

    plt.hist(counter)
    plt.title(f'Gene distribution for {folder}')
    plt.savefig(f'{folder}-{data_type}.jpg')
    plt.show()

    return None


def combined_data(
        data_sets: dict,
        data_info: dict
):
    """Method to get all data from all the datasets and save the combined set. """

    disease_genes = set()
    chemical_genes = set()
    diseases = set()
    chemicals = set()

    for el_type in data_sets:
        if el_type == 'diseases':
            file_name = 'disease_expression.json'
        else:
            file_name = 'chemical_expression.json'

        for folder_name in data_sets[el_type]:
            logger.warning(f'#####{folder_name}-{el_type}######')
            with open(os.path.join(DATA_DIR, folder_name, 'normalized', file_name)) as f:
                data_dict = json.load(f)

                for i in data_dict:
                    # Add chemicals and diseases in graph
                    if i in data_info[el_type]:
                        if el_type == 'diseases':
                            diseases.add(i)
                        else:
                            chemicals.add(i)

                        # Add genes for those chemicals
                        for gene_idx in data_dict[i]:
                            if el_type == 'diseases':
                                disease_genes.add(gene_idx)
                            else:
                                chemical_genes.add(gene_idx)

    logger.warning(f'Total disease genes - {len(disease_genes)}')
    logger.warning(f'Total chemical genes - {len(chemical_genes)}')
    logger.warning(f'Total diseases - {len(diseases)}')
    logger.warning(f'Total chemicals - {len(chemicals)}')

    if not os.path.exists(os.path.join(DATA_DIR, 'combined')):
        os.makedirs(os.path.join(DATA_DIR, 'combined'))

    with open(os.path.join(DATA_DIR, 'combined', 'dataset_genes.json'), 'w') as file:
        info = {
            'disease_gene': list(disease_genes),
            'chemical_gene': list(chemical_genes),
            'diseases': list(diseases),
            'chemicals': list(chemicals)
        }
        json.dump(info, file, ensure_ascii=False, indent=2)

    return None


def analyze_gene(
        data_sets: dict,
        data_info: dict
):
    """Overview of gene with respect to network.  """

    with open(os.path.join(DATA_DIR, 'combined', 'dataset_genes.json')) as file:
        info = json.load(file)

    for el_type in data_sets:
        if el_type == 'diseases':
            file_name = 'disease_expression.json'
        else:
            file_name = 'chemical_expression.json'

        for folder_name in data_sets[el_type]:
            logger.warning(f'#####{folder_name}-{el_type}######')
            d = 0
            c = 0
            u = 0

            with open(os.path.join(DATA_DIR, folder_name, 'normalized', file_name)) as f:
                data_dict = json.load(f)
                total = len(data_dict)
                for i in data_dict:
                    if i in data_info[el_type]:
                        if el_type == 'diseases':
                            if i in info['diseases']:
                                d += 1
                            else:
                                u += 1

                        if el_type == 'chemicals':
                            if i in info['chemicals']:
                                c += 1
                            else:
                                u += 1

            logger.warning(f'Found {d} diseases, {c} chemicals and {u} unmapped out of {total}')


def get_target_data(
        get_chemicals: bool,
        get_diseases: bool,
        custom_kg: bool,
):
    """Get specific data from combined data based on network. """

    with open(os.path.join(DATA_DIR, 'combined', 'dataset_genes.json'), 'r') as f:
        data_dict = json.load(f)

    with open(os.path.join(DATA_DIR, 'kg', 'normalized', 'data.json'), 'r') as f:
        if custom_kg:
            kg_data = json.load(f)['custom']
        else:
            kg_data = json.load(f)['openbio']

    filtered_dict = {
        'chemicals': [],
        'diseases': []
    }

    if get_chemicals:
        chemicals = data_dict['chemicals']

        for chemical in tqdm(chemicals, desc='Filtering chemicals'):
            if chemical in kg_data['chemicals']:
                filtered_dict['chemicals'].append(chemical)

    if get_diseases:
        chemicals = data_dict['diseases']

        for disease in tqdm(chemicals, desc='Filtering diseases'):
            if disease in kg_data['diseases']:
                filtered_dict['diseases'].append(disease)

    if custom_kg:
        g = 'customkg'
    else:
        g = 'openbiokg'

    with open(os.path.join(DATA_DIR, 'others', f'overlap_{g}.json'), 'w') as file:
        json.dump(filtered_dict, file, ensure_ascii=False, indent=2)


def main(
        do_gene_count: bool = False,
        get_data_list: bool = False,
        compare_with_network: bool = False
):
    try:
        with open(os.path.join(DATA_DIR, 'network', 'normalized', 'data.json')) as file:
            data_info = json.load(file)
    except FileNotFoundError:
        raise FileNotFoundError("Please get get_kg_data function first!")

    data_sets = {
        'diseases': ['creeds', 'geo', 'open_targets'],
        'chemicals': ['creeds', 'l1000']
    }

    if get_data_list:
        combined_data(data_sets=data_sets, data_info=data_info)

    if compare_with_network:
        analyze_gene(data_sets=data_sets, data_info=data_info)

    if do_gene_count:
        for el_type in data_sets:
            if el_type == 'diseases':
                file_name = 'disease_expression.json'
            else:
                file_name = 'chemical_expression.json'

            for folder_name in data_sets[el_type]:
                logger.warning(f'#####{folder_name}-{el_type}######')

                gene_count(
                    data_file=os.path.join(DATA_DIR, folder_name, 'normalized', file_name),
                    data_type=el_type,
                    data_info=data_info
                )


if __name__ == '__main__':
    get_target_data(get_chemicals=True, get_diseases=True, custom_kg=True)
