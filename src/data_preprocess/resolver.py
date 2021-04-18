# -*- coding: utf-8 -*-

"""
    Normalization of entities
"""

import pandas as pd
from typing import Optional
from ratelimit import limits

import requests

CLASS_NAME_MAP = {
    'protein': 'ncbigene',
    'pathology': 'pathologies',
    'others': 'others',
}

CLASS_XREF_MAP = {
    'protein': 'proteins',
    'pathology': 'pathologies',
}


def local_curies_to_names(class_name):
    df = pd.read_csv(
        f'../../data/xref/{CLASS_NAME_MAP[class_name]}.tsv',
        sep='\t',
        header=None
    )
    if df.shape[1] == 3:
        df.columns = ['curie', 'name_count', 'name']
    elif df.shape[1] == 2:
        df.columns = ['curie', 'name']
    else:
        raise ValueError('wrong number of columns')

    return dict(zip(df['curie'], df['name']))


def get_local_xrefs(class_name):
    df = pd.read_csv(f'../../data/xref/{CLASS_XREF_MAP[class_name]}.tsv', sep='\t', header=None)
    df.columns = ['source', 'target']
    return dict(zip(df['source'], df['target']))


# Get information from PubChem

def _extract_cid_from_response(response: requests.Response) -> Optional[str]:
    if response.status_code == 200:
        response_json = response.json()
        identifier_list = response_json.get("IdentifierList")
        if identifier_list:
            pubchem_compound_ids = identifier_list.get("CID", [])
            if pubchem_compound_ids and len(pubchem_compound_ids) > 0:
                # TODO add explanation of why there might be multiple and why its only taking the first
                return pubchem_compound_ids[0]
    return None


@limits(calls=5, period=1)
def _call_api_helper(url, method="GET", data=None) -> requests.Response:
    if method == "GET":
        return requests.get(url)
    else:
        return requests.post(
            url,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data=data,
        )


@limits(calls=5, period=1)
def _call_api(url, method="GET", data=None) -> Optional[str]:
    response = _call_api_helper(url=url, method=method, data=data)
    return _extract_cid_from_response(response)


def get_cid_by_smiles(smiles: str) -> Optional[str]:
    link = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
    return _call_api(link)
