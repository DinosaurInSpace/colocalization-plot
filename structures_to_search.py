# Contains structures to search against in SMARTS format

target_structures = {
    'Carboxylic_acid_or_conj_base': '[CX3](=O)[OX1H0-,OX2H1]',
    'CA_CB_aliphatic': '[CX4][CX3](=O)[OX1H0-,OX2H1]',
    'CA_CB_aromatic': '[c]:[cX3][CX3](=O)[OX1H0-,OX2H1]',
    'Hydroxyl': '[OX2H]',
    'Hydroxyl_in_alcohol': '[#6][OX2H]',
    'Hydroxyl_bonded_to_alkyl': '[OX2H][CX4]',
    'Phenol': '[OX2H][cX3]:[c]'
}
