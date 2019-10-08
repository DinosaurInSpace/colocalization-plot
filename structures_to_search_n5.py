"""
Contains structures to search against in SMARTS format
5x substrucutes from Sebastian Wolf et al BMC Bioinform.

Future to test:
140x neutral losses from Ma, Yan, et al Anal. Chem.
"""

#SMARTS
target_structures = {
    'Hydroxyl': '[OX2H]',
    'Cyano': '[CX2]#[NX1]',
    'Amino': '[NX3;H2,H1;!$(NC=O)]',
    'Aldehyde': '[CX3H1](=O)[#6]',
    'Carboxylic': '[CX3](=O)[OX1H0-,OX2H1]',
}

#Formula
target_loss_forumula = {
    'Hydroxyl': 'H2O',
    'Cyano': 'CN',
    'Amino': 'NH2',
    'Aldehyde': 'COH',
    'Carboxylic': 'CO2H',
}

#Exact mass Da, float
target_nl_mass = {
    'Hydroxyl': 18.0106,
    'Cyano': 27.0109,
    'Amino': 17.0266,
    'Aldehyde': 30.0106,
    'Carboxylic': 46.0055,
}