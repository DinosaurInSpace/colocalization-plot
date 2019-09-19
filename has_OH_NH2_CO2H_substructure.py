from rdkit import Chem
from __main__ import *
from contextlib import redirect_stdout
import io


def neutral_loss_finder(structure):
    # Can find arbitrary structures, as long as targets are in SMARTS.
    # https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html

    # Hydroxyl in alcohol
    hydroxy = Chem.MolFromSmarts('[CX4][OX2H]')

    # Primary or secondary amine, not amide.
    amine = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')  # Similiar issue, especially oxidation.

    # Carboxylic acid or conjugate base.
    carboxyl = Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]')

    targets = {"-OH": hydroxy, "-NH2": amine, "-CO2H": carboxyl}

    print([
        {"-OH": structure.HasSubstructMatch(targets["-OH"])},
        {"-NH2": structure.HasSubstructMatch(targets["-NH2"])},
        {"-CO2H": structure.HasSubstructMatch(targets["-CO2H"])}
        ])


def structure_converter(raw_structure, structure_type):
    # Checks for type of encoded structure

    if "INCHI" or "inchi" or "Inchi" in structure_type:
        structure = Chem.MolFromInchi(raw_structure)
        return structure

    elif "SMILES" or "smiles" or "Smiles" in structure_type:
        structure = Chem.MolFromSmiles(raw_structure)
        return structure

    elif "SMARTS" or "smarts" or "Smarts" in structure_type:
        structure = Chem.MolFromSmarts(raw_structure)
        return structure

    else:
        return "Error"


def stdout_to_var(structure):
    # Has SubstructMatch only delivers to standard out, instead capture as file.

    '''
    Instead of printing to file, print to variable!
    with redirect_stdout(open("temp_stdout_capture.txt", "a")):
        neutral_loss_finder(structure)
    '''

    with io.StringIO() as buf, redirect_stdout(buf):
        neutral_loss_finder(structure)
        output = buf.getvalue()
        return(output)