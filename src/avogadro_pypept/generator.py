#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-pypept/

import warnings

import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools

import pyPept.sequence as _pypept_sequence
from pyPept.sequence import SequenceConstants


def _get_monomer_info(path):
    # Replacement for pyPept.sequence.get_monomer_info that is compatible
    # with pandas >= 2.x. Upstream uses df.loc[idx, col] = list_value, which
    # newer pandas rejects with "Must have equal len keys and value when
    # setting with an iterable". Build each parsed column up-front and assign
    # the whole column at once with object dtype.
    df_group = PandasTools.LoadSDF(path)
    sep = SequenceConstants.csv_separator
    list_groups = ("m_Rgroups", "m_RgroupIdx", "m_attachmentPointIdx")
    for group in list_groups:
        parsed = []
        for raw in df_group[group]:
            parts = raw.split(sep)
            if group == "m_Rgroups":
                parsed.append([None if v == "None" else v for v in parts])
            else:
                parsed.append([None if v == "None" else int(v) for v in parts])
        df_group[group] = pd.Series(parsed, index=df_group.index, dtype=object)
    df_group = df_group.set_index("symbol")
    df_group = df_group.rename(columns={"ROMol": "m_romol"})
    return df_group


_pypept_sequence.get_monomer_info = _get_monomer_info

from pyPept.sequence import Sequence, correct_pdb_atoms  # noqa: E402
from pyPept.molecule import Molecule  # noqa: E402
from pyPept.converter import Converter  # noqa: E402
from pyPept.conformer import Conformer, SecStructPredictor  # noqa: E402


SS_SYMBOL = {"Coil": "C", "Helix": "H", "Sheet": "E"}


def _to_biln(sequence: str, fmt: str) -> str:
    fmt = fmt.lower()
    if fmt == "biln":
        return sequence
    if fmt == "helm":
        return Converter(helm=sequence).get_biln()
    if fmt == "fasta":
        # FASTA is a string of single-letter codes; BILN joins residues with '-'
        return "-".join(list(sequence.strip()))
    raise ValueError(f"Unknown peptide format: {fmt}")


def generate(opts):
    options = opts.get("options", {})
    fmt = options.get("format", "BILN")
    raw_sequence = options.get("sequence", "")
    ss_choice = options.get("secondary_structure", "Predict")

    biln = _to_biln(raw_sequence, fmt)

    seq = Sequence(biln)
    seq = correct_pdb_atoms(seq)
    mol = Molecule(seq)
    romol = mol.get_molecule(fmt="ROMol")

    fasta = Conformer.get_peptide(biln)
    if ss_choice == "Predict":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ss_value = SecStructPredictor.predict_active_ss(fasta)
    else:
        ss_value = SS_SYMBOL[ss_choice] * len(fasta)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        romol = Conformer.generate_conformer(romol, ss_value, generate_pdb=False)

    return Chem.MolToPDBBlock(romol)


def run(avo_input):
    return {
        "moleculeFormat": "pdb",
        "pdb": generate(avo_input),
    }
