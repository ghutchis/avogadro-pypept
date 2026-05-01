# avogadro-pypept

Avogadro 2 plugin to build peptides from BILN, HELM, or FASTA sequences using
[pyPept](https://github.com/Boehringer-Ingelheim/pyPept).

The plugin adds a **Build → Peptide…** menu item that prompts for:

- **Sequence Format** — `BILN`, `HELM`, or `FASTA`
- **Peptide Sequence** — e.g. `ac-D-T-H-F-E-I-A-am` (BILN) or `DTHFEIA` (FASTA)
- **Secondary Structure** — `Predict` (uses pyPept's `SecStructPredictor`),
  `Coil`, `Helix`, or `Sheet`

The generated 3D peptide is returned to Avogadro as a PDB block with residue
information preserved.

## Installation

This plugin installs pyPept directly from its GitHub repository:

```bash
pip install git+https://github.com/ghutchis/avogadro-pypept.git
```

Or with [pixi](https://pixi.sh):

```bash
pixi install
```

## Citation

Ochoa, R., Brown, J. B., Fox, T. *pyPept: a python library to generate
atomistic 2D and 3D representations of peptides.* J. Cheminform. **15**, 79
(2023). DOI: [10.1186/s13321-023-00748-2](https://doi.org/10.1186/s13321-023-00748-2)
