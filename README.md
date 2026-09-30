# Grids of chemical structures for slides and papers

Draw SMILES as a clean grid of 2D structures. **Slide** format uses RDKit colours. **Print** format follows the ChemDraw *ACS Document 1996* style (Arial labels, 14.4 pt bonds at most). Analogues sharing a scaffold are drawn in the same orientation.

```bash
pip install git+https://github.com/ersilia-os/molecule-grids.git          # add [app] for the web app
```

PNG, PDF and GIF export need the cairo library (`brew install cairo` or `apt install libcairo2`).

## Python

```python
from molecule_grids import draw_grid, read_molecules

grid = draw_grid(["CCO", "c1ccccc1O"], names=["Ethanol", "Phenol"], format="print", width=0.5)
grid.save("grid.svg")  # .svg, .png, .pdf or .gif

smiles, names = read_molecules("hits.csv")  # .csv with a smiles column, or .smi lines "SMILES name"
draw_grid(smiles, names, format="slide", group=True, number=True).save("hits.png")
```

`draw_grid` takes SMILES or RDKit molecules. The returned `Grid` displays inline in Jupyter and reports its size, `bond_pt`, `columns`, `max_columns` and `capacity`.

## Sizing

Sizes follow [stylia](https://github.com/ersilia-os/stylia): `width` is a fraction of the full format width (7.09 in for print, 13 in for slides).

- Each molecule sits in a square cell, so the height follows from the number of rows.
- All molecules share one bond length, set so the largest molecule fits its cell.
- Columns are limited so that atom labels stay at least 5 pt (print) or 8 pt (slide). Rows are added as needed.
- A figure holds up to 100 molecules; more raise `TooManyMolecules`. Save as `.gif` to get one frame per 100.

## Command line

| Command | Description |
|---|---|
| `molecule-grids draw INPUT -o OUT` | Draw a grid from a `.smi`/`.csv` file (`-f print\|slide`, `-w`, `-c`, `--number`, `--group`) |
| `molecule-grids examples -o FILE` | Write random global-health drugs as a `.smi` file |
| `molecule-grids app` | Launch the Streamlit app |

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
