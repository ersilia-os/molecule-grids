# Grids of chemical structures for slides and papers

Draw SMILES as a clean grid of 2D structures. **Slide** format uses RDKit colours. **Print** format follows the ChemDraw *ACS Document 1996* style (Arial labels, 14.4 pt bonds at most). Analogues sharing a scaffold are drawn in the same orientation.

```bash
pip install git+https://github.com/ersilia-os/molecule-grids.git          # add [app] for the web app
```

PNG, PDF and GIF export need the cairo library (`brew install cairo` or `apt install libcairo2`).

## Python

```python
from molecule_grids import draw_grid, read_molecules

grid = draw_grid(["CCO", "c1ccccc1O"], names=["Ethanol", "Phenol"], format="print", width=0.5, size="large")
grid.save("grid.svg")  # .svg, .png, .pdf or .gif

smiles, names = read_molecules("hits.csv")  # .csv with a smiles column, or .smi lines "SMILES name"
draw_grid(smiles, names, format="slide", group=True, number=True, frame=True).save("hits.png")
```

`draw_grid` takes SMILES or RDKit molecules. The returned `Grid` displays inline in Jupyter and reports its size, `bond_pt`, `columns` and `capacity`.

## Sizing

Sizes follow [stylia](https://github.com/ersilia-os/stylia): `width` is a fraction of the full format width (7.09 in for print, 13 in for slides).

- Slide figures are at most 7.31 in tall (a 16:9 slide); more molecules split into pages (GIF frames). Print figures grow as needed.
- Molecules have a fixed size: `size="small"`, `"medium"` (default) or `"large"`, i.e. 70, 100 or 140 % of the standard bond (14.4 pt in print, ACS; 20 pt on slides).
- Each molecule sits in a square cell. The width holds as many columns as fit; rows are added as needed.
- `squeeze` drops the square cells and reorders (and turns) molecules to fill the figure at the same bond length: `"grid"` (aligned, also `True`), `"rows"` (flowing like text) or `"free"` (interlocking outlines).
- A figure holds up to 100 molecules; more raise `TooManyMolecules`. Save as `.gif` to get one frame per 100.

## Command line

| Command | Description |
|---|---|
| `molecule-grids draw INPUT -o OUT` | Draw a grid from a `.smi`/`.csv` file (`-f print\|slide`, `-w`, `-s small\|medium\|large`, `--number`, `--group`, `--frame`, `--squeeze [grid\|rows\|free]`) |
| `molecule-grids examples -o FILE` | Write random global-health drugs as a `.smi` file |
| `molecule-grids app` | Launch the Streamlit app |

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
