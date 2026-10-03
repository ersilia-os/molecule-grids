# Grids of chemical structures for slides and papers

Draw SMILES as a clean grid of 2D structures, for slides or print. Two drawing styles: **medicinal chemist** (default; ChemDraw *ACS Document 1996*, black and white, Arial labels) and **computational** (RDKit colours). Analogues sharing a scaffold are drawn in the same orientation.

```bash
pip install git+https://github.com/ersilia-os/molecule-grids.git          # add [app] for the web app
```

On Linux servers, RDKit's drawing code needs three small X libraries: `apt install libxrender1 libxext6 libexpat1` (listed in `packages.txt` for Streamlit Community Cloud). Nothing else is needed: PNG and PDF are rendered with self-contained wheels (resvg, typst).

## Python

```python
from molecule_grids import draw_grid, read_molecules

grid = draw_grid(["CCO", "c1ccccc1O"], names=["Ethanol", "Phenol"], format="print", size="large")
grid.save("grid.svg")  # .svg, .png or .pdf; one file per page (grid_1.svg, ...) if needed

smiles, names = read_molecules("hits.csv")  # .csv with a smiles column, or .smi lines "SMILES name"
draw_grid(smiles, names, format="slide", style="computational", number=True).save("hits.png")
```

`draw_grid` takes SMILES or RDKit molecules. The returned `Grid` displays inline in Jupyter and reports its size (`width_mm`, `height_mm`, and `width`, the proportion used), `bond_pt`, `columns` and `capacity`.

## Sizing

Sizes follow [stylia](https://github.com/ersilia-os/stylia).

- **Width:** `width="auto"` (default) picks the proportion of the full width (180 mm print, two journal columns; 330 mm slide) whose figure is closest to 3:2 in print or 16:9 on slides. A number from 0.25 to 1 sets it yourself. With `fit=True` (default) the figure is never wider than its content needs.
- **Height:** at most one page, 247 mm in print (Nature) and 186 mm on slides (16:9), or less with `max_height` (a fraction of the page). More molecules go to further pages; nothing is cut.
- **Molecule size:** `size="small"`, `"medium"` (default) or `"large"` sets atom labels and captions to stylia's font sizes (print 5/6/8 pt, slides 8/10/13 pt); bonds are 1.44 times that, as in ChemDraw ACS.
- **Layout:** `squeeze` reorders (and turns) molecules to fill the figure at the same size: `"rows"` (default; flowing like text), `"grid"` (aligned columns) or `"free"` (interlocking outlines). `squeeze=False` uses aligned cells in input order, wide molecules spanning columns. `group=True` keeps analogues sharing a scaffold together in any layout.
- **Pages:** `grid.pages` holds one SVG per page, and `save` writes one numbered file per page.

## Command line

| Command | Description |
|---|---|
| `molecule-grids draw INPUT -o OUT` | Draw a grid from a `.smi`/`.csv` file (`-f print\|slide`, `--style medicinal\|computational`, `-w auto\|0.25–1`, `--max-height`, `-s small\|medium\|large`, `--squeeze rows\|grid\|free\|off`, `--group`, `--number`, `--no-fit`, `--no-frame`) |
| `molecule-grids examples -o FILE` | Write 20 to 50 random global-health drugs as a `.smi` file (`-n` to choose) |
| `molecule-grids app` | Launch the Streamlit app |

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
