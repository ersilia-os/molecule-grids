import click

from molecule_grids import FORMATS, draw_grid, read_molecules


@click.command()
@click.argument("input", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(dir_okay=False),
    help="Output .svg, .png or .pdf; one numbered file per page.",
)
@click.option(
    "-f",
    "--format",
    "fmt",
    type=click.Choice(sorted(FORMATS)),
    default="slide",
    show_default=True,
    help="Page: slide or print.",
)
@click.option(
    "--style",
    type=click.Choice(["medicinal", "computational"]),
    default="medicinal",
    show_default=True,
    help="Medicinal chemist (ChemDraw ACS 1996) or computational (RDKit colours).",
)
@click.option("-w", "--width", type=float, default=1.0, show_default=True, help="Fraction of the full format width.")
@click.option(
    "--max-height", type=float, default=1.0, show_default=True, help="Tallest figure, as a fraction of the page."
)
@click.option(
    "-s",
    "--size",
    type=click.Choice(["small", "medium", "large"]),
    default="medium",
    show_default=True,
    help="Molecule size: labels at stylia's font sizes (print 5/6/8 pt, slide 8/10/13 pt).",
)
@click.option("--number", is_flag=True, help="Caption with bold compound numbers 1, 2, 3...")
@click.option("--group", is_flag=True, help="Place analogues sharing a scaffold side by side.")
@click.option("--frame/--no-frame", default=True, show_default=True, help="Outline the whole figure.")
@click.option(
    "--squeeze",
    type=click.Choice(["rows", "grid", "free", "off"]),
    default="rows",
    show_default=True,
    help="Reorder and turn molecules to fill the figure; 'off' keeps square cells in input order.",
)
def draw(input, output, fmt, style, width, max_height, size, number, group, frame, squeeze):
    """Draw a grid from INPUT (.smi/.txt with 'SMILES name' lines, or .csv with a smiles column)."""
    smiles, names = read_molecules(input)
    if not smiles:
        raise click.ClickException(f"no valid molecules in {input}")
    try:
        grid = draw_grid(
            smiles,
            names,
            format=fmt,
            style=style,
            width=width,
            max_height=max_height,
            size=size,
            number=number,
            group=group,
            frame=frame,
            squeeze=False if squeeze == "off" else squeeze,
        )
        grid.save(output)
    except (ValueError, OSError) as e:
        raise click.ClickException(str(e)) from e
    click.echo(repr(grid))
