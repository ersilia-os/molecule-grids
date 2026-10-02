import click

from molecule_grids import FORMATS, TooManyMolecules, draw_grid, read_molecules


@click.command()
@click.argument("input", type=click.Path(exists=True, dir_okay=False))
@click.option("-o", "--output", required=True, type=click.Path(dir_okay=False), help="Output .svg, .png, .pdf or .gif.")
@click.option("-f", "--format", "fmt", type=click.Choice(sorted(FORMATS)), default="slide", show_default=True)
@click.option("-w", "--width", type=float, default=1.0, show_default=True, help="Fraction of the full format width.")
@click.option(
    "-s",
    "--size",
    type=click.Choice(["small", "medium", "large"]),
    default="medium",
    show_default=True,
    help="Molecule size: 70, 100 or 140 % of the standard bond.",
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
def draw(input, output, fmt, width, size, number, group, frame, squeeze):
    """Draw a grid from INPUT (.smi/.txt with 'SMILES name' lines, or .csv with a smiles column)."""
    smiles, names = read_molecules(input)
    if not smiles:
        raise click.ClickException(f"no valid molecules in {input}")
    try:
        grid = draw_grid(
            smiles,
            names,
            format=fmt,
            width=width,
            size=size,
            number=number,
            group=group,
            frame=frame,
            squeeze=False if squeeze == "off" else squeeze,
        )
        grid.save(output)
    except (TooManyMolecules, ValueError, OSError) as e:
        raise click.ClickException(str(e)) from e
    click.echo(repr(grid))
