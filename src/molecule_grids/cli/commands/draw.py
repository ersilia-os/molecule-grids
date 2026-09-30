import click

from molecule_grids import FORMATS, TooManyMolecules, draw_grid, read_molecules


@click.command()
@click.argument("input", type=click.Path(exists=True, dir_okay=False))
@click.option("-o", "--output", required=True, type=click.Path(dir_okay=False), help="Output .svg, .png, .pdf or .gif.")
@click.option("-f", "--format", "fmt", type=click.Choice(sorted(FORMATS)), default="slide", show_default=True)
@click.option("-w", "--width", type=float, default=1.0, show_default=True, help="Fraction of the full format width.")
@click.option("-c", "--columns", type=int, default=None, help="Number of columns [default: from width].")
@click.option("--number", is_flag=True, help="Caption with bold compound numbers 1, 2, 3...")
@click.option("--group", is_flag=True, help="Place analogues sharing a scaffold side by side.")
@click.option("--frame", is_flag=True, help="Outline the whole figure.")
@click.option("--squeeze", is_flag=True, help="Tight grid: reorder and shrink cells to fill the figure.")
def draw(input, output, fmt, width, columns, number, group, frame, squeeze):
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
            columns=columns,
            number=number,
            group=group,
            frame=frame,
            squeeze=squeeze,
        )
        grid.save(output)
    except (TooManyMolecules, ValueError, OSError) as e:
        raise click.ClickException(str(e)) from e
    click.echo(repr(grid))
