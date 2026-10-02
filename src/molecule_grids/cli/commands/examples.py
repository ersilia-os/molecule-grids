import random

import click

from molecule_grids.data import examples as ex


@click.command()
@click.option("-o", "--output", required=True, type=click.Path(dir_okay=False), help="Output .smi file.")
@click.option(
    "-n",
    "--n",
    "n",
    type=int,
    default=None,
    help=f"Number of molecules [default: random, {ex.SIZES[0]} to {ex.SIZES[1]}].",
)
@click.option("--seed", type=int, default=None, help="Random seed.")
def examples(output, n, seed):
    """Write random global-health drugs, including analogue families, as a .smi file."""
    try:
        picked = ex.shuffled(n, random.Random(seed))
    except ValueError as e:
        raise click.BadParameter(str(e), param_hint="-n") from e
    with open(output, "w") as f:
        f.write(ex.as_text(picked))
    click.echo(f"Wrote {len(picked)} molecules to {output}")
