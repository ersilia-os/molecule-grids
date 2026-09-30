import random

import click

from molecule_grids.data import examples as ex


@click.command()
@click.option("-o", "--output", required=True, type=click.Path(dir_okay=False), help="Output .smi file.")
@click.option("-n", "--n", "n", type=int, default=ex.SHUFFLE_SIZE, show_default=True, help="Number of molecules.")
@click.option("--seed", type=int, default=None, help="Random seed.")
def examples(output, n, seed):
    """Write random global-health drugs, including analogue families, as a .smi file."""
    pool = [m for fam in ex.FAMILIES for m in fam] + ex.SINGLES
    if not 1 <= n <= len(pool):
        raise click.BadParameter(f"must be between 1 and {len(pool)}", param_hint="-n")
    rng = random.Random(seed)
    picked = ex.shuffled(n, rng) if n >= 9 else rng.sample(pool, n)
    with open(output, "w") as f:
        f.write(ex.as_text(picked))
    click.echo(f"Wrote {n} molecules to {output}")
