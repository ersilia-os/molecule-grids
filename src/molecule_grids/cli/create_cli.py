"""Command-line interface: ``molecule-grids``."""

import click

from molecule_grids.cli.commands.app import app
from molecule_grids.cli.commands.draw import draw
from molecule_grids.cli.commands.examples import examples


@click.group()
@click.version_option(package_name="molecule-grids")
def cli():
    """Grids of 2D chemical structures for slides and manuscripts."""


cli.add_command(draw)
cli.add_command(examples)
cli.add_command(app)
