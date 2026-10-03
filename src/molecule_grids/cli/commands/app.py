import subprocess
import sys
from importlib import resources

import click


@click.command()
@click.option("--port", type=int, default=None, help="Port to serve on.")
def app(port):
    """Launch the Streamlit app (needs the [app] extra)."""
    try:
        import streamlit  # noqa: F401
    except ImportError:
        raise click.ClickException('Streamlit is not installed: pip install "molecule-grids[app]"') from None
    here = resources.files("molecule_grids") / "app"
    cmd = [sys.executable, "-m", "streamlit", "run", str(here / "streamlit_app.py")]
    if port:
        cmd += ["--server.port", str(port)]
    # Streamlit reads .streamlit/config.toml (the Ersilia theme) from the working directory.
    raise SystemExit(subprocess.call(cmd, cwd=str(here)))
