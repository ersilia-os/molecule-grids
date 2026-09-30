from click.testing import CliRunner

from molecule_grids.cli.create_cli import cli


def test_examples_then_draw(tmp_path):
    runner = CliRunner()
    smi, svg = tmp_path / "mols.smi", tmp_path / "grid.svg"
    result = runner.invoke(cli, ["examples", "-o", str(smi), "-n", "6", "--seed", "1"])
    assert result.exit_code == 0, result.output
    result = runner.invoke(cli, ["draw", str(smi), "-o", str(svg), "-f", "print", "-w", "1", "--group"])
    assert result.exit_code == 0, result.output
    assert svg.read_text().startswith("<?xml")


def test_draw_csv(tmp_path):
    csv, png = tmp_path / "mols.csv", tmp_path / "grid.png"
    csv.write_text("name,smiles\nEthanol,CCO\nBenzene,c1ccccc1\n")
    result = CliRunner().invoke(cli, ["draw", str(csv), "-o", str(png)])
    assert result.exit_code == 0, result.output
    assert png.stat().st_size > 0
