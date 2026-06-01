"""Tests for the command-line interface."""

import json

from simplemd.cli import main


def test_cli_runs_and_reports(capsys):
    exit_code = main(["--system", "LJ", "-N", "3", "--steps", "5", "--seed", "0"])
    assert exit_code == 0
    assert "trajectory shape" in capsys.readouterr().out


def test_cli_writes_output_file(tmp_path):
    out = tmp_path / "traj.json"
    exit_code = main(["--system", "LJ", "-N", "3", "--steps", "5", "--seed", "0", "-o", str(out)])
    assert exit_code == 0
    assert out.exists()
    data = json.loads(out.read_text())
    assert len(data) == 5  # 5 recorded frames
