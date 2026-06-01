"""Tests for the viewer's pure (non-UI) helpers.

The interactive pygame loop needs a display and is not unit-tested here; we
cover the trajectory loader, which is the part with logic worth testing.
"""

import json

import numpy as np

from simplemd.viz import load_trajectory


def test_load_trajectory_reads_json_array(tmp_path):
    traj = [[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], [[0.1, 0.0, 0.0], [1.1, 1.0, 1.0]]]
    path = tmp_path / "traj.json"
    path.write_text(json.dumps(traj))
    loaded = load_trajectory(str(path))
    assert isinstance(loaded, np.ndarray)
    np.testing.assert_array_equal(loaded, np.array(traj))
