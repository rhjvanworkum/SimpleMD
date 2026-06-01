# SimpleMD

A small, from-scratch molecular dynamics simulator: a **Lennard-Jones soft fluid** and a **TIP-4P water** model, with a clean Python API, a CLI, and an optional HTTP service.

## Install

SimpleMD uses [uv](https://docs.astral.sh/uv/) for reproducible environments.

```bash
git clone https://github.com/rhjvanworkum/SimpleMD
cd SimpleMD
uv sync                 # core install
uv sync --extra web     # + Flask/gunicorn for the HTTP API
uv sync --extra viz     # + pygame for the trajectory viewer
```

## Quickstart

### Python API

```python
from simplemd import run_simulation

# Returns a NumPy array of shape (n_frames, n_bodies, n_dim).
traj = run_simulation(
    system="LJ",            # "LJ" (Lennard-Jones) or "TIP4P" (water)
    integrator="leapfrog",  # "leapfrog", "verlet", or "pred-corr"
    forces="PBP",           # "PBP" (pair-by-pair)
    N=4,                    # N**3 unit cells
    density=0.8,
    temperature=1.0,
    delta_t=0.005,
    steps=500,
    seed=0,                 # optional: makes the run reproducible
)
print(traj.shape)
```

### Command line

```bash
uv run simplemd --system LJ -N 4 --steps 500 --seed 0 -o trajectory.json
```

### HTTP API

```bash
uv run --extra web flask --app simplemd.web run        # dev server
# or for production:
uv run --extra web gunicorn "simplemd.web:create_app()"
```

`POST /api/run_job` with a JSON body:

```json
{"system": "LJ", "integrator": "leapfrog", "forces": "PBP",
 "N": 3, "density": 0.8, "temperature": 1.0, "delta_t": 0.005,
 "steps": 100, "seed": 0}
```

returns `{"output": [...]}` — the trajectory as nested lists.

## Examples

Runnable, self-contained scripts in [`examples/`](examples/):

```bash
uv run examples/run_lennard_jones.py
uv run examples/run_tip4p_water.py
```

To replay a saved trajectory (needs the `viz` extra):

```bash
uv run --extra viz python -m simplemd.viz output/test_2d_10_cd.json
```

## Project structure

```
src/simplemd/
├── api.py            # run_simulation(), build_system(), registries (public API)
├── web.py            # Flask app factory (create_app) — POST /api/run_job
├── cli.py            # command-line entry point
├── config.py         # Settings (run parameters, incl. optional seed)
├── rng.py            # opt-in RNG seeding for reproducible runs
├── reporter.py       # records trajectories during a run
├── utils.py          # vector/quaternion math helpers
├── viz.py            # pygame trajectory viewer (optional)
├── systems/          # System base class + LJ fluid + TIP-4P water
├── integrators/      # leapfrog, velocity-verlet, predictor-corrector, thermostats
├── interactions/     # force evaluation (pair-by-pair, cell-division)
└── particles/        # Particle and Molecule data classes
tests/                # pytest suite (mirrors the package) + golden data
examples/             # runnable demonstrations
```

## Reproducibility

Simulations draw random initial velocities/orientations. Passing `seed=...`
(API/CLI) or `"seed"` (HTTP) makes a run fully reproducible; omitting it
preserves the original non-deterministic behavior.

## Notes / known limitations

- The `CD` (cell-division) force method and the `pred-corr` integrator **with
  the LJ model** are known to be broken in the original code (documented by
  tests in `tests/test_known_bugs.py`) and are not yet fixed.
- `pred-corr` works with the TIP-4P model.

See [CONTRIBUTING.md](CONTRIBUTING.md) for the development workflow.

## License

MIT — see [LICENSE](LICENSE).
