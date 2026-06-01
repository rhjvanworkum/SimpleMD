# Contributing to SimpleMD

Thanks for your interest! This project uses [uv](https://docs.astral.sh/uv/) for
environment management and standard Python tooling.

## Development setup

```bash
uv sync --extra web --extra viz   # install everything, including dev tools
```

The dev dependency group (pytest, ruff, mypy) is installed automatically by
`uv sync`.

## Running the checks

```bash
uv run pytest                 # tests (+ coverage config in pyproject)
uv run pytest --cov=simplemd  # with coverage report
uv run ruff check .           # lint
uv run ruff format .          # format (use --check in CI)
uv run mypy                   # type check
```

All four must pass before a PR is merged; CI runs them on every push and PR.

## Conventions

- **Preserve behavior.** The Lennard-Jones path is pinned to a tight tolerance
  by a characterization test (`tests/test_regression.py`) against golden data
  captured from a seeded run (bit-for-bit on a fixed machine; ~1e-9 across
  platforms). If you intend to change numerics, update the golden data
  deliberately and explain why.
- **Reproducibility.** Use `seed=...` in tests so runs are deterministic.
  Initialisation uses the global `numpy`/`random` generators on purpose (the
  `NPY002` lint is disabled) to keep seeding simple and stable.
- **Type hints** on public functions; mypy is lenient for now and tightened
  over time.
- Keep new dependencies minimal — prefer the standard library and NumPy.

## Physics regression tests

`tests/test_bugfixes.py` guards three behaviors that were previously broken and
have since been fixed: cell-division force evaluation (must match the
pair-by-pair reference), the predictor-corrector integrator for the LJ model,
and Newton's third law in the TIP-4P site forces. Keep these green when touching
the force or integrator code.
