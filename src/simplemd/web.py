"""Flask HTTP API.

Exposes ``POST /api/run_job``, preserving the original request/response
contract:

* request JSON: ``{N, density, temperature, delta_t, steps, system,
  integrator, forces}`` (optional ``seed``)
* response JSON: ``{"output": <trajectory as nested lists>}``

Use :func:`create_app` as an application factory (Flask / gunicorn)::

    gunicorn "simplemd.web:create_app()"
"""

from __future__ import annotations

import json

import numpy as np
from flask import Flask, request

from simplemd.api import run_simulation


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that serialises NumPy arrays as nested lists."""

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def create_app() -> Flask:
    """Create and configure the Flask application."""
    app = Flask(__name__)

    @app.route("/api/run_job", methods=["POST"])
    def run_job():
        data = request.json
        traj = run_simulation(
            system=data["system"],
            integrator=data["integrator"],
            forces=data["forces"],
            N=data["N"],
            density=data["density"],
            temperature=data["temperature"],
            delta_t=data["delta_t"],
            steps=data["steps"],
            seed=data.get("seed"),
        )
        return json.dumps({"output": traj}, cls=NumpyEncoder)

    return app


# Module-level app for `gunicorn simplemd.web:app` style invocation.
app = create_app()
