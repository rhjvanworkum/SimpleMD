"""Tests for the Flask HTTP API."""

import json

import pytest

from simplemd.web import create_app


@pytest.fixture
def client():
    app = create_app()
    app.config.update(TESTING=True)
    return app.test_client()


def test_run_job_returns_trajectory(client):
    payload = {
        "system": "LJ",
        "integrator": "leapfrog",
        "forces": "PBP",
        "N": 3,
        "density": 0.8,
        "temperature": 1.0,
        "delta_t": 0.005,
        "steps": 10,
        "seed": 0,
    }
    resp = client.post("/api/run_job", json=payload)
    assert resp.status_code == 200
    body = json.loads(resp.get_data(as_text=True))
    assert "output" in body
    # 10 frames, 27 particles, 3 dims
    assert len(body["output"]) == 10
    assert len(body["output"][0]) == 27
    assert len(body["output"][0][0]) == 3


def test_run_job_seed_is_reproducible(client):
    payload = {
        "system": "LJ",
        "integrator": "leapfrog",
        "forces": "PBP",
        "N": 3,
        "density": 0.8,
        "temperature": 1.0,
        "delta_t": 0.005,
        "steps": 5,
        "seed": 99,
    }
    a = json.loads(client.post("/api/run_job", json=payload).get_data(as_text=True))
    b = json.loads(client.post("/api/run_job", json=payload).get_data(as_text=True))
    assert a == b
