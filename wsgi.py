"""WSGI entry point for production servers (gunicorn / mod_wsgi).

Requires the package to be installed (``pip install .`` or ``uv sync``).
"""

from simplemd.web import app as application

__all__ = ["application"]
