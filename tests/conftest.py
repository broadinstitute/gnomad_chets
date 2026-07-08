"""Shared pytest configuration for the ``gnomad_chets.v4`` Hail test suite.

Provides a single session-scoped Hail context for every ``v4/test_*.py``
module (mirrors ``gnomad_methods/tests/conftest.py``). Tests build tiny
in-memory Tables with ``hl.Table.parallelize`` and run against the local
Spark backend, so the whole suite completes in seconds with no Dataproc.

Run with::

    pytest v4/ -v
"""

import hail as hl
import pytest


@pytest.fixture(scope="session", autouse=True)
def setup_hail():
    """Initialize Hail once per test session and stop it at the end."""
    if not hl.utils.java.Env.is_fully_initialized():
        hl.init(log="/dev/null", quiet=True, idempotent=True)
    yield
    hl.stop()
