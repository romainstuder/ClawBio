"""Shared test fixtures and path setup for stability-predictor tests."""

from __future__ import annotations

import sys
from pathlib import Path

# Make `stability_predictor_core` and the CLI module importable as a top-level
# package when running `pytest skills/stability-predictor/tests`.
SKILL_DIR = Path(__file__).resolve().parents[1]
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))
