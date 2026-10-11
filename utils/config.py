"""Config loading for the RNA fingerprinting pipeline (Python side).

Pure/side-effect-free: everything reads from config.conf at the repo root.
Scripts in pipeline/ import this instead of hard-coding paths or tunables.
Mirrors the INI consumed by utils/config.R so Python and R share one config.
"""
from __future__ import annotations

import configparser
from pathlib import Path

# config.conf lives one level up from this file: <repo>/utils/config.py -> <repo>/config.conf
_DEFAULT_CONF = Path(__file__).resolve().parents[1] / "config.conf"


def load_config(path: str | Path | None = None) -> configparser.ConfigParser:
    """Load the INI config. Raises if the file is missing."""
    conf_path = Path(path) if path else _DEFAULT_CONF
    if not conf_path.is_file():
        raise FileNotFoundError(f"config not found: {conf_path}")
    cfg = configparser.ConfigParser()
    cfg.read(conf_path)
    return cfg


def get_path(cfg: configparser.ConfigParser, section: str, key: str,
             required: bool = True) -> Path | None:
    """Fetch a path value; error if required and blank/missing."""
    val = cfg.get(section, key, fallback="").strip()
    if not val:
        if required:
            raise ValueError(f"[{section}] {key} is blank in config.conf — fill it in.")
        return None
    return Path(val)


def query_section(cfg: configparser.ConfigParser, cohort: str | None = None) -> str:
    """Return the [query.<cohort>] section name for the active cohort.

    Defaults to [run].cohort when `cohort` is None.
    """
    cohort = cohort or cfg.get("run", "cohort").strip()
    section = f"query.{cohort}"
    if not cfg.has_section(section):
        raise KeyError(f"no [{section}] section in config.conf")
    return section
