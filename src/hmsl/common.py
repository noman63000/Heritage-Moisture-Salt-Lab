"""Small shared utilities. All persisted scientific inputs are explicit and hashed."""
from __future__ import annotations
import hashlib, json, os, re
from pathlib import Path
from datetime import datetime, timezone

class InputError(ValueError):
    """A user-fixable input error, never a reason to substitute synthetic data."""

class Cancelled(RuntimeError):
    pass

def now() -> str:
    return datetime.now(timezone.utc).isoformat()

def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()

def digest(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, allow_nan=False).encode()).hexdigest()

def load_json(path):
    with open(path, encoding='utf-8') as f:
        return json.load(f)

def atomic_json(path, data):
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False), encoding='utf-8')
    os.replace(tmp, p)

def safe_path(root, rel):
    root = Path(root).resolve(); p = (root / str(rel)).resolve()
    if not p.is_relative_to(root):
        raise InputError('Path is outside the project folder.')
    return p

def safe_name(name):
    s = re.sub(r'[^a-zA-Z0-9_.-]+', '_', str(name)).strip('._')
    if not s or len(s) > 120:
        raise InputError('Use a short name containing letters, numbers, hyphens or underscores.')
    return s
