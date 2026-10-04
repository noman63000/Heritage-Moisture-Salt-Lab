"""Read-only resource telemetry and concurrency limits. No solver settings here."""
from __future__ import annotations
import ctypes
import os
import shutil
from pathlib import Path
from .safety import SafetyError

MAX_WORKERS = 21
GIB = 1024 ** 3

def resolve_workers(value, selected_count: int) -> int:
    if value == 'all':
        return max(1, min(MAX_WORKERS, selected_count))
    if isinstance(value, bool):
        raise SafetyError('Choose 1 to 21 simultaneous cases, or All selected.')
    try:
        n = int(value)
        if str(n) != str(value):
            raise ValueError()
    except (TypeError, ValueError, OverflowError):
        raise SafetyError('Choose a whole number from 1 to 21, or All selected.')
    if not 1 <= n <= MAX_WORKERS:
        raise SafetyError('Choose 1 to 21 simultaneous cases, or All selected.')
    return n

def memory_info() -> dict:
    """Available RAM, not just unused pages. None means telemetry unavailable."""
    try:
        if os.name == 'nt':
            class Status(ctypes.Structure):
                _fields_ = [('length', ctypes.c_ulong), ('load', ctypes.c_ulong)] + [
                    (n, ctypes.c_ulonglong) for n in ('total_phys', 'avail_phys',
                     'total_page', 'avail_page', 'total_virtual', 'avail_virtual', 'extended')]
            s = Status(); s.length = ctypes.sizeof(s)
            if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(s)):
                raise OSError('GlobalMemoryStatusEx failed')
            return {'total_bytes': int(s.total_phys), 'available_bytes': int(s.avail_phys)}
        p = Path('/proc/meminfo')
        if p.is_file():
            rows = {k: int(v.strip().split()[0]) * 1024 for k,v in
                    (line.split(':',1) for line in p.read_text().splitlines())}
            return {'total_bytes': rows['MemTotal'], 'available_bytes': rows.get('MemAvailable')}
    except (OSError, ValueError, KeyError, AttributeError):
        pass
    return {'total_bytes': None, 'available_bytes': None}

def resources(root) -> dict:
    mem = memory_info(); total = mem['total_bytes']; available = mem['available_bytes']
    disk = shutil.disk_usage(root).free
    reserve = max(2 * GIB, int((total or 0) * .05))
    critical = max(512 * 1024**2, int((total or 0) * .01))
    reasons = []
    if disk < 2 * GIB:
        reasons.append('Less than 2 GB free disk; waiting before starting another case.')
    if available is not None and available < reserve:
        reasons.append('Available RAM is below the reserve; waiting before starting another case.')
    emergencies = []
    if disk < 512 * 1024**2:
        emergencies.append('Critically low disk space (below 512 MB).')
    if available is not None and available < critical:
        emergencies.append('Critically low available RAM; requesting a safe pause.')
    return {'logical_cpus': os.cpu_count(), 'memory_total_GB': total/GIB if total else None,
            'memory_available_GB': available/GIB if available is not None else None,
            'memory_reserve_GB': reserve/GIB, 'free_disk_GB': disk/GIB,
            'can_launch': not reasons, 'launch_hold_reason': ' '.join(reasons),
            'critical': bool(emergencies), 'critical_reason': ' '.join(emergencies),
            'memory_monitor_available': available is not None}

class CasePause:
    """Picklable combination of a queue event and an individual case event."""
    def __init__(self, global_event, case_event):
        self.global_event = global_event; self.case_event = case_event
    def is_set(self):
        return self.global_event.is_set() or self.case_event.is_set()
    def set(self):
        # Worker resource emergencies apply to the entire queue.
        self.global_event.set()
