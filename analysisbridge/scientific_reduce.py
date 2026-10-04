#!/usr/bin/env python3
"""Streaming scientific reduction for large CSV/TSV and gzip-compressed tabular outputs.

The reducer intentionally uses only Python's standard library. It never expands a .csv.gz
into a giant temporary CSV. Instead it streams rows and writes compact derivatives that
are suitable for inspection and downstream analysis.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import heapq
import json
import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


META_CANDIDATES = [
    "run_id", "case_id", "simulation_id", "model", "scenario", "period", "calendar",
    "material", "site", "location", "experiment", "member", "variant", "status",
]
TIME_CANDIDATES = [
    "datetime", "timestamp", "time", "date_time", "datetime_utc", "timestamp_utc", "date",
]


def _safe_name(text: str, max_len: int = 120) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("._") or "table"
    return text[:max_len]


def _open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", errors="replace", newline="")
    return path.open("r", encoding="utf-8-sig", errors="replace", newline="")


def _delimiter_for(path: Path, sample: str) -> str:
    lower = path.name.lower()
    if lower.endswith(".tsv") or lower.endswith(".tsv.gz"):
        return "\t"
    try:
        return csv.Sniffer().sniff(sample, delimiters=",;\t|").delimiter
    except Exception:
        return ","


def _float_or_none(v: str) -> Optional[float]:
    if v is None:
        return None
    s = str(v).strip()
    if not s or s.lower() in {"nan", "na", "none", "null", "inf", "+inf", "-inf"}:
        return None
    try:
        x = float(s)
    except Exception:
        return None
    return x if math.isfinite(x) else None


def _parse_dt_value(s: str) -> Optional[datetime]:
    if not s:
        return None
    s = str(s).strip()
    if not s:
        return None
    if s.endswith("Z"):
        s = s[:-1] + "+00:00"
    try:
        return datetime.fromisoformat(s)
    except Exception:
        pass
    for fmt in (
        "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d",
        "%Y/%m/%d %H:%M:%S", "%Y/%m/%d %H:%M", "%Y/%m/%d",
        "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M", "%d/%m/%Y",
        "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M", "%m/%d/%Y",
    ):
        try:
            return datetime.strptime(s, fmt)
        except Exception:
            continue
    return None


def _season(month: int) -> str:
    if month in (12, 1, 2):
        return "DJF"
    if month in (3, 4, 5):
        return "MAM"
    if month in (6, 7, 8):
        return "JJA"
    return "SON"


def _percentile(vals: List[float], p: float) -> Optional[float]:
    if not vals:
        return None
    a = sorted(vals)
    if len(a) == 1:
        return a[0]
    pos = (len(a) - 1) * p
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return a[lo]
    f = pos - lo
    return a[lo] * (1 - f) + a[hi] * f


@dataclass
class Stat:
    n: int = 0
    mean: float = 0.0
    m2: float = 0.0
    min: float = math.inf
    max: float = -math.inf

    def add(self, x: float):
        self.n += 1
        d = x - self.mean
        self.mean += d / self.n
        self.m2 += d * (x - self.mean)
        if x < self.min:
            self.min = x
        if x > self.max:
            self.max = x

    def sd(self) -> Optional[float]:
        return math.sqrt(self.m2 / (self.n - 1)) if self.n > 1 else (0.0 if self.n == 1 else None)

    def as_dict(self) -> Dict[str, Optional[float]]:
        return {
            "n": self.n,
            "mean": self.mean if self.n else None,
            "sd": self.sd(),
            "min": self.min if self.n else None,
            "max": self.max if self.n else None,
        }


class Reservoir:
    """Small deterministic reservoir for approximate quantiles."""
    __slots__ = ("k", "values", "seen", "state")

    def __init__(self, k: int = 4096, seed: int = 1):
        self.k = k
        self.values: List[float] = []
        self.seen = 0
        self.state = seed & 0x7FFFFFFF or 1

    def add(self, x: float):
        self.seen += 1
        if len(self.values) < self.k:
            self.values.append(x)
            return
        self.state = (1103515245 * self.state + 12345) & 0x7FFFFFFF
        j = self.state % self.seen
        if j < self.k:
            self.values[j] = x


@dataclass
class DailyAgg:
    sums: Dict[str, float] = field(default_factory=lambda: defaultdict(float))
    counts: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    mins: Dict[str, float] = field(default_factory=dict)
    maxs: Dict[str, float] = field(default_factory=dict)

    def add(self, col: str, x: float):
        self.sums[col] += x
        self.counts[col] += 1
        self.mins[col] = x if col not in self.mins else min(self.mins[col], x)
        self.maxs[col] = x if col not in self.maxs else max(self.maxs[col], x)


class SpellTop:
    """Keep only the N longest/severest threshold spells."""
    def __init__(self, n_keep: int = 50):
        self.n_keep = n_keep
        self.heap = []
        self.counter = 0

    def push(self, rec: Dict):
        score = (float(rec.get("duration_rows", 0)), abs(float(rec.get("extreme_value", 0.0))))
        self.counter += 1
        item = (score, self.counter, rec)
        if len(self.heap) < self.n_keep:
            heapq.heappush(self.heap, item)
        elif score > self.heap[0][0]:
            heapq.heapreplace(self.heap, item)

    def rows(self) -> List[Dict]:
        return [x[2] for x in sorted(self.heap, key=lambda z: (z[0], z[1]), reverse=True)]


class ScientificReducer:
    def __init__(self, output_root: Path, log_cb=None, top_events: int = 20, max_key_metrics: int = 12):
        self.root = output_root / "SCIENTIFIC_REDUCTION"
        self.root.mkdir(parents=True, exist_ok=True)
        self.log_cb = log_cb or (lambda msg: None)
        self.top_events = max(5, int(top_events))
        self.max_key_metrics = max(4, int(max_key_metrics))

    def log(self, msg: str):
        self.log_cb(msg)

    @staticmethod
    def eligible(path: Path) -> bool:
        n = path.name.lower()
        return n.endswith(".csv.gz") or n.endswith(".tsv.gz")

    def _infer(self, headers: List[str], sample_rows: List[Dict[str, str]]) -> Dict:
        lower = {h.lower().strip(): h for h in headers}
        time_col = None
        for c in TIME_CANDIDATES:
            if c in lower:
                time_col = lower[c]
                break
        if time_col is None:
            for h in headers:
                hl = h.lower()
                if "time" in hl or "date" in hl:
                    time_col = h
                    break

        run_col = None
        for c in ("run_id", "case_id", "simulation_id", "run", "case"):
            if c in lower:
                run_col = lower[c]
                break

        meta_cols = []
        for c in META_CANDIDATES:
            if c in lower and lower[c] not in meta_cols:
                meta_cols.append(lower[c])
        if run_col in meta_cols:
            meta_cols.remove(run_col)

        ratios = {}
        for h in headers:
            if h == time_col or h in meta_cols:
                continue
            present = 0
            numeric = 0
            for r in sample_rows:
                v = (r.get(h) or "").strip()
                if not v:
                    continue
                present += 1
                if _float_or_none(v) is not None:
                    numeric += 1
            ratios[h] = numeric / present if present else 0.0
        axis_names = {"elapsed_hours", "elapsed_hour", "hour_index", "hours_since_start", "simulation_hour", "year", "month", "day", "hour", "step", "time_index", "index"}
        numeric_cols = [h for h in headers if ratios.get(h, 0) >= 0.80 and h.lower().strip() not in axis_names]

        def score(col: str) -> Tuple[int, int]:
            s = col.lower()
            points = 0
            if any(x in s for x in ("temp", "temperature", "t_mean", "t_max", "t_min")):
                points += 100
            if any(x in s for x in ("theta", "moist", "water", "humidity", "rh")):
                points += 90
            if any(x in s for x in ("sulf", "salt", "concentration", "conc")):
                points += 85
            if any(x in s for x in ("rain", "precip")):
                points += 80
            if any(x in s for x in ("ice", "freeze", "phase")):
                points += 75
            if any(x in s for x in ("evap", "seep", "flux")):
                points += 60
            return (points, -headers.index(col))

        key_numeric = sorted(numeric_cols, key=score, reverse=True)[: self.max_key_metrics]
        if len(numeric_cols) <= self.max_key_metrics:
            key_numeric = numeric_cols

        temp_cols = [c for c in key_numeric if any(x in c.lower() for x in ("temp", "temperature", "t_mean", "t_max", "t_min"))]
        moisture_cols = [c for c in key_numeric if any(x in c.lower() for x in ("theta", "moist", "water_content")) and not any(x in c.lower() for x in ("ice", "freeze"))]
        sulfate_cols = [c for c in key_numeric if any(x in c.lower() for x in ("sulf", "salt", "concentration", "conc"))]
        rain_cols = [c for c in key_numeric if any(x in c.lower() for x in ("rain", "precip"))]
        ice_cols = [c for c in key_numeric if any(x in c.lower() for x in ("ice", "freeze"))]
        elapsed_col = None
        for c in ("elapsed_hours", "elapsed_hour", "hour_index", "hours_since_start", "simulation_hour"):
            if c in lower:
                elapsed_col = lower[c]
                break

        return {
            "time_column": time_col,
            "elapsed_hours_column": elapsed_col,
            "run_column": run_col,
            "metadata_columns": meta_cols,
            "numeric_columns": numeric_cols,
            "key_numeric_columns": key_numeric,
            "numeric_detection_ratio": ratios,
            "temperature_columns": temp_cols,
            "moisture_columns": moisture_cols,
            "sulfate_or_salt_columns": sulfate_cols,
            "rain_or_precip_columns": rain_cols,
            "ice_or_freeze_columns": ice_cols,
        }

    def _run_id(self, row: Dict[str, str], infer: Dict) -> str:
        rc = infer["run_column"]
        if rc and str(row.get(rc, "")).strip():
            return str(row[rc]).strip()
        pieces = []
        for c in infer["metadata_columns"]:
            if c.lower() in {"model", "scenario", "period"} and str(row.get(c, "")).strip():
                pieces.append(str(row[c]).strip())
        return "|".join(pieces) if pieces else "ALL"

    def _row_dt(self, row: Dict[str, str], infer: Dict) -> Optional[datetime]:
        tc = infer["time_column"]
        if tc:
            dt = _parse_dt_value(row.get(tc, ""))
            if dt:
                return dt
        lower = {k.lower(): k for k in row}
        if all(x in lower for x in ("year", "month", "day")):
            try:
                hour = int(float(row.get(lower.get("hour", ""), 0) or 0)) if "hour" in lower else 0
                return datetime(int(float(row[lower["year"]])), int(float(row[lower["month"]])), int(float(row[lower["day"]])), hour)
            except Exception:
                pass
        ec = infer.get("elapsed_hours_column")
        if ec:
            elapsed = _float_or_none(row.get(ec, ""))
            period = str(row.get(lower.get("period", ""), "") or "") if "period" in lower else ""
            m = re.search(r"(?:19|20|21)\d{2}", period)
            if elapsed is not None and m:
                try:
                    return datetime(int(m.group(0)), 1, 1) + timedelta(hours=elapsed)
                except Exception:
                    pass
        return None

    def _threshold_defs(self, infer: Dict) -> List[Tuple[str, str, float]]:
        defs = []
        for c in infer["temperature_columns"]:
            cl = c.lower()
            if ("_c" in cl or "celsius" in cl or cl.endswith("c")) and "_k" not in cl and "kelvin" not in cl:
                defs += [(c, "lt", 0.0), (c, "gt", 40.0), (c, "gt", 50.0), (c, "gt", 60.0)]
        for c in infer["ice_or_freeze_columns"]:
            defs += [(c, "gt", 0.0)]
        for c in infer["rain_or_precip_columns"]:
            defs += [(c, "gt", 0.0), (c, "ge", 1.0), (c, "ge", 10.0)]
        return defs

    @staticmethod
    def _cond(x: float, op: str, threshold: float) -> bool:
        if op == "lt": return x < threshold
        if op == "le": return x <= threshold
        if op == "gt": return x > threshold
        return x >= threshold

    def reduce(self, path: Path, relative_path: str = "") -> Tuple[List[Path], Dict]:
        self.log(f"  Scientific reduction: {path.name}")
        base = _safe_name(Path(path.name).name.replace(".csv.gz", "").replace(".tsv.gz", ""))
        out_dir = self.root / base
        suffix = 1
        while out_dir.exists() and any(out_dir.iterdir()):
            suffix += 1
            out_dir = self.root / f"{base}_{suffix}"
        out_dir.mkdir(parents=True, exist_ok=True)

        # Sample and infer schema.
        with _open_text(path) as f:
            sample_text = f.read(65536)
        delim = _delimiter_for(path, sample_text)
        sample_rows: List[Dict[str, str]] = []
        headers: List[str] = []
        with _open_text(path) as f:
            reader = csv.DictReader(f, delimiter=delim)
            headers = reader.fieldnames or []
            for i, r in enumerate(reader):
                sample_rows.append(r)
                if i >= 199:
                    break
        if not headers:
            raise RuntimeError("No CSV header found in compressed table")
        infer = self._infer(headers, sample_rows)

        source_lower = path.name.lower()
        kind_prefix = "HOURLY" if "hour" in source_lower else ("DAILY" if "daily" in source_lower else ("PHASE" if any(x in source_lower for x in ("phase", "cold", "freeze")) else "TABLE"))
        schema_path = out_dir / f"{kind_prefix}_SCHEMA.json"
        schema_obj = {
            "source_relative_path": relative_path,
            "source_file_name": path.name,
            "compressed_size_bytes": path.stat().st_size,
            "delimiter": "TAB" if delim == "\t" else delim,
            "headers": headers,
            **infer,
            "notes": [
                "Numeric roles are inferred automatically from the first 200 non-header rows.",
                "Quantiles are approximate, based on deterministic bounded reservoirs.",
                "Temperature thresholds are expressed in the source column's units; they are intended for Celsius-like model outputs.",
                "Rain/precipitation thresholds are only emitted when a rain/precip column is detected; units are inherited from the source.",
                "If no timestamp/date is present but elapsed_hours and a YYYY-YYYY period are available, calendar time is reconstructed from January 1 of the first period year. This is an analysis convenience and may not exactly represent non-Gregorian calendars.",
            ],
        }
        schema_path.write_text(json.dumps(schema_obj, indent=2, ensure_ascii=False), encoding="utf-8")

        run_stats: Dict[str, Dict[str, Stat]] = defaultdict(lambda: defaultdict(Stat))
        reservoirs: Dict[str, Dict[str, Reservoir]] = defaultdict(dict)
        metadata: Dict[str, Dict[str, str]] = defaultdict(dict)
        month_stats: Dict[Tuple[str, str], Dict[str, Stat]] = defaultdict(lambda: defaultdict(Stat))
        season_stats: Dict[Tuple[str, str], Dict[str, Stat]] = defaultdict(lambda: defaultdict(Stat))
        top_hi: Dict[Tuple[str, str], List[Tuple[float, str]]] = defaultdict(list)
        top_lo: Dict[Tuple[str, str], List[Tuple[float, str]]] = defaultdict(list)
        threshold_counts = defaultdict(int)
        threshold_defs = self._threshold_defs(infer)
        spell_state: Dict[Tuple[str, str, str, float], Dict] = {}
        spell_top: Dict[Tuple[str, str, str, float], SpellTop] = defaultdict(lambda: SpellTop(50))

        key = infer["key_numeric_columns"]
        rain_set = set(infer["rain_or_precip_columns"])
        daily_path = out_dir / f"{kind_prefix}_DAILY_REDUCTION.csv"
        daily_headers = ["run_id"] + infer["metadata_columns"] + ["date"]
        for c in key:
            if c in rain_set:
                daily_headers += [f"{c}__sum", f"{c}__mean", f"{c}__max"]
            else:
                daily_headers += [f"{c}__mean", f"{c}__min", f"{c}__max"]

        current_daily: Dict[str, Tuple[str, DailyAgg]] = {}
        row_count = 0
        parsed_time_rows = 0
        first_dt: Dict[str, datetime] = {}
        last_dt: Dict[str, datetime] = {}
        run_rows = defaultdict(int)

        def flush_daily(writer, run: str):
            item = current_daily.get(run)
            if not item:
                return
            day, agg = item
            rec = {"run_id": run, "date": day}
            rec.update(metadata.get(run, {}))
            for c in key:
                n = agg.counts.get(c, 0)
                if c in rain_set:
                    rec[f"{c}__sum"] = agg.sums.get(c, 0.0) if n else ""
                    rec[f"{c}__mean"] = agg.sums.get(c, 0.0) / n if n else ""
                    rec[f"{c}__max"] = agg.maxs.get(c, "") if n else ""
                else:
                    rec[f"{c}__mean"] = agg.sums.get(c, 0.0) / n if n else ""
                    rec[f"{c}__min"] = agg.mins.get(c, "") if n else ""
                    rec[f"{c}__max"] = agg.maxs.get(c, "") if n else ""
            writer.writerow(rec)

        with daily_path.open("w", encoding="utf-8-sig", newline="") as df:
            dw = csv.DictWriter(df, fieldnames=daily_headers, extrasaction="ignore")
            dw.writeheader()
            with _open_text(path) as f:
                reader = csv.DictReader(f, delimiter=delim)
                for row in reader:
                    row_count += 1
                    run = self._run_id(row, infer)
                    run_rows[run] += 1
                    for c in infer["metadata_columns"]:
                        v = str(row.get(c, "") or "").strip()
                        if v and c not in metadata[run]:
                            metadata[run][c] = v
                    dt = self._row_dt(row, infer)
                    if dt:
                        parsed_time_rows += 1
                        first_dt[run] = dt if run not in first_dt else min(first_dt[run], dt)
                        last_dt[run] = dt if run not in last_dt else max(last_dt[run], dt)
                        day = dt.date().isoformat()
                        if run in current_daily and current_daily[run][0] != day:
                            flush_daily(dw, run)
                            current_daily[run] = (day, DailyAgg())
                        elif run not in current_daily:
                            current_daily[run] = (day, DailyAgg())

                    timestamp_text = dt.isoformat() if dt else str(row.get(infer["time_column"] or "", ""))
                    vals = {}
                    for c in key:
                        x = _float_or_none(row.get(c, ""))
                        if x is None:
                            continue
                        vals[c] = x
                        run_stats[run][c].add(x)
                        if c not in reservoirs[run]:
                            seed = int(hashlib.sha1(f"{run}|{c}".encode()).hexdigest()[:8], 16)
                            reservoirs[run][c] = Reservoir(4096, seed)
                        reservoirs[run][c].add(x)
                        if dt:
                            month_stats[(run, dt.strftime("%Y-%m"))][c].add(x)
                            season_stats[(run, f"{dt.year}-{_season(dt.month)}")][c].add(x)
                            current_daily[run][1].add(c, x)

                        hi = top_hi[(run, c)]
                        item = (x, timestamp_text)
                        if len(hi) < self.top_events:
                            heapq.heappush(hi, item)
                        elif x > hi[0][0]:
                            heapq.heapreplace(hi, item)
                        lo = top_lo[(run, c)]
                        item2 = (-x, timestamp_text)
                        if len(lo) < self.top_events:
                            heapq.heappush(lo, item2)
                        elif -x > lo[0][0]:
                            heapq.heapreplace(lo, item2)

                    for c, op, threshold in threshold_defs:
                        x = vals.get(c)
                        if x is None:
                            continue
                        k = (run, c, op, threshold)
                        active = self._cond(x, op, threshold)
                        if active:
                            threshold_counts[k] += 1
                        st = spell_state.get(k)
                        if active:
                            if st is None:
                                spell_state[k] = {"start": timestamp_text, "end": timestamp_text, "duration_rows": 1, "extreme_value": x, "last_dt": dt}
                            else:
                                gap_break = bool(dt and st.get("last_dt") and abs((dt - st["last_dt"]).total_seconds()) > 5400)
                                if gap_break:
                                    rec = {"run_id": run, "metric": c, "operator": op, "threshold": threshold, **{kk: vv for kk, vv in st.items() if kk != "last_dt"}}
                                    spell_top[k].push(rec)
                                    spell_state[k] = {"start": timestamp_text, "end": timestamp_text, "duration_rows": 1, "extreme_value": x, "last_dt": dt}
                                else:
                                    st["end"] = timestamp_text
                                    st["duration_rows"] += 1
                                    st["last_dt"] = dt
                                    if op in ("gt", "ge"):
                                        st["extreme_value"] = max(st["extreme_value"], x)
                                    else:
                                        st["extreme_value"] = min(st["extreme_value"], x)
                        elif st is not None:
                            rec = {"run_id": run, "metric": c, "operator": op, "threshold": threshold, **{kk: vv for kk, vv in st.items() if kk != "last_dt"}}
                            spell_top[k].push(rec)
                            spell_state.pop(k, None)

                    if row_count % 500000 == 0:
                        self.log(f"    {row_count:,} rows streamed from {path.name}")
            for run in list(current_daily):
                flush_daily(dw, run)

        for k, st in list(spell_state.items()):
            run, c, op, threshold = k
            rec = {"run_id": run, "metric": c, "operator": op, "threshold": threshold, **{kk: vv for kk, vv in st.items() if kk != "last_dt"}}
            spell_top[k].push(rec)

        outputs = [schema_path, daily_path]

        # Run-level stats and approximate quantiles.
        run_path = out_dir / f"{kind_prefix}_RUN_STATS.csv"
        run_fields = ["run_id"] + infer["metadata_columns"] + ["rows", "start_time", "end_time", "duration_days"]
        for c in key:
            run_fields += [f"{c}__n", f"{c}__mean", f"{c}__sd", f"{c}__min", f"{c}__p01", f"{c}__p05", f"{c}__p50", f"{c}__p95", f"{c}__p99", f"{c}__max"]
        with run_path.open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=run_fields)
            w.writeheader()
            for run in sorted(run_rows):
                rec = {"run_id": run, "rows": run_rows[run], **metadata.get(run, {})}
                if run in first_dt:
                    rec["start_time"] = first_dt[run].isoformat()
                    rec["end_time"] = last_dt[run].isoformat()
                    rec["duration_days"] = (last_dt[run] - first_dt[run]).total_seconds() / 86400.0
                for c in key:
                    st = run_stats[run].get(c)
                    if not st or not st.n:
                        continue
                    vals = reservoirs[run][c].values
                    rec[f"{c}__n"] = st.n
                    rec[f"{c}__mean"] = st.mean
                    rec[f"{c}__sd"] = st.sd()
                    rec[f"{c}__min"] = st.min
                    rec[f"{c}__p01"] = _percentile(vals, .01)
                    rec[f"{c}__p05"] = _percentile(vals, .05)
                    rec[f"{c}__p50"] = _percentile(vals, .50)
                    rec[f"{c}__p95"] = _percentile(vals, .95)
                    rec[f"{c}__p99"] = _percentile(vals, .99)
                    rec[f"{c}__max"] = st.max
                w.writerow(rec)
        outputs.append(run_path)

        # Monthly and seasonal summaries.
        def write_group_stats(target: Path, groups: Dict, period_name: str):
            fields = ["run_id"] + infer["metadata_columns"] + [period_name]
            for c in key:
                fields += [f"{c}__n", f"{c}__mean", f"{c}__min", f"{c}__max"]
            with target.open("w", encoding="utf-8-sig", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields)
                w.writeheader()
                for (run, per), stats in sorted(groups.items()):
                    rec = {"run_id": run, period_name: per, **metadata.get(run, {})}
                    for c in key:
                        st = stats.get(c)
                        if st and st.n:
                            rec[f"{c}__n"] = st.n
                            rec[f"{c}__mean"] = st.mean
                            rec[f"{c}__min"] = st.min
                            rec[f"{c}__max"] = st.max
                    w.writerow(rec)
            outputs.append(target)

        monthly_path = out_dir / f"{kind_prefix}_MONTHLY_STATS.csv"
        seasonal_path = out_dir / f"{kind_prefix}_SEASONAL_STATS.csv"
        write_group_stats(monthly_path, month_stats, "year_month")
        write_group_stats(seasonal_path, season_stats, "year_season")

        # Exact extreme rows.
        ext_path = out_dir / f"{kind_prefix}_EXTREME_EVENTS.csv"
        ext_fields = ["run_id"] + infer["metadata_columns"] + ["metric", "direction", "rank", "timestamp", "value"]
        with ext_path.open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=ext_fields)
            w.writeheader()
            for (run, c), heap in sorted(top_hi.items()):
                vals = sorted(heap, reverse=True)
                for rank, (x, ts) in enumerate(vals, 1):
                    w.writerow({"run_id": run, **metadata.get(run, {}), "metric": c, "direction": "maximum", "rank": rank, "timestamp": ts, "value": x})
            for (run, c), heap in sorted(top_lo.items()):
                vals = sorted([(-neg, ts) for neg, ts in heap])
                for rank, (x, ts) in enumerate(vals, 1):
                    w.writerow({"run_id": run, **metadata.get(run, {}), "metric": c, "direction": "minimum", "rank": rank, "timestamp": ts, "value": x})
        outputs.append(ext_path)

        threshold_path = out_dir / f"{kind_prefix}_THRESHOLD_COUNTS.csv"
        with threshold_path.open("w", encoding="utf-8-sig", newline="") as f:
            fields = ["run_id"] + infer["metadata_columns"] + ["metric", "operator", "threshold", "hours_or_rows", "fraction_of_run_rows"]
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for k, n in sorted(threshold_counts.items()):
                run, c, op, threshold = k
                w.writerow({"run_id": run, **metadata.get(run, {}), "metric": c, "operator": op, "threshold": threshold, "hours_or_rows": n, "fraction_of_run_rows": n / run_rows[run] if run_rows[run] else ""})
        outputs.append(threshold_path)

        spells_path = out_dir / f"{kind_prefix}_TOP_THRESHOLD_SPELLS.csv"
        with spells_path.open("w", encoding="utf-8-sig", newline="") as f:
            fields = ["run_id"] + infer["metadata_columns"] + ["metric", "operator", "threshold", "start", "end", "duration_rows", "extreme_value"]
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for k, store in sorted(spell_top.items()):
                for rec in store.rows():
                    out = {**rec, **metadata.get(rec["run_id"], {})}
                    w.writerow(out)
        outputs.append(spells_path)

        # Post-process the compact daily file for wetting/drying cycles, rain response, and duration-matched summaries.
        post = self._postprocess_daily(daily_path, infer, metadata, out_dir)
        outputs.extend(post)

        report = {
            "source_relative_path": relative_path,
            "source_file_name": path.name,
            "rows_streamed": row_count,
            "runs_detected": len(run_rows),
            "rows_with_parsed_time": parsed_time_rows,
            "time_parse_fraction": parsed_time_rows / row_count if row_count else 0.0,
            "key_numeric_columns": key,
            "scientific_output_files": [str(p.relative_to(self.root.parent)) for p in outputs],
            "status": "PASS" if row_count else "EMPTY",
        }
        report_path = out_dir / "SCIENTIFIC_REDUCTION_REPORT.json"
        report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        outputs.append(report_path)
        self.log(f"    Completed scientific reduction: {row_count:,} rows, {len(run_rows)} run(s)")
        return outputs, report

    def _postprocess_daily(self, daily_path: Path, infer: Dict, metadata: Dict[str, Dict[str, str]], out_dir: Path) -> List[Path]:
        by_run: Dict[str, List[Dict[str, str]]] = defaultdict(list)
        with daily_path.open("r", encoding="utf-8-sig", newline="") as f:
            for r in csv.DictReader(f):
                by_run[r["run_id"]].append(r)
        for run in by_run:
            by_run[run].sort(key=lambda r: r["date"])

        outputs = []
        moisture = infer["moisture_columns"][:2]
        sulfate = infer["sulfate_or_salt_columns"][:2]
        temp = infer["temperature_columns"][:2]
        rain = infer["rain_or_precip_columns"][:1]

        cycle_path = out_dir / "WETTING_DRYING_CYCLES.csv"
        with cycle_path.open("w", encoding="utf-8-sig", newline="") as f:
            fields = ["run_id"] + infer["metadata_columns"] + ["metric", "daily_points", "amplitude_threshold", "wetting_days", "drying_days", "reversal_cycles", "max_1d_increase", "max_1d_decrease"]
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for run, rows in sorted(by_run.items()):
                for c in moisture:
                    vals = []
                    for r in rows:
                        x = _float_or_none(r.get(f"{c}__mean", ""))
                        if x is not None:
                            vals.append(x)
                    if len(vals) < 2:
                        continue
                    rng = max(vals) - min(vals)
                    eps = max(1e-9, 0.02 * rng)
                    diffs = [vals[i] - vals[i-1] for i in range(1, len(vals))]
                    signs = [1 if d > eps else -1 if d < -eps else 0 for d in diffs]
                    nz = [s for s in signs if s]
                    reversals = sum(1 for a, b in zip(nz, nz[1:]) if a != b)
                    w.writerow({"run_id": run, **metadata.get(run, {}), "metric": c, "daily_points": len(vals), "amplitude_threshold": eps, "wetting_days": sum(1 for s in signs if s > 0), "drying_days": sum(1 for s in signs if s < 0), "reversal_cycles": reversals, "max_1d_increase": max(diffs), "max_1d_decrease": min(diffs)})
        outputs.append(cycle_path)

        rain_path = out_dir / "RAIN_RESPONSE_EVENTS.csv"
        with rain_path.open("w", encoding="utf-8-sig", newline="") as f:
            fields = ["run_id"] + infer["metadata_columns"] + ["rain_metric", "event_date", "rain_total", "response_metric", "baseline", "plus_1d", "plus_3d", "plus_7d", "delta_1d", "delta_3d", "delta_7d"]
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            if rain:
                rc = rain[0]
                response_cols = moisture + sulfate + temp
                for run, rows in sorted(by_run.items()):
                    candidates = []
                    for i, r in enumerate(rows):
                        p = _float_or_none(r.get(f"{rc}__sum", ""))
                        if p is not None and p > 0:
                            candidates.append((p, i))
                    for p, i in sorted(candidates, reverse=True)[:25]:
                        for c in response_cols:
                            base = _float_or_none(rows[i].get(f"{c}__mean", ""))
                            vals = {}
                            for lag in (1, 3, 7):
                                vals[lag] = _float_or_none(rows[i+lag].get(f"{c}__mean", "")) if i + lag < len(rows) else None
                            w.writerow({
                                "run_id": run, **metadata.get(run, {}), "rain_metric": rc, "event_date": rows[i]["date"], "rain_total": p,
                                "response_metric": c, "baseline": base,
                                "plus_1d": vals[1], "plus_3d": vals[3], "plus_7d": vals[7],
                                "delta_1d": vals[1] - base if vals[1] is not None and base is not None else "",
                                "delta_3d": vals[3] - base if vals[3] is not None and base is not None else "",
                                "delta_7d": vals[7] - base if vals[7] is not None and base is not None else "",
                            })
        outputs.append(rain_path)

        # Duration matching: historical/ERA baseline when present, otherwise shortest run.
        baseline_run = None
        for run in sorted(by_run):
            m = {k.lower(): str(v).lower() for k, v in metadata.get(run, {}).items()}
            if "historical" in m.get("scenario", "") or "era" in m.get("model", ""):
                baseline_run = run
                break
        if baseline_run is None and by_run:
            baseline_run = min(by_run, key=lambda r: len(by_run[r]))
        baseline_days = len(by_run.get(baseline_run, [])) if baseline_run else 0
        matched_metrics = (temp + moisture + sulfate)[:6]
        matched_path = out_dir / "MATCHED_DURATION_RUN_SUMMARY.csv"
        matched_rows = []
        fields = ["run_id"] + infer["metadata_columns"] + ["baseline_run", "matched_days"]
        for c in matched_metrics:
            fields += [f"{c}__mean", f"{c}__min", f"{c}__max"]
        with matched_path.open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for run, rows in sorted(by_run.items()):
                n = min(baseline_days, len(rows)) if baseline_days else len(rows)
                rec = {"run_id": run, **metadata.get(run, {}), "baseline_run": baseline_run or "", "matched_days": n}
                for c in matched_metrics:
                    mean_vals = [_float_or_none(r.get(f"{c}__mean", "")) for r in rows[:n]]
                    mean_vals = [x for x in mean_vals if x is not None]
                    min_vals = [_float_or_none(r.get(f"{c}__min", r.get(f"{c}__mean", ""))) for r in rows[:n]]
                    min_vals = [x for x in min_vals if x is not None]
                    max_vals = [_float_or_none(r.get(f"{c}__max", r.get(f"{c}__mean", ""))) for r in rows[:n]]
                    max_vals = [x for x in max_vals if x is not None]
                    if mean_vals:
                        rec[f"{c}__mean"] = sum(mean_vals) / len(mean_vals)
                    if min_vals:
                        rec[f"{c}__min"] = min(min_vals)
                    if max_vals:
                        rec[f"{c}__max"] = max(max_vals)
                w.writerow(rec)
                matched_rows.append(rec)
        outputs.append(matched_path)

        ensemble_path = out_dir / "MATCHED_DURATION_ENSEMBLE_SUMMARY.csv"
        group_cols = [c for c in infer["metadata_columns"] if c.lower() in {"scenario", "period"}]
        if not group_cols:
            group_cols = [c for c in infer["metadata_columns"] if c.lower() == "model"][:1]
        groups = defaultdict(list)
        for r in matched_rows:
            groups[tuple(str(r.get(c, "")) for c in group_cols)].append(r)
        efields = group_cols + ["model_count_or_run_count", "baseline_run", "matched_days"]
        for c in matched_metrics:
            efields += [f"{c}__ensemble_mean", f"{c}__ensemble_min", f"{c}__ensemble_max"]
        with ensemble_path.open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=efields)
            w.writeheader()
            for g, rows in sorted(groups.items()):
                rec = {c: g[i] for i, c in enumerate(group_cols)}
                rec["model_count_or_run_count"] = len(rows)
                rec["baseline_run"] = baseline_run or ""
                rec["matched_days"] = min([int(float(r.get("matched_days", 0) or 0)) for r in rows] or [0])
                for c in matched_metrics:
                    vals = [_float_or_none(r.get(f"{c}__mean", "")) for r in rows]
                    vals = [x for x in vals if x is not None]
                    if vals:
                        rec[f"{c}__ensemble_mean"] = sum(vals) / len(vals)
                        rec[f"{c}__ensemble_min"] = min(vals)
                        rec[f"{c}__ensemble_max"] = max(vals)
                w.writerow(rec)
        outputs.append(ensemble_path)

        change_path = out_dir / "MATCHED_DURATION_CHANGE_VS_BASELINE.csv"
        change_fields = ["run_id"] + infer["metadata_columns"] + ["baseline_run", "matched_days"]
        for c in matched_metrics:
            change_fields += [f"{c}__mean", f"{c}__mean_abs_change", f"{c}__mean_pct_change", f"{c}__max", f"{c}__max_abs_change", f"{c}__max_pct_change"]
        base_rec = next((r for r in matched_rows if r.get("run_id") == baseline_run), None)
        with change_path.open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=change_fields)
            w.writeheader()
            for rec0 in matched_rows:
                rec = {"run_id": rec0.get("run_id", ""), **metadata.get(rec0.get("run_id", ""), {}), "baseline_run": baseline_run or "", "matched_days": rec0.get("matched_days", "")}
                for c in matched_metrics:
                    for stat in ("mean", "max"):
                        v = _float_or_none(rec0.get(f"{c}__{stat}", ""))
                        b = _float_or_none(base_rec.get(f"{c}__{stat}", "")) if base_rec else None
                        rec[f"{c}__{stat}"] = v if v is not None else ""
                        if v is not None and b is not None:
                            rec[f"{c}__{stat}_abs_change"] = v - b
                            rec[f"{c}__{stat}_pct_change"] = ((v - b) / abs(b) * 100.0) if b != 0 else ""
                w.writerow(rec)
        outputs.append(change_path)
        return outputs
