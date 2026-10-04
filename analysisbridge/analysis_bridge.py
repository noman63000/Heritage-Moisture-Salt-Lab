#!/usr/bin/env python3
"""AnalysisBridge - local analysis-preparation tool for large folders and archives.

Runs fully on the user's computer. It recursively scans folders/archives, validates
inputs, hashes and inventories files, summarizes common structured formats, safely
extracts ZIP files into a workspace, and creates upload-friendly analysis packs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import queue
import shutil
import sqlite3
import sys
import threading
import time
import traceback
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from scientific_reduce import ScientificReducer

APP_NAME = "AnalysisBridge"
APP_VERSION = "1.2.0"
TEXT_EXTS = {".txt", ".log", ".md", ".ini", ".cfg", ".yaml", ".yml", ".xml"}
TABULAR_EXTS = {".csv", ".tsv"}
JSON_EXTS = {".json", ".jsonl", ".ndjson"}
ARCHIVE_EXTS = {".zip"}
DEFAULT_PACK_MB = 180
READ_CHUNK = 4 * 1024 * 1024
MAX_PREVIEW_BYTES = 512 * 1024
MAX_TEXT_PREVIEW_LINES = 80
MAX_CSV_SAMPLE_ROWS = 5000
MAX_JSON_ARRAY_SAMPLE = 5000


def human_bytes(n: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    x = float(n)
    for u in units:
        if x < 1024 or u == units[-1]:
            return f"{x:.2f} {u}"
        x /= 1024
    return f"{n} B"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(READ_CHUNK)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def safe_rel_name(text: str, max_len: int = 140) -> str:
    bad = '<>:"/\\|?*\x00'
    out = "".join("_" if c in bad else c for c in text).strip().strip(".")
    if not out:
        out = "item"
    return out[:max_len]


def guess_encoding(path: Path) -> str:
    raw = path.read_bytes()[:65536]
    for enc in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            raw.decode(enc)
            return enc
        except UnicodeDecodeError:
            pass
    return "latin-1"


def detect_delimiter(sample: str, ext: str) -> str:
    if ext == ".tsv":
        return "\t"
    try:
        return csv.Sniffer().sniff(sample, delimiters=",;\t|").delimiter
    except Exception:
        return ","


@dataclass
class FileRecord:
    source_root: str
    path: str
    relative_path: str
    extension: str
    size_bytes: int
    sha256: str
    category: str
    status: str = "ok"
    notes: str = ""
    duplicate_of: str = ""


class AnalysisBridgeEngine:
    def __init__(self, output_dir: Path, pack_mb: int = DEFAULT_PACK_MB, mode: str = "fast", log_cb=None, progress_cb=None, scientific: bool = True):
        self.output_dir = output_dir
        self.pack_limit = max(20, pack_mb) * 1024 * 1024
        self.mode = mode
        self.scientific = scientific
        self.log_cb = log_cb or (lambda msg: None)
        self.progress_cb = progress_cb or (lambda a, b, msg="": None)
        self.workspace = output_dir / "workspace"
        self.extracted = self.workspace / "extracted"
        self.summaries = self.output_dir / "summaries"
        self.previews = self.output_dir / "previews"
        self.upload_dir = self.output_dir / "CHATGPT_UPLOAD"
        self.records: List[FileRecord] = []
        self.issues: List[Dict[str, str]] = []
        self.hash_first: Dict[str, str] = {}
        self.summary_files: List[Path] = []
        self.scientific_files: List[Path] = []
        self.scientific_reports: List[Dict] = []
        self.reducer = None
        self.total_input_bytes = 0
        self.total_processed_bytes = 0

    def log(self, msg: str):
        self.log_cb(msg)

    def setup(self):
        if self.output_dir.exists():
            # Preserve existing output by creating a timestamped run folder.
            if any(self.output_dir.iterdir()):
                ts = time.strftime("%Y%m%d_%H%M%S")
                self.output_dir = self.output_dir.parent / f"{self.output_dir.name}_{ts}"
                self.workspace = self.output_dir / "workspace"
                self.extracted = self.workspace / "extracted"
                self.summaries = self.output_dir / "summaries"
                self.previews = self.output_dir / "previews"
                self.upload_dir = self.output_dir / "CHATGPT_UPLOAD"
        for p in (self.output_dir, self.workspace, self.extracted, self.summaries, self.previews, self.upload_dir):
            p.mkdir(parents=True, exist_ok=True)
        if self.scientific:
            self.reducer = ScientificReducer(self.output_dir, log_cb=self.log)

    def collect_paths(self, inputs: List[Path]) -> List[Tuple[Path, Path]]:
        items: List[Tuple[Path, Path]] = []
        seen = set()
        for inp in inputs:
            inp = inp.resolve()
            if not inp.exists():
                self.issues.append({"path": str(inp), "issue": "missing input"})
                continue
            root = inp if inp.is_dir() else inp.parent
            candidates = inp.rglob("*") if inp.is_dir() else [inp]
            for p in candidates:
                if not p.is_file():
                    continue
                try:
                    if p.resolve().is_relative_to(self.output_dir.resolve()):
                        continue
                except Exception:
                    pass
                key = str(p.resolve()).lower()
                if key in seen:
                    continue
                seen.add(key)
                items.append((root, p))
                try:
                    self.total_input_bytes += p.stat().st_size
                except OSError:
                    pass
        return items

    def category_for(self, path: Path) -> str:
        name = path.name.lower()
        if name.endswith(".csv.gz") or name.endswith(".tsv.gz"):
            return "compressed_tabular"
        ext = path.suffix.lower()
        if ext in TABULAR_EXTS:
            return "tabular"
        if ext in JSON_EXTS:
            return "json"
        if ext in TEXT_EXTS:
            return "text"
        if ext in ARCHIVE_EXTS:
            return "archive"
        return "other"

    def safe_extract_zip(self, zip_path: Path, dest: Path) -> List[Path]:
        extracted_files = []
        dest_resolved = dest.resolve()
        with zipfile.ZipFile(zip_path, "r") as zf:
            for info in zf.infolist():
                if info.is_dir():
                    continue
                target = dest / info.filename
                try:
                    target_resolved = target.resolve()
                except Exception:
                    self.issues.append({"path": str(zip_path), "issue": f"invalid member path: {info.filename}"})
                    continue
                if os.path.commonpath([str(dest_resolved), str(target_resolved)]) != str(dest_resolved):
                    self.issues.append({"path": str(zip_path), "issue": f"blocked unsafe ZIP member: {info.filename}"})
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(info, "r") as src, target.open("wb") as dst:
                    shutil.copyfileobj(src, dst, length=READ_CHUNK)
                extracted_files.append(target)
        return extracted_files

    def preview_text(self, path: Path, encoding: str) -> Optional[Path]:
        out = self.previews / (safe_rel_name(path.name) + ".preview.txt")
        try:
            lines = []
            with path.open("r", encoding=encoding, errors="replace") as f:
                for i, line in enumerate(f):
                    if i >= MAX_TEXT_PREVIEW_LINES:
                        break
                    lines.append(line.rstrip("\n\r"))
            out.write_text("\n".join(lines), encoding="utf-8")
            return out
        except Exception as e:
            self.issues.append({"path": str(path), "issue": f"preview failed: {e}"})
            return None

    def summarize_text(self, path: Path, rel: str) -> Dict:
        encoding = guess_encoding(path)
        line_count = 0
        nonempty = 0
        max_line_len = 0
        first_lines = []
        last_lines = []
        try:
            with path.open("r", encoding=encoding, errors="replace") as f:
                for line in f:
                    line_count += 1
                    txt = line.rstrip("\n\r")
                    if txt.strip():
                        nonempty += 1
                    max_line_len = max(max_line_len, len(txt))
                    if len(first_lines) < 10:
                        first_lines.append(txt[:500])
                    last_lines.append(txt[:500])
                    if len(last_lines) > 10:
                        last_lines.pop(0)
            self.preview_text(path, encoding)
            return {
                "type": "text",
                "relative_path": rel,
                "encoding": encoding,
                "lines": line_count,
                "nonempty_lines": nonempty,
                "max_line_length": max_line_len,
                "first_lines": first_lines,
                "last_lines": last_lines,
            }
        except Exception as e:
            return {"type": "text", "relative_path": rel, "error": str(e)}

    def summarize_csv(self, path: Path, rel: str) -> Dict:
        encoding = guess_encoding(path)
        ext = path.suffix.lower()
        rows = 0
        blank_rows = 0
        malformed = 0
        headers = []
        width_counts = Counter()
        samples = []
        delimiter = ","
        try:
            with path.open("r", encoding=encoding, errors="replace", newline="") as f:
                sample = f.read(65536)
                delimiter = detect_delimiter(sample, ext)
                f.seek(0)
                reader = csv.reader(f, delimiter=delimiter)
                for i, row in enumerate(reader):
                    if i == 0:
                        headers = [c.strip() for c in row]
                        width_counts[len(row)] += 1
                        continue
                    rows += 1
                    width_counts[len(row)] += 1
                    if not any(str(c).strip() for c in row):
                        blank_rows += 1
                    if headers and len(row) != len(headers):
                        malformed += 1
                    if len(samples) < 8:
                        samples.append(row[:30])
                    if self.mode == "fast" and rows >= MAX_CSV_SAMPLE_ROWS and path.stat().st_size > 250 * 1024 * 1024:
                        break
            # write a normalized small preview CSV
            preview_path = self.previews / (safe_rel_name(path.name) + ".preview.csv")
            with preview_path.open("w", encoding="utf-8", newline="") as pf:
                w = csv.writer(pf)
                if headers:
                    w.writerow(headers)
                w.writerows(samples)
            estimated = self.mode == "fast" and rows >= MAX_CSV_SAMPLE_ROWS and path.stat().st_size > 250 * 1024 * 1024
            return {
                "type": "tabular",
                "relative_path": rel,
                "encoding": encoding,
                "delimiter": "TAB" if delimiter == "\t" else delimiter,
                "headers": headers,
                "rows_scanned": rows,
                "row_count_is_partial": estimated,
                "blank_rows_scanned": blank_rows,
                "rows_with_unexpected_column_count": malformed,
                "column_count_distribution": dict(width_counts),
                "sample_rows": samples,
            }
        except Exception as e:
            return {"type": "tabular", "relative_path": rel, "error": str(e)}

    def summarize_json(self, path: Path, rel: str) -> Dict:
        encoding = guess_encoding(path)
        ext = path.suffix.lower()
        try:
            if ext in {".jsonl", ".ndjson"}:
                count = 0
                invalid = 0
                key_counts = Counter()
                sample = []
                with path.open("r", encoding=encoding, errors="replace") as f:
                    for line in f:
                        if not line.strip():
                            continue
                        count += 1
                        try:
                            obj = json.loads(line)
                            if isinstance(obj, dict):
                                key_counts.update(obj.keys())
                            if len(sample) < 5:
                                sample.append(obj)
                        except Exception:
                            invalid += 1
                        if self.mode == "fast" and count >= 100000 and path.stat().st_size > 250 * 1024 * 1024:
                            break
                return {
                    "type": "jsonl",
                    "relative_path": rel,
                    "encoding": encoding,
                    "records_scanned": count,
                    "invalid_records_scanned": invalid,
                    "common_keys": key_counts.most_common(100),
                    "sample": sample,
                    "record_count_is_partial": self.mode == "fast" and count >= 100000 and path.stat().st_size > 250 * 1024 * 1024,
                }
            if path.stat().st_size > 300 * 1024 * 1024 and self.mode == "fast":
                return {
                    "type": "json",
                    "relative_path": rel,
                    "encoding": encoding,
                    "note": "Large JSON not fully parsed in fast mode; use Deep mode for complete validation.",
                    "size_bytes": path.stat().st_size,
                }
            with path.open("r", encoding=encoding, errors="replace") as f:
                obj = json.load(f)
            if isinstance(obj, dict):
                return {
                    "type": "json",
                    "relative_path": rel,
                    "encoding": encoding,
                    "root_type": "object",
                    "top_level_keys": list(obj.keys())[:500],
                    "top_level_key_count": len(obj),
                }
            if isinstance(obj, list):
                key_counts = Counter()
                for x in obj[:MAX_JSON_ARRAY_SAMPLE]:
                    if isinstance(x, dict):
                        key_counts.update(x.keys())
                return {
                    "type": "json",
                    "relative_path": rel,
                    "encoding": encoding,
                    "root_type": "array",
                    "item_count": len(obj),
                    "sampled_common_keys": key_counts.most_common(100),
                }
            return {"type": "json", "relative_path": rel, "encoding": encoding, "root_type": type(obj).__name__}
        except Exception as e:
            return {"type": "json", "relative_path": rel, "error": str(e)}

    def write_summary(self, path: Path, rel: str, category: str) -> Optional[Path]:
        if category == "tabular":
            summary = self.summarize_csv(path, rel)
        elif category == "json":
            summary = self.summarize_json(path, rel)
        elif category == "text":
            summary = self.summarize_text(path, rel)
        else:
            return None
        fname = safe_rel_name(rel.replace(os.sep, "__")) + ".summary.json"
        out = self.summaries / fname
        out.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
        self.summary_files.append(out)
        return out

    def process_one(self, root: Path, path: Path, rel_override: Optional[str] = None, source_root_override: Optional[str] = None):
        try:
            size = path.stat().st_size
            digest = sha256_file(path)
            rel = rel_override or str(path.relative_to(root))
            category = self.category_for(path)
            rec = FileRecord(
                source_root=source_root_override or str(root),
                path=str(path),
                relative_path=rel,
                extension=path.suffix.lower(),
                size_bytes=size,
                sha256=digest,
                category=category,
            )
            if digest in self.hash_first:
                rec.duplicate_of = self.hash_first[digest]
                rec.status = "duplicate"
            else:
                self.hash_first[digest] = rel
            self.records.append(rec)
            if category in {"tabular", "json", "text"} and rec.status != "duplicate":
                self.write_summary(path, rel, category)
            elif category == "compressed_tabular" and rec.status != "duplicate":
                if self.scientific and self.reducer is not None:
                    try:
                        outputs, report = self.reducer.reduce(path, rel)
                        self.scientific_files.extend(outputs)
                        self.scientific_reports.append(report)
                        summary_out = self.summaries / (safe_rel_name(rel.replace(os.sep, "__")) + ".summary.json")
                        summary_out.write_text(json.dumps({
                            "type": "compressed_tabular",
                            "relative_path": rel,
                            "scientific_reduction": report,
                        }, ensure_ascii=False, indent=2), encoding="utf-8")
                        self.summary_files.append(summary_out)
                    except Exception as e:
                        self.issues.append({"path": str(path), "issue": f"scientific reduction failed: {e}"})
                else:
                    summary_out = self.summaries / (safe_rel_name(rel.replace(os.sep, "__")) + ".summary.json")
                    summary_out.write_text(json.dumps({
                        "type": "compressed_tabular",
                        "relative_path": rel,
                        "size_bytes": size,
                        "note": "Scientific reduction disabled for this run."
                    }, ensure_ascii=False, indent=2), encoding="utf-8")
                    self.summary_files.append(summary_out)
            elif category == "archive" and rec.status != "duplicate":
                self.process_zip(path, rel)
            self.total_processed_bytes += size
        except Exception as e:
            self.issues.append({"path": str(path), "issue": f"processing failed: {e}"})

    def process_zip(self, zip_path: Path, rel: str):
        dest = self.extracted / safe_rel_name(Path(rel).stem + "_" + hashlib.md5(str(zip_path).encode()).hexdigest()[:8])
        try:
            files = self.safe_extract_zip(zip_path, dest)
            self.log(f"  Extracted {len(files)} files from {Path(rel).name}")
            for child in files:
                child_rel = f"{rel}::/{child.relative_to(dest)}"
                # Avoid infinite recursion from deliberately nested archives by allowing two levels.
                nested_depth = child_rel.count("::/")
                if child.suffix.lower() == ".zip" and nested_depth > 2:
                    self.issues.append({"path": child_rel, "issue": "nested ZIP depth limit reached"})
                    continue
                self.process_one(dest, child, rel_override=child_rel, source_root_override=str(zip_path))
        except zipfile.BadZipFile:
            self.issues.append({"path": str(zip_path), "issue": "invalid/corrupt ZIP"})
        except Exception as e:
            self.issues.append({"path": str(zip_path), "issue": f"ZIP extraction failed: {e}"})

    def write_inventory(self):
        inv = self.output_dir / "INVENTORY.csv"
        fields = list(FileRecord.__dataclass_fields__.keys())
        with inv.open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in self.records:
                w.writerow(asdict(r))
        issues_path = self.output_dir / "ISSUES.csv"
        with issues_path.open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["path", "issue"])
            w.writeheader()
            w.writerows(self.issues)
        return inv, issues_path

    def write_sqlite(self):
        db = self.output_dir / "analysis_catalog.sqlite"
        con = sqlite3.connect(db)
        try:
            con.execute("DROP TABLE IF EXISTS files")
            con.execute("DROP TABLE IF EXISTS issues")
            con.execute("""CREATE TABLE files (
                source_root TEXT, path TEXT, relative_path TEXT, extension TEXT,
                size_bytes INTEGER, sha256 TEXT, category TEXT, status TEXT,
                notes TEXT, duplicate_of TEXT
            )""")
            con.executemany(
                "INSERT INTO files VALUES (?,?,?,?,?,?,?,?,?,?)",
                [(r.source_root, r.path, r.relative_path, r.extension, r.size_bytes, r.sha256, r.category, r.status, r.notes, r.duplicate_of) for r in self.records],
            )
            con.execute("CREATE INDEX idx_files_sha ON files(sha256)")
            con.execute("CREATE INDEX idx_files_cat ON files(category)")
            con.execute("CREATE TABLE issues (path TEXT, issue TEXT)")
            con.executemany("INSERT INTO issues VALUES (?,?)", [(i["path"], i["issue"]) for i in self.issues])
            con.commit()
        finally:
            con.close()
        return db

    def run_fingerprint(self) -> str:
        h = hashlib.sha256()
        for r in sorted(self.records, key=lambda x: (x.relative_path, x.sha256)):
            h.update(r.relative_path.encode("utf-8", errors="replace"))
            h.update(b"\0")
            h.update(r.sha256.encode("ascii"))
            h.update(b"\n")
        return h.hexdigest()

    def make_master_summary(self):
        cats = Counter(r.category for r in self.records)
        exts = Counter(r.extension or "(none)" for r in self.records)
        dups = sum(1 for r in self.records if r.status == "duplicate")
        total_bytes = sum(r.size_bytes for r in self.records)
        summary = {
            "app": APP_NAME,
            "version": APP_VERSION,
            "generated_at_local": time.strftime("%Y-%m-%d %H:%M:%S"),
            "mode": self.mode,
            "scientific_reduction_enabled": self.scientific,
            "scientific_tables_reduced": len(self.scientific_reports),
            "scientific_rows_streamed": sum(int(r.get("rows_streamed", 0)) for r in self.scientific_reports),
            "run_fingerprint_sha256": self.run_fingerprint(),
            "files_catalogued_including_extracted_archive_members": len(self.records),
            "total_catalogued_bytes": total_bytes,
            "total_catalogued_human": human_bytes(total_bytes),
            "duplicates": dups,
            "issues": len(self.issues),
            "by_category": dict(cats),
            "top_extensions": exts.most_common(50),
        }
        p = self.output_dir / "MASTER_SUMMARY.json"
        p.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        return p, summary

    def pack_candidates(self, core_files: List[Path]) -> List[Path]:
        cands = []
        cands.extend(core_files)
        cands.extend(sorted(self.summary_files))
        cands.extend(sorted(self.previews.glob("*")))
        cands.extend(sorted(set(self.scientific_files)))
        # Include original small files where useful, but avoid duplicates and huge binaries.
        for r in self.records:
            if r.status == "duplicate":
                continue
            if r.category in {"tabular", "json", "text", "compressed_tabular"} and r.size_bytes <= 10 * 1024 * 1024:
                p = Path(r.path)
                if p.exists() and p.is_file():
                    cands.append(p)
        # unique by resolved path
        out, seen = [], set()
        for p in cands:
            try:
                key = str(p.resolve())
            except Exception:
                key = str(p)
            if key not in seen and p.exists() and p.is_file():
                seen.add(key)
                out.append(p)
        return out

    def create_packs(self, candidates: List[Path]) -> List[Path]:
        packs = []
        pack_index = 1
        current_size = 0
        zf = None
        current_path = None

        def open_new():
            nonlocal pack_index, current_size, zf, current_path
            current_path = self.upload_dir / f"ANALYSIS_PACK_{pack_index:03d}.zip"
            zf = zipfile.ZipFile(current_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6)
            current_size = 0
            packs.append(current_path)
            pack_index += 1

        def close_current():
            nonlocal zf
            if zf is not None:
                zf.close()
                zf = None

        open_new()
        for p in candidates:
            try:
                size = p.stat().st_size
                if current_size > 0 and current_size + size > self.pack_limit:
                    close_current()
                    open_new()
                if p.is_relative_to(self.output_dir):
                    arc = str(p.relative_to(self.output_dir))
                else:
                    arc = "source_small_files/" + safe_rel_name(p.name)
                zf.write(p, arcname=arc)
                current_size += size
            except Exception as e:
                self.issues.append({"path": str(p), "issue": f"packaging failed: {e}"})
        close_current()
        return packs

    def write_readme(self, packs: List[Path], summary: Dict):
        text = f"{APP_NAME} {APP_VERSION} - Analysis Package\n\n"
        text += "PURPOSE\n-------\n"
        text += "This folder was created automatically from raw folders/ZIP files. No manual extraction is needed.\n"
        text += "Upload ANALYSIS_INDEX.zip first, then ANALYSIS_PACK_001.zip onward in numeric order.\n\n"
        text += "KEY FILES\n---------\n"
        text += "MASTER_SUMMARY.json  - overall run statistics\n"
        text += "INVENTORY.csv        - one row per catalogued file/archive member with SHA-256\n"
        text += "ISSUES.csv           - corrupt, unsafe, unreadable or skipped items\n"
        text += "analysis_catalog.sqlite - searchable local catalogue\n"
        text += "summaries/            - machine-readable summaries of CSV/JSON/text/log files\n"
        text += "previews/             - small human-readable samples\n"
        text += "SCIENTIFIC_REDUCTION/ - compact derivatives streamed from large .csv.gz files\n\n"
        text += f"RUN STATS\n---------\nFiles catalogued: {summary['files_catalogued_including_extracted_archive_members']}\n"
        text += f"Duplicates: {summary['duplicates']}\nIssues: {summary['issues']}\nMode: {summary['mode']}\n"
        text += f"Scientific tables reduced: {summary.get('scientific_tables_reduced', 0)}\n"
        text += f"Scientific rows streamed: {summary.get('scientific_rows_streamed', 0):,}\n"
        text += "Upload packs are in CHATGPT_UPLOAD and are numbered in the order they should be uploaded.\n\n"
        text += "IMPORTANT\n---------\n"
        text += "Original input files are never modified. ZIP extraction is isolated in the run workspace.\n"
        text += "Keep the full run folder locally; only CHATGPT_UPLOAD packs normally need to be uploaded.\n"
        p = self.output_dir / "README_FIRST.txt"
        p.write_text(text, encoding="utf-8")
        return p

    def write_final_audit(self, summary: Dict) -> Path:
        audit = {
            "app": APP_NAME,
            "version": APP_VERSION,
            "run_fingerprint_sha256": summary["run_fingerprint_sha256"],
            "catalogued_files": len(self.records),
            "unique_files": sum(1 for r in self.records if r.status != "duplicate"),
            "duplicates": sum(1 for r in self.records if r.status == "duplicate"),
            "issues": len(self.issues),
            "status": "PASS" if not self.issues else "PASS_WITH_ISSUES",
            "original_inputs_modified": False,
            "unsafe_zip_members_blocked": sum(1 for i in self.issues if "unsafe ZIP member" in i.get("issue", "")),
            "created_at_local": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        p = self.output_dir / "FINAL_AUDIT.json"
        p.write_text(json.dumps(audit, indent=2), encoding="utf-8")
        return p

    def write_upload_index(self, packs: List[Path], core_files: List[Path]) -> Path:
        manifest = self.output_dir / "UPLOAD_MANIFEST.csv"
        with manifest.open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f)
            w.writerow(["upload_order", "file_name", "size_bytes", "sha256"])
            for i, pack in enumerate(packs, start=1):
                w.writerow([i + 1, pack.name, pack.stat().st_size, sha256_file(pack)])

        index_zip = self.upload_dir / "ANALYSIS_INDEX.zip"
        with zipfile.ZipFile(index_zip, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
            for p in core_files + [manifest]:
                if p.exists() and p.is_file():
                    zf.write(p, arcname=p.name)
        # Insert the index itself as upload order 1 after its bytes are final.
        rows = []
        with manifest.open("r", encoding="utf-8-sig", newline="") as f:
            rows = list(csv.reader(f))
        with manifest.open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f)
            w.writerow(rows[0])
            w.writerow([1, index_zip.name, index_zip.stat().st_size, sha256_file(index_zip)])
            for row in rows[1:]:
                w.writerow(row)
        return index_zip

    def run(self, inputs: List[Path]) -> Path:
        self.setup()
        self.log(f"{APP_NAME} {APP_VERSION} started")
        items = self.collect_paths(inputs)
        if not items:
            raise RuntimeError("No readable input files found.")
        self.log(f"Found {len(items)} top-level files ({human_bytes(self.total_input_bytes)})")
        total = len(items)
        for idx, (root, path) in enumerate(items, start=1):
            self.progress_cb(idx - 1, total, f"Processing {path.name}")
            self.log(f"[{idx}/{total}] {path}")
            self.process_one(root, path)
        self.progress_cb(total, total, "Building catalogue and upload packs")
        inv, issues = self.write_inventory()
        db = self.write_sqlite()
        master, summary = self.make_master_summary()
        audit = self.write_final_audit(summary)
        # Create the handoff instructions before packaging so the packs contain the final core metadata.
        readme = self.write_readme([], summary)
        handoff = self.output_dir / "COPY_THIS_WITH_UPLOAD.txt"
        handoff.write_text(
            "These files were prepared locally by AnalysisBridge.\n"
            "Please verify FINAL_AUDIT.json and MASTER_SUMMARY.json first, then use INVENTORY.csv as the source-of-truth file map.\n"
            "Treat SHA-256 values in INVENTORY.csv as file identity and relative_path as provenance.\n"
            "If multiple ANALYSIS_PACK files are supplied, combine them logically as one run sharing the same run fingerprint.\n"
            "Use ISSUES.csv before drawing conclusions from missing or unreadable content.\n"
            "If SCIENTIFIC_REDUCTION is present, prefer its run/month/season/extreme/event/matched-duration tables over raw hourly uploads.\n",
            encoding="utf-8",
        )
        core = [master, audit, inv, issues, db, readme, handoff]
        candidates = self.pack_candidates(core)
        packs = self.create_packs(candidates)
        index_zip = self.write_upload_index(packs, core)
        self.log(f"Done. Created index + {len(packs)} data pack(s) in {self.upload_dir}")
        self.progress_cb(total, total, "Complete")
        return self.output_dir


def launch_gui():
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    class App(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title(f"{APP_NAME} {APP_VERSION}")
            self.geometry("920x650")
            self.minsize(820, 560)
            self.inputs: List[str] = []
            self.q = queue.Queue()
            self.build_ui()
            self.after(100, self.poll_queue)

        def build_ui(self):
            pad = {"padx": 10, "pady": 6}
            title = ttk.Label(self, text="AnalysisBridge", font=("Segoe UI", 20, "bold"))
            title.pack(anchor="w", padx=14, pady=(12, 0))
            ttk.Label(self, text="Raw folders / ZIPs -> verified catalogue + streamed scientific reductions + upload-ready packs").pack(anchor="w", padx=14, pady=(0, 8))

            top = ttk.Frame(self)
            top.pack(fill="x", **pad)
            ttk.Button(top, text="Add Folder", command=self.add_folder).pack(side="left", padx=(0, 6))
            ttk.Button(top, text="Add Files / ZIPs", command=self.add_files).pack(side="left", padx=6)
            ttk.Button(top, text="Clear", command=self.clear_inputs).pack(side="left", padx=6)

            self.listbox = tk.Listbox(self, height=8, selectmode="extended")
            self.listbox.pack(fill="x", padx=14, pady=4)

            opts = ttk.LabelFrame(self, text="Output settings")
            opts.pack(fill="x", padx=14, pady=8)
            ttk.Label(opts, text="Output folder:").grid(row=0, column=0, sticky="w", **pad)
            self.out_var = tk.StringVar(value=str(Path.home() / "Desktop" / "AnalysisBridge_Output"))
            ttk.Entry(opts, textvariable=self.out_var).grid(row=0, column=1, sticky="ew", **pad)
            ttk.Button(opts, text="Browse", command=self.choose_output).grid(row=0, column=2, **pad)
            ttk.Label(opts, text="Pack size (MB):").grid(row=1, column=0, sticky="w", **pad)
            self.pack_var = tk.IntVar(value=DEFAULT_PACK_MB)
            ttk.Spinbox(opts, from_=20, to=450, textvariable=self.pack_var, width=10).grid(row=1, column=1, sticky="w", **pad)
            ttk.Label(opts, text="Mode:").grid(row=1, column=1, sticky="e", padx=(0, 120), pady=6)
            self.mode_var = tk.StringVar(value="fast")
            ttk.Combobox(opts, textvariable=self.mode_var, values=["fast", "deep"], state="readonly", width=12).grid(row=1, column=1, sticky="e", **pad)
            self.science_var = tk.BooleanVar(value=True)
            ttk.Checkbutton(opts, text="Scientific reduction of large .csv.gz outputs (recommended)", variable=self.science_var).grid(row=2, column=1, sticky="w", **pad)
            opts.columnconfigure(1, weight=1)

            self.start_btn = ttk.Button(self, text="Build Analysis Packs", command=self.start)
            self.start_btn.pack(anchor="w", padx=14, pady=(4, 8))
            self.progress = ttk.Progressbar(self, mode="determinate")
            self.progress.pack(fill="x", padx=14, pady=4)
            self.status_var = tk.StringVar(value="Ready")
            ttk.Label(self, textvariable=self.status_var).pack(anchor="w", padx=14)

            logbox = ttk.LabelFrame(self, text="Activity")
            logbox.pack(fill="both", expand=True, padx=14, pady=10)
            self.log_text = tk.Text(logbox, height=12, wrap="word")
            self.log_text.pack(fill="both", expand=True, padx=6, pady=6)
            self.log_text.configure(state="disabled")

        def add_folder(self):
            p = filedialog.askdirectory()
            if p:
                self.add_input(p)

        def add_files(self):
            paths = filedialog.askopenfilenames(title="Select files or ZIP archives")
            for p in paths:
                self.add_input(p)

        def add_input(self, p):
            if p not in self.inputs:
                self.inputs.append(p)
                self.listbox.insert("end", p)

        def clear_inputs(self):
            self.inputs.clear()
            self.listbox.delete(0, "end")

        def choose_output(self):
            p = filedialog.askdirectory()
            if p:
                self.out_var.set(str(Path(p) / "AnalysisBridge_Output"))

        def append_log(self, msg):
            self.log_text.configure(state="normal")
            self.log_text.insert("end", msg + "\n")
            self.log_text.see("end")
            self.log_text.configure(state="disabled")

        def start(self):
            if not self.inputs:
                messagebox.showwarning("No input", "Add at least one folder, file, or ZIP archive.")
                return
            self.start_btn.configure(state="disabled")
            self.progress["value"] = 0
            self.status_var.set("Starting...")
            out = Path(self.out_var.get()).expanduser()
            mode = self.mode_var.get()
            pack_mb = int(self.pack_var.get())
            scientific = bool(self.science_var.get())
            t = threading.Thread(target=self.worker, args=(list(self.inputs), out, pack_mb, mode, scientific), daemon=True)
            t.start()

        def worker(self, inputs, out, pack_mb, mode, scientific):
            try:
                eng = AnalysisBridgeEngine(
                    out,
                    pack_mb=pack_mb,
                    mode=mode,
                    log_cb=lambda m: self.q.put(("log", m)),
                    progress_cb=lambda a, b, m="": self.q.put(("progress", a, b, m)),
                    scientific=scientific,
                )
                result = eng.run([Path(p) for p in inputs])
                self.q.put(("done", str(result)))
            except Exception:
                self.q.put(("error", traceback.format_exc()))

        def poll_queue(self):
            try:
                while True:
                    item = self.q.get_nowait()
                    kind = item[0]
                    if kind == "log":
                        self.append_log(item[1])
                    elif kind == "progress":
                        _, a, b, msg = item
                        self.progress["maximum"] = max(1, b)
                        self.progress["value"] = a
                        self.status_var.set(msg)
                    elif kind == "done":
                        self.start_btn.configure(state="normal")
                        self.status_var.set("Complete")
                        result = item[1]
                        messagebox.showinfo("Complete", f"Analysis packs created successfully.\n\nUpload ANALYSIS_INDEX.zip first from:\n{result}\\CHATGPT_UPLOAD")
                    elif kind == "error":
                        self.start_btn.configure(state="normal")
                        self.status_var.set("Error")
                        self.append_log(item[1])
                        messagebox.showerror("Error", "The run failed. See Activity for details.")
            except queue.Empty:
                pass
            self.after(100, self.poll_queue)

    App().mainloop()


def cli_main(argv=None):
    p = argparse.ArgumentParser(description="Build upload-friendly analysis packs from raw folders/files/ZIPs.")
    p.add_argument("inputs", nargs="+", help="Input files or folders")
    p.add_argument("-o", "--output", default="AnalysisBridge_Output", help="Output directory")
    p.add_argument("--pack-mb", type=int, default=DEFAULT_PACK_MB, help="Approximate maximum source-size per upload pack")
    p.add_argument("--mode", choices=["fast", "deep"], default="fast")
    p.add_argument("--no-science", action="store_true", help="Disable streamed scientific reduction of compressed tabular outputs")
    args = p.parse_args(argv)
    eng = AnalysisBridgeEngine(Path(args.output), args.pack_mb, args.mode, log_cb=print, scientific=not args.no_science)
    out = eng.run([Path(x) for x in args.inputs])
    print(f"\nOutput: {out}")
    print(f"Upload folder: {out / 'CHATGPT_UPLOAD'}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) == 1:
        launch_gui()
    else:
        raise SystemExit(cli_main())
