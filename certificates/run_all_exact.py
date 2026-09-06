#!/usr/bin/env python3
"""Run Paper 2's exact finite-verifier suite."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import platform
import subprocess
import sys
import time
from pathlib import Path


if not __debug__:
    raise SystemExit("exact verifiers require assertions; do not use python -O")


HERE = Path(__file__).resolve().parent
COMMANDS = (
    ("odd-dihedral q6", "verify_odd_dihedral.py"),
    ("high-charge q6", "verify_high_charge_q6.py"),
    ("high-charge q8 atlas", "verify_high_charge_q8_atlas.py"),
    ("high-charge q8 normalized", "verify_high_charge_q8_normalized.py"),
)


def source_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run every Paper 2 exact verifier, continue after dependency or "
            "verification failures, and optionally save a complete transcript."
        )
    )
    parser.add_argument(
        "--transcript",
        type=Path,
        help="write the complete UTF-8 run transcript to this path",
    )
    parser.add_argument(
        "--q6-compiler",
        help="compiler executable passed through to the high-charge q6 wrapper",
    )
    parser.add_argument(
        "--q6-z3",
        help="pinned Z3 executable passed through to the high-charge q6 wrapper",
    )
    parser.add_argument("--q6-audit-dir", type=Path,
                        help="retain exact q6 nonlinear inputs and solver responses")
    args = parser.parse_args()

    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    transcript: list[str] = []

    def emit(message: str = "") -> None:
        print(message, flush=True)
        transcript.append(message)

    emit("PAPER 2 EXACT-VERIFIER TRANSCRIPT v1")
    emit(f"started_utc: {dt.datetime.now(dt.timezone.utc).isoformat()}")
    emit(f"platform: {platform.platform()}")
    emit(f"python: {sys.version.splitlines()[0]}")
    emit(f"executable: {sys.executable}")
    emit(f"certificate_directory: {HERE}")
    emit(f"sha256 {Path(__file__).name}: {source_digest(Path(__file__))}")

    failures = []
    for label, script, *arguments in COMMANDS:
        if script == "verify_high_charge_q6.py":
            if args.q6_compiler:
                arguments.extend(("--compiler", args.q6_compiler))
            if args.q6_z3:
                arguments.extend(("--z3", args.q6_z3))
            if args.q6_audit_dir:
                arguments.extend(("--audit-dir", str(args.q6_audit_dir.resolve())))
        command = [sys.executable, str(HERE / script), *arguments]
        emit(f"\n=== {label} ===")
        emit(f"sha256 {script}: {source_digest(HERE / script)}")
        emit("command: " + subprocess.list2cmdline(command))
        started = time.monotonic()
        process = subprocess.run(
            command,
            cwd=HERE,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )
        if process.stdout:
            for line in process.stdout.rstrip().splitlines():
                emit(line)
        if process.stderr:
            for line in process.stderr.rstrip().splitlines():
                emit("STDERR: " + line)
        runtime = time.monotonic() - started
        emit(f"result: exit={process.returncode}, runtime={runtime:.3f}s")
        if process.returncode != 0:
            failures.append(label)

    if failures:
        emit("\nPAPER 2 EXACT-VERIFIER SUITE: FAIL")
        emit("failed_entries: " + ", ".join(failures))
    else:
        emit("\nPAPER 2 EXACT-VERIFIER SUITE: PASS")
    emit(f"finished_utc: {dt.datetime.now(dt.timezone.utc).isoformat()}")

    if args.transcript is not None:
        transcript_path = args.transcript
        if not transcript_path.is_absolute():
            transcript_path = transcript_path.resolve()
        transcript_path.parent.mkdir(parents=True, exist_ok=True)
        transcript_path.write_text("\n".join(transcript) + "\n", encoding="utf-8")
        print(f"transcript_written: {transcript_path}", flush=True)

    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
