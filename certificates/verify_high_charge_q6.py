#!/usr/bin/env python3
"""Replay the exact arbitrary-complex q6 high-charge base certificate."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

if not __debug__:
    raise SystemExit("exact verifier requires assertions; do not use python -O")


HERE = Path(__file__).resolve().parent
WRAPPER_SOURCE = Path(__file__).resolve()
ORBIT_SOURCE = HERE / "verify_high_charge_q6_orbits.cpp"
ALGEBRA_VERIFIER = HERE / "verify_high_charge_q6_algebra.py"
EXPECTED_Z3_VERSION = "4.14.0"


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def resolve_executable(requested: str | None, candidates: tuple[str, ...]) -> str | None:
    if requested:
        resolved = shutil.which(requested)
        if resolved:
            return resolved
        path = Path(requested).expanduser().resolve()
        return str(path) if path.is_file() else None
    return next((path for name in candidates if (path := shutil.which(name))), None)


def compiler_command(executable: str) -> list[str]:
    """Return a GCC/Clang-compatible compiler command prefix."""

    if Path(executable).name.lower() in {"zig", "zig.exe"}:
        return [executable, "c++"]
    return [executable]


def version_line(command: list[str], argument: str = "--version") -> str:
    process = subprocess.run(
        [*command, argument],
        text=True,
        capture_output=True,
        timeout=10,
        check=False,
    )
    output = (process.stdout or process.stderr).strip()
    if process.returncode != 0 or not output:
        raise SystemExit(
            "could not obtain version from " + subprocess.list2cmdline(command)
        )
    return output.splitlines()[0]


def source_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Replay the exact arbitrary-complex q6 certificate."
    )
    parser.add_argument(
        "--compiler",
        help="path or command name for clang++, c++, g++, or zig",
    )
    parser.add_argument(
        "--z3",
        help="path or command name for the pinned Z3 4.14.0 executable",
    )
    parser.add_argument("--audit-dir", type=Path,
                        help="retain the orbit handshake and exact nonlinear query inputs")
    args = parser.parse_args()

    compiler = resolve_executable(args.compiler, ("clang++", "c++", "g++", "zig"))
    z3 = resolve_executable(args.z3, ("z3",))
    missing = []
    if compiler is None:
        missing.append("a GCC- or Clang-compatible C++17 compiler")
    if z3 is None:
        missing.append(f"Z3 {EXPECTED_Z3_VERSION}")
    if missing:
        raise SystemExit(
            "missing q6 dependencies: "
            + "; ".join(missing)
            + ". Install them on PATH or pass --compiler/--z3 with explicit paths."
        )

    compiler_prefix = compiler_command(compiler)
    compiler_version = version_line(compiler_prefix)
    z3_version = version_line([z3])
    match = re.search(r"\b(\d+\.\d+\.\d+)\b", z3_version)
    if match is None or match.group(1) != EXPECTED_Z3_VERSION:
        found = match.group(1) if match else "unknown"
        raise SystemExit(
            f"q6 certificate pins Z3 {EXPECTED_Z3_VERSION}; found {found}: "
            f"{z3_version}"
        )

    print(
        "compiler: "
        + compiler_version
        + " ("
        + subprocess.list2cmdline(compiler_prefix)
        + ")",
        flush=True,
    )
    print(f"z3: {z3_version} ({z3})", flush=True)
    for source in (WRAPPER_SOURCE, ORBIT_SOURCE, ALGEBRA_VERIFIER):
        print(f"sha256 {source.name}: {source_digest(source)}", flush=True)

    with tempfile.TemporaryDirectory(prefix="high-charge-q6-") as directory:
        executable = Path(directory) / (
            "verify_high_charge_q6_orbits.exe"
            if os.name == "nt"
            else "verify_high_charge_q6_orbits"
        )
        orbit_handshake = Path(directory) / "q6-orbits-v1.txt"
        run([
            *compiler_prefix,
            "-O3",
            "-std=c++17",
            str(ORBIT_SOURCE),
            "-o",
            str(executable),
        ])
        run([str(executable), "--emit-orbits", str(orbit_handshake)])
        audit_arguments = []
        if args.audit_dir is not None:
            audit_directory = args.audit_dir.resolve()
            audit_directory.mkdir(parents=True, exist_ok=True)
            if (audit_directory / "queries.json").exists():
                raise SystemExit("audit directory already contains a run; use a new directory")
            shutil.copyfile(orbit_handshake, audit_directory / "q6-orbits-v1.txt")
            audit_arguments = ["--audit-dir", str(audit_directory)]
        run([
            sys.executable,
            str(ALGEBRA_VERIFIER),
            "--cpp-orbits",
            str(orbit_handshake),
            "--z3",
            z3,
            *audit_arguments,
        ])


if __name__ == "__main__":
    main()
