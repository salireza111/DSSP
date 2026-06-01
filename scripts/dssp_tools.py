#!/usr/bin/env python3
"""Developer utilities for JS DSSP output, local mkdssp output, and comparison."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
JS_MODULE_URL = (REPO_ROOT / "src" / "dssp.js").as_uri()

STRUCTURE_NAMES = {
    "H": "alpha helix",
    "B": "isolated beta bridge",
    "E": "extended beta strand",
    "G": "3-10 helix",
    "I": "pi helix",
    "P": "poly-proline II helix",
    "T": "hydrogen-bonded turn",
    "S": "bend",
    "-": "loop",
}

JS_DRIVER = r"""
import { readFileSync } from "node:fs";

const {
  assignDSSP,
  filterPDBByChain,
  toCSV,
  toDsspJSON
} = await import(process.env.JS_DSSP_MODULE_URL);

const pdbPath = process.env.PDB_PATH;
const outputFormat = process.env.OUTPUT_FORMAT || "json";
const chain = process.env.PDB_CHAIN || "";
const asaSamples = Number.parseInt(process.env.ASA_SAMPLES || "401", 10);
const hbondCutoff = Number.parseFloat(process.env.HBOND_CUTOFF || "-0.5");

let pdbText = readFileSync(pdbPath, "utf8");
if (chain) pdbText = filterPDBByChain(pdbText, chain);

const result = assignDSSP(pdbText, { asaSamples, hbondCutoff });

if (outputFormat === "csv") {
  process.stdout.write(`${toCSV(result)}\n`);
} else {
  process.stdout.write(`${JSON.stringify(toDsspJSON(result), null, 2)}\n`);
}
"""


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    try:
        if args.command == "js":
            output = run_js_dssp(args.pdb, args)
            write_output(output, args.output)
            return 0
        if args.command == "real":
            parsed = run_real_dssp(args.pdb, args)
            output = format_real_output(parsed, args.format)
            write_output(output, args.output)
            return 0
        if args.command == "benchmark":
            js_result = json.loads(run_js_dssp(args.pdb, args, output_format="json"))
            real_result = run_real_dssp(args.pdb, args)
            report = compare_results(js_result, real_result, args.pdb, args)
            output = json.dumps(report, indent=2) + "\n" if args.json else format_benchmark_report(report)
            write_output(output, args.output)
            return 0
    except Exception as exc:  # noqa: BLE001 - command-line tool should report cleanly.
        print(f"error: {exc}", file=sys.stderr)
        return 1

    parser.print_help(sys.stderr)
    return 2


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export or benchmark this project's JS DSSP-compatible output against local mkdssp/dssp."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    js_parser = subparsers.add_parser("js", help="export DSSP output from src/dssp.js")
    add_pdb_argument(js_parser)
    add_js_options(js_parser)
    add_export_options(js_parser)

    real_parser = subparsers.add_parser("real", help="export DSSP output from installed mkdssp/dssp")
    add_pdb_argument(real_parser)
    add_real_options(real_parser)
    add_export_options(real_parser)

    benchmark_parser = subparsers.add_parser("benchmark", help="compare JS output against installed mkdssp/dssp")
    add_pdb_argument(benchmark_parser)
    add_js_options(benchmark_parser)
    add_real_options(benchmark_parser)
    benchmark_parser.add_argument("--json", action="store_true", help="write the benchmark report as JSON")
    benchmark_parser.add_argument("-o", "--output", help="write the report to this file instead of stdout")

    return parser


def add_pdb_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("pdb", help="input legacy PDB file")
    parser.add_argument("--chain", help="optional one-character chain ID to isolate before calculation")


def add_js_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--node-bin", default="node", help="Node.js executable for running src/dssp.js")
    parser.add_argument("--asa-samples", type=int, default=401, help="surface points per atom for JS ASA")
    parser.add_argument("--hbond-cutoff", type=float, default=-0.5, help="JS H-bond cutoff in kcal/mol")


def add_real_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--dssp-bin",
        help="mkdssp/dssp executable; defaults to DSSP_BIN, MKDSSP_BIN, mkdssp, then dssp",
    )


def add_export_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--format", choices=("json", "csv"), default="json", help="export format")
    parser.add_argument("-o", "--output", help="write output to this file instead of stdout")


def run_js_dssp(pdb_path: str, args: argparse.Namespace, output_format: str | None = None) -> str:
    node_bin = shutil.which(args.node_bin) if os.path.basename(args.node_bin) == args.node_bin else args.node_bin
    if not node_bin:
        raise RuntimeError("Node.js was not found. Install node or pass --node-bin.")

    env = os.environ.copy()
    env.update(
        {
            "JS_DSSP_MODULE_URL": JS_MODULE_URL,
            "PDB_PATH": str(Path(pdb_path).resolve()),
            "OUTPUT_FORMAT": output_format or args.format,
            "PDB_CHAIN": args.chain or "",
            "ASA_SAMPLES": str(args.asa_samples),
            "HBOND_CUTOFF": str(args.hbond_cutoff),
        }
    )
    completed = subprocess.run(
        [node_bin, "--input-type=module", "-e", JS_DRIVER],
        cwd=REPO_ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(clean_process_error("JS DSSP failed", completed))
    return completed.stdout


def run_real_dssp(pdb_path: str, args: argparse.Namespace) -> dict:
    dssp_bin = resolve_dssp_binary(args.dssp_bin)
    pdb = Path(pdb_path).resolve()
    if not pdb.exists():
        raise RuntimeError(f"PDB file does not exist: {pdb}")

    with tempfile.TemporaryDirectory(prefix="dssp-tools-") as tmp:
        tmp_dir = Path(tmp)
        input_path = pdb
        if args.chain:
            filtered = filter_pdb_by_chain(pdb.read_text(encoding="utf-8", errors="replace"), args.chain)
            input_path = tmp_dir / "chain-filtered.pdb"
            input_path.write_text(filtered, encoding="utf-8")

        dssp_path = tmp_dir / "output.dssp"
        run_dssp_command(dssp_bin, input_path, dssp_path)
        parsed = parse_classic_dssp(dssp_path.read_text(encoding="utf-8", errors="replace"))

    parsed["dssp_bin"] = dssp_bin
    parsed["source"] = "real"
    return parsed


def resolve_dssp_binary(explicit: str | None) -> str:
    candidates = [explicit, os.environ.get("DSSP_BIN"), os.environ.get("MKDSSP_BIN"), "mkdssp", "dssp"]
    for candidate in candidates:
        if not candidate:
            continue
        if os.path.sep in candidate:
            path = Path(candidate).expanduser()
            if path.exists():
                return str(path)
        found = shutil.which(candidate)
        if found:
            return found
    raise RuntimeError("No mkdssp/dssp executable found. Install mkdssp or pass --dssp-bin.")


def run_dssp_command(dssp_bin: str, input_path: Path, output_path: Path) -> None:
    attempts = [
        [dssp_bin, "--output-format", "dssp", str(input_path), str(output_path)],
        [dssp_bin, str(input_path), str(output_path)],
        [dssp_bin, "-i", str(input_path), "-o", str(output_path)],
    ]
    errors = []

    for command in attempts:
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        if output_path.exists() and output_path.stat().st_size > 0:
            return
        if "#  RESIDUE AA STRUCTURE" in completed.stdout:
            output_path.write_text(completed.stdout, encoding="utf-8")
            return
        errors.append(clean_process_error(" ".join(command), completed))

    raise RuntimeError("Could not run DSSP successfully:\n" + "\n\n".join(errors))


def filter_pdb_by_chain(pdb_text: str, chain: str) -> str:
    normalized = chain.strip().upper()
    if len(normalized) != 1:
        raise RuntimeError("Legacy PDB chain IDs must be one character.")

    kept_atoms = 0
    kept_lines = []
    chain_records = {"ATOM", "HETATM", "ANISOU", "TER"}
    for line in pdb_text.splitlines():
        record = line[:6].strip()
        if record in {"CONECT", "MASTER"}:
            continue
        if record in chain_records:
            line_chain = (line[21:22].strip() or "_").upper()
            if line_chain != normalized:
                continue
            if record == "ATOM":
                kept_atoms += 1
        kept_lines.append(line)

    if kept_atoms == 0:
        raise RuntimeError(f"No ATOM records were found for chain {normalized}.")
    return "\n".join(kept_lines) + "\n"


def parse_classic_dssp(text: str) -> dict:
    rows = []
    summary = {}
    in_rows = False

    for line in text.splitlines():
        if "TOTAL NUMBER OF RESIDUES" in line:
            values = line.split()
            summary["residueCount"] = int_or_none(values[0] if values else "")
            summary["chainCount"] = int_or_none(values[1] if len(values) > 1 else "")
        elif "ACCESSIBLE SURFACE OF PROTEIN" in line:
            summary["asaTotal"] = float_or_none(line.split()[0] if line.split() else "")
        elif "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)" in line:
            summary["hbondCount"] = int_or_none(line.split()[0] if line.split() else "")
        elif "#  RESIDUE AA STRUCTURE" in line:
            in_rows = True
            continue

        if not in_rows or not line.strip():
            continue

        dssp_index = int_or_none(line[0:5])
        if dssp_index is None:
            continue

        aa = line[13:14].strip()
        if aa == "!":
            continue

        structure = line[16:17].strip() or "-"
        nh_o_1 = parse_bond_field(line[39:50])
        nh_o_2 = parse_bond_field(line[61:72])
        o_hn_1 = parse_bond_field(line[50:61])
        o_hn_2 = parse_bond_field(line[72:83])
        row = {
            "dssp_index": dssp_index,
            "chain": line[11:12].strip() or "_",
            "residue_number": line[5:10].strip(),
            "insertion_code": line[10:11].strip(),
            "aa": aa,
            "structure": structure,
            "structure_name": STRUCTURE_NAMES.get(structure, "unknown"),
            "bridge_1": int_or_none(line[25:29]),
            "bridge_2": int_or_none(line[29:33]),
            "asa": float_or_none(line[34:38]),
            "rsa": None,
            "phi": float_or_none(line[103:109]),
            "psi": float_or_none(line[109:115]),
            "kappa": float_or_none(line[91:97]),
            "alpha": float_or_none(line[97:103]),
            "tco": float_or_none(line[83:91]),
            "x_ca": float_or_none(line[115:122]),
            "y_ca": float_or_none(line[122:129]),
            "z_ca": float_or_none(line[129:136]),
            "nh_o_1_offset": bond_value(nh_o_1, "offset"),
            "nh_o_1_energy": bond_value(nh_o_1, "energy"),
            "nh_o_2_offset": bond_value(nh_o_2, "offset"),
            "nh_o_2_energy": bond_value(nh_o_2, "energy"),
            "o_hn_1_offset": bond_value(o_hn_1, "offset"),
            "o_hn_1_energy": bond_value(o_hn_1, "energy"),
            "o_hn_2_offset": bond_value(o_hn_2, "offset"),
            "o_hn_2_energy": bond_value(o_hn_2, "energy"),
        }
        rows.append(row)

    if "residueCount" not in summary:
        summary["residueCount"] = len(rows)
    if "chainCount" not in summary:
        summary["chainCount"] = len({row["chain"] for row in rows})
    if "hbondCount" not in summary:
        summary["hbondCount"] = None

    summary["counts"] = count_structures(rows)
    if "asaTotal" not in summary:
        summary["asaTotal"] = sum(row["asa"] for row in rows if is_number(row["asa"]))

    return {
        "source": "real",
        "summary": summary,
        "sequence": "".join(row["aa"] for row in rows),
        "secondary": "".join(row["structure"] for row in rows),
        "residues": rows,
    }


def parse_bond_field(text: str) -> dict | None:
    stripped = text.strip()
    if not stripped or "," not in stripped:
        return None
    offset_text, energy_text = stripped.split(",", 1)
    return {
        "offset": int_or_none(offset_text),
        "energy": float_or_none(energy_text),
    }


def bond_value(bond: dict | None, key: str) -> int | float | None:
    if not bond:
        return None
    return bond.get(key)


def format_real_output(parsed: dict, output_format: str) -> str:
    if output_format == "json":
        return json.dumps(json_result(parsed), indent=2) + "\n"
    return rows_to_csv(parsed["residues"]) + "\n"


def json_result(parsed: dict) -> dict:
    summary = parsed.get("summary", {})
    return {
        "source": parsed.get("source"),
        "dssp_bin": parsed.get("dssp_bin"),
        "summary": {
            "residueCount": summary.get("residueCount"),
            "chainCount": summary.get("chainCount"),
            "hbondCount": summary.get("hbondCount"),
            "asaTotal": round_or_none(summary.get("asaTotal"), 1),
            "meanRsa": round_or_none(summary.get("meanRsa"), 4),
            "counts": summary.get("counts", {}),
        },
        "sequence": parsed.get("sequence", ""),
        "secondary": parsed.get("secondary", ""),
        "residues": parsed.get("residues", []),
    }


def rows_to_csv(rows: list[dict]) -> str:
    columns = [
        "dssp_index",
        "chain",
        "residue_number",
        "insertion_code",
        "aa",
        "structure",
        "structure_name",
        "bridge_1",
        "bridge_2",
        "asa",
        "rsa",
        "phi",
        "psi",
        "kappa",
        "alpha",
        "tco",
        "x_ca",
        "y_ca",
        "z_ca",
        "nh_o_1_offset",
        "nh_o_1_energy",
        "nh_o_2_offset",
        "nh_o_2_energy",
        "o_hn_1_offset",
        "o_hn_1_energy",
        "o_hn_2_offset",
        "o_hn_2_energy",
    ]
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=columns)
    writer.writeheader()
    for row in rows:
        writer.writerow({key: csv_value(key, row.get(key)) for key in columns})
    return out.getvalue().rstrip("\n")


def csv_value(column: str, value: object) -> object:
    if value is None:
        return ""
    decimals = {
        "asa": 1,
        "rsa": 4,
        "phi": 1,
        "psi": 1,
        "kappa": 1,
        "alpha": 1,
        "tco": 3,
        "x_ca": 1,
        "y_ca": 1,
        "z_ca": 1,
        "nh_o_1_energy": 1,
        "nh_o_2_energy": 1,
        "o_hn_1_energy": 1,
        "o_hn_2_energy": 1,
    }.get(column)
    if decimals is not None and is_number(value):
        return f"{float(value):.{decimals}f}"
    return value


def compare_results(js_result: dict, real_result: dict, pdb_path: str, args: argparse.Namespace) -> dict:
    js_rows = js_result["residues"]
    real_rows = real_result["residues"]
    real_by_key = {residue_key(row): row for row in real_rows}

    matched_pairs = []
    unmatched_js = []
    for js_row in js_rows:
        key = residue_key(js_row)
        real_row = real_by_key.get(key)
        if real_row:
            matched_pairs.append((js_row, real_row))
        elif len(unmatched_js) < 20:
            unmatched_js.append(display_key(js_row))

    js_keys = {residue_key(row) for row in js_rows}
    unmatched_real = [display_key(row) for row in real_rows if residue_key(row) not in js_keys][:20]

    ss_matches = 0
    ss_differences = []
    aa_differences = []
    totals = {key: 0.0 for key in ("asa", "phi", "psi", "kappa", "alpha", "tco")}
    counts = {key: 0 for key in totals}

    for js_row, real_row in matched_pairs:
        js_ss = js_row.get("structure") or "-"
        real_ss = real_row.get("structure") or "-"
        if js_ss == real_ss:
            ss_matches += 1
        elif len(ss_differences) < 25:
            ss_differences.append(
                {
                    "residue": display_key(real_row),
                    "aa": real_row.get("aa"),
                    "js": js_ss,
                    "real": real_ss,
                }
            )

        if (js_row.get("aa") or "").upper() != (real_row.get("aa") or "").upper() and len(aa_differences) < 25:
            aa_differences.append(
                {
                    "residue": display_key(real_row),
                    "js": js_row.get("aa"),
                    "real": real_row.get("aa"),
                }
            )

        for metric in totals:
            add_metric_difference(totals, counts, metric, js_row.get(metric), real_row.get(metric))

    compared = len(matched_pairs)
    real_hbond_count = real_result["summary"].get("hbondCount")
    js_hbond_count = js_result["summary"].get("hbondCount")
    real_asa_total = real_result["summary"].get("asaTotal")
    js_asa_total = js_result["summary"].get("asaTotal")

    return {
        "pdb": str(Path(pdb_path).resolve()),
        "chain": args.chain or None,
        "dssp_bin": real_result.get("dssp_bin"),
        "js": {
            "residue_count": len(js_rows),
            "hbond_count": js_hbond_count,
            "asa_total": js_asa_total,
            "structure_counts": count_structures(js_rows),
        },
        "real": {
            "residue_count": len(real_rows),
            "hbond_count": real_hbond_count,
            "asa_total": real_asa_total,
            "structure_counts": count_structures(real_rows),
        },
        "matched_residues": compared,
        "secondary_structure": {
            "matches": ss_matches,
            "agreement": ss_matches / compared if compared else None,
            "first_differences": ss_differences,
        },
        "mean_absolute_difference": {
            metric: totals[metric] / counts[metric] if counts[metric] else None for metric in totals
        },
        "deltas": {
            "hbond_count": numeric_delta(js_hbond_count, real_hbond_count),
            "asa_total": numeric_delta(js_asa_total, real_asa_total),
        },
        "aa_differences": aa_differences,
        "unmatched_js": unmatched_js,
        "unmatched_real": unmatched_real,
    }


def add_metric_difference(totals: dict, counts: dict, metric: str, left: object, right: object) -> None:
    if not is_number(left) or not is_number(right):
        return
    if metric != "asa" and (abs(float(left)) >= 359.9 or abs(float(right)) >= 359.9):
        return
    totals[metric] += abs(float(left) - float(right))
    counts[metric] += 1


def format_benchmark_report(report: dict) -> str:
    agreement = report["secondary_structure"]["agreement"]
    agreement_text = "n/a" if agreement is None else f"{agreement * 100:.1f}%"
    lines = [
        f"Benchmark PDB: {report['pdb']}",
        f"DSSP binary: {report['dssp_bin']}",
        f"Chain: {report['chain'] or 'all'}",
        (
            "Residues: "
            f"JS {report['js']['residue_count']}, "
            f"real {report['real']['residue_count']}, "
            f"matched {report['matched_residues']}"
        ),
        (
            "Secondary-structure agreement: "
            f"{agreement_text} "
            f"({report['secondary_structure']['matches']}/{report['matched_residues']})"
        ),
        (
            "H-bonds: "
            f"JS {format_optional(report['js']['hbond_count'])}, "
            f"real {format_optional(report['real']['hbond_count'])}, "
            f"delta {format_optional(report['deltas']['hbond_count'])}"
        ),
        (
            "ASA total: "
            f"JS {format_optional(report['js']['asa_total'], 1)}, "
            f"real {format_optional(report['real']['asa_total'], 1)}, "
            f"delta {format_optional(report['deltas']['asa_total'], 1)}"
        ),
        "Mean absolute differences:",
    ]

    for metric, value in report["mean_absolute_difference"].items():
        decimals = 3 if metric == "tco" else 2
        lines.append(f"  {metric}: {format_optional(value, decimals)}")

    lines.append("Secondary-structure counts:")
    codes = sorted(set(report["js"]["structure_counts"]) | set(report["real"]["structure_counts"]))
    for code in codes:
        lines.append(
            f"  {code}: JS {report['js']['structure_counts'].get(code, 0)}, "
            f"real {report['real']['structure_counts'].get(code, 0)}"
        )

    differences = report["secondary_structure"]["first_differences"]
    if differences:
        lines.append("First SS differences:")
        for diff in differences:
            lines.append(f"  {diff['residue']} {diff['aa']}: JS {diff['js']}, real {diff['real']}")

    if report["aa_differences"]:
        lines.append("First AA differences:")
        for diff in report["aa_differences"]:
            lines.append(f"  {diff['residue']}: JS {diff['js']}, real {diff['real']}")

    if report["unmatched_js"]:
        lines.append("Unmatched JS residues: " + ", ".join(report["unmatched_js"]))
    if report["unmatched_real"]:
        lines.append("Unmatched real DSSP residues: " + ", ".join(report["unmatched_real"]))

    return "\n".join(lines) + "\n"


def count_structures(rows: list[dict]) -> dict:
    counts = {code: 0 for code in STRUCTURE_NAMES}
    for row in rows:
        code = row.get("structure") or "-"
        counts[code] = counts.get(code, 0) + 1
    return counts


def residue_key(row: dict) -> tuple[str, str, str]:
    return (
        str(row.get("chain") or "_"),
        str(row.get("residue_number") or ""),
        str(row.get("insertion_code") or ""),
    )


def display_key(row: dict) -> str:
    insertion = row.get("insertion_code") or ""
    return f"{row.get('chain') or '_'}:{row.get('residue_number') or '?'}{insertion}"


def numeric_delta(left: object, right: object) -> float | None:
    if not is_number(left) or not is_number(right):
        return None
    return float(left) - float(right)


def is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def round_or_none(value: object, digits: int) -> float | None:
    if not is_number(value):
        return None
    return round(float(value), digits)


def int_or_none(text: object) -> int | None:
    try:
        stripped = str(text).strip()
        return int(stripped) if stripped else None
    except ValueError:
        return None


def float_or_none(text: object) -> float | None:
    try:
        stripped = str(text).strip().replace(",", ".")
        value = float(stripped) if stripped else math.nan
        return value if math.isfinite(value) else None
    except ValueError:
        return None


def format_optional(value: object, decimals: int = 0) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, int) or decimals == 0:
        return str(round(float(value))) if isinstance(value, float) else str(value)
    if is_number(value):
        return f"{float(value):.{decimals}f}"
    return str(value)


def clean_process_error(prefix: str, completed: subprocess.CompletedProcess) -> str:
    details = [f"{prefix} (exit {completed.returncode})"]
    if completed.stderr.strip():
        details.append(completed.stderr.strip())
    if completed.stdout.strip():
        details.append(completed.stdout.strip()[:1000])
    return "\n".join(details)


def write_output(text: str, output_path: str | None) -> None:
    if output_path:
        Path(output_path).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    raise SystemExit(main())
