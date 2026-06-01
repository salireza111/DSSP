import { execFileSync } from "node:child_process";
import { mkdtempSync, readFileSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { basename, join, resolve } from "node:path";
import { assignDSSP } from "../src/dssp.js";

const pdbPath = resolve(process.argv[2] ?? "assets/1YIO.pdb");
const mkdsspPath = process.env.MKDSSP_BIN || process.argv[3] || "mkdssp";
const tempDir = mkdtempSync(join(tmpdir(), "dssp-validate-"));
const dsspPath = join(tempDir, `${basename(pdbPath)}.dssp`);

try {
  execFileSync(mkdsspPath, ["--output-format", "dssp", pdbPath, dsspPath], { stdio: ["ignore", "ignore", "pipe"] });
  const pdbText = readFileSync(pdbPath, "utf8");
  const jsResult = assignDSSP(pdbText);
  const reference = parseClassicDssp(readFileSync(dsspPath, "utf8"));
  const report = compareResults(jsResult.residues, reference);

  printReport(pdbPath, mkdsspPath, report);
  if (report.residueCountCompared === 0) {
    process.exitCode = 2;
  }
} finally {
  rmSync(tempDir, { force: true, recursive: true });
}

function parseClassicDssp(text) {
  const rows = [];
  let inRows = false;
  for (const line of text.split(/\r?\n/)) {
    if (line.includes("#  RESIDUE AA STRUCTURE")) {
      inRows = true;
      continue;
    }
    if (!inRows || !line.trim()) continue;
    const dsspIndex = Number.parseInt(line.slice(0, 5), 10);
    if (!Number.isFinite(dsspIndex)) continue;

    rows.push({
      dsspIndex,
      resSeq: line.slice(5, 10).trim(),
      chainId: line.slice(11, 12).trim() || "_",
      aa: line.slice(13, 14).trim(),
      ss: line.slice(16, 17).trim() || "-",
      bridge1: numberOrNull(line.slice(25, 29)),
      bridge2: numberOrNull(line.slice(29, 33)),
      asa: numberOrNull(line.slice(34, 38)),
      tco: numberOrNull(line.slice(83, 91)),
      kappa: numberOrNull(line.slice(91, 97)),
      alpha: numberOrNull(line.slice(97, 103)),
      phi: numberOrNull(line.slice(103, 109)),
      psi: numberOrNull(line.slice(109, 115))
    });
  }
  return rows;
}

function compareResults(jsRows, referenceRows) {
  const count = Math.min(jsRows.length, referenceRows.length);
  const differences = [];
  const totals = {
    asa: 0,
    phi: 0,
    psi: 0,
    kappa: 0,
    alpha: 0,
    tco: 0
  };
  const counts = Object.fromEntries(Object.keys(totals).map((key) => [key, 0]));
  let ssMatches = 0;

  for (let i = 0; i < count; i += 1) {
    const js = jsRows[i];
    const ref = referenceRows[i];
    const jsSs = js.ss || "-";
    if (jsSs === ref.ss) ssMatches += 1;
    else if (differences.length < 20) {
      differences.push({
        index: i + 1,
        residue: `${ref.chainId}:${ref.resSeq}`,
        aa: ref.aa,
        js: jsSs,
        mkdssp: ref.ss
      });
    }

    addDiff(totals, counts, "asa", js.asa, ref.asa);
    addDiff(totals, counts, "phi", js.phi, ref.phi);
    addDiff(totals, counts, "psi", js.psi, ref.psi);
    addDiff(totals, counts, "kappa", js.kappa, ref.kappa);
    addDiff(totals, counts, "alpha", js.alpha, ref.alpha);
    addDiff(totals, counts, "tco", js.tco, ref.tco);
  }

  return {
    residueCountJs: jsRows.length,
    residueCountMkdssp: referenceRows.length,
    residueCountCompared: count,
    ssMatches,
    ssAgreement: count ? ssMatches / count : 0,
    meanAbsoluteDifference: Object.fromEntries(Object.entries(totals).map(([key, value]) => [
      key,
      counts[key] ? value / counts[key] : null
    ])),
    firstSecondaryStructureDifferences: differences
  };
}

function addDiff(totals, counts, key, left, right) {
  if (!Number.isFinite(left) || !Number.isFinite(right)) return;
  if (Math.abs(right) >= 359.9 && key !== "asa") return;
  totals[key] += Math.abs(left - right);
  counts[key] += 1;
}

function printReport(pdbPath, mkdsspPath, report) {
  const percent = (100 * report.ssAgreement).toFixed(1);
  console.log(`Validated ${pdbPath} against ${mkdsspPath}`);
  console.log(`Residues compared: ${report.residueCountCompared} (JS ${report.residueCountJs}, mkdssp ${report.residueCountMkdssp})`);
  console.log(`Secondary-structure agreement: ${percent}% (${report.ssMatches}/${report.residueCountCompared})`);
  console.log("Mean absolute differences:");
  for (const [metric, value] of Object.entries(report.meanAbsoluteDifference)) {
    console.log(`  ${metric}: ${value == null ? "n/a" : value.toFixed(metric === "tco" ? 3 : 2)}`);
  }
  if (report.firstSecondaryStructureDifferences.length > 0) {
    console.log("First SS differences:");
    for (const diff of report.firstSecondaryStructureDifferences) {
      console.log(`  ${diff.index} ${diff.residue} ${diff.aa}: JS ${diff.js}, mkdssp ${diff.mkdssp}`);
    }
  }
}

function numberOrNull(text) {
  const value = Number.parseFloat(String(text).replace(",", ".").trim());
  return Number.isFinite(value) ? value : null;
}
