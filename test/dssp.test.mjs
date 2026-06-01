import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import {
  DSSP_COLUMNS,
  assignDSSP,
  calculateSolventAccessibility,
  filterPDBByChain,
  hydrogenBondEnergy,
  parsePDB,
  toCSV,
  toDsspJSON,
  toDsspRows
} from "../src/dssp.js";

const miniPdb = `ATOM      1  N   ALA A   1       0.000   1.200   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.200   1.700   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.100   0.600   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.800  -0.600   0.000  1.00 20.00           O
ATOM      5  N   GLY A   2       3.300   0.900   0.000  1.00 20.00           N
ATOM      6  CA  GLY A   2       4.200  -0.100   0.000  1.00 20.00           C
ATOM      7  C   GLY A   2       5.600   0.400   0.000  1.00 20.00           C
ATOM      8  O   GLY A   2       5.900   1.600   0.000  1.00 20.00           O
END`;

test("parsePDB groups backbone atoms into residues", () => {
  const residues = parsePDB(miniPdb);
  assert.equal(residues.length, 2);
  assert.equal(residues[0].aa, "A");
  assert.equal(residues[1].aa, "G");
  assert.ok(residues[0].atoms.CA);
  assert.ok(residues[1].atoms.O);
});

test("filterPDBByChain keeps only the requested legacy PDB chain", () => {
  const chainA = miniPdb.replace(/\nEND$/, "");
  const chainB = chainA.replaceAll(" A   ", " B   ");
  const mixedChainPdb = `${chainA}\n${chainB}\nCONECT    1    2\nEND`;
  const filtered = filterPDBByChain(mixedChainPdb, "b");
  const residues = parsePDB(filtered);

  assert.equal(residues.length, 2);
  assert.deepEqual([...new Set(residues.map((residue) => residue.chainId))], ["B"]);
  assert.doesNotMatch(filtered, /^CONECT/m);
});

test("filterPDBByChain reports a missing chain", () => {
  assert.throws(() => filterPDBByChain(miniPdb, "Z"), /chain Z/);
});

test("assignDSSP returns sequence, secondary string, and inferred hydrogens", () => {
  const result = assignDSSP(miniPdb);
  assert.equal(result.sequence, "AG");
  assert.equal(result.secondary.length, 2);
  assert.ok(result.residues[1].atoms.H.inferred);
  assert.equal(result.summary.residueCount, 2);
  assert.ok(result.summary.asaTotal > 0);
  assert.ok(result.summary.meanRsa > 0);
});

test("calculateSolventAccessibility adds residue ASA and RSA values", () => {
  const residues = parsePDB(miniPdb);
  calculateSolventAccessibility(residues, { asaSamples: 48 });

  assert.ok(residues[0].asa > 0);
  assert.ok(residues[0].rsa > 0);
  assert.ok(residues[0].atomRecords.some((atom) => atom.asa > 0));
});

test("hydrogenBondEnergy follows the DSSP electrostatic sign convention", () => {
  const acceptor = {
    atoms: {
      C: { x: 0, y: 0, z: 0 },
      O: { x: 1.23, y: 0, z: 0 }
    }
  };
  const donor = {
    atoms: {
      N: { x: 2.85, y: 0, z: 0 },
      H: { x: 1.85, y: 0, z: 0 }
    }
  };
  assert.ok(hydrogenBondEnergy(acceptor, donor) < -0.5);
  assert.equal(hydrogenBondEnergy(acceptor, { ...donor, isProline: true }), Number.POSITIVE_INFINITY);
});

test("toCSV includes residue rows and hydrogen-bond columns", () => {
  const csv = toCSV(assignDSSP(miniPdb));
  const [header] = csv.split("\n");
  assert.equal(header, DSSP_COLUMNS.join(","));
  assert.match(csv, /1,A,1,,A/);
  assert.match(csv, /x_ca,y_ca,z_ca/);
});

test("toDsspJSON uses the same residue schema as CSV", () => {
  const result = assignDSSP(miniPdb);
  const json = toDsspJSON(result);
  const rows = toDsspRows(result);

  assert.equal(json.source, "js");
  assert.equal(json.dssp_bin, null);
  assert.deepEqual(Object.keys(json.residues[0]), DSSP_COLUMNS);
  assert.deepEqual(json.residues, rows);
  assert.equal(json.residues[0].nh_o_1_offset, 0);
  assert.equal(json.residues[0].nh_o_1_energy, 0);
});

test("assignDSSP identifies helices and sheets in the bundled 1YIO demo", () => {
  const demoPdb = readFileSync(new URL("../assets/1YIO.pdb", import.meta.url), "utf8");
  const result = assignDSSP(demoPdb);
  assert.equal(result.residues.length, 198);
  assert.ok(result.summary.counts.H > 0);
  assert.ok(result.summary.counts.E > 0);
  assert.ok(result.summary.asaTotal > 0);
});

function syntheticHelixPdb() {
  const residues = ["ALA", "GLU", "LEU", "LYS", "GLN", "ALA", "ARG", "TYR"];
  const nitrogen = residues.map((_, i) => {
    const theta = i * 100 * Math.PI / 180;
    return { x: 3 * Math.cos(theta), y: 3 * Math.sin(theta), z: i * 1.5 };
  });
  const oxygens = residues.map((_, i) => {
    const target = nitrogen[i + 4] ?? nitrogen[i];
    const direction = normalize({ x: Math.cos(i), y: Math.sin(i), z: 0.2 });
    return {
      x: target.x + 2.8 * direction.x,
      y: target.y + 2.8 * direction.y,
      z: target.z + 2.8 * direction.z
    };
  });

  let serial = 1;
  const lines = ["HEADER    SYNTHETIC ALPHA HELIX TEST"];
  for (let i = 0; i < residues.length; i += 1) {
    const n = nitrogen[i];
    const next = nitrogen[i + 1] ?? { x: n.x + 1, y: n.y, z: n.z };
    const chainDirection = normalize({ x: next.x - n.x, y: next.y - n.y, z: next.z - n.z });
    const c = {
      x: next.x - 1.2 * chainDirection.x,
      y: next.y - 1.2 * chainDirection.y,
      z: next.z - 1.2 * chainDirection.z
    };
    const ca = { x: (n.x + c.x) / 2, y: (n.y + c.y) / 2, z: (n.z + c.z) / 2 + 0.7 };
    let h = { x: n.x + 0.8, y: n.y, z: n.z };
    if (i >= 4) {
      const towardAcceptor = normalize({
        x: oxygens[i - 4].x - n.x,
        y: oxygens[i - 4].y - n.y,
        z: oxygens[i - 4].z - n.z
      });
      h = {
        x: n.x + 1.01 * towardAcceptor.x,
        y: n.y + 1.01 * towardAcceptor.y,
        z: n.z + 1.01 * towardAcceptor.z
      };
    }

    lines.push(atomLine(serial++, "N", residues[i], i + 1, n, "N"));
    lines.push(atomLine(serial++, "H", residues[i], i + 1, h, "H"));
    lines.push(atomLine(serial++, "CA", residues[i], i + 1, ca, "C"));
    lines.push(atomLine(serial++, "C", residues[i], i + 1, c, "C"));
    lines.push(atomLine(serial++, "O", residues[i], i + 1, oxygens[i], "O"));
  }
  lines.push("END");
  return lines.join("\n");
}

function atomLine(serial, atomName, resName, resSeq, point, element) {
  return `${"ATOM".padEnd(6)}${String(serial).padStart(5)} ${atomName.padStart(4)} ${resName.padStart(3)} A${String(resSeq).padStart(4)}    ${point.x.toFixed(3).padStart(8)}${point.y.toFixed(3).padStart(8)}${point.z.toFixed(3).padStart(8)}  1.00 20.00          ${element.padStart(2)}`;
}

function normalize(vector) {
  const length = Math.hypot(vector.x, vector.y, vector.z) || 1;
  return { x: vector.x / length, y: vector.y / length, z: vector.z / length };
}
