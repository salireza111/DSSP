const HBOND_Q1Q2 = 0.084;
const HBOND_F = 332;
const DEFAULT_HBOND_CUTOFF = -0.5;
const PEPTIDE_BOND_MAX = 2.5;
const MIN_DISTANCE = 0.5;
const MIN_HBOND_ENERGY = -9.9;
const MAX_CA_DISTANCE = 9.0;
const WATER_PROBE_RADIUS = 1.4;
const DEFAULT_ASA_SAMPLES = 401;
const CHAIN_FILTER_RECORDS = new Set(["ATOM", "HETATM", "ANISOU", "TER"]);

const DSSP_ATOM_RADII = {
  N: 1.65,
  CA: 1.87,
  C: 1.76,
  O: 1.4,
  SIDE: 1.8
};

const MAX_ASA = {
  A: 129,
  R: 274,
  N: 195,
  D: 193,
  C: 167,
  Q: 225,
  E: 223,
  G: 104,
  H: 224,
  I: 197,
  L: 201,
  K: 236,
  M: 224,
  F: 240,
  P: 159,
  S: 155,
  T: 172,
  W: 285,
  Y: 263,
  V: 174,
  U: 167,
  O: 236
};

const sphereSampleCache = new Map();

const AA3_TO_1 = {
  ALA: "A",
  ARG: "R",
  ASN: "N",
  ASP: "D",
  CYS: "C",
  GLN: "Q",
  GLU: "E",
  GLY: "G",
  HIS: "H",
  ILE: "I",
  LEU: "L",
  LYS: "K",
  MET: "M",
  PHE: "F",
  PRO: "P",
  SER: "S",
  THR: "T",
  TRP: "W",
  TYR: "Y",
  VAL: "V",
  SEC: "U",
  PYL: "O",
  ASX: "B",
  GLX: "Z",
  MSE: "M"
};

const STRUCTURE_NAMES = {
  H: "alpha helix",
  B: "isolated beta bridge",
  E: "extended beta strand",
  G: "3-10 helix",
  I: "pi helix",
  P: "poly-proline II helix",
  T: "hydrogen-bonded turn",
  S: "bend",
  "-": "loop"
};

export const DSSP_COLUMNS = [
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
  "o_hn_2_energy"
];

const DSSP_DECIMALS = {
  asa: 1,
  rsa: 4,
  phi: 1,
  psi: 1,
  kappa: 1,
  alpha: 1,
  tco: 3,
  x_ca: 1,
  y_ca: 1,
  z_ca: 1,
  nh_o_1_energy: 1,
  nh_o_2_energy: 1,
  o_hn_1_energy: 1,
  o_hn_2_energy: 1
};

export function filterPDBByChain(pdbText, chainId) {
  const requestedChain = normalizeChainId(chainId);
  const text = String(pdbText ?? "");
  if (!requestedChain) return text;

  let keptAtoms = 0;
  const filteredLines = text.split(/\r?\n/).filter((line) => {
    const record = line.slice(0, 6).trim();
    if (record === "CONECT" || record === "MASTER") return false;
    if (!CHAIN_FILTER_RECORDS.has(record)) return true;

    const lineChain = (line.slice(21, 22).trim() || "_").toUpperCase();
    if (lineChain !== requestedChain) return false;
    if (record === "ATOM") keptAtoms += 1;
    return true;
  });

  if (keptAtoms === 0) {
    throw new Error(`No ATOM records were found for chain ${requestedChain}.`);
  }
  return filteredLines.join("\n");
}

function normalizeChainId(chainId) {
  const normalized = String(chainId ?? "").trim().toUpperCase();
  if (!normalized) return "";
  if (normalized.length !== 1) {
    throw new Error("Legacy PDB chain IDs are one character.");
  }
  return normalized;
}

export function parsePDB(pdbText, options = {}) {
  const includeHetatm = options.includeHetatm ?? false;
  const keepAltLocs = new Set(["", "A", "1"]);
  const residues = [];
  const residueByKey = new Map();
  let inFirstModel = false;
  let sawModel = false;

  const lines = pdbText.split(/\r?\n/);
  for (const line of lines) {
    const record = line.slice(0, 6).trim();
    if (record === "MODEL") {
      if (sawModel) break;
      sawModel = true;
      inFirstModel = true;
      continue;
    }
    if (record === "ENDMDL" && sawModel) break;
    if (sawModel && !inFirstModel) continue;
    if (record !== "ATOM" && !(includeHetatm && record === "HETATM")) continue;

    const atomName = line.slice(12, 16).trim();
    const altLoc = line.slice(16, 17).trim();
    if (!keepAltLocs.has(altLoc)) continue;

    const resName = line.slice(17, 20).trim().toUpperCase();
    const chainId = line.slice(21, 22).trim() || "_";
    const resSeq = line.slice(22, 26).trim();
    const iCode = line.slice(26, 27).trim();
    const x = Number.parseFloat(line.slice(30, 38));
    const y = Number.parseFloat(line.slice(38, 46));
    const z = Number.parseFloat(line.slice(46, 54));
    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
      continue;
    }

    const occupancy = Number.parseFloat(line.slice(54, 60));
    const key = `${chainId}\u0000${resSeq}\u0000${iCode}\u0000${resName}`;
    let residue = residueByKey.get(key);
    if (!residue) {
      residue = {
        index: residues.length,
        dsspIndex: residues.length + 1,
        chainId,
        resSeq,
        iCode,
        resName,
        aa: AA3_TO_1[resName] ?? "X",
        atoms: {},
        atomRecords: []
      };
      residues.push(residue);
      residueByKey.set(key, residue);
    }

    const atom = {
      name: atomName,
      altLoc,
      x,
      y,
      z,
      occupancy: Number.isFinite(occupancy) ? occupancy : 0,
      element: inferElement(line, atomName)
    };
    residue.atomRecords.push(atom);
    keepBestAtom(residue.atoms, normalizedAtomName(atomName), atom);
  }

  return annotateBackbone(residues);
}

export function assignDSSP(input, options = {}) {
  const residues = Array.isArray(input) ? cloneResidues(input) : parsePDB(String(input ?? ""), options);
  annotateBackbone(residues);
  assignHydrogens(residues);
  calculateAngles(residues);
  calculateSolventAccessibility(residues, options);

  const hbondCutoff = options.hbondCutoff ?? DEFAULT_HBOND_CUTOFF;
  const candidatePairs = residuePairsWithinCaDistance(residues);
  const { energyMatrix, hbonds } = calculateHBondEnergies(residues, candidatePairs, hbondCutoff);
  const hbondMap = buildHbondMapFromBestColumns(residues, hbondCutoff);

  assignBetaSheets(residues, hbondMap, candidatePairs);
  assignHelicesTurnsAndBends(residues, hbondMap, options);
  assignPPHelices(residues, options.minPolyProlineStretchLength ?? 3);

  return {
    residues,
    hbonds,
    energyMatrix,
    sequence: residues.map((residue) => residue.aa).join(""),
    secondary: residues.map((residue) => residue.ss).join(""),
    summary: summarize(residues, hbonds),
    options: { hbondCutoff }
  };
}

export function hydrogenBondEnergy(acceptorResidue, donorResidue) {
  if (donorResidue?.isProline) return Number.POSITIVE_INFINITY;
  const o = acceptorResidue?.atoms?.O;
  const c = acceptorResidue?.atoms?.C;
  const n = donorResidue?.atoms?.N;
  const h = donorResidue?.atoms?.H;
  if (!o || !c || !n || !h) return Number.POSITIVE_INFINITY;

  const dON = distance(o, n);
  const dCH = distance(c, h);
  const dOH = distance(o, h);
  const dCN = distance(c, n);
  if ([dON, dCH, dOH, dCN].some((value) => value < MIN_DISTANCE)) {
    return MIN_HBOND_ENERGY;
  }

  const energy = HBOND_Q1Q2 * HBOND_F * (1 / dON + 1 / dCH - 1 / dOH - 1 / dCN);
  return Math.max(MIN_HBOND_ENERGY, Math.round(energy * 1000) / 1000);
}

export function calculateSolventAccessibility(residues, options = {}) {
  const probeRadius = options.probeRadius ?? WATER_PROBE_RADIUS;
  const sampleCount = Math.max(24, Math.trunc(options.asaSamples ?? DEFAULT_ASA_SAMPLES));
  const sampleDirections = sphereSampleDirections(sampleCount);

  for (const residue of residues) {
    residue.asa = 0;
    residue.rsa = null;
    for (const atom of residue.atomRecords ?? []) {
      atom.asa = 0;
    }
  }

  const atoms = [];
  for (const residue of residues) {
    for (const atom of residue.atomRecords ?? []) {
      const element = normalizedElement(atom);
      if (!element || element === "H" || element === "D") continue;
      const vdwRadius = dsspAtomRadius(atom);
      atoms.push({
        atom,
        residue,
        radius: vdwRadius + probeRadius
      });
    }
  }

  if (atoms.length === 0) return residues;

  const maxExpandedRadius = Math.max(...atoms.map((entry) => entry.radius));
  const cellSize = maxExpandedRadius * 2;
  const grid = buildAtomGrid(atoms, cellSize);

  for (const entry of atoms) {
    const neighbors = atomNeighbors(entry, grid, cellSize);
    let exposed = 0;
    for (const direction of sampleDirections) {
      const point = {
        x: entry.atom.x + direction.x * entry.radius,
        y: entry.atom.y + direction.y * entry.radius,
        z: entry.atom.z + direction.z * entry.radius
      };
      if (!isBuriedPoint(point, neighbors)) exposed += 1;
    }

    const atomAsa = 4 * Math.PI * entry.radius * entry.radius * exposed / sampleCount;
    entry.atom.asa = atomAsa;
    entry.residue.asa += atomAsa;
  }

  for (const residue of residues) {
    const maxAsa = MAX_ASA[residue.aa];
    residue.rsa = maxAsa ? residue.asa / maxAsa : null;
  }

  return residues;
}

export function toCSV(result) {
  const rows = [DSSP_COLUMNS];
  for (const residue of toDsspRows(result)) {
    rows.push(DSSP_COLUMNS.map((column) => formatDsspCell(column, residue[column])));
  }

  return rows.map((row) => row.map(csvCell).join(",")).join("\n");
}

export function toDsspJSON(result, options = {}) {
  return {
    source: options.source ?? "js",
    dssp_bin: options.dsspBin ?? null,
    summary: dsspSummary(result.summary),
    sequence: result.sequence,
    secondary: result.secondary,
    residues: toDsspRows(result)
  };
}

export function toDsspRows(result) {
  const rows = [];
  for (const residue of result.residues) {
    const ca = residue.atoms?.CA;
    const nhO1 = residue.donorHBonds?.[0];
    const nhO2 = residue.donorHBonds?.[1];
    const oHn1 = residue.acceptorHBonds?.[0];
    const oHn2 = residue.acceptorHBonds?.[1];

    rows.push([
      residue.dsspIndex,
      residue.chainId,
      residue.resSeq,
      residue.iCode,
      residue.aa,
      residue.ss,
      STRUCTURE_NAMES[residue.ss] ?? "",
      residue.bridgePartners?.[0] ?? 0,
      residue.bridgePartners?.[1] ?? 0,
      roundNumber(residue.asa, 1),
      roundNumber(residue.rsa, 4),
      dsspAngleValue(residue.phi),
      dsspAngleValue(residue.psi),
      dsspAngleValue(residue.kappa),
      dsspAngleValue(residue.alpha),
      roundNumber(residue.tco, 3) ?? 0,
      roundNumber(ca?.x, 1),
      roundNumber(ca?.y, 1),
      roundNumber(ca?.z, 1),
      formatBondOffset(nhO1, residue),
      formatBondEnergy(nhO1),
      formatBondOffset(nhO2, residue),
      formatBondEnergy(nhO2),
      formatBondOffset(oHn1, residue),
      formatBondEnergy(oHn1),
      formatBondOffset(oHn2, residue),
      formatBondEnergy(oHn2)
    ]);
  }

  return rows.map((row) => Object.fromEntries(DSSP_COLUMNS.map((column, index) => [column, row[index]])));
}

export function formatDsspCell(column, value) {
  if (value == null) return "";
  const digits = DSSP_DECIMALS[column];
  if (digits != null && Number.isFinite(value)) return value.toFixed(digits);
  return value;
}

export function structureName(code) {
  return STRUCTURE_NAMES[code] ?? "unknown";
}

function cloneResidues(residues) {
  return residues.map((residue, index) => ({
    ...residue,
    index,
    dsspIndex: index + 1,
    atoms: { ...residue.atoms },
    atomRecords: residue.atomRecords ? [...residue.atomRecords] : []
  }));
}

function inferElement(line, atomName) {
  const fromColumn = line.slice(76, 78).trim();
  if (fromColumn) return fromColumn.toUpperCase();
  const stripped = atomName.replace(/[0-9]/g, "").trim();
  return stripped ? stripped[0].toUpperCase() : "";
}

function normalizedAtomName(atomName) {
  const name = atomName.trim().toUpperCase();
  if (name === "O1" || name === "OT1") return "O";
  if (name === "H" || name === "HN" || name === "HT1" || name === "1H") return "H";
  return name;
}

function keepBestAtom(atoms, atomName, atom) {
  if (!["N", "CA", "C", "O", "H"].includes(atomName)) return;
  const current = atoms[atomName];
  if (!current || atom.occupancy > current.occupancy || current.altLoc) {
    atoms[atomName] = atom;
  }
}

function annotateBackbone(residues) {
  for (let i = 0; i < residues.length; i += 1) {
    const residue = residues[i];
    residue.index = i;
    residue.dsspIndex = i + 1;
    residue.hasBackbone = Boolean(residue.atoms.N && residue.atoms.CA && residue.atoms.C && residue.atoms.O);
    residue.isProline = residue.resName === "PRO" || residue.aa === "P";
    residue.prev = null;
    residue.next = null;
    residue.chainBreakBefore = i === 0;
    residue.chainBreakAfter = i === residues.length - 1;
    residue.ss = "-";
    residue.preliminary = "";
    residue.isBend = false;
    residue.sheet = 0;
    residue.strand = 0;
    residue.betaPartners = [];
    residue.bridgePartners = [];
    residue.donorHBonds = [];
    residue.acceptorHBonds = [];
    residue.helixFlags = {
      G: "none",
      H: "none",
      I: "none",
      P: "none"
    };
  }

  for (let i = 1; i < residues.length; i += 1) {
    const prev = residues[i - 1];
    const residue = residues[i];
    const continuous =
      prev.chainId === residue.chainId &&
      prev.atoms.C &&
      residue.atoms.N &&
      distance(prev.atoms.C, residue.atoms.N) <= PEPTIDE_BOND_MAX;
    if (continuous) {
      prev.next = residue;
      residue.prev = prev;
      prev.chainBreakAfter = false;
      residue.chainBreakBefore = false;
    }
  }

  return residues;
}

function assignHydrogens(residues) {
  for (const residue of residues) {
    if (residue.isProline || !residue.atoms.N || !residue.prev?.atoms?.C || !residue.prev?.atoms?.O) {
      continue;
    }
    const n = residue.atoms.N;
    const cPrev = residue.prev.atoms.C;
    const oPrev = residue.prev.atoms.O;
    const direction = normalize(sub(cPrev, oPrev));
    if (!direction) continue;
    residue.atoms.H = {
      name: "H",
      inferred: true,
      x: n.x + 1.01 * direction.x,
      y: n.y + 1.01 * direction.y,
      z: n.z + 1.01 * direction.z,
      occupancy: 0,
      element: "H"
    };
  }
}

function calculateAngles(residues) {
  for (let i = 0; i < residues.length; i += 1) {
    const residue = residues[i];
    residue.phi = null;
    residue.psi = null;
    residue.kappa = null;
    residue.alpha = null;
    residue.tco = null;

    if (residue.prev?.atoms?.C && residue.atoms.N && residue.atoms.CA && residue.atoms.C) {
      residue.phi = dihedral(residue.prev.atoms.C, residue.atoms.N, residue.atoms.CA, residue.atoms.C);
    }
    if (residue.atoms.N && residue.atoms.CA && residue.atoms.C && residue.next?.atoms?.N) {
      residue.psi = dihedral(residue.atoms.N, residue.atoms.CA, residue.atoms.C, residue.next.atoms.N);
    }
    if (residue.prev?.atoms?.O && residue.prev?.atoms?.C && residue.atoms.O && residue.atoms.C) {
      const prevCO = sub(residue.prev.atoms.O, residue.prev.atoms.C);
      const thisCO = sub(residue.atoms.O, residue.atoms.C);
      residue.tco = dot(normalize(prevCO), normalize(thisCO));
    }
    if (i >= 2 && i + 2 < residues.length) {
      const a = residues[i - 2].atoms.CA;
      const b = residue.atoms.CA;
      const c = residues[i + 2].atoms.CA;
      if (a && b && c && sameContinuousSegment(residues, i - 2, i + 2)) {
        residue.kappa = angleFromVectors(sub(b, a), sub(c, b));
      }
    }
    if (i >= 1 && i + 2 < residues.length) {
      const a = residues[i - 1].atoms.CA;
      const b = residue.atoms.CA;
      const c = residues[i + 1].atoms.CA;
      const d = residues[i + 2].atoms.CA;
      if (a && b && c && d && sameContinuousSegment(residues, i - 1, i + 2)) {
        residue.alpha = dihedral(a, b, c, d);
      }
    }
  }
}

function residuePairsWithinCaDistance(residues) {
  const pairs = [];
  const maxDistanceSq = MAX_CA_DISTANCE * MAX_CA_DISTANCE;
  for (let i = 0; i + 1 < residues.length; i += 1) {
    const caI = residues[i].atoms.CA;
    if (!caI) continue;
    for (let j = i + 1; j < residues.length; j += 1) {
      const caJ = residues[j].atoms.CA;
      if (!caJ) continue;
      if (distanceSq(caI, caJ) <= maxDistanceSq) {
        pairs.push([i, j]);
      }
    }
  }
  return pairs;
}

function calculateHBondEnergies(residues, candidatePairs, cutoff) {
  const energyMatrix = Array.from({ length: residues.length }, () => Array(residues.length).fill(Number.POSITIVE_INFINITY));
  const hbonds = [];
  for (const residue of residues) {
    residue.donorHBonds = [];
    residue.acceptorHBonds = [];
  }

  for (const [i, j] of candidatePairs) {
    registerHBond(residues, energyMatrix, hbonds, i, j, cutoff);
    if (shouldCalculateReverseHBond(i, j)) {
      registerHBond(residues, energyMatrix, hbonds, j, i, cutoff);
    }
  }

  hbonds.sort((a, b) => a.energy - b.energy);
  return { energyMatrix, hbonds };
}

function shouldCalculateReverseHBond(i, j) {
  // Matches PDB-REDO DSSP: calculate i->j for each near i<j pair, skip only adjacent reverse j->i.
  return j !== i + 1;
}

function registerHBond(residues, energyMatrix, hbonds, donor, acceptor, cutoff) {
  const energy = hydrogenBondEnergy(residues[acceptor], residues[donor]);
  energyMatrix[acceptor][donor] = energy;
  const bond = {
    acceptor,
    donor,
    energy,
    offsetFromDonor: residues[acceptor].dsspIndex - residues[donor].dsspIndex,
    offsetFromAcceptor: residues[donor].dsspIndex - residues[acceptor].dsspIndex
  };

  insertBestBond(residues[donor].donorHBonds, bond);
  insertBestBond(residues[acceptor].acceptorHBonds, bond);
  if (energy < cutoff) hbonds.push(bond);
}

function insertBestBond(list, bond) {
  if (!Number.isFinite(bond.energy)) return;
  list.push(bond);
  list.sort((a, b) => a.energy - b.energy);
  if (list.length > 2) list.length = 2;
}

function buildHbondMapFromBestColumns(residues, cutoff) {
  const map = Array.from({ length: residues.length }, () => Array(residues.length).fill(false));
  for (const residue of residues) {
    for (const bond of residue.donorHBonds) {
      if (bond.energy < cutoff) {
        map[bond.acceptor][bond.donor] = true;
      }
    }
  }
  return map;
}

function assignBetaSheets(residues, hbondMap, candidatePairs) {
  const bridges = [];

  for (const [i, j] of candidatePairs) {
    const type = testBridge(residues, hbondMap, i, j);
    if (type === "none") continue;

    let found = false;
    for (const bridge of bridges) {
      if (bridge.type !== type || i !== bridge.i.at(-1) + 1) continue;

      if (type === "parallel" && bridge.j.at(-1) + 1 === j) {
        bridge.i.push(i);
        bridge.j.push(j);
        found = true;
        break;
      }

      if (type === "antiparallel" && bridge.j[0] - 1 === j) {
        bridge.i.push(i);
        bridge.j.unshift(j);
        found = true;
        break;
      }
    }

    if (!found) {
      bridges.push({
        type,
        i: [i],
        j: [j],
        chainI: residues[i].chainId,
        chainJ: residues[j].chainId,
        sheet: 0,
        ladder: 0,
        links: []
      });
    }
  }

  bridges.sort((a, b) => a.chainI.localeCompare(b.chainI) || a.i[0] - b.i[0]);
  extendBridgeLadders(residues, bridges);
  assignSheetsAndLadders(bridges);
  applyBridgeAssignments(residues, bridges);
  assignStrands(residues);

  return bridges;
}

function testBridge(residues, hbondMap, i, j) {
  const a = residues[i]?.prev;
  const b = residues[i];
  const c = residues[i]?.next;
  const d = residues[j]?.prev;
  const e = residues[j];
  const f = residues[j]?.next;
  if (!a || !b || !c || !d || !e || !f) return "none";
  if (!sameContinuousSegment(residues, a.index, c.index) || !sameContinuousSegment(residues, d.index, f.index)) return "none";

  const parallel =
    (donatesTo(hbondMap, c, e) && donatesTo(hbondMap, e, a)) ||
    (donatesTo(hbondMap, f, b) && donatesTo(hbondMap, b, d));
  if (parallel) return "parallel";

  const antiparallel =
    (donatesTo(hbondMap, c, d) && donatesTo(hbondMap, f, a)) ||
    (donatesTo(hbondMap, e, b) && donatesTo(hbondMap, b, e));
  return antiparallel ? "antiparallel" : "none";
}

function donatesTo(hbondMap, donor, acceptor) {
  return Boolean(donor && acceptor && hbondMap[acceptor.index]?.[donor.index]);
}

function extendBridgeLadders(residues, bridges) {
  for (let i = 0; i < bridges.length; i += 1) {
    for (let j = i + 1; j < bridges.length; j += 1) {
      const a = bridges[i];
      const b = bridges[j];
      const ibi = a.i[0];
      const iei = a.i.at(-1);
      const jbi = a.j[0];
      const jei = a.j.at(-1);
      const ibj = b.i[0];
      const iej = b.i.at(-1);
      const jbj = b.j[0];
      const jej = b.j.at(-1);

      if (
        a.type !== b.type ||
        !sameContinuousSegment(residues, Math.min(ibi, ibj), Math.max(iei, iej)) ||
        !sameContinuousSegment(residues, Math.min(jbi, jbj), Math.max(jei, jej)) ||
        ibj - iei >= 6 ||
        (iei >= ibj && ibi <= iej)
      ) {
        continue;
      }

      const bulge = a.type === "parallel"
        ? ((jbj - jei < 6 && ibj - iei < 3) || jbj - jei < 3)
        : ((jbi - jej < 6 && ibj - iei < 3) || jbi - jej < 3);

      if (bulge) {
        a.i.push(...b.i);
        if (a.type === "parallel") a.j.push(...b.j);
        else a.j.unshift(...b.j);
        bridges.splice(j, 1);
        j -= 1;
      }
    }
  }
}

function assignSheetsAndLadders(bridges) {
  const ladderSet = new Set(bridges);
  let sheet = 1;
  let ladder = 0;

  while (ladderSet.size > 0) {
    const first = ladderSet.values().next().value;
    const sheetSet = new Set([first]);
    ladderSet.delete(first);

    let changed = true;
    while (changed) {
      changed = false;
      for (const a of [...sheetSet]) {
        for (const b of [...ladderSet]) {
          if (bridgesLinked(a, b)) {
            sheetSet.add(b);
            ladderSet.delete(b);
            changed = true;
            break;
          }
        }
        if (changed) break;
      }
    }

    for (const bridge of sheetSet) {
      bridge.ladder = ladder;
      bridge.sheet = sheet;
      bridge.links = [...sheetSet];
      ladder += 1;
    }
    sheet += 1;
  }
}

function bridgesLinked(a, b) {
  return arraysOverlap(a.i, b.i) || arraysOverlap(a.i, b.j) || arraysOverlap(a.j, b.i) || arraysOverlap(a.j, b.j);
}

function arraysOverlap(a, b) {
  const set = new Set(a);
  return b.some((value) => set.has(value));
}

function applyBridgeAssignments(residues, bridges) {
  for (const bridge of bridges) {
    const betaI = bridge.i.some((index) => residues[index].betaPartners[0]) ? 1 : 0;
    const betaJ = bridge.j.some((index) => residues[index].betaPartners[0]) ? 1 : 0;
    const ss = bridge.i.length > 1 ? "E" : "B";

    if (bridge.type === "parallel") {
      for (let k = 0; k < bridge.i.length; k += 1) {
        setBetaPartner(residues[bridge.i[k]], betaI, residues[bridge.j[k]], bridge, true);
        setBetaPartner(residues[bridge.j[k]], betaJ, residues[bridge.i[k]], bridge, true);
      }
    } else {
      const reversedJ = [...bridge.j].reverse();
      const reversedI = [...bridge.i].reverse();
      for (let k = 0; k < bridge.i.length; k += 1) {
        setBetaPartner(residues[bridge.i[k]], betaI, residues[reversedJ[k]], bridge, false);
        setBetaPartner(residues[bridge.j[k]], betaJ, residues[reversedI[k]], bridge, false);
      }
    }

    for (let i = bridge.i[0]; i <= bridge.i.at(-1); i += 1) {
      if (residues[i].ss !== "E") residues[i].ss = ss;
      residues[i].sheet = bridge.sheet;
    }
    for (let i = bridge.j[0]; i <= bridge.j.at(-1); i += 1) {
      if (residues[i].ss !== "E") residues[i].ss = ss;
      residues[i].sheet = bridge.sheet;
    }
  }

  for (const residue of residues) {
    residue.bridgePartners = residue.betaPartners
      .filter(Boolean)
      .map((partner) => partner.residue.dsspIndex);
  }
}

function setBetaPartner(residue, slot, partnerResidue, bridge, parallel) {
  residue.betaPartners[slot] = {
    residue: partnerResidue,
    ladder: bridge.ladder,
    sheet: bridge.sheet,
    parallel
  };
}

function assignStrands(residues) {
  let strand = 0;
  const maxSheet = Math.max(0, ...residues.map((residue) => residue.sheet || 0));
  for (let sheet = 1; sheet <= maxSheet; sheet += 1) {
    let previous = null;
    for (const residue of residues) {
      if (residue.sheet !== sheet) continue;
      if (!previous || residue.index !== previous.index + 1) strand += 1;
      residue.strand = strand;
      previous = residue;
    }
  }
}

function assignHelicesTurnsAndBends(residues, hbondMap, options = {}) {
  for (const [code, stride] of [["G", 3], ["H", 4], ["I", 5]]) {
    for (let i = 0; i + stride < residues.length; i += 1) {
      if (sameContinuousSegment(residues, i, i + stride) && donatesTo(hbondMap, residues[i + stride], residues[i])) {
        setHelixFlag(residues[i + stride], code, "end");
        for (let j = i + 1; j < i + stride; j += 1) {
          if (residues[j].helixFlags[code] === "none") setHelixFlag(residues[j], code, "middle");
        }
        setHelixFlag(residues[i], code, residues[i].helixFlags[code] === "end" ? "startAndEnd" : "start");
      }
    }
  }

  for (const residue of residues) {
    residue.isBend = Number.isFinite(residue.kappa) && residue.kappa > 70;
  }

  for (let i = 1; i + 4 < residues.length; i += 1) {
    if (isHelixStart(residues[i], "H") && isHelixStart(residues[i - 1], "H")) {
      for (let j = i; j <= i + 3; j += 1) residues[j].ss = "H";
    }
  }

  for (let i = 1; i + 3 < residues.length; i += 1) {
    if (isHelixStart(residues[i], "G") && isHelixStart(residues[i - 1], "G")) {
      let empty = true;
      for (let j = i; empty && j <= i + 2; j += 1) {
        empty = residues[j].ss === "-" || residues[j].ss === "G";
      }
      if (empty) {
        for (let j = i; j <= i + 2; j += 1) residues[j].ss = "G";
      }
    }
  }

  const preferPiHelices = options.preferPiHelices ?? true;
  for (let i = 1; i + 5 < residues.length; i += 1) {
    if (isHelixStart(residues[i], "I") && isHelixStart(residues[i - 1], "I")) {
      let empty = true;
      for (let j = i; empty && j <= i + 4; j += 1) {
        empty = residues[j].ss === "-" || residues[j].ss === "I" || (preferPiHelices && residues[j].ss === "H");
      }
      if (empty) {
        for (let j = i; j <= i + 4; j += 1) residues[j].ss = "I";
      }
    }
  }

  for (let i = 1; i + 1 < residues.length; i += 1) {
    if (residues[i].ss !== "-") continue;
    let isTurn = false;
    for (const [code, stride] of [["G", 3], ["H", 4], ["I", 5]]) {
      for (let k = 1; k < stride && !isTurn; k += 1) {
        isTurn = i >= k && isHelixStart(residues[i - k], code);
      }
    }
    if (isTurn) residues[i].ss = "T";
    else if (residues[i].isBend) residues[i].ss = "S";
  }
}

function setHelixFlag(residue, code, flag) {
  residue.helixFlags[code] = flag;
}

function isHelixStart(residue, code) {
  return residue?.helixFlags?.[code] === "start" || residue?.helixFlags?.[code] === "startAndEnd";
}

function assignPPHelices(residues, stretchLength) {
  const epsilon = 29;
  const phiMin = -75 - epsilon;
  const phiMax = -75 + epsilon;
  const psiMin = 145 - epsilon;
  const psiMax = 145 + epsilon;

  for (let i = 1; i + stretchLength < residues.length; i += 1) {
    let matches = true;
    for (let j = 0; j < stretchLength; j += 1) {
      const residue = residues[i + j];
      matches =
        matches &&
        residue.phi >= phiMin &&
        residue.phi <= phiMax &&
        residue.psi >= psiMin &&
        residue.psi <= psiMax;
    }
    if (!matches) continue;

    if (stretchLength === 2) {
      setHelixFlag(residues[i], "P", residues[i].helixFlags.P === "end" ? "middle" : "start");
      setHelixFlag(residues[i + 1], "P", "end");
    } else {
      setHelixFlag(residues[i], "P", residues[i].helixFlags.P === "end" ? "startAndEnd" : "start");
      for (let j = 1; j < stretchLength - 1; j += 1) setHelixFlag(residues[i + j], "P", "middle");
      setHelixFlag(residues[i + stretchLength - 1], "P", "end");
    }

    for (let j = 0; j < stretchLength; j += 1) {
      if (residues[i + j].ss === "-") residues[i + j].ss = "P";
    }
  }
}

function sphereSampleDirections(count) {
  const cached = sphereSampleCache.get(count);
  if (cached) return cached;

  const directions = [];
  const n = Math.floor((count - 1) / 2);
  const points = 2 * n + 1;
  const goldenRatio = (1 + Math.sqrt(5)) / 2;
  for (let i = -n; i <= n; i += 1) {
    const latitude = Math.asin(2 * i / points);
    const longitude = positiveModulo(i, goldenRatio) * 2 * Math.PI / goldenRatio;
    directions.push({
      x: Math.sin(longitude) * Math.cos(latitude),
      y: Math.cos(longitude) * Math.cos(latitude),
      z: Math.sin(latitude)
    });
  }

  sphereSampleCache.set(count, directions);
  return directions;
}

function buildAtomGrid(atoms, cellSize) {
  const grid = new Map();
  for (const entry of atoms) {
    const key = cellKeyForPoint(entry.atom, cellSize);
    const bucket = grid.get(key);
    if (bucket) bucket.push(entry);
    else grid.set(key, [entry]);
  }
  return grid;
}

function atomNeighbors(entry, grid, cellSize) {
  const neighbors = [];
  const cx = Math.floor(entry.atom.x / cellSize);
  const cy = Math.floor(entry.atom.y / cellSize);
  const cz = Math.floor(entry.atom.z / cellSize);

  for (let dx = -1; dx <= 1; dx += 1) {
    for (let dy = -1; dy <= 1; dy += 1) {
      for (let dz = -1; dz <= 1; dz += 1) {
        const bucket = grid.get(cellKey(cx + dx, cy + dy, cz + dz));
        if (!bucket) continue;
        for (const candidate of bucket) {
          if (candidate === entry) continue;
          const maxContact = entry.radius + candidate.radius;
          if (distance(entry.atom, candidate.atom) < maxContact) {
            neighbors.push(candidate);
          }
        }
      }
    }
  }

  return neighbors;
}

function isBuriedPoint(point, neighbors) {
  for (const neighbor of neighbors) {
    if (distance(point, neighbor.atom) < neighbor.radius) return true;
  }
  return false;
}

function cellKeyForPoint(point, cellSize) {
  return cellKey(
    Math.floor(point.x / cellSize),
    Math.floor(point.y / cellSize),
    Math.floor(point.z / cellSize)
  );
}

function cellKey(x, y, z) {
  return `${x},${y},${z}`;
}

function normalizedElement(atom) {
  const element = String(atom.element || "").trim().toUpperCase();
  if (element) return element;
  const stripped = String(atom.name || "").replace(/[0-9]/g, "").trim().toUpperCase();
  if (!stripped) return "";
  if (stripped === "CL" || stripped === "BR") return stripped;
  return stripped[0];
}

function dsspAtomRadius(atom) {
  const name = normalizedAtomName(atom.name);
  if (name === "N") return DSSP_ATOM_RADII.N;
  if (name === "CA") return DSSP_ATOM_RADII.CA;
  if (name === "C") return DSSP_ATOM_RADII.C;
  if (name === "O") return DSSP_ATOM_RADII.O;
  return DSSP_ATOM_RADII.SIDE;
}

function positiveModulo(value, divisor) {
  return ((value % divisor) + divisor) % divisor;
}

function summarize(residues, hbonds) {
  const counts = Object.fromEntries(Object.keys(STRUCTURE_NAMES).map((code) => [code, 0]));
  let asaTotal = 0;
  let rsaTotal = 0;
  let rsaCount = 0;
  for (const residue of residues) {
    counts[residue.ss] = (counts[residue.ss] ?? 0) + 1;
    if (Number.isFinite(residue.asa)) asaTotal += residue.asa;
    if (Number.isFinite(residue.rsa)) {
      rsaTotal += residue.rsa;
      rsaCount += 1;
    }
  }
  return {
    residueCount: residues.length,
    chainCount: new Set(residues.map((residue) => residue.chainId)).size,
    hbondCount: hbonds.length,
    asaTotal,
    meanRsa: rsaCount ? rsaTotal / rsaCount : null,
    counts
  };
}

function sameContinuousSegment(residues, start, end) {
  if (start < 0 || end >= residues.length || start > end) return false;
  for (let i = start + 1; i <= end; i += 1) {
    if (residues[i].prev !== residues[i - 1]) return false;
  }
  return true;
}

function distance(a, b) {
  return Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
}

function distanceSq(a, b) {
  return (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2;
}

function sub(a, b) {
  return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
}

function add(a, b) {
  return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z };
}

function scale(a, factor) {
  return { x: a.x * factor, y: a.y * factor, z: a.z * factor };
}

function dot(a, b) {
  if (!a || !b) return Number.NaN;
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

function cross(a, b) {
  return {
    x: a.y * b.z - a.z * b.y,
    y: a.z * b.x - a.x * b.z,
    z: a.x * b.y - a.y * b.x
  };
}

function norm(a) {
  return Math.hypot(a.x, a.y, a.z);
}

function normalize(a) {
  const length = norm(a);
  if (!Number.isFinite(length) || length < MIN_DISTANCE) return null;
  return scale(a, 1 / length);
}

function angleFromVectors(a, b) {
  const aUnit = normalize(a);
  const bUnit = normalize(b);
  if (!aUnit || !bUnit) return null;
  const cosine = clamp(dot(aUnit, bUnit), -1, 1);
  return radiansToDegrees(Math.acos(cosine));
}

function dihedral(a, b, c, d) {
  const v12 = sub(a, b);
  const v43 = sub(d, c);
  const z = sub(b, c);
  const p = cross(z, v12);
  const x = cross(z, v43);
  const y = cross(z, x);
  const u = dot(x, x);
  const v = dot(y, y);
  if (u <= 0 || v <= 0) return 360;

  const xComponent = dot(p, x) / Math.sqrt(u);
  const yComponent = dot(p, y) / Math.sqrt(v);
  if (xComponent === 0 && yComponent === 0) return 360;
  return radiansToDegrees(Math.atan2(yComponent, xComponent));
}

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

function radiansToDegrees(radians) {
  return radians * 180 / Math.PI;
}

function dsspSummary(summary = {}) {
  return {
    residueCount: summary.residueCount ?? 0,
    chainCount: summary.chainCount ?? 0,
    hbondCount: summary.hbondCount ?? 0,
    asaTotal: roundNumber(summary.asaTotal, 1),
    meanRsa: roundNumber(summary.meanRsa, 4),
    counts: { ...(summary.counts ?? {}) }
  };
}

function roundNumber(value, digits) {
  if (!Number.isFinite(value)) return null;
  const factor = 10 ** digits;
  return Math.round(value * factor) / factor;
}

function dsspAngleValue(value) {
  return Number.isFinite(value) ? roundNumber(value, 1) : 360;
}

function formatBondOffset(bond, residue) {
  if (!bond) return 0;
  return bond.donor === residue.index ? bond.offsetFromDonor : bond.offsetFromAcceptor;
}

function formatBondEnergy(bond) {
  return bond ? roundNumber(bond.energy, 1) : 0;
}

function csvCell(value) {
  const text = String(value ?? "");
  if (/[",\n]/.test(text)) {
    return `"${text.replaceAll('"', '""')}"`;
  }
  return text;
}

export const constants = {
  HBOND_Q1Q2,
  HBOND_F,
  DEFAULT_HBOND_CUTOFF,
  PEPTIDE_BOND_MAX,
  MIN_HBOND_ENERGY,
  MAX_CA_DISTANCE,
  WATER_PROBE_RADIUS,
  DEFAULT_ASA_SAMPLES,
  DSSP_ATOM_RADII,
  MAX_ASA
};
