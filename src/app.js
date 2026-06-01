import {
  DSSP_COLUMNS,
  assignDSSP,
  filterPDBByChain,
  formatDsspCell,
  parsePDB,
  structureName,
  toCSV,
  toDsspJSON,
  toDsspRows
} from "./dssp.js";

const fileInput = document.querySelector("#pdb-file");
const pdbText = document.querySelector("#pdb-text");
const runButton = document.querySelector("#run-analysis");
const clearButton = document.querySelector("#clear-input");
const demoButton = document.querySelector("#load-demo");
const fetchPdbForm = document.querySelector("#pdb-fetch-form");
const pdbCodeInput = document.querySelector("#pdb-code");
const pdbChainInput = document.querySelector("#pdb-chain");
const fetchPdbButton = document.querySelector("#fetch-pdb");
const cutoffInput = document.querySelector("#hbond-cutoff");
const statusEl = document.querySelector("#status");
const summaryEl = document.querySelector("#summary");
const sequenceEl = document.querySelector("#sequence");
const tableBody = document.querySelector("#residue-table tbody");
const exportCsvButton = document.querySelector("#export-csv");
const exportJsonButton = document.querySelector("#export-json");
const dropZone = document.querySelector("#drop-zone");
const residueViewer = document.querySelector("#residue-viewer");
const viewerSurface = document.querySelector("#viewer-mol");
const viewerTitle = document.querySelector("#viewer-title");
const viewerStructure = document.querySelector("#viewer-structure");
const viewerCaption = document.querySelector("#viewer-caption");
const viewerRibbonToggle = document.querySelector("#viewer-ribbon-toggle");
const viewerAngleToggle = document.querySelector("#viewer-angle-toggle");
const viewerHbondToggle = document.querySelector("#viewer-hbond-toggle");

let latestResult = null;
let viewerResidueIndex = null;
let viewerHideTimer = null;
let molViewer = null;
let viewerPinned = false;
let viewerRibbonMode = false;
let viewerAngleMode = true;
let viewerHbondMode = false;

const DEMO_PDB_URL = "./assets/1YIO.pdb";
const RCSB_PDB_URL = "https://files.rcsb.org/download/";

fileInput.addEventListener("change", async () => {
  const file = fileInput.files?.[0];
  if (!file) return;
  pdbText.value = await file.text();
  setStatus(`Loaded ${file.name}. Ready to calculate.`, "ok");
});

dropZone.addEventListener("dragover", (event) => {
  event.preventDefault();
  dropZone.classList.add("is-dragging");
});

dropZone.addEventListener("dragleave", () => {
  dropZone.classList.remove("is-dragging");
});

dropZone.addEventListener("drop", async (event) => {
  event.preventDefault();
  dropZone.classList.remove("is-dragging");
  const file = event.dataTransfer.files?.[0];
  if (!file) return;
  pdbText.value = await file.text();
  setStatus(`Loaded ${file.name}. Ready to calculate.`, "ok");
});

runButton.addEventListener("click", runAnalysis);

tableBody.addEventListener("pointerover", handleResiduePointerOver);
tableBody.addEventListener("pointermove", handleResiduePointerMove);
tableBody.addEventListener("pointerout", scheduleHideViewer);
tableBody.addEventListener("focusin", handleResidueFocusIn);
tableBody.addEventListener("focusout", scheduleHideViewer);
tableBody.addEventListener("click", handleResidueClick);
tableBody.addEventListener("keydown", handleResidueKeydown);
sequenceEl.addEventListener("pointerover", handleResiduePointerOver);
sequenceEl.addEventListener("pointermove", handleResiduePointerMove);
sequenceEl.addEventListener("pointerout", scheduleHideViewer);
sequenceEl.addEventListener("focusin", handleResidueFocusIn);
sequenceEl.addEventListener("focusout", scheduleHideViewer);
sequenceEl.addEventListener("click", handleResidueClick);
sequenceEl.addEventListener("keydown", handleResidueKeydown);
residueViewer.addEventListener("pointerenter", cancelViewerHide);
residueViewer.addEventListener("pointerleave", scheduleHideViewer);
viewerRibbonToggle.addEventListener("change", () => {
  viewerRibbonMode = viewerRibbonToggle.checked;
  renderResidueViewer();
});
viewerAngleToggle.addEventListener("click", () => {
  viewerAngleMode = !viewerAngleMode;
  syncViewerButtons();
  renderResidueViewer();
});
viewerHbondToggle.addEventListener("click", () => {
  viewerHbondMode = !viewerHbondMode;
  syncViewerButtons();
  refreshViewerCaption();
  renderResidueViewer();
});
document.addEventListener("pointerdown", handleDocumentPointerDown);
window.addEventListener("resize", () => {
  if (molViewer) {
    molViewer.resize();
    molViewer.render();
  }
});

clearButton.addEventListener("click", () => {
  pdbText.value = "";
  fileInput.value = "";
  latestResult = null;
  renderEmpty();
  setStatus("Paste a PDB file or drop one into the page.", "idle");
});

demoButton.addEventListener("click", loadDemoPdb);
fetchPdbForm.addEventListener("submit", fetchPdbById);

exportCsvButton.addEventListener("click", () => {
  if (!latestResult) return;
  downloadText("dssp-results.csv", toCSV(latestResult), "text/csv");
});

exportJsonButton.addEventListener("click", () => {
  if (!latestResult) return;
  downloadText("dssp-results.json", JSON.stringify(toDsspJSON(latestResult), null, 2), "application/json");
});

function runAnalysis() {
  const text = pdbText.value.trim();
  if (!text) {
    setStatus("No PDB text yet. Add a structure first.", "error");
    return;
  }

  try {
    const residues = parsePDB(text);
    if (residues.length === 0) {
      throw new Error("No ATOM records were found in this file.");
    }
    const cutoff = Number.parseFloat(cutoffInput.value);
    latestResult = assignDSSP(residues, {
      hbondCutoff: Number.isFinite(cutoff) ? cutoff : -0.5
    });
    renderResult(latestResult);
    setStatus(resultStatus("Calculated", latestResult), "ok");
  } catch (error) {
    latestResult = null;
    renderEmpty();
    setStatus(error.message, "error");
  }
}

async function loadDemoPdb() {
  setStatus("Loading 1YIO demo...", "idle");
  try {
    const response = await fetch(DEMO_PDB_URL);
    if (!response.ok) {
      throw new Error("Demo PDB could not be loaded.");
    }
    pdbText.value = await response.text();
    runAnalysis();
  } catch (error) {
    const message = error instanceof Error ? error.message : "Demo PDB could not be loaded.";
    setStatus(`${message} Use a local static server or choose the file manually.`, "error");
  }
}

async function fetchPdbById(event) {
  event.preventDefault();
  const pdbId = pdbCodeInput.value.trim().toUpperCase();
  const chainId = pdbChainInput.value.trim().toUpperCase();
  if (!/^[A-Z0-9]{4}$/.test(pdbId)) {
    setStatus("Enter a 4-character PDB ID, for example 1YIO.", "error");
    pdbCodeInput.focus();
    return;
  }
  if (chainId && !/^[A-Z0-9_]$/.test(chainId)) {
    setStatus("Enter one chain ID character, for example A, or leave it blank.", "error");
    pdbChainInput.focus();
    return;
  }

  pdbCodeInput.value = pdbId;
  pdbChainInput.value = chainId;
  setFetchBusy(true);
  const targetLabel = chainId ? `${pdbId} chain ${chainId}` : pdbId;
  setStatus(`Fetching ${targetLabel} from RCSB PDB...`, "idle");
  try {
    const response = await fetch(`${RCSB_PDB_URL}${encodeURIComponent(pdbId)}.pdb`);
    if (!response.ok) {
      const reason = response.status === 404
        ? `${pdbId} was not found as a legacy PDB file.`
        : `RCSB PDB returned HTTP ${response.status} for ${pdbId}.`;
      throw new Error(reason);
    }

    const pdb = await response.text();
    pdbText.value = chainId ? filterPDBByChain(pdb, chainId) : pdb;
    fileInput.value = "";
    runAnalysis();
    if (latestResult) {
      setStatus(resultStatus(`Fetched ${targetLabel}. Calculated`, latestResult), "ok");
    }
  } catch (error) {
    const message = error instanceof Error ? error.message : "Could not fetch that PDB file.";
    setStatus(`${message} Check the ID, internet access, or paste the PDB text manually.`, "error");
  } finally {
    setFetchBusy(false);
  }
}

function setFetchBusy(isBusy) {
  fetchPdbButton.disabled = isBusy;
  pdbCodeInput.disabled = isBusy;
  pdbChainInput.disabled = isBusy;
  fetchPdbButton.textContent = isBusy ? "Fetching..." : "Fetch PDB";
}

function renderResult(result) {
  exportCsvButton.disabled = false;
  exportJsonButton.disabled = false;
  renderSummary(result);
  renderSequence(result);
  renderRows(result);
}

function renderSummary(result) {
  const cards = [
    ["Residues", result.summary.residueCount],
    ["Chains", result.summary.chainCount],
    ["H-bonds", result.summary.hbondCount],
    ["Helix", countCodes(result.summary.counts, ["H", "G", "I", "P"])],
    ["Sheet", countCodes(result.summary.counts, ["E", "B"])],
    ["Loop", result.summary.counts["-"] ?? 0],
    ["Mean RSA", formatPercent(result.summary.meanRsa)]
  ];

  summaryEl.replaceChildren(
    ...cards.map(([label, value]) => {
      const card = document.createElement("div");
      card.className = "metric";
      const number = document.createElement("strong");
      number.textContent = value;
      const caption = document.createElement("span");
      caption.textContent = label;
      card.append(number, caption);
      return card;
    })
  );
}

function renderSequence(result) {
  const fragment = document.createDocumentFragment();
  for (const residue of result.residues) {
    const span = document.createElement("span");
    span.className = `residue-code ss-${residue.ss === "-" ? "loop" : residue.ss.toLowerCase()}`;
    span.dataset.residueIndex = String(residue.index);
    span.tabIndex = 0;
    span.setAttribute("role", "button");
    span.setAttribute("aria-label", `${residue.aa} ${residue.chainId}:${residue.resSeq}${residue.iCode || ""} ${structureName(residue.ss)}`);
    span.title = `${residue.dsspIndex} ${residue.chainId}:${residue.resSeq}${residue.iCode || ""} ${residue.resName} - ${structureName(residue.ss)}`;
    span.textContent = residue.aa;
    fragment.append(span);
  }
  sequenceEl.replaceChildren(fragment);
}

function renderRows(result) {
  const dsspRows = toDsspRows(result);
  const rows = result.residues.map((residue, index) => {
    const dsspRow = dsspRows[index];
    const tr = document.createElement("tr");
    tr.dataset.residueIndex = String(residue.index);
    tr.tabIndex = 0;
    const cells = DSSP_COLUMNS.map((column) => formatDsspCell(column, dsspRow[column]));
    for (const cell of cells) {
      const td = document.createElement("td");
      td.textContent = cell;
      tr.append(td);
    }
    return tr;
  });
  tableBody.replaceChildren(...rows);
}

function renderEmpty() {
  hideResidueViewer();
  viewerPinned = false;
  exportCsvButton.disabled = true;
  exportJsonButton.disabled = true;
  summaryEl.replaceChildren();
  sequenceEl.replaceChildren();
  tableBody.replaceChildren();
}

function countCodes(counts, codes) {
  return codes.reduce((total, code) => total + (counts[code] ?? 0), 0);
}

function resultStatus(prefix, result) {
  return `${prefix} ${result.residues.length} residues, ${result.hbonds.length} backbone H-bonds, and ${formatSummaryAsa(result.summary.asaTotal)} A^2 ASA.`;
}

function formatAngle(value) {
  return Number.isFinite(value) ? value.toFixed(1) : "";
}

function formatSummaryAsa(value) {
  return Number.isFinite(value) ? value.toFixed(0) : "";
}

function formatPercent(value) {
  return Number.isFinite(value) ? `${(value * 100).toFixed(0)}%` : "";
}

function handleResiduePointerOver(event) {
  if (viewerPinned) return;
  const target = residueTarget(event.target);
  if (!target) return;
  cancelViewerHide();
  showResidueViewer(Number(target.dataset.residueIndex), event, "hover");
}

function handleResiduePointerMove(event) {
  if (residueViewer.hidden || viewerPinned) return;
  positionResidueViewer(event, "hover");
}

function handleResidueFocusIn(event) {
  if (viewerPinned) return;
  const target = residueTarget(event.target);
  if (!target) return;
  const rect = target.getBoundingClientRect();
  showResidueViewer(Number(target.dataset.residueIndex), {
    clientX: rect.right,
    clientY: rect.top + rect.height / 2,
    target
  }, "focus");
}

function handleResidueClick(event) {
  const target = residueTarget(event.target);
  if (!target) return;
  event.stopPropagation();
  cancelViewerHide();
  viewerPinned = true;
  showResidueViewer(Number(target.dataset.residueIndex), event, "pinned");
}

function handleResidueKeydown(event) {
  if (event.key !== "Enter" && event.key !== " ") return;
  const target = residueTarget(event.target);
  if (!target) return;
  event.preventDefault();
  viewerPinned = true;
  const rect = target.getBoundingClientRect();
  showResidueViewer(Number(target.dataset.residueIndex), {
    clientX: rect.right,
    clientY: rect.top + rect.height / 2,
    target
  }, "pinned");
}

function residueTarget(target) {
  const element = target instanceof Element ? target.closest("[data-residue-index]") : null;
  if (!element || (!tableBody.contains(element) && !sequenceEl.contains(element))) return null;
  return element;
}

function cancelViewerHide() {
  window.clearTimeout(viewerHideTimer);
}

function scheduleHideViewer(event) {
  if (viewerPinned) return;
  if (event?.relatedTarget instanceof Node && residueViewer.contains(event.relatedTarget)) return;
  cancelViewerHide();
  viewerHideTimer = window.setTimeout(hideResidueViewer, 140);
}

function hideResidueViewer() {
  viewerResidueIndex = null;
  residueViewer.hidden = true;
  viewerPinned = false;
  markViewerFocus(null);
}

function showResidueViewer(index, point, mode = "hover") {
  if (!latestResult?.residues?.[index]) return;
  viewerResidueIndex = index;
  const residue = latestResult.residues[index];
  viewerTitle.textContent = `${residue.aa} ${residue.chainId}:${residue.resSeq}${residue.iCode || ""} (${residue.resName})`;
  viewerStructure.textContent = structureName(residue.ss);
  viewerStructure.className = `viewer-structure ss-${residue.ss === "-" ? "loop" : residue.ss.toLowerCase()}`;
  viewerRibbonToggle.checked = viewerRibbonMode;
  syncViewerButtons();
  refreshViewerCaption();
  residueViewer.hidden = false;
  positionResidueViewer(point, mode);
  markViewerFocus(index);
  renderResidueViewer();
}

function positionResidueViewer(point, mode = "hover") {
  const margin = 12;
  const rect = residueViewer.getBoundingClientRect();
  const target = point.target instanceof Element ? point.target : residueTarget(document.elementFromPoint(point.clientX, point.clientY));
  const targetRect = target?.getBoundingClientRect();
  const tableRect = document.querySelector(".table-wrap")?.getBoundingClientRect();
  const preferRightOfTable = mode === "hover" && target && tableBody.contains(target) && tableRect;
  let left = preferRightOfTable ? tableRect.right + margin : point.clientX + 16;
  let top = preferRightOfTable && targetRect ? targetRect.top : point.clientY + 16;

  if (left + rect.width + margin > window.innerWidth) {
    left = preferRightOfTable ? window.innerWidth - rect.width - margin : point.clientX - rect.width - 16;
  }
  if (top + rect.height + margin > window.innerHeight) {
    top = window.innerHeight - rect.height - margin;
  }
  residueViewer.style.left = `${Math.max(margin, left)}px`;
  residueViewer.style.top = `${Math.max(margin, top)}px`;
}

function handleDocumentPointerDown(event) {
  if (!viewerPinned) return;
  if (residueViewer.contains(event.target) || residueTarget(event.target)) return;
  hideResidueViewer();
}

function syncViewerButtons() {
  viewerAngleToggle.setAttribute("aria-pressed", String(viewerAngleMode));
  viewerHbondToggle.setAttribute("aria-pressed", String(viewerHbondMode));
}

function refreshViewerCaption() {
  const residue = latestResult?.residues?.[viewerResidueIndex];
  if (!residue) return;
  viewerCaption.replaceChildren(...viewerCaptionParts(residue));
}

function markViewerFocus(index) {
  document.querySelectorAll(".is-viewer-focus").forEach((element) => {
    element.classList.remove("is-viewer-focus");
  });
  if (index == null) return;
  document.querySelectorAll(`[data-residue-index="${index}"]`).forEach((element) => {
    element.classList.add("is-viewer-focus");
  });
}

function bestBondText(residue) {
  const donor = residue.donorHBonds?.[0];
  const acceptor = residue.acceptorHBonds?.[0];
  if (donor) return `best N-H to O ${donor.offsetFromDonor}:${donor.energy.toFixed(2)} kcal/mol`;
  if (acceptor) {
    const prefix = acceptor.offsetFromAcceptor > 0 ? "+" : "";
    return `best O to H-N ${prefix}${acceptor.offsetFromAcceptor}:${acceptor.energy.toFixed(2)} kcal/mol`;
  }
  return "no backbone H-bond at current cutoff";
}

function viewerCaptionParts(residue) {
  const parts = [
    captionPill(`Residue ${residue.dsspIndex}`),
    captionPill(`Phi ${formatAngle(residue.phi)}`),
    captionPill(`Psi ${formatAngle(residue.psi)}`),
    captionPill(`DSSP ${residue.ss}`),
    captionPill(`RSA ${formatPercent(residue.rsa)}`)
  ];
  if (viewerHbondMode) {
    parts.push(captionPill(`H-bonds ${hbondGuidesForResidue(residue).length}`, "wide"));
  } else {
    parts.push(captionPill(bestBondText(residue), "wide"));
  }
  return parts;
}

function captionPill(text, variant = "") {
  const pill = document.createElement("span");
  pill.className = variant ? `caption-pill ${variant}` : "caption-pill";
  pill.textContent = text;
  return pill;
}

function renderResidueViewer() {
  if (!latestResult || viewerResidueIndex == null) return;
  if (!window.$3Dmol) {
    renderViewerMessage("3Dmol.js did not load.");
    return;
  }

  const residue = latestResult.residues[viewerResidueIndex];
  viewerSurface.replaceChildren();
  molViewer = window.$3Dmol.createViewer(viewerSurface, { backgroundColor: "white" });
  molViewer.addModel(pdbText.value, "pdb");
  molViewer.setStyle({}, { cartoon: { color: "spectrum", opacity: 0.55 } });

  if (viewerRibbonMode) {
    styleResidueRibbon(molViewer, residue, "#16a34a");
  } else {
    styleResidue(molViewer, residue, "#dc2626", 0.32);
    if (residue.prev) {
      styleResidue(molViewer, residue.prev, "#2563eb", 0.22);
    }
    if (residue.next) {
      styleResidue(molViewer, residue.next, "#15803d", 0.22);
    }
  }

  ["N", "CA", "C"].forEach((atomName) => {
    const atom = residue.atoms[atomName];
    if (!atom) return;
    molViewer.addSphere({ center: atom, radius: 0.28, color: "#facc15" });
    molViewer.addLabel(atomName, {
      position: atom,
      backgroundColor: "white",
      backgroundOpacity: 0.86,
      fontColor: "black",
      fontSize: 11
    });
  });
  if (viewerAngleMode && residue.prev?.atoms?.C && residue.atoms.C) {
    addAngleGuide(molViewer, residue.prev.atoms.C, residue.atoms.C, "#16a34a", `Phi ${formatAngle(residue.phi)}`);
  }
  if (viewerAngleMode && residue.next?.atoms?.N && residue.atoms.N) {
    addAngleGuide(molViewer, residue.atoms.N, residue.next.atoms.N, "#f97316", `Psi ${formatAngle(residue.psi)}`);
  }
  if (viewerHbondMode) {
    addDsspHbondGuides(molViewer, residue);
  }

  molViewer.zoomTo(selectionForResidue(residue));
  molViewer.zoom(1.15);
  molViewer.render();
  molViewer.resize();
}

function renderViewerMessage(message) {
  molViewer = null;
  viewerSurface.replaceChildren();
  const empty = document.createElement("div");
  empty.className = "viewer-empty";
  empty.textContent = message;
  viewerSurface.append(empty);
}

function styleResidue(viewer, residue, color, radius) {
  viewer.setStyle(selectionForResidue(residue), { stick: { color, radius } });
}

function styleResidueRibbon(viewer, residue, color) {
  viewer.setStyle(selectionForResidue(residue), {
    cartoon: { color, opacity: 1 }
  });
}

function addAngleGuide(viewer, start, end, color, label) {
  viewer.addCylinder({
    start,
    end,
    radius: 0.1,
    color,
    fromCap: 1,
    toCap: 1
  });
  viewer.addLabel(label, {
    position: midpoint(start, end),
    backgroundColor: "white",
    backgroundOpacity: 0.86,
    fontColor: "black",
    fontSize: 11
  });
}

function addDsspHbondGuides(viewer, residue) {
  for (const guide of hbondGuidesForResidue(residue)) {
    addDashedGuide(viewer, guide.start, guide.end, "#facc15");
    viewer.addLabel(guide.label, {
      position: midpoint(guide.start, guide.end),
      backgroundColor: "white",
      backgroundOpacity: 0.88,
      fontColor: "black",
      fontSize: 10
    });
  }
}

function hbondGuidesForResidue(residue) {
  const guides = [];
  for (const bond of residue.donorHBonds ?? []) {
    const acceptor = latestResult?.residues?.[bond.acceptor];
    const start = atomPoint(residue.atoms.H) ?? atomPoint(residue.atoms.N);
    const end = atomPoint(acceptor?.atoms?.O);
    if (start && end) {
      guides.push({
        start,
        end,
        label: `N-H->O ${bond.offsetFromDonor}:${bond.energy.toFixed(2)}`
      });
    }
  }
  for (const bond of residue.acceptorHBonds ?? []) {
    const donor = latestResult?.residues?.[bond.donor];
    const start = atomPoint(donor?.atoms?.H) ?? atomPoint(donor?.atoms?.N);
    const end = atomPoint(residue.atoms.O);
    const prefix = bond.offsetFromAcceptor > 0 ? "+" : "";
    if (start && end) {
      guides.push({
        start,
        end,
        label: `O->H-N ${prefix}${bond.offsetFromAcceptor}:${bond.energy.toFixed(2)}`
      });
    }
  }
  return guides.slice(0, 4);
}

function addDashedGuide(viewer, start, end, color) {
  const segments = 9;
  for (let i = 0; i < segments; i += 2) {
    viewer.addCylinder({
      start: interpolatePoint(start, end, i / segments),
      end: interpolatePoint(start, end, (i + 1) / segments),
      radius: 0.055,
      color,
      fromCap: 1,
      toCap: 1
    });
  }
}

function interpolatePoint(start, end, fraction) {
  return {
    x: start.x + (end.x - start.x) * fraction,
    y: start.y + (end.y - start.y) * fraction,
    z: start.z + (end.z - start.z) * fraction
  };
}

function selectionForResidue(residue) {
  const residueNumber = Number.parseInt(residue.resSeq, 10);
  const selection = Number.isFinite(residueNumber) ? { resi: residueNumber } : { resi: residue.resSeq };
  if (residue.chainId && residue.chainId !== "_") {
    selection.chain = residue.chainId;
  }
  return selection;
}

function atomPoint(atom) {
  if (!atom) return null;
  return { x: atom.x, y: atom.y, z: atom.z };
}

function midpoint(a, b) {
  return {
    x: (a.x + b.x) / 2,
    y: (a.y + b.y) / 2,
    z: (a.z + b.z) / 2
  };
}

function setStatus(message, type) {
  statusEl.textContent = message;
  statusEl.dataset.type = type;
}

function downloadText(filename, text, mimeType) {
  const blob = new Blob([text], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.append(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

renderEmpty();
