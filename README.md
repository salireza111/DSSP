## Simple DSSP
#### DSSP-compatible JavaScript implementation
A static, browser-side JavaScript implementation of DSSP-compatible secondary-structure and accessibility assignment for legacy PDB files.

Open `index.html` through a static server, paste or drop a PDB file, or load the bundled `1YIO.pdb` demo, and the app calculates:

- backbone hydrogen-bond energies with the Kabsch-Sander electrostatic formula
- inferred amide hydrogens when PDB hydrogens are missing
- H/G/I/P helix, E/B beta bridge or strand, T turn, S bend, and loop labels
- phi, psi, kappa, alpha, TCO, bridge partners, and best H-bonds
- per-residue solvent-accessible surface area (ASA) and relative solvent accessibility (RSA) with DSSP-style radii, a 1.4 Å water probe, and 401 Fibonacci surface points
- CSV and JSON exports
- direct fetch by 4-character PDB ID from the RCSB PDB legacy PDB download endpoint, with optional chain filtering
- 3Dmol.js residue popover when hovering, clicking, or focusing a result row/sequence tile, with remembered ribbon, phi/psi, and DSSP H-bond overlay toggles plus click-to-pin behavior

![ScreenShot](./Simple%20DSSP.png)

## Run

```bash
npm run serve
```

Then open `http://127.0.0.1:4173`.

The app is plain HTML/CSS/JS. Any static web host can serve it.

## Validate Against mkdssp

If `mkdssp` is installed, compare the JavaScript output against the official binary:

```bash
npm run validate:mkdssp
node scripts/validate-mkdssp.mjs path/to/file.pdb
```

Set `MKDSSP_BIN=/path/to/mkdssp` if the executable is not on `PATH`.

## Underhood DSSP Tools

For benchmarking and machine-readable exports, use the Python utility:

```bash
python3 scripts/dssp_tools.py benchmark assets/1YIO.pdb
```

It runs `src/dssp.js`, runs the installed `mkdssp`/`dssp`, matches residues by chain and residue number, then reports secondary-structure agreement, H-bond count delta, ASA total delta, and mean absolute differences for ASA, phi, psi, kappa, alpha, and TCO.

The main-page table, JS CSV/JSON exports, and real-DSSP CSV/JSON exports use the same DSSP-style residue field order:
`dssp_index, chain, residue_number, insertion_code, aa, structure, structure_name, bridge_1, bridge_2, asa, rsa, phi, psi, kappa, alpha, tco, x_ca, y_ca, z_ca`, followed by split H-bond offset and energy columns.

Useful commands:

```bash
# Benchmark JS against real DSSP, optionally for one chain.
python3 scripts/dssp_tools.py benchmark path/to/file.pdb --chain A

# Save the benchmark report as JSON.
python3 scripts/dssp_tools.py benchmark path/to/file.pdb --json -o benchmark.json

# Export this project's JS DSSP output.
python3 scripts/dssp_tools.py js path/to/file.pdb --format json -o js-dssp.json
python3 scripts/dssp_tools.py js path/to/file.pdb --format csv -o js-dssp.csv

# Export real local DSSP output parsed from mkdssp/dssp.
python3 scripts/dssp_tools.py real path/to/file.pdb --format json -o real-dssp.json
python3 scripts/dssp_tools.py real path/to/file.pdb --format csv -o real-dssp.csv
```

The tool looks for `DSSP_BIN`, `MKDSSP_BIN`, `mkdssp`, then `dssp`. You can also pass the binary directly:

```bash
python3 scripts/dssp_tools.py benchmark path/to/file.pdb --dssp-bin /path/to/mkdssp
```

For slower but denser JS ASA sampling:

```bash
python3 scripts/dssp_tools.py benchmark path/to/file.pdb --asa-samples 960
```

The same utility is available through npm:

```bash
npm run benchmark:dssp
npm run dssp -- js assets/1YIO.pdb --format csv -o js-dssp.csv
```

## Use As A Module

```js
import { assignDSSP } from "./src/dssp.js";

const result = assignDSSP(pdbText);
console.log(result.sequence);
console.log(result.secondary);
console.log(result.residues);
```

## Notes

This is a JavaScript port of the core DSSP assignment rules for static browser use. It now mirrors key `mkdssp` behavior more closely: DSSP carbonyl-vector amide hydrogen placement, top-two H-bond columns, 9 Å C-alpha neighbor filtering, beta bridge/ladder/bulge/sheet grouping, helix priority including pi-helix preference, PPII `P`, and DSSP-style ASA sampling.

Important limitations:

- Exact parity is not guaranteed for every structure; use `npm run validate:mkdssp` to compare against the same `mkdssp` version you want to target.
- It currently reads legacy PDB, not mmCIF, so DSSP 4 mmCIF categories and some metadata-derived features are not reproduced.
- Disulfide annotations, complete DSSP/mmCIF output categories, and automated `mkdssp` comparison tests are not implemented yet.
- RSA depends on the maximum-ASA reference table used here.


Useful references:

- [PDB-REDO DSSP documentation](https://pdb-redo.eu/dssp/about)
- [DSSP 4 source repository](https://github.com/PDB-REDO/dssp)
- [Kabsch and Sander DSSP paper DOI](https://doi.org/10.1002/bip.360221211)
