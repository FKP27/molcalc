#!/bin/bash

cd molcalc/static/

mkdir -p rdkit
cd rdkit

RDKIT_DIST_BASE="https://unpkg.com/@rdkit/rdkit/dist"

wget "$RDKIT_DIST_BASE/RDKit_minimal.js" -O rdkit.js
wget "$RDKIT_DIST_BASE/RDKit_minimal.wasm" -O RDKit_minimal.wasm

if ! grep -q "Auto-init RDKit" rdkit.js; then
  cat <<'EOF' >> rdkit.js

// Auto-init RDKit and expose a global for existing code.
if (typeof initRDKitModule === "function") {
  initRDKitModule().then(function(RDKit) {
    window.RDKit = RDKit;
  });
}
EOF
fi
