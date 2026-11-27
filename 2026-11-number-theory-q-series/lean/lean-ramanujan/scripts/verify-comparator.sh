#!/usr/bin/env bash
# Runs leanprover/comparator on comparator.json (RamanujanStatement vs RamanujanSolution, with the
# NanoDa kernel replay), using the tool revisions that zeta23 pins for Lean v4.33.0-rc2.  Modelled
# on zeta23/scripts/verify-comparator.sh (Apache-2.0, Anthropic PBC).  Comparator is pointed at
# the Landrun binary directly: at this pin it already passes the `--` delimiter itself, and
# zeta23's landrun-wrapper.sh rejects that.
# Requires Linux with Landlock, git, go, cargo, lake.  Run from anywhere, after `lake build`.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cache=${COMPARATOR_CACHE:-"$root/.lake/comparator-tools"}
mkdir -p "$cache/bin"
comparator_commit=75c730e9b6ef5c2c3b334fad7c3d51fe20624c88
lean4export_commit=9fb131bb100eb32ccf6836f14e4f8328d13b6792
landrun_commit=811cfff51ceaf3d9843708aa6d22e9b84ccac8b4
nanoda_commit=68d5ca9db226849b41a6fff59d796ff19d0a8840
checkout_exact() {
  if [ ! -d "$2/.git" ]; then git clone --filter=blob:none "$1" "$2"; fi
  git -C "$2" fetch --depth 1 origin "$3"
  git -C "$2" checkout --detach "$3"
}
checkout_exact https://github.com/leanprover/lean4export.git "$cache/lean4export" "$lean4export_commit"
checkout_exact https://github.com/leanprover/comparator.git "$cache/comparator" "$comparator_commit"
checkout_exact https://github.com/robsimmons/nanoda_lib.git "$cache/nanoda" "$nanoda_commit"
GOBIN="$cache/bin" go install "github.com/zouuup/landrun/cmd/landrun@$landrun_commit"
(cd "$cache/comparator" && lake build comparator)
(cd "$cache/lean4export" && lake build lean4export)
(cd "$cache/nanoda" && cargo build --release --locked)
cd "$root"
COMPARATOR_LEAN4EXPORT="$cache/lean4export/.lake/build/bin/lean4export" \
COMPARATOR_NANODA="$cache/nanoda/target/release/nanoda_bin" \
COMPARATOR_LANDRUN="$cache/bin/landrun" \
  lake env "$cache/comparator/.lake/build/bin/comparator" "$root/comparator.json"
