#!/bin/bash
# sha256 of all outputs and data files of this directory (run from monotonicity/hkmn23/)
{ ls *.out *.txt *.json 2>/dev/null; find raw -type f \( -name "*.json" -o -name "*.jsonl" -o -name "*.log" -o -name "*.time" \) | sort; } | xargs sha256sum > SHA256SUMS
wc -l SHA256SUMS
