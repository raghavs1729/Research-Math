#!/bin/bash
set -x
cd "$(dirname "$0")"
if [ ! -x "$HOME/.elan/bin/lake" ]; then
  curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -o /tmp/elan-init.sh && sh /tmp/elan-init.sh -y --default-toolchain none
fi
export PATH="$HOME/.elan/bin:$PATH"
lake update 2>&1 | tail -5
lake exe cache get 2>&1 | tail -5
lake build 2>&1 | tail -5
