#!/bin/sh
# replay.sh -- AUDITOR'S from-scratch re-verification of monotonicity/cert/.
#
#   sh replay.sh                 everything, into a fresh directory replay-<timestamp>/
#   OUT=dir STAGES="0 1 5" sh replay.sh      selected stages, outputs in dir/
#   WORKERS=2 (default; at most 3)           worker processes for stage 4
#   MUTATION=1                               also run the mutation tests (stage 8)
#
# Nothing under monotonicity/cert/ or prototype/ is written (python runs with -B;
# every output goes to $OUT).  The script stops at the first failing check.
#
# Stages and measured run times (python-flint 0.9.0, mpmath 1.4.1, Python 3.13;
# wall-clock seconds of ONE process on the shared, heavily loaded audit machine,
# 2026-11-24 -- an idle core is roughly 1.5x faster):
#   0  sha256 of zero lists / certify_arb.py / audited code / results        <1 s
#   1  python-flint semantics relied upon                                    3 s
#   2  zero lists complete and simple: certify_arb.py --K 0 (hypothesis H3)  138 s
#   3  selftest; aux certificates recomputed and required to be byte-identical
#      to cert/aux_results.txt; auditor's independent aux check; tail sanity  5 s
#   4  both priority ranges re-certified FROM SCRATCH with the audited code
#      (317 + 296 cells): 1560 s + 1322 s with ONE worker (as run for the audit);
#      about 25 min wall with the default WORKERS=2
#   5  result files: exact tiling, all PASS, hi < 0, replay == audited        <1 s
#   6  auditor's independent enclosure of delta' on EVERY audited cell
#      (613 cells, 4 sub-cells each); must be negative and meet the audited one   361 s
#   7  plain mpmath evaluation of delta' at 20 points per modulus (in-cell test),
#      certified density grid, finite differences of delta                    375 s
#   8  (MUTATION=1) mutation tests, artificial positive family, driver tests,
#      enclosure-property tests of the unmutated code                         ~700 s
# Total without stage 8: about 63 min on one core, about 40 min wall with WORKERS=2.
set -e
cd "$(dirname "$0")"
AUD=$(pwd)
CERT=$AUD/../cert
EXP=$AUD/../../prototype/explore
OUT=${OUT:-replay-$(date +%Y%m%d-%H%M%S)}
STAGES=${STAGES:-"0 1 2 3 4 5 6 7"}
W=${WORKERS:-2}
[ "${MUTATION:-0}" = 1 ] && STAGES="$STAGES 8"
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
export AUDIT_OUT="$OUT"
PY="python3 -B"
has() { case " $STAGES " in *" $1 "*) return 0;; *) return 1;; esac; }
say() { echo; echo "=== stage $* ==="; }

if has 0; then say "0: integrity (SHA256SUMS)"
  (cd "$AUD/../.." && sha256sum -c monotonicity/cert-audit/SHA256SUMS)
fi
if has 1; then say "1: python-flint semantics"
  $PY flint_semantics.py > "$OUT/flint_semantics.log"; tail -1 "$OUT/flint_semantics.log"
fi
if has 2; then say "2: zero lists complete and simple (certify_arb.py --K 0)"
  $PY "$EXP/certify_arb.py" --K 0 > "$OUT/certify_arb_K0.log"
  grep -E "ordinates certified|winding ball|RESULT" "$OUT/certify_arb_K0.log"
  grep -q "chi4: 511/511 ordinates certified" "$OUT/certify_arb_K0.log"
  grep -q "chi3: 537/537 ordinates certified" "$OUT/certify_arb_K0.log"
  grep -q "RESULT: all certifications PASSED" "$OUT/certify_arb_K0.log"
fi
if has 3; then say "3: selftest, aux certificates, independent aux check"
  $PY "$CERT/certify_derivative.py" selftest > "$OUT/selftest.log"; tail -1 "$OUT/selftest.log"
  $PY "$CERT/certify_derivative.py" aux --out "$OUT/aux_results.txt" > /dev/null
  cmp "$OUT/aux_results.txt" "$CERT/aux_results.txt" && echo "aux_results.txt reproduced byte for byte"
  $PY aux_check.py > "$OUT/aux_check.log"; tail -1 "$OUT/aux_check.log"
  $PY tail_sanity.py > "$OUT/tail_sanity.log"; tail -1 "$OUT/tail_sanity.log"
fi
if has 4; then say "4: priority ranges from scratch ($W workers)"
  for q in 4 3; do
    if [ -e "$OUT/results_q$q.txt" ] && [ "${RESUME:-0}" != 1 ]; then
      echo "refusing to resume $OUT/results_q$q.txt (the driver would skip its cells); use a fresh OUT or RESUME=1"; exit 1
    fi
  done
  $PY "$CERT/certify_derivative.py" range --q 4 --lo 3.6 --hi 17.4 --ratio 1.005 --tol 0.02 --workers "$W" --out "$OUT/results_q4.txt"
  $PY "$CERT/certify_derivative.py" range --q 3 --lo 4.4 --hi 19.2 --ratio 1.005 --tol 0.02 --workers "$W" --out "$OUT/results_q3.txt"
fi
if has 5; then say "5: result files (tiling, PASS, hi<0, replay == audited)"
  AUDIT_REPLAY="$OUT" $PY verify_results.py --require-replay > "$OUT/verify_results.log" || { cat "$OUT/verify_results.log"; exit 1; }
  cat "$OUT/verify_results.log"
fi
if has 6; then say "6: independent enclosures on every audited cell"
  for q in 4 3; do
    if [ -e "$OUT/indep_q$q.txt" ] && [ "${RESUME:-0}" != 1 ]; then echo "refusing to resume $OUT/indep_q$q.txt"; exit 1; fi
    $PY compare.py indep --q $q --nsub 4 | tee "$OUT/compare_indep_q$q.log"
    $PY compare.py summary --q $q > /dev/null      # non-zero exit unless all cells meet and are negative
  done
fi
if has 7; then say "7: plain mpmath values, finite differences, certified density grid"
  for q in 4 3; do
    if [ $q = 4 ]; then LO=3.6; HI=17.4; FD="3.7 9 17.3"; else LO=4.4; HI=19.2; FD="4.5 10 19.1"; fi
    if [ ! -e "$OUT/plain_q$q.txt" ]; then $PY plain_eval.py --q $q --npoints 20 --lo $LO --hi $HI --out "$OUT/plain_q$q.txt"; fi
    $PY compare.py plain --q $q > "$OUT/compare_plain_q$q.log"; tail -1 "$OUT/compare_plain_q$q.log"
    $PY compare.py grid --q $q > "$OUT/compare_grid_q$q.log"; tail -1 "$OUT/compare_grid_q$q.log"
    $PY plain_eval.py --q $q --fd $FD > "$OUT/plain_fd_q$q.txt"; cat "$OUT/plain_fd_q$q.txt"
  done
fi
if has 8; then say "8: mutation tests"
  cmp mutants/certify_derivative_copy.py "$CERT/certify_derivative.py"
  export MUT_TMP="$OUT/mutants_tmp"; mkdir -p "$MUT_TMP"
  $PY mutants/mutate.py all "$OUT/mutation_results.json" > "$OUT/mutation_all.log"; cat "$OUT/mutation_all.log"
  $PY mutants/mutate.py art > "$OUT/artificial_family.log"; cat "$OUT/artificial_family.log"
  $PY mutants/mutate.py driver > "$OUT/driver_checks.log"; tail -1 "$OUT/driver_checks.log"
  $PY mutants/proptest.py baseline > "$OUT/proptest_baseline.log"; cat "$OUT/proptest_baseline.log"
  # (mutate.py art / driver and proptest.py also refresh their .json summaries in mutants/)
fi
echo
echo "replay finished: stages [$STAGES] OK, outputs in $OUT"
