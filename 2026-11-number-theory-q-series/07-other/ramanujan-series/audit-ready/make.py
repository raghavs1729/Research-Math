"""Generate audit-ready/T11.md and T12.md from RESULTS.md (run from monotonicity/ramanujan-series):
self-contained statement + conventions + proof + certificate, so that an independent audit or a Lean
formalization can start from one file."""
s = open('RESULTS.md').read()


def section(a, b):
    i0 = s.index(a)
    i1 = s.index(b, i0 + 10)
    return s[i0:i1].rstrip().rstrip('-').rstrip()


conv = section("## 0. Conventions and the two general facts used", "## Theorem 1.")
HEAD = """# {name} — self-contained statement, certificate and proof (audit-ready)

Generated from `RESULTS.md` (section 0 and {name}) by `audit-ready/make.py`, so that an independent audit or a Lean
formalization can start from this file alone. Machine-readable certificate: `{gp}` (reproduced at the end).
Exact check script: `{script}` (output `{out}`).

Status: {status}

"""
STATUS = {
    'T11': "proof read line by line by two independent audits (no error found). Lean: certificate identities only "
           "(`lean/RamanujanSeries/Sun1pi2.lean`); Part B (Laplace method) is a written proof, not machine-checked. "
           "Novelty: not found in any source read; canonical data absent from the Omega / Phi / 2H2 families "
           "(RESULTS.md, last section); Chu-Zhang 2014 and Chu 2021 unread.",
    'T13': "NOT yet independently audited; not formalized. Step 1 of the proof is quoted from Au (arXiv:2312.14051, "
           "Example I) and step 3 uses a published harmonic-number identity (Au; Wei). Novelty: an observation of Cohen and "
           "Guillera (2021) and eq. (2) of D. Chen's MathOverflow question 486961 (2025), where Guillera indicated the proof "
           "route; NOT a new statement — only the written-out proof is ours.",
    'T12': "NOT yet independently audited; not formalized. Novelty: a conjecture of Guillera (2018); not found proved in "
           "any source read; in blocks of three terms it is an integer Dougall pattern (see the Remark), so Chu-Zhang "
           "2014 and Chu 2021 (unread) are a residual risk.",
}
parts = {'T11': ('Theorem 11', section("## Theorem 11.", "## Theorem 12 (added after the audit)."), 'certificates/T11.gp',
                 'proofs/verify_T11_sun_1pi2_sqrt2.py'),
         'T12': ('Theorem 12', section("## Theorem 12 (added after the audit).", "## Theorem 13 (added 2026-11-25)."), 'certificates/T12.gp',
                 'proofs/verify_T12_guillera_shift_48.py'),
         'T13': ('Theorem 13', section("## Theorem 13 (added 2026-11-25).", "## What resisted, and why"), 'certificates/T13.gp',
                 'proofs/verify_T13_cullen_half_shift.py')}
for key, (name, body, gp, script) in parts.items():
    txt = (HEAD.format(name=name, gp=gp, script=script, out=script.replace('.py', '.out.txt'), status=STATUS[key])
           + conv + "\n\n---\n\n" + body + "\n\n---\n\n## Certificate file `" + gp + "` (PARI/GP; `sg` = sigma)\n\n```\n"
           + open(gp).read() + "```\n\nCheck: `gp -q < " + gp + "` prints `1` on every line; `python3 " + script
           + "` ends with `ALL CHECKS PASSED`.\n")
    open('audit-ready/%s.md' % key, 'w').write(txt)
    print(key, len(txt.splitlines()), 'lines')
