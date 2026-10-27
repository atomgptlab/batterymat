"""Build the fixed-cell olivine Li0 inputs (ISIF=2 at the relaxed Li16 cell) from the olivine Li16 results.

Usage: python build_olivine_fixcell.py   (run in extra_2026-10-26/, after olivine-{LFP,LMP}/Li16/results/ exist)
Positions = olivine-<tag>/Li16/results/CONTCAR with all 16 Li removed; the cell is kept. INCAR, KPOINTS and
POTCAR_spec = olivine-<tag>/Li0 with ISIF 3 -> 2, so the pair Li16 (ISIF 3) / Li0-fixcell (ISIF 2) follows the
protocol of the benchmark chains (cell relaxed at x = 1 only).
"""
import shutil
from pathlib import Path
import numpy as np

for tag in ['LFP', 'LMP']:
    src = Path(f'olivine-{tag}/Li16/results/CONTCAR'); ref = Path(f'olivine-{tag}/Li0/POSCAR')
    out = Path(f'olivine-{tag}/Li0-fixcell'); out.mkdir(exist_ok=True)
    c = src.read_text().splitlines(); r = ref.read_text().splitlines()
    els, cnt = c[5].split(), [int(x) for x in c[6].split()]
    assert els[0] == 'Li' and cnt[0] == 16 and els[1:] == r[5].split() and cnt[1:] == [int(x) for x in r[6].split()]
    assert c[7].strip().lower().startswith('d')
    pos = c[8:8 + sum(cnt)][16:]
    # same atoms in the same order as the Li0 input (so its MAGMOM applies): compare against the unrelaxed positions
    p0 = np.array([[float(v) for v in l.split()[:3]] for l in r[8:8 + sum(cnt[1:])]])
    p1 = np.array([[float(v) for v in l.split()[:3]] for l in pos])
    dd = (p1 - p0 + 0.5) % 1.0 - 0.5
    print(tag, 'max displacement from unrelaxed Li0 input (fractional) %.4f' % np.abs(dd).max())
    assert np.abs(dd).max() < 0.05
    head = [f'olivine {"LiFePO4" if tag == "LFP" else "LiMnPO4"} x=0 at the relaxed x=1 cell: olivine-{tag}/Li16/results/CONTCAR with all Li removed, for ISIF=2']
    (out / 'POSCAR').write_text('\n'.join(head + c[1:5] + ['  ' + '  '.join(els[1:]), '  ' + '  '.join(str(x) for x in cnt[1:]), 'Direct'] + pos) + '\n')
    inc = Path(f'olivine-{tag}/Li0/INCAR').read_text()
    assert inc.count('ISIF = 3') == 1
    (out / 'INCAR').write_text(inc.replace('ISIF = 3', 'ISIF = 2'))
    for f in ['KPOINTS', 'POTCAR_spec']: shutil.copy(f'olivine-{tag}/Li0/{f}', out / f)
    print(out, 'written')
