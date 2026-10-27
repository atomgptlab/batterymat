"""Build olivine LiMnPO4 / LiFePO4 endpoint inputs (Li16 and Li0, 2x2x1, ISIF=3) from JVASP-91261.

Usage: BATTERYMAT_DFT3D_CACHE=<dft_3d.pkl> python build_olivine.py <out_dir> <template INCAR>
The template is LFP-ISIF3/Li0/INCAR (current input template with ISIF 3, NSW 300 and a cold start);
only the LDAU tags and MAGMOM are replaced. LiFePO4 is built by putting Fe on the Mn sites.
"""
import os, sys, pickle, math, itertools
import numpy as np, pandas as pd
from pathlib import Path
from pymatgen.core import Structure, Lattice

CACHE = os.environ['BATTERYMAT_DFT3D_CACHE']
OUT = Path(sys.argv[1])
TEMPLATE = Path(sys.argv[2])  # 6g Li0 INCAR (current template, ISIF 3, NSW 300, cold start)

d = pickle.load(open(CACHE, 'rb')); df = pd.DataFrame(d) if not isinstance(d, pd.DataFrame) else d
r = df[df['jid'] == 'JVASP-91261'].iloc[0]; a = r['atoms']
prim = Structure(Lattice(a['lattice_mat']), a['elements'], a['coords'], coords_are_cartesian=a['cartesian'])
print('JVASP-91261', prim.composition.reduced_formula, prim.get_space_group_info(0.1), [round(x, 3) for x in prim.lattice.abc], [round(x, 2) for x in prim.lattice.angles])
dim = [max(1, math.ceil(7.0 / n)) for n in prim.lattice.abc]
print('dim', dim)
sup = prim * dim

U = {'Mn': 3.90, 'Fe': 5.30}; M = {'Mn': 4, 'Fe': 5}; PAW = {'Li': 'Li_sv', 'Mn': 'Mn_pv', 'Fe': 'Fe_pv', 'P': 'P', 'O': 'O'}

# AFM: two-color the nearest-neighbour TM graph (corner-sharing octahedra within the bc layers)
tm = [i for i, s in enumerate(sup) if s.specie.symbol == 'Mn']
D = np.array([[sup.get_distance(i, j) if i != j else 99 for j in tm] for i in tm])
dmin = D.min(); nn = D < dmin + 0.3
print('TM-TM nearest %.3f A, neighbours per TM %s, next shell %.3f A' % (dmin, sorted(set(nn.sum(1))), np.sort(D[~nn & (D < 98)])[0]))
color = {}
for start in range(len(tm)):
    if start in color: continue
    color[start] = 1; stack = [start]
    while stack:
        u = stack.pop()
        for v in np.where(nn[u])[0]:
            if v not in color: color[v] = -color[u]; stack.append(v)
            elif color[v] == color[u]: sys.exit('NN graph not bipartite')
# relative sign of the layers: choose the one that makes the next shell most antiparallel
comps = []; seen = set()
for s0 in range(len(tm)):
    if s0 in seen: continue
    comp = {s0}; stack = [s0]
    while stack:
        u = stack.pop()
        for v in np.where(nn[u])[0]:
            if v not in comp: comp.add(v); stack.append(v)
    seen |= comp; comps.append(sorted(comp))
print('in-layer components', len(comps), [len(c) for c in comps])
shell2 = (D > dmin + 0.3) & (D < np.sort(D[~nn & (D < 98)])[0] + 0.3)
best = None
for flips in itertools.product([1, -1], repeat=len(comps) - 1):
    f = (1,) + flips; col = dict(color)
    for c, s in zip(comps, f):
        for k in c: col[k] = color[k] * s
    anti = sum(1 for i in range(len(tm)) for j in range(i + 1, len(tm)) if shell2[i, j] and col[i] != col[j])
    tot = sum(1 for i in range(len(tm)) for j in range(i + 1, len(tm)) if shell2[i, j])
    if best is None or anti > best[0]: best = (anti, tot, col)
print('second shell antiparallel %d of %d' % best[:2])
col = best[2]
assert all(col[i] != col[j] for i in range(len(tm)) for j in range(len(tm)) if nn[i, j])
spin = {tm[k]: col[k] for k in range(len(tm))}
assert sum(spin.values()) == 0

def write(tmel, with_li, path):
    path.mkdir(parents=True, exist_ok=True)
    keep = [i for i, s in enumerate(sup) if with_li or s.specie.symbol != 'Li']
    def key(i):
        el = sup[i].specie.symbol; order = {'Li': 0, 'Mn': 1, 'P': 2, 'O': 3}[el]
        return (order, -spin.get(i, 0), tuple(np.round(sup[i].frac_coords % 1.0, 6)))
    idx = sorted(keep, key=key)
    els = ['Li', tmel, 'P', 'O'] if with_li else [tmel, 'P', 'O']
    cnt = [sum(1 for i in idx if sup[i].specie.symbol == ('Mn' if e == tmel else e)) for e in els]
    form = ''.join(f'{e}{c}' for e, c in zip(els, cnt))
    lines = [f'olivine {"LiFePO4" if tmel=="Fe" else "LiMnPO4"} {"x=1" if with_li else "x=0"}: JVASP-91261 (mp-18997) 2x2x1' + (', Fe substituted for Mn' if tmel == 'Fe' else '') + (', all Li removed' if not with_li else '') + f', {form}, for ISIF=3', '1.0']
    for v in sup.lattice.matrix: lines.append('  ' + ''.join(f'{x:20.14f}' for x in v))
    lines.append('  ' + '  '.join(els)); lines.append('  ' + '  '.join(str(c) for c in cnt)); lines.append('Direct')
    for i in idx: lines.append(''.join(f'{x:20.14f}' for x in sup[i].frac_coords % 1.0))
    (path / 'POSCAR').write_text('\n'.join(lines) + '\n')
    nup = sum(1 for i in idx if spin.get(i, 0) == 1); ndn = sum(1 for i in idx if spin.get(i, 0) == -1)
    nonli = cnt[0] if with_li else 0
    mag = (f'{nonli}*0 ' if with_li else '') + f'{nup}*{M[tmel]} {ndn}*-{M[tmel]} {cnt[-2] + cnt[-1]}*0'
    if with_li: ldaul, ldauu, ldauj = '-1 2 -1 -1', f'0.00 {U[tmel]:.2f} 0.00 0.00', '0.00 0.00 0.00 0.00'
    else: ldaul, ldauu, ldauj = '2 -1 -1', f'{U[tmel]:.2f} 0.00 0.00', '0.00 0.00 0.00'
    inc = []
    for ln in TEMPLATE.read_text().splitlines():
        k = ln.split('=')[0].strip()
        ln = {'LDAUL': f'LDAUL = {ldaul}', 'LDAUU': f'LDAUU = {ldauu}', 'LDAUJ': f'LDAUJ = {ldauj}', 'MAGMOM': f'MAGMOM = {mag}'}.get(k, ln)
        inc.append(ln)
    (path / 'INCAR').write_text('\n'.join(inc) + '\n')
    (path / 'KPOINTS').write_text('Automatic mesh\n0\nGamma\n' + ' '.join(str(max(1, round(3 / n))) for n in dim) + '\n0 0 0\n')
    (path / 'POTCAR_spec').write_text('\n'.join(PAW[e] for e in els) + '\n')
    print(path, form, '| MAGMOM', mag)

for tmel, tag in [('Mn', 'LMP'), ('Fe', 'LFP')]:
    write(tmel, True, OUT / f'olivine-{tag}' / 'Li16')
    write(tmel, False, OUT / f'olivine-{tag}' / 'Li0')
