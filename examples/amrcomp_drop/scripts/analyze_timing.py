#!/usr/bin/env python3
"""Analyze timing monitor output for load imbalance and communication overhead."""

# Parse the timing file
lines = open('monitor/timing').readlines()
header = lines[0].split()
data = []
for line in lines[1:]:
    vals = line.split()
    if len(vals) == len(header):
        data.append([float(v) for v in vals])

# Skip timestep 0 (init) and first 10 (startup transient)
data = [row for row in data if row[0] >= 11]
n = len(data)
cols = {h: i for i, h in enumerate(header)}

def mean(v): return sum(v)/len(v) if v else 0
def std(v):
    m = mean(v)
    return (sum((x-m)**2 for x in v)/(len(v)-1))**0.5 if len(v)>1 else 0

print(f'Analyzing {n} timesteps (skipping first 10 + init)')
print(f'Columns: {", ".join(header)}')
print()

# Phase timing
print('='*80)
print('PHASE TIMING (seconds per timestep, max across ranks = wall-clock)')
print('='*80)
for p in ['dQdt_max', 'plic_max', 'relax_max', 'visc_max']:
    v = [row[cols[p]] for row in data]
    print(f'{p:>16s}: mean={mean(v):.4f}  std={std(v):.4f}  min={min(v):.4f}  max={max(v):.4f}')

# dQdt breakdown
print()
print('='*80)
print('dQdt BREAKDOWN (seconds per timestep)')
print('='*80)
dqdt = [row[cols['dQdt_max']] for row in data]
for p in ['prim_max', 'sl_max', 'fv_max', 'div_max']:
    v = [row[cols[p]] for row in data]
    frac = [vi/di*100 for vi,di in zip(v, dqdt) if di > 0]
    print(f'{p:>16s}: mean={mean(v):.4f}  ({mean(frac):.1f}% of dQdt)')
comm = [row[cols['dQdt_max']] - row[cols['prim_max']] - row[cols['sl_max']] - row[cols['fv_max']] - row[cols['div_max']] for row in data]
frac = [c/d*100 for c,d in zip(comm, dqdt) if d > 0]
print(f"{'comm_overhead':>16s}: mean={mean(comm):.4f}  ({mean(frac):.1f}% of dQdt)")

# PLIC breakdown
print()
print('='*80)
print('PLIC BREAKDOWN (seconds per timestep)')
print('='*80)
pt = [row[cols['plic_max']] for row in data]
pn = [row[cols['plicnet_max']] for row in data]
pp = [row[cols['polygon_max']] for row in data]
pc = [t-n-p for t,n,p in zip(pt, pn, pp)]
for label, v in [('plicnet_max', pn), ('polygon_max', pp), ('fill_plic_comm', pc)]:
    frac = [vi/ti*100 for vi,ti in zip(v, pt) if ti > 0]
    print(f'{label:>16s}: mean={mean(v):.4f}  ({mean(frac):.1f}% of plic)')

# Load imbalance
print()
print('='*80)
print('LOAD IMBALANCE  (max-min)/max  (0=perfect, 1=worst)')
print('='*80)
for label, mx, mn in [('prim','prim_max','prim_min'), ('sl','sl_max','sl_min'),
                       ('fv','fv_max','fv_min'), ('div','div_max','div_min'),
                       ('plicnet','plicnet_max','plicnet_min'), ('polygon','polygon_max','polygon_min')]:
    imb = [(row[cols[mx]]-row[cols[mn]])/row[cols[mx]] if row[cols[mx]]>0 else 0 for row in data]
    print(f'{label:>16s}: mean={mean(imb):.3f}  std={std(imb):.3f}  max={max(imb):.3f}')

# Load distribution diagnostics
if 'cells_max' in cols:
    print()
    print('='*80)
    print('LOAD DISTRIBUTION (finest-level cells per rank)')
    print('='*80)
    cm = [row[cols['cells_max']] for row in data]
    cn = [row[cols['cells_min']] for row in data]
    mm = [row[cols['mixed_max']] for row in data]
    mn = [row[cols['mixed_min']] for row in data]
    cell_imb = [(mx-mn)/mx if mx>0 else 0 for mx,mn in zip(cm,cn)]
    mix_imb  = [(mx-mn)/mx if mx>0 else 0 for mx,mn in zip(mm,mn)]
    print(f'  cells_max  : mean={mean(cm):.0f}')
    print(f'  cells_min  : mean={mean(cn):.0f}')
    print(f'  cell imbal : mean={mean(cell_imb):.3f}  (max-min)/max')
    print(f'  mixed_max  : mean={mean(mm):.0f}')
    print(f'  mixed_min  : mean={mean(mn):.0f}')
    print(f'  mixed imbal: mean={mean(mix_imb):.3f}  (max-min)/max')
    print(f'  mixed ratio: mean max/min = {mean([mx/mn if mn>0 else float("inf") for mx,mn in zip(mm,mn)]):.1f}x')
    # Correlation: does mixed imbalance predict SL imbalance?
    sl_imb = [(row[cols['sl_max']]-row[cols['sl_min']])/row[cols['sl_max']] if row[cols['sl_max']]>0 else 0 for row in data]
    print(f'  SL imbal   : mean={mean(sl_imb):.3f}')
    print(f'  Correlation: mixed_imb vs sl_imb')
    if len(sl_imb) > 1:
        mx_m = mean(mix_imb); sl_m = mean(sl_imb)
        cov = sum((a-mx_m)*(b-sl_m) for a,b in zip(mix_imb, sl_imb))/(len(sl_imb)-1)
        r = cov/(std(mix_imb)*std(sl_imb)) if std(mix_imb)>0 and std(sl_imb)>0 else 0
        print(f'    Pearson r = {r:.3f}')

# Total time budget
print()
print('='*80)
print('TOTAL TIME BUDGET PER TIMESTEP')
print('='*80)
total = [row[cols['dQdt_max']]+row[cols['plic_max']]+row[cols['relax_max']]+row[cols['visc_max']] for row in data]
for label, key in [('dQdt','dQdt_max'), ('plic','plic_max'), ('relax','relax_max'), ('visc','visc_max')]:
    v = [row[cols[key]] for row in data]
    frac = [vi/ti*100 for vi,ti in zip(v, total) if ti > 0]
    print(f'{label:>16s}: mean={mean(v):.4f}s  ({mean(frac):.1f}%)')
print(f"{'TOTAL':>16s}: mean={mean(total):.4f}s")
