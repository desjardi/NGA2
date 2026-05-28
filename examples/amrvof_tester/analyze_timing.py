#!/usr/bin/env python3
"""Analyze amrvof timing data from monitor/timing."""
import numpy as np
import matplotlib.pyplot as plt
import sys, os

# Read data
fname = sys.argv[1] if len(sys.argv) > 1 else 'monitor/timing'
data = np.loadtxt(fname, skiprows=1)
if data.ndim == 1:
    data = data.reshape(1, -1)

# Column mapping
step    = data[:, 0].astype(int)
time    = data[:, 1]
nm_max  = data[:, 2].astype(int)
nm_min  = data[:, 3].astype(int)
adv_max = data[:, 4];  adv_min = data[:, 5]
sl_max  = data[:, 6];  sl_min  = data[:, 7]
pl_max  = data[:, 8];  pl_min  = data[:, 9]
pn_max  = data[:, 10]; pn_min  = data[:, 11]
pg_max  = data[:, 12]; pg_min  = data[:, 13]

# ── Summary table ──
print(f"\n{'='*60}")
print(f"  AMRVOF TIMING ANALYSIS  ({len(step)} steps, t=[{time[0]:.3f}, {time[-1]:.3f}])")
print(f"{'='*60}\n")
print(f"  Comm-free timers (the ones that matter):")
print(f"  {'':10s} {'avg_max':>9s} {'avg_min':>9s} {'ratio':>6s} {'imbal':>9s}")
for name, mx, mn in [('SL',sl_max,sl_min),('PLICnet',pn_max,pn_min),('Polygon',pg_max,pg_min)]:
    am, an = mx.mean(), mn.mean()
    r = f"{am/an:.2f}x" if an > 0 else "n/a"
    print(f"  {name:10s} {am:9.5f} {an:9.5f} {r:>6s} {(am-an):9.5f}")

print(f"\n  Timers with blocking comms (imbalance hidden):")
print(f"  {'':10s} {'avg_max':>9s} {'avg_min':>9s}")
for name, mx, mn in [('advance',adv_max,adv_min),('PLIC',pl_max,pl_min)]:
    print(f"  {name:10s} {mx.mean():9.5f} {mn.mean():9.5f}")

print(f"\n  advance_vof breakdown (% of advance_max):")
sl_pct = sl_max.mean()/adv_max.mean()*100
print(f"    SL:      {sl_pct:5.1f}%  (remap + fluxes + update)")
print(f"    Other:   {100-sl_pct:5.1f}%  (setup, fill, comms)")
print(f"\n  build_plic breakdown (% of plic_max):")
pn_pct = pn_max.mean()/pl_max.mean()*100 if pl_max.mean() > 0 else 0
pg_pct = pg_max.mean()/pl_max.mean()*100 if pl_max.mean() > 0 else 0
print(f"    PLICnet: {pn_pct:5.1f}%")
print(f"    Polygon: {pg_pct:5.1f}%")
print(f"    Other:   {100-pn_pct-pg_pct:5.1f}%  (setup, comms)")

print(f"\n  Load distribution (mixed cells at finest level):")
print(f"    avg max: {nm_max.mean():.0f}   avg min: {nm_min.mean():.0f}   ratio: {nm_max.mean()/nm_min.mean():.2f}x")
print(f"    SL cost/cell ratio: {sl_max.mean()/sl_min.mean():.2f}x  vs  cell count ratio: {nm_max.mean()/nm_min.mean():.2f}x")
if sl_max.mean()/sl_min.mean() > nm_max.mean()/nm_min.mean() * 1.1:
    print(f"    ⚠ SL imbalance exceeds cell count imbalance → per-cell cost varies significantly")

# ── Plots ──
fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

# Panel 1: SL timing with min/max band
ax = axes[0]
ax.fill_between(time, sl_min*1e3, sl_max*1e3, alpha=0.3, color='C0', label='SL (min–max)')
ax.plot(time, sl_max*1e3, 'C0-', lw=0.5)
ax.plot(time, sl_min*1e3, 'C0-', lw=0.5)
ax.fill_between(time, pn_min*1e3, pn_max*1e3, alpha=0.3, color='C1', label='PLICnet (min–max)')
ax.plot(time, pn_max*1e3, 'C1-', lw=0.5)
ax.fill_between(time, pg_min*1e3, pg_max*1e3, alpha=0.3, color='C2', label='Polygon (min–max)')
ax.set_ylabel('Wall time (ms)')
ax.set_title('Comm-free component timing')
ax.legend(loc='upper left', fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 2: SL imbalance
ax = axes[1]
sl_imbal = (sl_max - sl_min) * 1e3
ax.plot(time, sl_imbal, 'C3-', lw=0.8, label='SL imbalance (max−min)')
ax.axhline(sl_imbal.mean(), color='C3', ls='--', lw=0.5, label=f'avg = {sl_imbal.mean():.1f} ms')
ax.set_ylabel('Imbalance (ms)')
ax.set_title('SL load imbalance')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 3: Mixed cell counts
ax = axes[2]
ax.fill_between(time, nm_min, nm_max, alpha=0.3, color='C4')
ax.plot(time, nm_max, 'C4-', lw=0.5, label='nmixed max')
ax.plot(time, nm_min, 'C4--', lw=0.5, label='nmixed min')
ax.set_ylabel('Mixed cells per rank')
ax.set_xlabel('Time')
ax.set_title('Mixed cell distribution across ranks')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
outfile = os.path.splitext(fname)[0] + '_analysis.png'
plt.savefig(outfile, dpi=150)
print(f"\n  Plot saved to: {outfile}")
plt.show()
