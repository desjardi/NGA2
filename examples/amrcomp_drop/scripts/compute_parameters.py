#!/usr/bin/env python3
"""
Convert known ratios (ambient / pre-shock) to input-file ratios
(liquid / post-shock gas) for both viscosity and thermal conductivity.

Viscosity
---------
  Code (line 203):  mu_l  = visc_ratio / Re
  Code (line 201):  mu_g  = (1+Suth_T)*T^Suth_n / (Re*(T+Suth_T))
  At T2=1:          mu_g2 = 1/Re
  => visc_ratio = muL1/muG2

  Conversion:  muL1/muG2  = (muL1/muG1) * (muG1/muG2)

Thermal conductivity
--------------------
  Code (line 210):  k_g = GammaG * CvG * mu_g / Pr
  Code (line 212):  k_l = diff_ratio * GammaG * CvG / (Re * Pr)
  At T2=1:          k_g2 = GammaG * CvG / (Re * Pr)
  => diff_ratio = kL1/kG2

  Since k_g = GammaG*CvG*mu_g/Pr  and  GammaG, CvG, Pr are constants:
      kG1/kG2 = muG1/muG2   (same Sutherland correction!)

  Conversion:  kL1/kG2  = (kL1/kG1) * (kG1/kG2)
                         = (kL1/kG1) * (muG1/muG2)

Usage: python3 compute_visc_diff_ratio.py input
"""
import re, sys, math

def read_input(filepath):
    """Parse the NGA2 input file into a dict."""
    params = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            m = re.match(r'^(.+?):\s+(.+)$', line)
            if m:
                params[m.group(1).strip()] = m.group(2).strip()
    return params

def main():
    # --- Read input file ---
    infile = sys.argv[1] if len(sys.argv) > 1 else 'input'
    params = read_input(infile)

    # --- Prompt user for necessary quantities and physical ratios (with defaults) ---
    print("\n  Please provide gammas:")
    GammaG_input = input("  GammaG (default=1.3): ").strip()
    GammaG = float(GammaG_input) if GammaG_input else 1.3
    GammaL_input = input("  GammaL (default=6.12): ").strip()
    GammaL = float(GammaL_input) if GammaL_input else 6.12

    print("\n  Please provide the shock Mach number:")
    Ms_input = input("  Shock Mach number Ms (default=5.12): ").strip()
    Ms = float(Ms_input) if Ms_input else 5.12

    print("\n  Please provide the following ambient ratios (no input corresponds to default ambient ratios for water/air):")
    density_input = input("  rhoL1/rhoG1 (Ambient density ratio, default=816.0): ").strip()
    rhoL1_over_rhoG1 = float(density_input) if density_input else 816.0

    sound_speed_input = input("  cL1/cG1 (Ambient sound speed ratio, default=4.3): ").strip()
    sound_speed_ratio = float(sound_speed_input) if sound_speed_input else 4.3

    visc_input = input("  muL1/muG1 (Ambient viscosity ratio, default=58.3): ").strip()
    muL1_over_muG1 = float(visc_input) if visc_input else 58.3

    cond_input = input("  kL1/kG1  (Ambient thermal conductivity ratio, default=23.1): ").strip()
    kL1_over_kG1 = float(cond_input) if cond_input else 23.1

    exp_input = input("  Sutherland exponent (default=1.5): ").strip()
    Suth_n = float(exp_input) if exp_input else 1.5

    temp_input = input("  Sutherland temperature (default=0.07286): ").strip()
    Suth_T = float(temp_input) if temp_input else 0.07286

    # --- Replicate Rankine-Hugoniot from simulation.f90 ---
    rhoG2 = 1.0
    rhoG1 = ((GammaG - 1.0) * Ms**2 + 2.0) / ((GammaG + 1.0) * Ms**2)
    
    pG1 = 0.25 * rhoG1 / GammaG * ((GammaG + 1.0) * Ms / (Ms**2 - 1.0))**2
    pG2 = pG1 * (2.0 * GammaG / (GammaG + 1.0) * (Ms**2 - 1.0) + 1.0)

    u1 = Ms * math.sqrt(GammaG * pG1 / rhoG1)
    u2 = u1 * rhoG1 / rhoG2

    u2=abs(u2-u1)
    M2=u2/math.sqrt(GammaG*pG2/rhoG2)
    u1=0
    M1=u1/math.sqrt(GammaG*pG1/rhoG1)  

    # Calculate liquid state properties
    rhoL = rhoL1_over_rhoG1 * rhoG1
    PinfL = pG1 * (rhoL1_over_rhoG1 * sound_speed_ratio**2 * GammaG / GammaL - 1.0)
    ML=u2 / math.sqrt(GammaL * (pG1 + PinfL) / rhoL)

    # CvG from T2=1 normalization
    CvG = pG2 / (rhoG2 * (GammaG - 1.0))

    # Pre-shock gas temperature (PinfG=0 for ideal gas)
    TG1 = pG1 / (CvG * rhoG1 * (GammaG - 1.0))

    # Sutherland ratio: muG1/muG2 (at T2=1, the denominator evaluates to 1)
    muG1_over_muG2 = ((1.0 + Suth_T)/(TG1 + Suth_T)) * (TG1**Suth_n)

    # Since k_g = GammaG*CvG*mu_g/Pr, and GammaG, CvG, Pr are constant:
    #   kG1/kG2 = muG1/muG2
    kG1_over_kG2 = muG1_over_muG2

    # Convert viscosity:  muL1/muG2 = (muL1/muG1) * (muG1/muG2)
    muL1_over_muG2 = muL1_over_muG1 * muG1_over_muG2

    # Convert diffusivity: kL1/kG2 = (kL1/kG1) * (kG1/kG2)
    kL1_over_kG2 = kL1_over_kG1 * kG1_over_kG2

    # --- Report ---
    print(f"{'='*60}")
    print(f"  Ratio conversion:  ambient (pre-shock, 1)  ->  code (post-shock, 2)")
    print(f"{'='*60}")
    print(f"  GammaG              = {GammaG}")
    print(f"  GammaL              = {GammaL}")
    print(f"  Shock Mach Ms       = {Ms}")
    print(f"  Post-shock Mach M2  = {M2}")
    print(f"  Liquid Mach ML      = {ML}")
    print(f"  Suth_n              = {Suth_n}")
    print(f"  Suth_T              = {Suth_T}")
    print(f"{'─'*60}")
    print(f"  rhoG1 (nondim)      = {rhoG1:.10f}")
    print(f"  pG1   (nondim)      = {pG1:.10f}")
    print(f"  TG1   (nondim)      = {TG1:.10f}")
    print(f"  muG1/muG2 = kG1/kG2 = {muG1_over_muG2:.10f}")
    print(f"{'─'*60}")
    print(f"  Viscosity:")
    print(f"    Input:  muL1/muG1 (known) = {muL1_over_muG1}")
    print(f"    Output: muL1/muG2 (for code) = {muL1_over_muG2:.6f}")
    print(f"  Thermal conductivity:")
    print(f"    Input:  kL1/kG1  (known)  = {kL1_over_kG1}")
    print(f"    Output: kL1/kG2  (for code)  = {kL1_over_kG2:.6f}")
    print(f"{'='*60}")
    print(f"\n  Set in your input file:")
    print(f"    GammaG:            {GammaG}")
    print(f"    GammaL:            {GammaL}")
    print(f"    Gas Mach number:   {M2:.6f}")
    print(f"    Density ratio:     {rhoL:.6f}")
    print(f"    Liquid Mach number:{ML:.6f}")
    print(f"    Viscosity ratio:   {muL1_over_muG2:.6f}")
    print(f"    Diffusivity ratio: {kL1_over_kG2:.6f}")
    print(f"    Sutherland exponent: {Suth_n}")
    print(f"    Sutherland temperature: {Suth_T}")

if __name__ == '__main__':
    main()
