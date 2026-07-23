# Course Notes: Ds1 -> Ds gamma in LCSR

## Step 2: OPE bookkeeping before traces

The correlator is organized as

```text
Pi_OPE = Pi_hard + Pi_soft .
```

`Pi_hard` is the perturbative photon-emission sector.  The photon is attached to
either the strange-quark line or the heavy-quark line through one electromagnetic
insertion.  In the base analysis both quark lines use free massive propagators,

```text
S_s(k) = i (kslash + m_s)/(k^2 - m_s^2),
S_Q(k) = i (kslash + m_Q)/(k^2 - m_Q^2).
```

All explicit powers of `m_s` are kept.  We do not insert the usual vacuum
condensate expansion of the quark propagator in this base OPE.

`Pi_soft` is the nonperturbative photon sector.  Here the relevant light-quark
bilinear is replaced by photon distribution-amplitude matrix elements,

```text
<gamma(q)| sbar(x) Gamma s(0) |0>.
```

Condensate parameters that appear in the standard photon DA normalizations, such
as `chi <sbar s>`, are retained as photon-wavefunction inputs.  They are not
ordinary condensate corrections to a propagator and should not be double counted
in `Pi_hard`.

The two axial-channel interpolating currents are

```text
J_A,mu = sbar gamma_mu gamma5 Q,
J_B,mu = i sbar sigma_{mu alpha} P^alpha gamma5 Q / (m_Q + m_s),
```

where `P = p + q` is the momentum carried by the initial `Ds1` channel.  This
choice is important: `J_B` interpolates the axial meson, so it must contract with
the axial-channel momentum rather than the final `Ds` momentum.

The first Mathematica artifact is
`scripts/step2_current_building_blocks.wl`.  It loads FeynCalc, defines the A,
B, and pseudoscalar vertices, and performs minimal Dirac-algebra checks before
we build the full hard-photon traces.

The second Mathematica artifact is
`scripts/step2_hard_photon_numerators.wl`.  It computes numerator-level traces
for the hard photon emission sector:

```text
N_hard = e_Q N_heavy + e_s N_strange.
```

At this stage these are only Dirac numerators.  Denominators, color factors,
loop integration, double dispersion, and Borel transformation are not yet
attached.  This separation is intentional: it lets us verify the current
structure, signs, and `m_s` dependence before doing the analytic integration.

## LaTeX notes

The main pedagogical notes now live in
`notes/Ds1_to_Ds_gamma_LCSR_notes.tex`, with the compiled PDF at
`notes/Ds1_to_Ds_gamma_LCSR_notes.pdf`.  The notes are written as a living
derivation for later paper writing, while this Markdown file remains a compact
work log.
