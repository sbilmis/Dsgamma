# Photon normalization and scales — discussion, 25 October 2026

Status: source/code review complete; the author agrees to the investigation
plan recorded in `../../project.org`, but explicitly requires fresh-session
approval before calculations. No numerical evolution, matching calculation,
ensemble, or production change was performed. This issue is separate from
the settled tensor-current normalization typo in the earlier mixing paper.

## What the referee asks

Specific comment 2 of `/Users/sbilmis/Downloads/dv14652_report_1_a.pdf`
requests (i) the evolution from the lattice normalization at 2 GeV to 1 GeV
and the quark-mass/current scales, and (ii) compatibility of the finite-mass
lattice subtraction with the photon DA and the retained perturbative terms.
The second issue cannot be answered by documenting multiplicative running alone.

## Verified sources

- Bacchio et al., [arXiv:2404.10644](https://arxiv.org/pdf/2404.10644),
  Eqs. (17)–(18): the finite quantity is
  `Hbar_q = (1 - m_q d/dm_q) H_q`. This removes the UV divergence and also
  modifies finite mass terms; the chiral limit is preserved. The subtraction
  acts on the full correlator. Its finite nonperturbative effects cannot be
  presumed to be an entirely perturbative conversion. Their earlier lattice
  references [20,21] are Bali et al. 1209.6015 and 2004.08778.
- Ball–Braun–Kivel, [hep-ph/0207307](https://arxiv.org/pdf/hep-ph/0207307),
  Eq. (2.20): the finite-mass tensor matrix element has a perturbative
  logarithmic term in addition to the nonperturbative photon DA. The finite
  constant requires matching. Eqs. (3.4)–(3.5) give the homogeneous LO
  running of `f_perp = chi <qbar q>` as
  `[alpha_s(mu)/alpha_s(mu0)]^[4/(33-2*n_f)]` for a fixed flavor number.
  This formula alone does not establish finite-mass subtraction matching.
- Janowski–Pullin–Zwicky,
  [arXiv:2106.13616](https://arxiv.org/pdf/2106.13616), Sec. 4.2.1,
  Eqs. (4.11)–(4.20): a finite-mass DA prescription is constructed together
  with the perturbative contribution. The authors explicitly distinguish
  their susceptibility from the earlier lattice definitions in their
  refs. [65,66], which are the same Bali papers cited above. This supplies a
  concrete matching framework to examine, not an already validated conversion
  for our two-current transition kernels. Do not copy its normalization into
  production without translating the definitions and hard/soft separation.

PDFs and extracted text are saved under `references/`.

## Verified local implementation

- `scripts/lattice_photon_normalization_comparison.py:48–54` stores
  `FPERP_S_2GEV=-0.0510`, error `0.0018`, and the explicit estimate
  `TENSOR_RUNNING_2_TO_1=1.08`. It does not compute alpha_s or document
  threshold matching. The resulting stored normalization is -55.08 MeV.
  The factor is not newly calculated or verified in this discussion.
- Lines 84–87 set `chi=fperp/ss`, preserving the sampled product `chi*ss`.
  That algebraic reparameterization does not demonstrate scheme matching.
- `scripts/rohrwild_transition_exact.py:330–376` keeps finite masses in the
  hard contribution and inserts `ss*chi*phi_gamma` in twist two. No explicit
  conversion for the lattice finite-mass subtraction is visible in this
  assembly. This observation is not proof that the existing result is wrong.
- Table I of the author's current `r1-2.pdf` lists mc(mc), mb(mb), ms(2 GeV)
  separately from condensates/photon inputs at 1 GeV. The mass/current scale
  prescription must be stated and checked consistently; do not describe all
  current inputs as already evolved to a common scale.

## Agreed investigation, pending execution approval

1. Specify the mass, tensor-current and DA schemes/scales; implement a
   reproducible RG evolution with alpha_s inputs and flavor thresholds.
2. Establish the finite-ms definition used in our factorization, using the
   above primary sources. Identify the conversion to the lattice subtraction
   or what additional mass-dependent information it requires. A single
   published normalization need not determine that conversion.
3. Adopt a matched input/prescription only after this check; propagate any
   remaining matching uncertainty separately from the lattice error.
4. Then assess its absolute effect on the basis and suppressed amplitudes in
   the isolated rerun. No smallness conclusion follows from ms alone.

Do not apply the mass-derivative subtraction to the full heavy-meson sum
rule, reinterpret divergent unsubtracted lattice values as continuum inputs,
or infer a width correction from the size of the lattice subtraction.
No new photon input or hard/soft prescription has been selected here.
