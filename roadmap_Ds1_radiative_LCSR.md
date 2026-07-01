# Roadmap: Ds1(2460)→Ds γ and Ds1(2536)→Ds γ in QCD Light-Cone Sum Rules
### (with the c̄s → b̄s / Bs1 counterparts)

**Method:** Light-Cone Sum Rules (LCSR) in the electromagnetic background-field
formulation, with photon distribution amplitudes (DAs) for the soft photon
emission. This is the Rohrwild / Aliev–Özpineci machinery, extended to the two
1⁺ interpolating currents and their physical mixing.

---

## 0. Where this sits in the literature (done — see chat)

| Reference | States | Method | Overlap / gap |
|---|---|---|---|
| Colangelo, De Fazio, Özpineci 2005 (hep-ph/0505195) | Ds1(2460) only | LCSR, single axial current | **Closest precedent.** No 2536, no two-current mixing |
| Rohrwild 2007 (0708.1405) | D*→Dγ | LCSR + photon DAs, background field | **Template** for OPE/DA machinery |
| Aliev et al., PRD 54 (1996) 857 | B*,D*→Bγ,Dγ | LCSR + photon wavefunction | Template, cross-check limit |
| Wang 2006 (hep-ph/0612225) | Ds1(2460) | LCSR via φ(1020) VMD | Indirect; different mechanism |
| Yang 2007 (0705.0692) | 1³P1 / 1¹P1 axial DAs | conformal expansion + SR | Supplies current structure & mixing |
| Bondar–Milstein 2025 (PRD 111 114019) | both Ds1, Bs1, Bc1 | relativistic potential model | **Numerical benchmark** for our Γ |
| arXiv 2510.17477 (2025) | both, molecular | LCSR magnetic moments | Different observable & picture |

**Novelty:** unified LCSR for *both* physical 1⁺ states in the mixed
(γμγ₅ ⊕ iσ·p γ₅) current basis, m_s ≠ 0, photon DAs, + Bs sector.

---

## 1. Physical framework & definitions

### 1.1 Interpolating currents (flavor / spin eigenbasis)
Two independent 1⁺ currents (Q = c or b, q = s):

- **A-type (axial-vector):**  J^{(0)}_{Aμ} = q̄ γ_μ γ₅ Q
- **B-type (tensor-axial):**  J^{(0)}_{Bμ} = i q̄ σ_{μα} P^α γ₅ Q / (m_Q + m_s),
  with `P = p + q` the momentum carried by the initial `Ds1` channel.

The 1/(m_Q+m_s) fixes the mass dimension of J_B to equal that of J_A (both
dimension 3). *[Confirmed against user's Eq. (3); the denominator is required.]*

### 1.2 Physical states via mixing (user's Eq. 2)
- Ds1(2460):  J_{A_i μ} =  sinθ · J^{(0)}_{Aμ} + cosθ · J^{(0)}_{Bμ}
- Ds1(2536):  J_{A'_i μ} = cosθ · J^{(0)}_{Aμ} − sinθ · J^{(0)}_{Bμ}

θ is the ³P₁–¹P₁ (equivalently j_ℓ=1/2 ↔ 3/2) mixing angle. In the heavy-quark
limit θ → ideal (≈ 35.3°); we keep it as an input and study sensitivity.

### 1.3 Decay constants (couplings of currents to states)
- ⟨0| J^{(0)}_{Aμ} |Ds1(p,η)⟩ = f_A m_{Ds1} η_μ
- ⟨0| J^{(0)}_{Bμ} |Ds1(p,η)⟩ = f_B m_{Ds1} η_μ  (with the 1/(m_Q+m_s) absorbed)
- ⟨0| J_{Ds} |Ds(p)⟩ = f_{Ds} m_{Ds}² / (m_Q + m_s),  J_{Ds} = q̄ i γ₅ Q
  (pseudoscalar Ds interpolating current)

Values from Pullin–Zwicky (f_{H1}, f^T_{H1}) and Yang; f_{Ds} from PDG/lattice.

### 1.4 Transition matrix element & width (E1 transition, 1⁺ → 0⁻ γ)
The gauge-invariant amplitude (real photon, q²=0, ε_γ·q=0):

  ⟨Ds(p) γ(q,ε)| j^{em}_μ |Ds1(p+q, η)⟩
     = e · G · [ (η·ε*_γ)(p·q) − (η·q)(p·ε*_γ) ]

G is the form factor (target of the sum rule) at q²=0. This is a **P→S electric
dipole (E1)** transition (Ds1 is p-wave, Ds is s-wave), hence the non-ε_{μναβ}
structure — contrast D*(1⁻)→D(0⁻)γ which is ε_{μναβ} (M1). *Cross-check to
build in.*

Decay width (E1 scaling):
  Γ = (α_em / 3) · G² · |q_γ|³ · (kinematic factor to be derived cleanly),
  |q_γ| = (m²_{Ds1} − m²_{Ds}) / (2 m_{Ds1}).

The exact prefactor is fixed in Step 5 by matching the amplitude to the
polarization sum.

---

## 2. Correlation function & OPE (FeynCalc-heavy)

### 2.1 The correlator in the photon background field
  Π_μ(p,q) = i ∫ d⁴x e^{i p·x} ⟨γ(q)| T{ J_{Ds1,μ}(x) J_{Ds}^†(0) } |0⟩

with J_{Ds1,μ} the *mixed* current (Step 1.2). Linear response to the external
photon field F_{αβ} = i(ε_α q_β − ε_β q_α).

### 2.2 Wick contractions → two topologies
1. **Perturbative / hard photon emission:** photon attaches to the c(b) line or
   the s line through one external electromagnetic-field insertion. This sector
   uses free **massive** quark propagators and is done in momentum space. Keep
   all explicit powers of `m_s`; do not add ordinary vacuum-condensate
   corrections to the quark propagators in the base analysis.
2. **Soft photon emission:** photon couples to the s q̄ pair at long distance →
   replaced by **photon DA** matrix elements ⟨γ| s̄ Γ s |0⟩ (twist 2,3,4;
   2- and 3-particle). Condensate parameters that appear in photon-DA
   normalizations, such as χ⟨s̄s⟩, are kept as part of the photon wavefunction,
   not as condensate-expanded propagator terms.

### 2.3 Propagators and photon insertion
- **Base hard-emission propagators:**
  `S_s(k)=i(k̸+m_s)/(k²-m_s²)` and `S_Q(k)=i(k̸+m_Q)/(k²-m_Q²)`, with one
  electromagnetic insertion on either the strange or heavy line. The strange
  quark is never treated as massless unless we explicitly take an `m_s → 0`
  cross-check.
- **No ordinary condensate-expanded propagators in the base OPE:** omit terms
  such as `−⟨s̄s⟩/12`, mixed-condensate terms, and gluon-condensate propagator
  corrections at this stage. These may be added later only as separately labeled
  higher-power corrections.
- **Soft-emission replacement:** when the photon is soft/nonperturbative,
  replace the relevant light-quark bilinear by photon DA matrix elements rather
  than by a vacuum-condensate-expanded light propagator.
- **Representation:** use momentum space for the hard-emission traces; use the
  natural light-cone coordinate representation for photon DAs, then Fourier/Borel
  transform to the same invariant amplitude.

Reference set for the hard insertion and photon DAs: Ball–Braun–Kivel photon
DAs; Rohrwild; Colangelo–De Fazio–Özpineci; background-field LCSR conventions
from Balitsky–Braun / Colangelo–Khodjamirian.

### 2.4 Traces & simplification → FeynCalc
Dirac traces, σ-matrix algebra, contraction with F_{αβ}, projection onto the E1
Lorentz structure. Output: invariant amplitude Π^{OPE}(p², (p+q)²) per structure,
per current (A, B), including the A–B interference needed for mixing.

---

## 3. Photon distribution amplitudes (from Rohrwild / Ball–Braun–Kivel)

Tabulate the complete set with numerical parameters at μ = 1 GeV:
- **twist-2:** φ_γ(u), magnetic susceptibility χ (≈ 3.15 GeV⁻² [BBK] / 2.85 [Rohrwild])
- **twist-3:** ψ^v, ψ^a
- **twist-4:** h_γ, 𝔸, and 3-particle 𝒮, 𝒯_i, …
- **3-particle params:** f_{3γ}, ω^V_γ, ω^A_γ, κ, κ⁺, ζ_i

Each enters the OPE convolution ∫du. Encoded as a reusable module (numbers +
functional forms) for both the D and B sectors.

---

## 4. Hadronic (phenomenological) side

Double dispersion relation in p² (Ds channel) and (p+q)² (Ds1 channel):

  Π^{had} = [ f_{Ds1} f_{Ds} m_{Ds1} m²_{Ds}/(m_Q+m_s) · G ]
            / [ (m²_{Ds1} − (p+q)²)(m²_{Ds} − p²) ] + (continuum) + (excited)

Ground-state pole isolates G. Continuum & excited states modeled by
quark-hadron duality above thresholds s₀ (Ds1 channel), s₀' (Ds channel).

---

## 5. Sum rule construction & extraction of G

1. Equate Π^{OPE} = Π^{had}.
2. **Double Borel transform** (p² → M₁², (p+q)² → M₂²). Since m_{Ds1} ≈ m_{Ds}
   are not equal, keep M₁²,M₂² distinct or use M² with u₀ = M₁²/(M₁²+M₂²);
   continuum subtraction accordingly.
3. Solve for G_A and G_B (one sum rule per current), including interference.
4. **Assemble physical form factors via mixing:**
   - G[Ds1(2460)] = sinθ · G_A + cosθ · G_B
   - G[Ds1(2536)] = cosθ · G_A − sinθ · G_B
5. Insert into Γ (Step 1.4).

---

## 6. Numerical analysis, stability, results

- Inputs: quark masses (m_c, m_b, m_s), condensates ⟨q̄q⟩, ⟨s̄s⟩, m₀², ⟨G²⟩,
  χ & photon-DA params, decay constants f_{Ds1,A}, f_{Ds1,B}, f_{Ds}, θ.
- **Working-window analysis:** Borel M² stability plots; continuum-threshold
  s₀ dependence; pole vs continuum fraction; convergence of the twist expansion.
- Deliver Γ[Ds1(2460)→Dsγ], Γ[Ds1(2536)→Dsγ] with uncertainty budget.
- **Benchmark** against Bondar–Milstein Table I and Colangelo et al.

## 7. Bs (b̄s) counterparts

Repeat Steps 2–6 with Q = b: Bs1(5778)/Bs1(?)→Bs γ. Reuse all symbolic output
(same currents, same DAs), swap heavy-quark inputs, re-tune Borel windows.

## 8. Cross-checks (applied continuously, not just at the end)
- **Mass dimension** of every intermediate expression.
- **m_s → 0** limit: SU(3) structure, χ⟨q̄q⟩ term behavior.
- **Heavy-quark limit** m_Q → ∞: correct 1/m_Q scaling of G and f's; match HQET.
- **Gauge invariance / transversality** of the amplitude structure.
- **D*→Dγ reduction** as an independent sanity limit of the machinery.
- Independent re-derivation of Dirac traces in a second symbolic engine.

---

## Software & division of labor
- **FeynCalc (user's Mathematica):** Dirac/σ traces, contractions — I will
  generate ready-to-run `.wl`/notebook cells.
- **Independent symbolic cross-check:** I re-run the traces in a second engine
  (sympy / custom gamma-algebra) to catch sign/convention errors.
- **Package-X / FeynHelpers:** if the perturbative photon-emission diagram needs
  one-loop scalar integrals.
- **Python (me):** photon-DA numerics, Borel/continuum analysis, stability
  plots, width computation, uncertainty propagation, all figures & tables.

### Execution workflow (settled)
Mathematica/FeynCalc runs locally through the direct kernel path:
`/Applications/Wolfram.app/Contents/MacOS/WolframKernel`. Plain `wolframscript`
does not reliably find the kernel, so scripts are run as
`WolframKernel -noprompt -script path/to/script.wl`.

1. I generate and run fully-commented `.wl` FeynCalc scripts.
2. Outputs are written to `outputs/` for review and later projection/integration.
3. I independently reproduce selected traces in `sympy` / custom gamma algebra
   as the second-engine cross-check.

### Code style requirement (course-oriented)
**Every code file — `.wl` and `.py` — is written with pedagogical comments:**
each block states what it computes, why, the physics/identity being used, and
what a student should check. The goal is that a master's student can read the
code top-to-bottom and follow the derivation without external notes. A parallel
LaTeX note file, `notes/Ds1_to_Ds_gamma_LCSR_notes.tex`, grows with each step
and is compiled to PDF for later paper-writing use. `course_notes_Ds.md` remains
a compact work log.

### Step-by-step discipline
Finalize the **Ds (charm) sector completely** — both Ds1(2460) and Ds1(2536) —
and get user sign-off before starting the Bs (bottom) sector. One step at a
time, each producing its artifact + course-notes section, with a verification
pause.

## Settled Step-2 choices
1. **OPE split:** include both hard photon emission and soft photon DA terms.
   Hard emission uses free massive quark propagators with one EM insertion; soft
   emission uses photon DA matrix elements.
2. **No condensate-expanded propagators in the base analysis:** ordinary vacuum
   condensate terms are omitted from the quark propagators. Condensate
   normalizations inside photon DAs are retained as photon-DA inputs.
3. **Representation:** momentum space for hard-emission traces; light-cone
   coordinate representation for photon DAs, followed by Fourier/Borel matching.
4. **Strange-quark mass:** keep explicit `m_s` terms throughout, including
   possible `m_s²` contributions.
5. **Tensor-current momentum:** `J_B` contracts with `P^α=(p+q)^α`, the momentum
   carried by the initial `Ds1` channel.
6. **Staging:** Stage 1 contains hard photon emission plus two-particle photon
   DAs. Stage 2 adds the background-gluon propagator insertion needed for
   three-particle photon DAs.

## Open decisions to settle before Step 2 numerics (need user input)
1. **Which Ds1(2536) partner mass/width** and mixing-angle value(s) to adopt
   (ideal 35.3° vs fitted); scan or fix?
2. **Twist truncation:** through twist-4 + 3-particle (full, Rohrwild-level) or
   start at twist-3 and add twist-4 as a correction?
