# Step 1 — Conventions, Currents, Decay Constants
### Ds1 → Ds γ radiative transitions in LCSR (charm sector)

This file fixes every convention the rest of the calculation depends on. Any
FeynCalc or sympy script must use exactly these choices, or signs and factors
of i will not match between engines.

---

## 1. Metric, Dirac, Levi-Civita, Fourier

| Object | Convention |
|---|---|
| Metric | g_{μν} = diag(+1, −1, −1, −1) (Bjorken–Drell / mostly-minus) |
| Dirac algebra | {γ_μ, γ_ν} = 2 g_{μν}; γ5 = i γ⁰γ¹γ²γ³; γ5² = 1 |
| σ_{μν} | σ_{μν} = (i/2)[γ_μ, γ_ν] |
| Slash | a̸ = γ_μ a^μ |
| Trace norm | Tr[1] = 4; Tr[γ_μγ_ν] = 4g_{μν}; Tr[γ5 γ_μγ_νγ_αγ_β] = −4i ε_{μναβ} |
| Levi-Civita | ε^{0123} = +1, ε_{0123} = −1 (⇒ the −4i above) |
| Fourier | f(x) = ∫ d⁴k/(2π)⁴ e^{−ik·x} f̃(k); propagator carries e^{−ik·x} |
| Charge | e > 0, electron charge = −e; quark charge Q_q e with Q_u=+2/3, Q_d=Q_s=−1/3, Q_c=+2/3, Q_b=−1/3 |
| α_em | α = e²/(4π) = 1/137.036 (real photon, q²=0) |

**FeynCalc note:** these are FeynCalc's defaults *except* verify
`$LeviCivitaSign = -1` to match ε_{0123} = −1. State it explicitly in each script.

---

## 2. Interpolating currents

Heavy quark Q ∈ {c, b}; light quark q = s. External momenta: the 1⁺ meson
carries P = p + q (p = final Ds momentum, q = photon momentum).

### 2.1 Flavor/spin eigenbasis (the two "bare" currents)
```
J_A^(0)_μ (x) = q̄(x) γ_μ γ5 Q(x)                         (axial-vector)
J_B^(0)_μ (x) = i q̄(x) σ_{μα} P^α γ5 Q(x) / (m_Q + m_s)  (tensor-axial)
```
where `P = p + q` is the momentum carried by the initial `Ds1` channel. The
tensor current interpolates the axial state, so its momentum contraction is with
the `Ds1` momentum, not the final `Ds` momentum.

**Mass-dimension check (crucial):**
- q̄ … Q has dimension 3 for a bilinear with no derivative.
- J_A: γ_μγ5 carries no dimension ⇒ dim(J_A) = 3. ✓
- J_B *without* the denominator: σ_{μα} P^α carries one power of momentum ⇒
  dim = 4. ✗ (mismatch with J_A)
- Dividing by (m_Q + m_s) (dimension 1) restores dim(J_B) = 3. ✓

This is why the user's Eq. (3) needs the (m_Q+m_s) in the denominator: it makes
the two currents dimensionally homogeneous so they can be linearly combined in
Eq. (2) and share a single mixing angle.

### 2.2 Physical states via mixing (user's Eq. 2)
```
Ds1(2460):  J_μ = sinθ · J_A^(0)_μ + cosθ · J_B^(0)_μ
Ds1(2536):  J_μ = cosθ · J_A^(0)_μ − sinθ · J_B^(0)_μ
```
θ is the ³P₁–¹P₁ (equivalently j_ℓ = 1/2 ↔ 3/2) mixing angle. Because the two
currents are orthogonal combinations, the physical form factors are the SAME
orthogonal rotation of the two "bare" form factors G_A, G_B (see §5).

### 2.3 Final-state pseudoscalar current
```
J_Ds (x) = q̄(x) i γ5 Q(x)
```
with ⟨0| J_Ds |Ds(p)⟩ = f_Ds m_Ds² / (m_Q + m_s).

---

## 3. Decay-constant definitions

```
⟨0| J_A^(0)_μ |Ds1(P, η)⟩ = f_A m_{Ds1} η_μ
⟨0| J_B^(0)_μ |Ds1(P, η)⟩ = f_B m_{Ds1} η_μ
⟨0| J_Ds       |Ds(p)⟩     = f_Ds m_Ds² / (m_Q + m_s)
```
η_μ = polarization vector of the 1⁺ state, with η·P = 0 and Σ η_μ η*_ν =
−g_{μν} + P_μ P_ν / m_{Ds1}². Numerical values and sources: `inputs_table.csv`.

---

## 4. Transition amplitude & Lorentz structure (E1)

The electromagnetic transition matrix element (real photon, ε_γ·q = 0):
```
⟨Ds(p) γ(q, ε_γ)| j^{em}_μ |Ds1(P, η)⟩
      = e · G · [ (ε_γ*·η)(p·q) − (ε_γ*·p)(η·q) ]
```
This is exactly Colangelo–De Fazio–Özpineci Eq. (4.1) (their coupling g₁ = our
G_A when only J_A is used). It is the **electric-dipole (E1)** structure
appropriate to a p-wave(1⁺) → s-wave(0⁻) transition — note it has NO ε_{μναβ},
in contrast to the magnetic-dipole (M1) D*(1⁻)→D(0⁻)γ, which does. This
sign/structure difference is a built-in cross-check.

**Gauge invariance:** replacing ε_γ → q gives G·[(q·η)(p·q) − (q·p)(η·q)] =
G(p·q)[(q·η) − (η·q)] ... — vanishes upon using p·q and the antisymmetric
combination; verified in course notes. The structure is transverse by
construction.

---

## 5. Width formula (derived, not quoted)

Averaging over the 3 initial polarizations and summing over the 2 photon
polarizations (full derivation in `course_notes_Ds.md` §5):
```
Σ_pol |M|² = 2 e² G² (p·q)²,     p·q = m_{Ds1} |q_γ|
⟨|M|²⟩ = (1/3) Σ|M|²
Γ = |q_γ| / (8π m_{Ds1}²) · ⟨|M|²⟩
```
which simplifies to the compact, dimensionally-correct result
```
┌─────────────────────────────────────────────┐
│   Γ[1⁺ → 0⁻ γ] = ( α / 3 ) · G² · |q_γ|³     │
└─────────────────────────────────────────────┘
   |q_γ| = ( m_{Ds1}² − m_Ds² ) / ( 2 m_{Ds1} )
```
G has mass dimension −1 (as in Colangelo, GeV⁻¹), so G²|q_γ|³ has dimension
+1 = energy. ✓

**Numerical validation:** with Colangelo's |g₁| ∈ [0.29, 0.37] GeV⁻¹ and
|q_γ| = 442 MeV for Ds1(2460), this formula gives Γ ∈ [17.7, 28.8] keV,
reproducing their published (19–29) keV. Formula + normalization confirmed.

---

## 6. Benchmarks to hit (charm sector)

| Quantity | Value | Source |
|---|---|---|
| Γ[Ds1(2460)→Dsγ] | 297 keV | Bondar–Milstein 2025 (potential model) |
| Γ[Ds1(2536)→Dsγ] | 12 keV | Bondar–Milstein 2025 |
| Γ[Ds1(2460)→Dsγ] | (19±29→ i.e. large) — via g₁ | Colangelo 2005 LCSR (single current) |
| Exp bound Γ[2536→Dsγ] | < 8 keV | PDG (from BM Table) |
| Exp bound Γ[2460→Dsγ] | (branching, small) | PDG |

Note the strong model dependence in the literature — the 2536→Dsγ width spans
6–61 keV across predictions. Our LCSR value with mixing is the deliverable.

---

## 7. Open normalization to settle in Step 3–4
- Exact relation between f_B (our J_B decay constant) and the Pullin–Zwicky
  f^T_{H1}: their tensor current is q̄ σ_{αβ} Q; ours contracts one index with
  P^α/(m_Q+m_s). The proportionality (a factor ~ m_{Ds1}/(m_Q+m_s)) will be
  pinned when we write the J_B two-point sum rule.
