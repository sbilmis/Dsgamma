# Referee replies and calculation TODO

## Canonical project.org handoff — latest 25 October 2026

The author agrees to the proposed investigation and requests a populated
`project.org` based on the personal `tproj` snippet. It is now the canonical
task queue, with one NEXT: ask for explicit execution approval in the fresh
Codex session before starting calculations. This TODO remains the history
and evidence index. No physics calculation is authorized merely by the
creation of the plan; no rerun or figure regeneration has started.

The project contains the central a=1/s_eff/eta/r prescription, technical
prerequisites, shared-input MC plan, current manuscript figures 1--4,
additional signed Pi_A/Pi_B/Pi_2 panels, response-only f2/f1, plot-style details,
and acceptance criteria. The current manuscript and referee documents are
copied under `recalculation_a1_20261025/baseline/reference_documents/` with
hash provenance. Read `../project.org` first, summarize the concrete plan,
ask approval, and proceed after the author answers. Production and
final_manuscript remain read-only. This latest plan supersedes the earlier
"further agreement" hold: the next required action is the explicit fresh-
session approval, not a repetition of the already discussed choices.

The author subsequently requested a local commit. Its scope is the handoff
and supporting review record; pre-existing production edits are excluded.
Caches/build products and the four large provisional raw .npz scans remain
local, with their scripts and summary diagnostics preserved in the checkpoint.
No push or physics execution was requested by this commit instruction.

## Photon-normalization discussion — 25 October 2026

The author asks to discuss specific comment 2 before the full rerun.
Read `referee_a_photon_normalization_20261025/README.md`. Source/code review
distinguishes RG evolution from finite-ms subtraction matching. The current
factor 1.08 is a stored LO estimate, not an implemented alpha_s evolution.
BBK Eq. (2.20) and Janowski–Pullin–Zwicky Sec. 4.2.1 explicitly discuss the
finite-mass hard/soft separation; the latter distinguishes its definition
from the earlier lattice subtraction used by the 2024 lattice reference.
This establishes a matching question, not a numerical correction or proof
that the existing predictions are wrong. No matching or numerical evolution
has been calculated, and no replacement input has been adopted. The explicit
hold on calculations remains in force. Production and final_manuscript
were not edited.

## Intercept choice; explicit hold on starting calculations — 25 October 2026

After the geometry explanation the author says to use the recommended eta
(typed "n"), and asks about the possible size of the change, unequal Borel
parameters, corresponding DA arguments, and readiness. The author explicitly
instructs: do not start calculations before agreement on these points.
`recalculation_a1_20261025/decisions.json` now records eta=m_A^2/m_P^2 and
an execution hold. The Borel ratio r is a separate choice, intended centrally
as m_A^2/m_P^2, giving u0=1/(1+r) in our manuscript convention. Read the
new README opening. This discussion only updates records; no calculations
or new figures were run. The previous 5.18-keV diagnostic must not be used
as a forecast of the full a=1 continuum result.

## Agreed a=1 and effective threshold; rerun preparation — 25 October 2026

The author now agrees to the linear domain a=1 and sigma_0=s_eff=s_0 from
Table II, and requests that the full recalculation and regenerated figures
stay in review, separate from previous results. The current turn asks what
must be discussed before that run. `recalculation_a1_20261025/` contains
the preparation, decisions, 409 protected-file hashes and 92 source copies.
No physics recalculation has been run in that directory yet.

The remaining intercept-ratio choice was asked asynchronously: recommend
eta=tilde_s_A/tilde_s_P=m_A^2/m_P^2, keeping it fixed during Borel-ratio
variation. The author replied that they could not decide and requested an
explanation. No intercept ratio is selected; a=1 and sigma_0=s_0 are settled.
A geometric illustration in `recalculation_a1_20261025/discussion/` shows
how fixed s_eff fixes the diagonal crossing while eta changes the slope.
The domain ratio eta is independent of the Borel ratio r. The paper itself
uses mass-based intercepts with equal Borel parameters numerically. No
new OPE or width was evaluated for this illustration.
The README separates this modelling choice from required technical checks:
general-Borel hard/soft kernels and DA conventions, finite-ms photon input
matching/scales, pole separation, leading-d5 residue correction and shared
QCD-input sampling, DA/EOM uncertainty design, then paired benchmarks and
the full ensemble. All four current manuscript figures and the relevant
response-only figures/tables are listed for both sectors.

## Pullin--Zwicky continuum formula discussion — later on 25 October 2026

The author decided to keep the Borel parameters unequal/general in the next
analysis, but explicitly requested agreement on the continuum prescription
first. The published JHEP09(2021)023 paper is accessible and now saved under
`referee_a_continuum_planning_20261025/references/`. Read the new final section
of that directory's README. C.4 and D.10--D.15 supply the relevant linear-domain
framework with boundary terms; D.18 is a symmetric-limit control that the
current uniform-soft subtraction does not reproduce for all F_n,0 kernels.
The raw uncut n=2 maps agree after exchanging the paper's hadronic labels.
Paper q is the axial/vector momentum, not our real-photon momentum.

Exact independent test kernels validate C.4 for the linear a=1 domain through
n=4. A separate a=2 curved-domain test exposes an apparent missing curvature
term in the printed general-a expression; this does not affect the a=1 check.
The proposed working choice is a=1 with sigma matched initially to existing
Table II s0, and a mass-based ratio of domain intercepts for discussion.
The domain choice is not yet adopted; no new physics numerics were performed.

The author also reminded us that the working windows were already selected.
The final transition windows and two-point acceptance criteria were verified
against saved production records. Start from these existing windows and
reassess them after formula changes; do not describe the selection as absent.

## Continuum discussion before further recalculation — 25 October 2026

The author asked whether manuscript predictions must change, requested
perturbative/total diagnostics, and asked to discuss continuum subtraction
and agree on a recalculation roadmap before doing further width calculations.
Read `review/referee_a_continuum_planning_20261025/README.md`. This turn
extracts existing component ratios and evaluates hard spectral-retention
diagnostics only; no new continuum variant, ensemble, patch or manuscript edit
is adopted. Do not treat this planning request as permission to skip discussion
and regenerate the production results immediately.

The signed pert/total ratio is NOT a pole fraction. For Ds2536 at its midpoint
it is 19.379 (symmetric) or 3.344 (mass-based diagnostic). Its basis A/B ratios
are 2.396/1.549 or 1.872/1.469. All four channels are tabulated in the new
directory. The hard absolute spectral retention at the midpoints is A/B:
Ds2460 0.694/0.771, Ds2536 0.726/0.799, Bs5750 0.513/0.781,
Bs5830 0.531/0.798. Bottom A falls to 0.396/0.413 at window edges. These
are perturbative proxies, not positive full-transition pole probabilities.

The referee's continuum comment is substantive: specify the two-variable
duality domain and derive each subtraction with surface terms. Colangelo
hep-ph/0505195 Eq. (4.4) and Rohrwild 0708.1405 Eq. (20) display different
treatments. Neither citation alone establishes the present rule. Older review
notes claiming the all-soft replacement was proven correct established an
implementation choice, not the missing duality derivation; do not repeat that
unqualified claim. The proposed roadmap groups convention/DA/momentum fixes,
continuum derivation, input matching and pole separation, window selection,
correlated uncertainties, and then synchronized manuscript/reply updates.
If a validated calculation confirms a material numerical or uncertainty change,
the manuscript must reflect it; the provisional 5.18-keV point is not yet a
replacement for the manuscript's 0.55-keV MC median.

## Unequal-Borel diagnostic completed — later on 25 October 2026

The author requested the actual calculation in a new folder. The isolated
package is `review/referee_a_borel_asymmetry_20261025/`; read its README before
using the numerical table. It scans all four channels, reproduces all 372
saved symmetric Borel points, and preserves all production/manuscript hashes.
At fixed reference inputs and Table II midpoints, r=m_A^2/m_P^2 changes widths
(keV) from 47.2989 to 50.2675, 0.182716 to 5.183620, 12.0623 to 11.6803,
and 8.66890 to 8.12156 for Ds2460, Ds2536, Bs5750, Bs5830 respectively.
These are fixed-point diagnostic values, NOT the MC medians, NOT a new
uncertainty extraction, and NOT yet final referee evidence.

The primary branch uses the common heavy-line fraction 1/(1+r) in the stored
off-shell kernels, and keeps the existing diagonal hard spectral and soft
continuum models. The printed manuscript uses u0 for the two-particle terms
but 1-u0 for the three-particle terms; the older derivation calls the common
fraction a0. Reflecting only the three-particle argument gives a very different
Ds2536 mass-ratio width, 0.053050 keV. This is a convention diagnostic only,
not a second physical prediction. Printed versus numerical antisymmetric DA
argument orientation and the B(u) coefficient also need reconciliation. No
DA signs or parameters were changed. Do not select the smaller result merely
because it supports stability, or claim the asymmetry effect is negligible.
The general hard double dispersion, continuum domain and individual-channel
admissibility remain unvalidated outside the preserved model. The table and
plot are for author review; a final reply paragraph has not been adopted.

## Referee A major comment (ii): three-panel diagnostic — 25 October 2026

The author subsequently accepted the Pi_2 plot and proposed the
mass-proportional Borel choice M1^2/M2^2=m_A^2/m_P^2. Its use was verified in
Rohrwild, hep-ph/0701085v3, Eq. (49); it is an auxiliary prescription, not
an exact relation. Using manuscript conventions, u0=1/(1+r), so the higher
charm-strange channel has r=1.65879 and u0=0.37611. The kinematic ratios for
all four channels and the source are recorded in the major-ii README.
The author wants this numerical comparison in the referee response only.
No unequal-Borel widths have yet been evaluated; a table of paired symmetric,
mass-proportional, and ratio-scan results is proposed. New scans must use
consistent general-Borel kernels, continuum treatment, and pole exponents.

Latest steering: the author proposed using Pi_2=cos(theta)Pi_A-sin(theta)Pi_B
in panel (c), with the manuscript's plot style. A revised review-only figure
is [ds2536_signed_invariants.png](referee_a_major_ii/ds2536_signed_invariants.png),
with SVG, LaTeX caption, midpoint table, and a separate plotting script.
It uses the unchanged saved data; the original g2 figure is retained.
All panels now have units GeV^3 and match the manuscript's sans-serif fonts,
boxed axes, inward ticks, and core colours. The suggested next presentation
is a sensitivity table covering both sectors and a focused unequal-Borel
scan plot. Ranges and DA constraints still need establishing; these new
sensitivity calculations have not been performed. See that directory's
README for the Borel-parameter convention and twist-four constraint check
to resolve before a consistent parameter scan. No production inputs changed.

The author moved through the two major comments of
`/Users/sbilmis/Downloads/dv14652_report_1_a.pdf` and explicitly requested
the proposed three-panel plot. The first diagnostic is now saved in
[referee_a_major_ii/](referee_a_major_ii/README.md): Pi_A, Pi_B, and the
signed perturbative/twist decomposition of g2 for Ds1(2536) -> Ds gamma.
It uses the preserved manuscript Borel-scan inputs, theta=26.6 degrees,
M^2=3--4.5 GeV^2, and s0=9--10 GeV^2. Central curves use s0=9.5;
bands are threshold variation only. The combined twist-four curve includes
all retained two- and three-particle terms, also tabulated separately.

All 93 saved Borel-scan points and 42 saved central numerical entries are
reproduced exactly. At the midpoint, g2 is the sum -0.470228 (perturbative),
+0.370088 (twist 2), +0.027027 (twist 3), and +0.048849 (combined twist 4)
GeV^-1, giving -0.024265 GeV^-1. Thus the projected higher-twist terms
are comparable to or larger than the residual. This is a sensitivity
observation, not a determination of omitted terms or proof of OPE failure.
Only this channel has been plotted so far. The other-channel diagnostics,
DA/asymmetric-Borel scans, continuum assessment, and a current uncertainty
budget remain pending. The older Stage-2 budget uses a different baseline
and must not be presented as the current manuscript's budget.

For major comment (i), the author cited arXiv:1007.0814 and preferred a
conceptual explanation of the zero mixed-correlator construction. Proposed
wording was discussed. The normalization typo is settled, but no new
present-OPE off-diagonal check or pole-contamination bound was performed;
do not mark that numerical request resolved. No production correction,
baseline replacement, or manuscript modification has been adopted.

## Final response-format review — later on 24 October 2026

The author clarified that the requested final review concerns the response
letter, not the manuscript. The supplied r1-2.pdf is a 25-page manuscript.
A targeted check of the same Downloads folder found the newer actual reply:
`/Users/sbilmis/Downloads/ref1_reply-1.pdf`, created at 16:47:23 on 24 October.
All four pages of that reply have been read and visually inspected. It
supersedes ref1_reply.pdf for the current response review. Its headings use
"Referee B"; preserve that if it matches the editor's designation.

The author opted for the short Comment 3 response: thank the referee,
describe the added individual-constant manuscript figure, and supply the
ratio plot in the response. No long weak-current derivation or ratio
benchmark values are included in this version. Do not automatically
replace it with the longer proposed response fragment.

Final presentation suggestions (not yet applied):

- Add a document title and a brief general acknowledgement before Comment 1.
- Page 2 contains only the closing paragraph of Comment 1. Allow Comment 2
  to follow it instead of forcing a separate page; retain the readable font.
- Label the response plot Figure R1 to distinguish it from manuscript
  Figure 1; use automatic label references.
- Comment 1 still says "throughout the sampled regions" although the short
  reply presents only midpoints. Prefer "at these representative points".
  The sentence about Eq. (46) can say "mixing angles adopted from our earlier
  study" rather than "resulting mixing angles".
- Expand the short ratio caption to specify the two sectors, fixed
  representative inputs with recomputation of both residues at each angle,
  and the meanings of the vertical reference lines. This preserves the
  short response without implying a new uncertainty prediction.
- A separate conclusion is unnecessary for this reply. No PDF, manuscript,
  or numerical source/result was modified during the format review.

## Latest response decisions — 24 October 2026

- The author selected the shorter Comment 1 response: identify the added
  method equations and normalization footnote, then report the original
  conditional reconstruction's midpoint shifts +0.185/+0.057 degrees.
  Retain the phrase "using a reconstruction of the earlier sum-rule
  expressions" and the central-input-check qualification. Use
  `-(0.240)^3 GeV^3` for the current condensate. The proposed minor edit
  narrows "throughout the sampled regions" to "at these representative
  points" when only midpoints are presented. This choice does not adopt
  the completed-current-OPE alternative in the earlier proposal file,
  validate the original program, or authorize production/manuscript changes.
- Comment 2 says the vector-channel study is nearing completion and will
  include the suggested ratios. The redundant final sentence about a
  separate publication can be deleted. No Conclusion addition is planned.
- The author has now moved to Comment 3. A fresh proposed fragment is
  [comment3_response_proposal_20261024.tex](referee1_revision/comment3_response_proposal_20261024.tex).
  It retains individual f1/f2 curves in manuscript Figure 1 and puts the
  existing f2/f1 figure in the reply only, with representative values
  0.418 (Ds) and 0.166 (Bs), explicitly without new uncertainty estimates.
  It includes the conditional weak-current relation needed for the referee's
  nonleptonic application; no branching-fraction calculation or promise to
  put that application in the vector-radiative paper is made.
- The last supplied manuscript still needs its Figure 1/2 image swap fixed
  and its Figure 3 discussion/caption changed from "central values" to fixed
  representative QCD/auxiliary inputs with angle-dependent residues. The
  numerical Figure 3 scan already includes that dependence. The fresh reply
  does not claim these wording fixes have already been implemented.
- The older combined `referee1_revision/referee1_reply.tex` and PDF are
  historical drafts with obsolete Comment 1/2 wording. Do not reuse them
  as the author's latest accepted response. The live reply is edited by the
  author. All saved figures still use the preserved original numerical baseline.

## Current-document review — 24 October 2026

The author supplied `/Users/sbilmis/Downloads/r1-1.pdf` and
`/Users/sbilmis/Downloads/ref1_reply.pdf` and asked to finalize the response
one comment at a time. These are the current reference documents; their
hashes are respectively
`3af0acef95de252d215cac59b071b047587be700a8101bf41977cf7e7182b271` and
`0f160492548abce15cf810917b2ad81f92340ac7dec7b5e1bd59736d2cc5159d`.
The earlier page locations below remain historical.

Comment 1 review:

- [x] Verify the method insertion: manuscript pp. 11--12, Eqs. (39)--(43),
  including the current-normalization footnote. Rotation and inverse-relation
  signs are consistent. Published angle inputs appear in Eq. (46), p. 13.
- [ ] Add an opening sentence in the reply stating exactly where the
  manuscript was expanded; keep the existing method explanation if desired.
- [ ] Replace the outdated author-review note and select the numerical reply
  prescription explicitly. The current reply p. 2 contains the old conditional
  reconstruction, +0.185/+0.057 degrees. The completed current-OPE check is
  +0.175/+0.045 degrees with grid ranges 0.116--0.331 and 0.036--0.060 degrees.
  It must be identified as an independent LO + quark + leading-mixed-condensate
  sensitivity check, with gluon-condensate and radiative terms omitted, not
  as a rerun of the original full calculation. Proposed complete wording is
  [comment1_response_proposal_20261024.tex](referee1_revision/comment1_response_proposal_20261024.tex).
  This is for discussion, not a silent adoption of the correction.
- [ ] Bring manuscript Appendix C, Eqs. (C2)--(C3), and the numerical
  implementation/results into agreement with the saved leading-d5 correction
  before final submission. The supplied PDF still has the original light-only
  d5 terms. No new calculation was run during this document review.

Immediate assembly correction: Figures 1 and 2 have exchanged image content.
Manuscript p. 14 shows the four-panel Borel plot under the new decay-constant
caption; p. 15 shows the two-panel decay-constant plot under the Borel caption.
Exchange the image arguments, preserving the intended captions and labels:
Figure 1 = f1/f2 versus theta; Figure 2 = form factors versus Borel mass.
This is visually confirmed, not an extraction artifact, and is not yet fixed.

Later-item observations are queued without finalizing them: the new response
for Comment 2 is present on reply p. 3; Comment 3 has only the quoted comment
on p. 4. The new decay-constant discussion/caption is inserted, but the old
"central values" wording for the width-angle scan persists on manuscript
pp. 15--16. Address this when reaching Comment 3. Do not claim all requested
changes are implemented or that the response is submission-ready.

No manuscript, uploaded PDF, production source, or numerical baseline has
been edited. Only review notes and a proposed response fragment were saved.

Latest discussion, 2026-10-24: the author returned to Referee 1 item 2 and
asked to describe the vector-final-state study as nearing completion. The
author prefers to keep this information in the referee response without a
Conclusion addition. The proposed response should omit the claim that a
statement has already been added to the Conclusion. This does not adopt the
pending numerical correction or resume ensemble regeneration.

Updated 2026-10-19. The proposed revision package is in [referee1_revision/](referee1_revision/README.md), with individual decay-constant curves for both sectors, a response-only f2/f1 plot, paste-ready LaTeX, exact insertion instructions using r1.pdf, and a reply draft. The live manuscript is in Prism and will be edited by the author. A conditional condensate-only sensitivity calculation is now complete in [referee1_item1/](referee1_item1/README.md); its connection to the original implementation remains unresolved, so the final item 1 response is still provisional.

## Handoff for the next session — 20 October 2026

The author left on 19 October and explicitly requested notes and continuation
tomorrow. Work is paused until the author returns. No scheduled or background
work was requested.

Where we stopped:

- The latest discussion concerned the mixed-condensate coefficient in Referee 1
  item 1. The independent tree-level leading-d5 derivation found a missing
  heavy-line background-gluon contribution in the current AB and BB entries;
  its AA entry vanishes. This explains why the pure axial-current check could
  pass. The normalization typo is already settled and must not be reopened.
- The derivation, seven symbolic controls, 24 covariant trace checks, and
  proposed Python/Wolfram implementation checks are complete. Read
  [the resolution report](mixed_condensate_resolution/README.md),
  [the LaTeX note](mixed_condensate_resolution/dimension5_correction.tex), and
  [the unapplied patch](mixed_condensate_resolution/proposed_implementation.patch).
- At the four representative preferred-angle points, the correction changes
  f by less than 0.2% and the widths by less than 0.35% with transition
  numerators held fixed. These are not updated ensemble predictions. Threshold
  fits, accepted ensembles, uncertainty summaries, and figures have not been
  regenerated with the completed d5 matrix.
- Keep the two condensate comparisons distinct: +0.185164/+0.057096 degrees
  (Ds/Bs) belong to the printed-appendix reconstruction; +0.174981/+0.045126
  degrees belong to the completed current LO+d3+d5 matrix. Neither is a verified
  rerun of the original 2025 program. The leading coefficient is derived;
  the original full-OPE implementation and uncertainty prescription remain
  separate provenance questions. The published angles remain external inputs.
- No production source, baseline result, or manuscript file was replaced.
  The proposed patch remains unapplied. All earlier proposed f(theta) and
  f2/f1 figures still reflect the original production OPE.

Resume in this order:

1. Discuss the saved coefficient finding and the exact scope of the correction
   with the author. Do not repeat the derivation or original-file search merely
   because a new session has started.
2. After agreement on adopting the correction, regenerate the affected
   two-point ensembles and dependent results in a separate review directory,
   including acceptance criteria and any threshold refits. Compare against
   the preserved baseline; keep the charm shared-input sampling issue separate.
3. Establish which validated OPE prescription will support the response-only
   condensate comparison, then finalize the quantitative wording for item 1.
   Do not present the old conditional numbers as newly verified.
4. Continue Referee 1 item by item. Preserve the agreed presentation scope:
   individual f1(theta), f2(theta) for both sectors in the manuscript; f2/f1
   in the referee response only. Update affected figures if the correction is
   adopted, and give the author insertion instructions for Prism.

Keep `final_manuscript/` read-only. Saving this handoff is not approval to
replace manuscript or production results. There is no running calculation
to resume or terminate.

## Current status

Latest follow-up (19 October 2026): the author asked to resolve the
mixed-condensate coefficient issue. The independent derivation in
[mixed_condensate_resolution/](mixed_condensate_resolution/README.md) identifies
a missing heavy-line background-gluon contribution in the current AB and BB
entries; its AA entry is zero. The leading dimension-five coefficient is now
derived and independently checked. The current code and the earlier printed
expression are both incomplete/different in this limit; selecting a factor of
one or one-half would not give the complete result. An isolated implementation,
fixed-point impact calculation, and unapplied consistent source patch are saved.
No production outputs or manuscript files have been replaced. The original
larger-OPE calculation and its uncertainty procedure are separate open questions.

Preceding calculation request (19 October 2026): the author asked to perform the
condensate-value comparison for the response only, expecting a small effect.
The completed conditional reconstruction gives shifts +0.185164 degrees (Ds)
and +0.057096 degrees (Bs) at the earlier window midpoints. The respective
21x21-grid ranges are +0.118680 to +0.368437 and +0.044416 to +0.079634 degrees.
An independent Wolfram transcription agrees with Python numerically. These
are smaller than the quoted errors across the sampled windows, but do not
resolve the earlier audit's original-implementation/OPE questions. See
[the comparison and provisional response wording](referee1_item1/README.md).

Latest discussion decision (19 October 2026): retain the published
$\theta_{D_s}=26.6^\circ\pm0.6^\circ$ and
$\theta_{B_s}=38.5^\circ\pm0.1^\circ$ as external inputs and return to
Referee 1, one item at a time. Further searching for the original angle code
is deferred for now. This is a choice of inputs, not a verified rerun or a
numerical resolution of the updated-input question. Keep the open audit and
shared-input charm MC tasks below separate from drafting the agreed revision
text. No manuscript insertion or numerical replacement is implied.

- [x] Complete an initial correctness/provenance audit and archive confirmed
  inactive figure bundles; see [the audit report](correctness_audit_20261019/REPORT.md).
- [x] Independently derive the complete leading mixed-condensate matrix and
  identify the omitted heavy-line contribution. Verify with symbolic control
  channels and covariant Dirac/color traces.
- [x] Quantify its fixed-point impact and prepare an unapplied consistent patch
  for the Python source, both Wolfram checks, and notebook builder.
- [ ] Review and adopt the correction, regenerate accepted ensembles and affected
  figures in a separate output directory, and compare the complete uncertainty
  results. The fixed-point changes do not substitute for this regeneration.
- [ ] Establish the original larger-OPE numerical prescription, including the
  printed label/sign mismatch and any higher-order terms, before describing a
  condensate check as a rerun of that calculation. The leading-d5 coefficient
  itself is no longer unresolved.
- [ ] Discuss and adopt a consistent shared-input charm uncertainty sampler.
  The exact replay and paired variant are saved separately; baseline outputs
  and all manuscript files remain unchanged.

- [x] Complete a paired central-condensate check in the explicitly marked appendix reconstruction, including a working-window grid and independent numerical verification.
- [ ] Discuss this conditional result and resolve its original-implementation qualification before finalizing the quantitative item 1 response.
- [x] Prepare the item 1 method clarification and the item 2 follow-up-study response.
- [ ] Verify which proposed text the user has inserted into Prism.
- [x] Prepare individual $f_1(\theta)$ and $f_2(\theta)$ curves for both charm-strange and bottom-strange sectors for the manuscript.
- [x] Prepare $f_2(\theta)/f_1(\theta)$ curves for the referee reply only, covering both sectors consistently. The user corrected an interim request for the inverse: show only f2/f1.
- [x] Explain the referee's factorization motivation in the response draft and author notes, distinguishing weak-current constants from rotated-current residues. No full branching-fraction calculation was performed.
- [x] Recompute the residue scan, verify agreement with existing Figure 2 at all angles, and prepare caption/prose clarification. The existing width figure needs no fixed-residue correction.
- [x] Assemble a Referee 1 reply draft with the response-only ratio figure. Item 1's unverified input-sensitivity conclusion is explicitly pending. All proposed changes remain outside final_manuscript.

### 19 October 2026 calculation and implementation check

Follow-up on the supplied Bc code folder: see
[the source-location note](correctness_audit_20261019/bc_folder_followup/README.md).
The Bc perturbative matrix agrees algebraically with the current one up to
a common factor after flavor-role matching. The original Ds/Bs numerical
implementation remains unidentified. The nearby earlier-paper project was
found, but its available computational files are plotting routines and its
raw curves correspond to an older draft. Do not use them as the final baseline.
The user recalls Mathematica and Python but cannot recall the filename.

The figure driver recomputed both residues at all 161 original angles in each sector. All 1,288 comparisons against both state rows in the saved Figure 2 CSV agree to within 1.11e-16 GeV. Preferred-angle scan ratios are 0.4175195 (Ds, 26.6 degrees) and 0.1662805 (Bs, 38.5 degrees). These are fixed-representative-input values, not marginalized ratio predictions; no uncertainty bands are shown. Full provenance and numerical data are in the revision package.

The author suggested that the original angle code should be in scripts/. Searches and code inspection found external published-angle inputs and a separate LO+d3+d5 matrix diagnostic. A limited check using earlier inputs gave midpoint diagnostic angles 30.04965 degrees (Ds) and 39.43370 degrees (Bs). This has **not** reproduced the original study: equivalence of its OPE and numerical prescription has not been established. It does not show that the published angles are incorrect. The user asked for clarification, and this distinction was explained. Do not report these diagnostic numbers in the referee response or use them to claim that the input update is negligible. See [input_check_status.md](referee1_revision/input_check_status.md).

The remaining sections below preserve the earlier discussion and detailed checklist; where an old deferred status appears, this dated update takes precedence.

## Sources and working rules

- Use the author-supplied [r1.pdf](/Users/sbilmis/Downloads/r1.pdf) for the present manuscript discussion. Check for a newer supplied version when resuming.
- Referee 1's report was supplied in the conversation; its requests are summarized below.
- Referee 2's report is [dv14652_report_1_a.pdf](/Users/sbilmis/Downloads/dv14652_report_1_a.pdf).
- The earlier mixing-angle paper is [Phys. Rev. D 111, 014005 (2025)](https://open.metu.edu.tr/bitstream/handle/11511/113265/PhysRevD.111.014005.pdf), DOI 10.1103/PhysRevD.111.014005, arXiv:2408.08014.
- Keep [final_manuscript/](../final_manuscript/) read-only, including build artifacts. The user will implement text themselves unless they explicitly approve edits.
- Preserve existing scripts and outputs. Save later calculations in a dedicated directory outside final_manuscript.
- Do not mark proposed manuscript text as implemented without checking an updated version.

## Established normalization fact

The author checked the earlier calculation and confirmed that it included the tensor-current factor $1/(m_Q+m_s)$, exactly as in the present heavy-strange analysis. The factor was accidentally omitted from the earlier printed definition. Do not reopen this question or apply an angle conversion on that account. See [AGENTS.md](../AGENTS.md).

This clarification does not itself determine the effect of updated numerical inputs or settle Referee 2's separate current-state-overlap question.

## Referee 1 item 1

Request: briefly explain the earlier mixing-angle method and its results, taking account of changes in inputs, particularly the quark condensate.

### Method explanation

- [x] Discussed the explanation with the user.
- [x] Supplied LaTeX describing the normalized currents, diagonalization condition, inverse relations, and cancellation of the difference between physical diagonal correlators.
- [x] Located the insertion in r1.pdf: page 11, after the paragraph following Eqs. (37)-(38) ending with the reference to Appendix C, before "Substituting the residues obtained from the two-point sum rules...".
- [ ] Verify implementation when the user supplies an updated manuscript.

Use hatted two-point amplitudes and the current rotation

$$
R(\theta)=
\begin{pmatrix}\sin\theta&\cos\theta\\
\cos\theta&-\sin\theta\end{pmatrix},
\qquad
\tan2\theta=-\frac{2\widehat\Pi_{AB}}
{\widehat\Pi_{AA}-\widehat\Pi_{BB}}.
$$

For $\Delta=\widehat\Pi_{11}-\widehat\Pi_{22}$, the inverse relations are
$\widehat\Pi_{AA}-\widehat\Pi_{BB}=-\cos(2\theta)\Delta$ and
$2\widehat\Pi_{AB}=\sin(2\theta)\Delta$. Preserve these signs when adapting the separate $B_c$ paper. State diagonalization as the condition used within the adopted two-state method.

### Decisions about the input check

The user expected a small effect because the input uncertainty ranges overlap, the angle depends on a ratio, and perturbative terms dominate in the existing analysis. The conditional calculation now supplies actual shifts, as recorded above. It supports smaller-than-quoted-error shifts on the sampled grid, not an unqualified negligible-effect claim or a verified rerun of the original study.

The manuscript's $\theta=90^\circ$ comparison supports the pure axial-current decay-constant limit. It does not directly test the finite-angle sensitivity involving the $AB$ and $BB$ entries.

The user wants this quantitative comparison in the referee reply only, without adding a numerical study to the manuscript. Discuss the actual magnitude and the reconstruction limitation before adopting response wording or changing inputs/predictions.

### Inputs for the deferred comparison

The earlier values are from Table II of the 2025 paper; present values are from Table I on page 12 of r1.pdf.

| Input | Earlier paper | Present manuscript | Action |
| --- | --- | --- | --- |
| $\langle\bar qq\rangle(1\,\mathrm{GeV})$ | $-(1.65\pm0.15)\times10^{-2}\,\mathrm{GeV}^3$ | $-(0.240\pm0.010)^3\,\mathrm{GeV}^3$ | Compare central values and distinguish uncertainty parametrizations |
| $\kappa_s=\langle\bar ss\rangle/\langle\bar qq\rangle$ | $0.80\pm0.20$ | $0.80\pm0.10$ | Central ratio unchanged; error narrowed |
| $m_0^2$ | $0.80\pm0.20\,\mathrm{GeV}^2$ | $0.80\pm0.20\,\mathrm{GeV}^2$ | Keep unchanged |
| $m_s$ | $0.126\,\mathrm{GeV}$ at $1\,\mathrm{GeV}$ | $0.093\pm0.011\,\mathrm{GeV}$ at $2\,\mathrm{GeV}$ | Compare at a common scale |
| $\overline m_c(\overline m_c)$ | $1.27\pm0.02\,\mathrm{GeV}$ | $1.27\pm0.02\,\mathrm{GeV}$ | Unchanged |
| $\overline m_b(\overline m_b)$ | $4.18^{+0.03}_{-0.02}\,\mathrm{GeV}$ | $4.18\pm0.03\,\mathrm{GeV}$ | Central value unchanged; note error convention |

The central light-condensate magnitude changes from 0.0165 to 0.013824 GeV cubed, a 16.22 percent decrease. The quoted intervals overlap: earlier [0.0150, 0.0180], present [0.012167, 0.015625] GeV cubed. The present interval does not cover the full earlier interval.

Varying both present condensate inputs over their quoted intervals gives a rectangular envelope for $|\langle\bar ss\rangle|$ of [0.0085169, 0.0140625] GeV cubed, containing the earlier central value 0.0132 GeV cubed. This envelope is not a propagated one-standard-deviation interval.

Compare angle shifts with the published uncertainties:
$\theta_{D_s}=26.6^\circ\pm0.6^\circ$ and
$\theta_{B_s}=38.5^\circ\pm0.1^\circ$.

### Deferred calculation checklist

- [ ] Locate the authors' earlier mixing-angle implementation and establish its OPE, input scales, Borel definition, continuum prescription, and angle branch.
- [ ] Establish the earlier-input baseline. Published mixing windows are $M^2=3$-6 GeV squared, $s_0=10\pm1$ GeV squared for $D_s$, and $M^2=9$-13 GeV squared, $s_0=43\pm1$ GeV squared for $B_s$ (earlier Table III). Distinguish these from the present transition and residue windows.
- [ ] Change the condensate normalization alone first. Change the mixed condensate consistently through $\langle\bar s g_s\sigma Gs\rangle=m_0^2\langle\bar ss\rangle$.
- [ ] Compare remaining updated inputs at consistent scales, holding the OPE and auxiliary-parameter prescription fixed. Track changed uncertainty assignments separately from central shifts.
- [ ] Evaluate old and new inputs at identical representative points and working-region boundaries for both sectors. Keep the physical angle branch continuous.
- [ ] Save old/new angles, their differences, and the shift divided by the corresponding quoted angle uncertainty. Report representative shifts and their range.
- [ ] Record signed perturbative and condensate contributions to $AB$ and $AA-BB$. Dominance in each diagonal entry alone does not establish dominance in their difference.
- [ ] Begin with this limited check. Extend the grid or use matched uncertainty sampling only if the findings or the intended uncertainty claim require it. A central-input check alone does not validate a full uncertainty distribution.
- [ ] If the effect is small, prepare a compact numerical table and a reply sentence with the actual shifts. Keep the numerical study in the reply as requested.
- [ ] If shifts are significant or the baseline fails, report the result for discussion. Do not assume negligible sensitivity or silently replace the published angles.
- [ ] Discuss the results before finalizing the item 1 reply.

### Implementation pointers

These files are inspection starting points, not verified replacements for the original mixing-angle calculation:

- [twopoint_ds1_matrix_sumrule.py](../scripts/twopoint_ds1_matrix_sumrule.py) explicitly uses an LO plus dimension-three and dimension-five OPE.
- [twopoint_bs1_matrix_sumrule.py](../scripts/twopoint_bs1_matrix_sumrule.py) continues that OPE to bottom, uses the published angle externally for physical residues, and labels its matrix angle as a diagnostic.
- [mixing_angle_inputs.py](../scripts/mixing_angle_inputs.py) stores the external published angles.
- [Ds1_AA_AB_BB_two_point_sumrule.nb](../notebooks/Ds1_AA_AB_BB_two_point_sumrule.nb) and the Mathematica scripts may help trace the implementation.

The earlier paper lists additional condensate inputs. Verify the terms actually retained in its implementation. Do not attribute an OPE-truncation difference to changed input values. Avoid existing main routines that overwrite established outputs.

Suggested later directory: review/referee1_item1/. Save the new driver, input provenance, pointwise CSV, and a brief summary there.

## Referee 1 item 2

Request: consider transitions to vector $D_s^*$ and $B_s^*$ final states, acknowledge the experimental usefulness of ratios to the pseudoscalar channels, and discuss possible mixing-angle sensitivity, especially for $D_{s1}(2536)$. The referee explicitly recognizes that the calculation may warrant a separate paper.

Latest wording preference (2026-10-24): the separate study is nearing
completion. Keep the follow-up statement in the response; no Conclusion
addition is planned. The older manuscript-addition suggestion below is
superseded by this preference. Continue to describe the suggested ratios as
planned for that study, without claiming they have already been completed.

- [x] Author confirmed on 2026-10-16 that another paper is being prepared for these calculations and will present the ratios following the referee's suggestion.
- [ ] Finalize a courteous reply acknowledging the usefulness of the suggestion and explaining the follow-up work.
- [ ] Discuss a brief qualitative manuscript sentence about the experimental value of ratios and their potential mixing sensitivity. No insertion is yet confirmed.
- [ ] Keep the vector-channel calculation and quantitative ratios assigned to the follow-up paper; no such calculation is requested in this revision session.

Suggested reply for discussion:

> We thank the referee for this useful suggestion. We agree that ratios of radiative widths into pseudoscalar and vector final states are experimentally more accessible than absolute widths, particularly for the bottom-strange mesons, and may provide additional sensitivity to the mixing angle. We are currently preparing a separate study of the radiative transitions into the vector states $D_s^*$ and $B_s^*$. Following the referee's suggestion, that study will also present the ratios of the radiative widths into pseudoscalar and vector final states and examine their dependence on the mixing angle, including the $D_{s1}(2536)$ channel.

Suggested qualitative manuscript addition for discussion:

> Ratios of radiative widths into pseudoscalar and vector final states are experimentally more accessible than absolute widths, particularly for bottom-strange mesons, and may provide complementary information on the mixing angle. The vector-final-state transitions and these ratios will be investigated in a separate study currently in preparation.

Do not claim that a manuscript sentence has been added until the user implements it. Do not claim calculated ratios or demonstrated mixing sensitivity before the follow-up analysis establishes them.

## Referee 1 item 3

Request: present the angle dependence of $f_2/f_1$ and discuss its connection to factorized $B$ decays, with arXiv:2411.03399 as experimental motivation.

- [x] Began discussion of an additional decay-constant and ratio figure; user requested verification of whether Figure 2 had fixed residues.
- [x] Verified that the current Figure 2 builder reprojects both two-point residues at each angle and passes them into the radiative form factors. See [build_figure2_data.py](../figures/figure2_mixing_angle_widths/build_figure2_data.py), especially project_residues and build_rows.
- [x] Confirmed that the saved [Figure 2 CSV](../figures/figure2_mixing_angle_widths/figure2_mixing_angle_widths.csv) contains 161 distinct values of each residue in each sector. For example, charm $f_1,f_2$ are 403.89, 171.95 MeV at 25 degrees and 397.45, 203.26 MeV at 45 degrees. These are representative scan values, not new marginalized results.
- [x] Confirmed that the local final_manuscript, arxiv, and Python figure PDF copies have identical hashes; visually matched that figure to Figure 2 on page 14 of r1.pdf.
- [ ] Propose a caption clarification: underlying QCD and auxiliary inputs are held fixed while $f_1(\theta)$ and $f_2(\theta)$ are recomputed. The builder selects a representative accepted two-point sample rather than simply fixing the residues at the preferred angle.
- [x] User specified the presentation scope: provide $f_1(\theta)$ and $f_2(\theta)$ for both charm-strange and bottom-strange sectors. Do not add a decay-constant ratio plot to the present manuscript.
- [ ] Prepare the individual-constant figure after completing discussion of the calculation details. Proposed layout: charm in the left panel and bottom in the right, with both residues in each panel. No new plot has yet been generated.
- [ ] If showing uncertainty bands, use shared sampled QCD inputs for both residues at each fixed angle and check the sum-rule working criteria over the displayed range.
- [x] Latest user decision: include the $f_2(\theta)/f_1(\theta)$ plot in the referee reply only, not in the manuscript. This supersedes the earlier plan to defer the ratio plot entirely to follow-up work.
- [ ] Generate the response-only ratio figure from the same inputs as the individual-constant curves. If uncertainties are shown, form the ratio sample by sample using shared inputs.
- [ ] Revisit why the referee requests the ratio. The user said this motivation remains unclear. The report explicitly connects it to a factorized $B$-decay branching-fraction ratio; clarify that connection and distinguish the weak-current constants from the rotated-current residues before interpreting the plot.
- [ ] Distinguish the decay-constant ratio in item 3 from the radiative-width ratios promised in item 2; explicitly state the follow-up scope when drafting the replies.
- [ ] Explain the current definitions as needed for the reply. A full numerical factorized branching-ratio application remains outside the agreed present scope.
- [ ] Resume calculation details, figure preparation, proposed modifications, and the reply in the next session; all are currently deferred.

Further Figure 2 verification on 2026-10-16: at 45 degrees the stored $D_{s1}(2536)$ width is 27.4139 keV with $f_2=203.259$ MeV. Freezing $f_2$ to the quoted 168 MeV, with the same transition numerator, would give 40.1286 keV. For $B_{s1}(5830)$ at 25 degrees, the stored width is 4.05817 keV with $f_2=167.095$ MeV; freezing it to the quoted 89 MeV would give 14.3047 keV. The PDF displays the variable-residue curves. These comparisons use the existing CSV and the $1/f_2^2$ width scaling; no new OPE calculation was run. The available committed Figure 2 builder also reprojects at each angle.

Under the manuscript's exclusive-current-coupling assumption, the inverse current rotation gives $J^A=\sin\theta J^{(1)}+\cos\theta J^{(2)}$, hence weak axial-current constants $f_{A,1}=\sin\theta f_1$ and $f_{A,2}=\cos\theta f_2$. This is a conditional algebraic relation, not a validation of state isolation. A factorized $B$-decay application must use the appropriate weak-current matrix elements and include phase space and $B\to D$ form factors. Do not directly identify the rotated-current residue ratio squared with the branching-fraction ratio.

The LHCb reference supplied by the referee is [arXiv:2411.03399](https://arxiv.org/abs/2411.03399), studying $D_{s1}(2460)^+\to D_s^+\pi^+\pi^-$ in reconstructed $B$ decays. Use the paper for the experimental motivation; verify specific branching-fraction claims before quoting them.

## Referee 2 queue

These items await discussion; do not execute them during the deferred Referee 1 input check.

| Report item | Work to discuss |
| --- | --- |
| Main i | State isolation, residue matrix, off-diagonal correlator, possible two-pole treatment. The printed normalization omission is already clarified |
| Main ii | Signed contributions, reliability of the suppressed amplitude, omitted terms relative to the residual, DA shape, twist-four, Borel asymmetry, radiative effects, uncertainty budget |
| Specific 1 | Continuum subtraction, duality region, boundary terms, and impact on the suppressed channel |
| Specific 2 | Photon-normalization evolution, mass/current scales, lattice finite-mass subtraction, DA and perturbative/soft separation |
| Specific 3 | Nonlocal DA definitions, symmetric-point vector-DA contribution, completeness of the selected projection |
| Specific 4 | Sampling distributions, two-point windows/thresholds, shared-input correlations, residue recomputation with angle |
| Specific 5 | Momentum assignments, Fourier phases, and light-cone expansion regime |
| Specific 6 | Primed integration arguments, denominators, and momentum-fraction notation |
| Specific 7 | Pure axial transition-kernel comparison with earlier LCSR using common inputs |

## Prompt for a later Codex session

Use this prompt when resuming the deferred calculation:

~~~text
Resume Referee 1 item 1 for the Ds/Bs radiative-decay paper in
/Users/sbilmis/Dsgamma/papers/Ds_Bs_gamma.

Read AGENTS.md and review/referee_reply_todo.md first. Perform only the
deferred input-parameter sensitivity check. Keep final_manuscript/
read-only and preserve existing outputs. Use the supplied
/Users/sbilmis/Downloads/r1.pdf unless the author supplies a newer version.

Established fact: the calculations underlying Phys. Rev. D 111, 014005
(2025), arXiv:2408.08014, included 1/(m_Q+m_s) in the tensor current.
Its printed omission was a typo. Do not reopen this or convert the angle
on that account.

Locate the original mixing-angle implementation and establish its
earlier-input baseline. The local two-point scripts use a specified
truncated OPE and must not automatically be treated as the original
published-angle calculation.

Compare the old inputs with Table I of r1.pdf, beginning with
<qq>(1 GeV) = -0.0165 GeV^3 versus -(0.240 GeV)^3. Track the old/new
errors on <ss>/<qq>, 0.20/0.10 around 0.80, and change the mixed
condensate consistently. Compare strange masses only with their quoted
scales treated consistently. The TODO records the inputs and original
mixing Borel/threshold windows.

Hold the OPE, normalization, Borel definition, and continuum prescription
fixed in paired evaluations. Compute old/new angles at identical
representative points and working-region boundaries in both sectors.
Report shifts relative to 0.6 degrees for Ds and 0.1 degrees for Bs.
Inspect signed perturbative and condensate pieces of AB and AA-BB to
test whether the ratio and perturbative dominance suppress the effect.
Maintain the same physical angle branch.

Start with a limited check; broaden it only if findings or uncertainty
claims require it. If the old baseline cannot be reproduced, explain
the mismatch rather than presenting a different calculation as an
input-only comparison. Do not assume that the effect is negligible.
The theta=90-degree decay-constant agreement supports the axial limit;
it is not a direct finite-angle sensitivity test.

Save a reproducible driver, input provenance, pointwise results, and
a concise summary under review/referee1_item1/. Report the findings
for discussion before drafting the final reply. The user prefers to
keep a negligible sensitivity result in the referee reply rather than
add a numerical study to the manuscript. If it is significant, discuss
the implications before changing angles or radiative predictions.
Do not proceed to the other referee items as part of this calculation.
~~~
