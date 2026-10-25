# Project instructions and established facts

## Canonical project.org handoff and fresh-session approval — latest 2026-10-25

The author agreed to the proposed investigation and asked for `project.org`
using the personal `tproj` snippet, including the calculation/replotting plan
and the manuscript figure styles. The author explicitly requires a fresh
Codex session to ASK FOR APPROVAL BEFORE STARTING CALCULATIONS.

Start with `project.org`, which is now the canonical task queue, and its one
NEXT task. Briefly summarize the scope, ask for the explicit review-only
execution approval, and wait for the answer. Creating the plan is not approval
to execute it. Do not run derivations, numerical physics checks, scans,
ensembles, or result plots before that approval. Read-only source review and
discussion remain permitted. After approval, record it in the project's Work
Log, release the execution hold in `review/recalculation_a1_20261025/decisions.json`,
and proceed through the technical prerequisites and task sequence without
repeatedly asking for the same permission.

The agreed future central prescription is a=1, sigma_0=s_eff=Table II s0,
eta=m_A^2/m_P^2 and, independently, r=M1^2/M2^2=m_A^2/m_P^2 (M1 axial).
Keep eta fixed during r scans. Photon matching/scales are an agreed
investigation, not a solved matching problem or an adopted replacement input.
This latest entry supersedes conflicting older dated planning/hold wording;
it does not retroactively validate any numerical diagnostic.

All new physics code/data/figures go under `review/recalculation_a1_20261025/`.
Preserve production and `final_manuscript/` as read-only. The project records
the current figures 1--4, extra signed Pi_A/Pi_B/Pi_2 panels, response-only
f2/f1, exact style references, and scientific completion criteria. Current
manuscript/report/reply PDFs and the tproj snippet are copied into the review
baseline's `reference_documents/` with hashes for a portable handoff.
The author subsequently requested a local Git checkpoint of this handoff and
supporting review material. Existing production edits remain outside that
checkpoint. This commit request does not authorize calculations, a push,
cross-machine synchronization or a background job. The execution hold remains.

## Explicit discussion hold after intercept choice, later on 2026-10-25

The author now requests the recommended eta=m_A^2/m_P^2 (typed "n") and asks
about expected result changes, unequal M1/M2 and DA evaluation, and readiness.
The author explicitly says DO NOT START CALCULATIONS BEFORE AGREEMENT on
these points. Only read saved evidence and update discussion/decision notes
in this turn; no new derivations, numerical checks, ensembles or figures.
The intended central Borel ratio is separately r=m_A^2/m_P^2 with M1 axial,
M2 pseudoscalar and u0=1/(1+r). Do not conflate eta with r or present an
argument-only substitution as the completed continuum treatment.
The new review decisions record the intercept choice and the execution hold.
The photon-input, pole-separation and kernel checks remain prerequisites
before any full result/figure regeneration can be called validated.

## Linear-domain agreement and isolated rerun preparation, 2026-10-25

The author explicitly agreed to a=1 and sigma_0=s_0 from Table II (also
called s_eff), and required all recalculations/figures to remain in review,
separate from previous results. The author asked what remains to discuss
before the full rerun. Preparation is in `review/recalculation_a1_20261025/`:
409 baseline hashes, 92 source snapshots, agreed decisions, and a checklist.
No new ensemble or result figure has been generated. The author answered
the intercept-ratio question by requesting an explanation and more discussion;
eta=tilde_s_A/tilde_s_P remains unselected. The recommendation is mass-based,
but do not treat it as agreed. A geometric explanation/illustration is saved
in the new directory's `discussion/`; it is not a physics recalculation.
Do not re-ask about a=1 or the identification sigma_0=s_0, which are settled.
Before final numbers, complete the kernel/convention and boundary checks,
photon finite-ms normalization matching, quantitative current/pole separation,
and the documented residue/shared-input corrections. Read the new README.

## Pullin--Zwicky continuum discussion, later on 2026-10-25

The author chose a general/unequal-Borel analysis for the forthcoming rerun,
but asked to discuss and agree on the JHEP09(2021)023 continuum formulas first.
Read the newest TODO and the final section of
`review/referee_a_continuum_planning_20261025/README.md`. The paper is accessible
and saved locally. Its full boundary prescription is not the uniform soft
factor currently used in production. Linear-domain C.4 checks through n=4
pass; an explicit curved-domain test identifies an apparent missing curvature
term in the printed general-a formula. Do not apply the general-a formula
uncritically or describe arbitrary-domain verification as complete.
The proposed a=1 domain and threshold-intercept assignment remain under
discussion; no new width calculation or production change is adopted.
Previously selected windows were verified; use them as the baseline to
reassess after changes, rather than implying that no selection was done.

## Continuum discussion and planning, later on 2026-10-25

The latest author request is to discuss continuum subtraction and establish
the recalculation roadmap BEFORE further width recalculations. See
`review/referee_a_continuum_planning_20261025/README.md` and the newest TODO.
Only existing perturbative/total ratios and hard spectral-retention diagnostics
were evaluated in that turn. No continuum variant, production patch, new
ensemble or manuscript change was adopted. Signed pert/total is not a pole
fraction; the absolute hard spectral retention is only a perturbative proxy.
The common soft continuum factor is an unvalidated prescription for this
correlator, despite categorical wording in some older local review notes.
Derive the double-spectral domain and boundary terms before selecting the
final subtraction rule. A confirmed material shift must enter the manuscript,
but the provisional 5.18-keV fixed point must not replace the published-input
MC median before those checks. Continue the discussion, not automatic large
recalculation, unless the author requests the next derivation/calculation.

## Unequal-Borel referee diagnostic, 2026-10-25

The author authorized a separate numerical comparison of r=M1^2/M2^2=1
with r=m_A^2/m_P^2, for the referee response only. Its code, data, checks,
table and plot are in `review/referee_a_borel_asymmetry_20261025/`.
Read that README and the newest TODO entry before using the results.
The fixed-input preserved-model scan is complete, but it is not a validated
general unequal-Borel extraction. The suppressed charm midpoint changes
0.182716 -> 5.183620 keV using the stored off-shell common-fraction convention;
these are not the manuscript MC medians. Do not call the effect negligible.
The manuscript/derivation fraction and antisymmetric-DA conventions require
reconciliation; the alternative argument diagnostic must not be substituted
as a preferred physical answer. The hard diagonal spectral and soft continuum
models were retained, not rederived. All production files and the manuscript
remain unchanged. No new uncertainty ensemble or manuscript edits adopted.

## Session handoff — paused at the author's request, 2026-10-19

The author asked to save notes, stop for now, and continue tomorrow
(2026-10-20). Resume only when the author returns; no background calculation
or scheduled action was requested. Start with the handoff at the top of
`review/referee_reply_todo.md` and the mixed-condensate resolution report.
The leading-d5 derivation and isolated checks are complete. The proposed
production patch is unapplied, full ensembles have not been regenerated,
and the item 1 quantitative reply remains provisional. Discuss the recorded
finding and next steps before adopting revised numerical results. Preserve
`final_manuscript/` as read-only and retain the existing production baseline.

On 2026-10-24 the author returned to Referee 1 item 2, asking to say that the
separate vector-final-state study is nearing completion and preferring no
Conclusion addition. Keep this statement in the referee reply and omit any
claim that the Conclusion has been changed. This narrow resumption does not
adopt the pending OPE correction or authorize replacing numerical results.

Later on 2026-10-24 the author supplied new current documents:
`/Users/sbilmis/Downloads/r1-1.pdf` (manuscript) and
`/Users/sbilmis/Downloads/ref1_reply.pdf` (reply), and requested finalization
one item at a time. These supersede r1.pdf and the older reply for placement
and wording review. Begin with Comment 1. The method insertion is now present
on manuscript pp. 11--12, Eqs. (39)--(43), including the normalization
footnote; the angle inputs are Eq. (46), p. 13. The current reply p. 2 still
contains the outdated author-review note and the earlier conditional
0.185/0.057-degree results. A qualified completed-current-OPE alternative is
saved as `review/referee1_revision/comment1_response_proposal_20261024.tex`;
it is proposed wording, not an adopted original-program rerun or submission
approval. The manuscript Appendix C, Eqs. (C2)--(C3), still has the
uncompleted mixed-condensate terms. Production correction/regeneration remains
pending. Also, Figures 1 and 2 have exchanged images: p. 14 has a Borel plot
under the decay-constant caption, and p. 15 has the decay-constant plot under
the Borel caption. Correct the image paths while preserving the intended
caption/label numbering. Do not report this as fixed until a new PDF verifies it.

## Tensor-current normalization in the earlier mixing-angle paper

The user confirmed on 2026-10-16 that the calculations underlying
Phys. Rev. D 111, 014005 (2025), DOI: 10.1103/PhysRevD.111.014005,
actually included the same tensor-current normalization factor
`1/(m_Q + m_s)` used in the present radiative-decay paper, with `Q = c, b`
for the heavy-strange systems. The factor was accidentally omitted from
the printed current definition in the earlier manuscript. The authors
checked their calculations and verified that the factor was present.

Treat this as an established author-confirmed fact in future discussions.
Do not repeatedly ask whether the earlier calculation included the factor,
or infer a numerical normalization mismatch from the printed omission.
No angle conversion is needed solely because of this typographical omission.
When discussing the published definition, distinguish the printed typo from
the correctly normalized calculation. This clarification does not by itself
settle separate questions about current-state overlaps or updated QCD inputs.

## Referee-response workflow

- Latest file for the final response-format review is
  `/Users/sbilmis/Downloads/ref1_reply-1.pdf` (4 pages, 24 October 2026,
  16:47:23), found after the author clarified that they meant the response
  letter. The r1-2.pdf supplied immediately before that clarification is
  the manuscript, not the letter. The author chose the short Comment 3
  text and included the ratio figure in the reply; do not automatically
  restore the longer weak-current discussion. Read the latest format-review
  notes at the top of `review/referee_reply_todo.md`. The letter currently
  calls this reviewer "Referee B"; do not relabel it without checking the
  intended editor designation. No edits to the supplied PDFs were made.
- Latest response choices on 2026-10-24: the author selected Comment 1 wording
  retaining the +0.185/+0.057-degree conditional reconstruction midpoint
  shifts, explicitly described as a reconstruction and central-input check.
  Do not substitute the +0.175/+0.045 alternative without discussion. This
  is not approval of an original-program validation or numerical OPE changes.
  Comment 2 states that the vector-channel study is nearing completion; no
  Conclusion addition. The current discussion has moved to Comment 3; use
  `review/referee1_revision/comment3_response_proposal_20261024.tex`, with
  individual f1/f2 curves in manuscript Figure 1 and the ratio in the reply.
  The older combined local reply remains a historical draft. Read the latest
  dated decisions in the TODO before using it.
- Keep `final_manuscript/` read-only. Do not edit its contents or generate
  build artifacts there unless the user explicitly agrees to and approves
  the proposed changes.
- Work through Referee 1's comments one item at a time: discuss the item
  with the user first, then perform the analysis, then draft the reply for
  that item. Do not advance through all items in one batch.
- Keep any analysis files and proposed manuscript or response text outside
  `final_manuscript/` until the user approves manuscript edits.
- Read `review/referee_reply_todo.md` when resuming referee-response work.
  It records decisions, deferred calculations, the report queue, and a
  prompt for a later Codex session. As of 2026-10-16, the user has deferred
  the remaining Referee 1 calculations, figure preparation, proposed
  modifications, and replies to the next session. Resume when requested.
- For this referee discussion, the user designated
  `/Users/sbilmis/Downloads/r1.pdf` as the manuscript reference. Use its
  notation and insertion locations unless a newer version is supplied;
  do not substitute the earlier local LaTeX snapshot.
- For Referee 1 item 2, the author states that a separate paper is being
  prepared for vector-final-state radiative decays and will include the
  ratios suggested by the referee. Do not treat this statement as a
  claim that those ratios have already been calculated.
- The current Figure 2 builder and saved data already recompute both
  two-point residues at every plotted angle. This was checked against
  the displayed figure in `r1.pdf` on 2026-10-16. See the TODO for evidence;
  do not assume Figure 2 used fixed residues solely from its caption.
- For Referee 1 item 3, the user wants individual f1(theta) and f2(theta)
  plots for both charm-strange and bottom-strange sectors, with no
  decay-constant ratio plot in the present manuscript. The latest decision
  is to include f2(theta)/f1(theta) plots in the referee reply only; this
  supersedes the earlier plan to defer the ratio plot entirely. Clarify
  the motivation involving factorized B decays when resuming. Distinguish
  this ratio from the radiative-width ratios discussed under item 2.
- On 2026-10-19 the user resumed preparation of the figures and Referee 1
  reply. The proposed revision package is in review/referee1_revision/.
  It contains both-sector f1/f2 individual curves, the response-only f2/f1
  ratio, insertion snippets keyed to r1.pdf, and a reply draft. The user
  corrected an interim f1/f2 request: only f2/f1 is wanted as the ratio.
  Prism remains the live manuscript, to be edited by the author.
- The input-sensitivity conclusion for item 1 remains unverified. The local
  scripts and notebook import the published angles externally; their
  truncated-matrix angle diagnostic has not been established as the original
  2025 implementation. An earlier-input diagnostic was recorded in
  review/referee1_revision/input_check_status.md, but its numbers are not
  a reproduction of the original study and do not establish that the
  published angles are incorrect. Keep those diagnostic numbers out of the
  referee reply. The author asked for clarification of this distinction.
- Later on 2026-10-19, the user chose to retain the published mixing angles
  as inputs and return to the referee discussion. Continue with Referee 1
  item by item; defer further original-code searching for now. This working
  choice does not validate the reconstructed diagnostic or establish that
  the effect of updated QCD inputs is negligible. Keep the input-sensitivity
  comparison and the separate shared-input charm MC issue recorded as open;
  do not silently replace numerical results or mark those checks complete.
- The user then explicitly requested the condensate-only quantitative check
  for the referee reply, not the manuscript. It is saved in
  `review/referee1_item1/`: the marked appendix reconstruction gives midpoint
  shifts +0.185164 degrees (Ds) and +0.057096 degrees (Bs), and 21x21-grid
  ranges +0.118680 to +0.368437 and +0.044416 to +0.079634 degrees.
  Python and an independent high-precision Wolfram transcription agree.
  These shifts are below the quoted errors on the sampled grid, but the
  Bs shift is not negligible relative to its narrow error. The original
  implementation remains unidentified. The subsequent leading-d5 derivation
  below supersedes the unresolved-coefficient status, but does not validate
  the printed-appendix reconstruction; do not present it as a verified
  original rerun. Read that directory's README before using its draft.

## Mixed-condensate derivation, later on 2026-10-19

- The user explicitly asked to resolve the coefficient issue. Read
  `review/mixed_condensate_resolution/README.md`. Independent symbolic and
  covariant Dirac/color calculations identify a missing heavy-line
  background-gluon contribution in the current two-point code. With
  C5=<sbar g sigma.G s>, m=mQ, d=mQ+ms, K=M2, its Borel matrix is
  C5 exp(-m^2/K)(1-m^2/K) [[0,1/(12d)],[1/(12d),m/(6d^2)]].
- This resolves the tree-level leading dimension-five coefficient, at ms=0
  in the Wilson coefficient with the current normalization retained. The
  AA addition vanishes, explaining why the pure axial check did not detect
  it. The complete ms=0 AA-BB coefficient is
  -C5 exp(-m^2/K)[m/(3K)+1/(6m)]; neither the existing light-only coefficient
  nor the earlier printed factor-two expression is the complete result.
- All seven symbolic controls and 24 covariant trace checks pass. At the
  four representative preferred-angle residue points, f changes by less
  than 0.2%, and widths by less than 0.35% at fixed transition numerators.
  These are fixed-point comparisons, not updated uncertainties or refits.
- The isolated implementation and an unapplied four-source patch are saved
  in that review directory. Production sources, accepted samples, and
  manuscript files remain unchanged. Discuss before adopting and regenerating
  numerical results; keep the charm sampling issue separate.
- In the completed current LO+d3+d5 OPE the old/new condensate shifts are
  +0.174981 degrees (Ds) and +0.045126 degrees (Bs) at old-window midpoints.
  These are not replacements for the +0.185164/+0.057096-degree values in
  the different, printed-appendix reconstruction. Original-program provenance,
  the larger OPE, and its uncertainty procedure remain separate open issues.
  Do not again describe the leading-d5 calculation itself as unresolved.

## Correctness audit and organization, 2026-10-19

- The user asked to verify the discrepancy and organize inactive work into
  an archive. Read `review/correctness_audit_20261019/REPORT.md` and
  `CURRENT_WORKFLOW.md` when resuming.
- Independent cut-loop traces agree with the current LO perturbative matrix.
  The printed 2025 appendix has an apparent label/sign interchange in its
  leading terms. A marked full-formula reconstruction using that assignment
  gives 26.73607 degrees (Ds) and 38.32237 degrees (Bs) at old-window midpoints.
  This is NOT a verified rerun of the original code or uncertainty analysis.
  The initially observed factor-two mixed-condensate difference was followed
  up in the independent derivation documented above. Do not promote these
  checks into a claim that either published or current full results are
  proven correct or incorrect.
- The preferred charm MC was reproduced exactly. Its transition mc, ms and
  strange condensate are sampled independently of the selected residue's
  inputs, unlike the preferred bottom implementation. A separate paired
  shared-input variant is saved in the audit; it has not replaced scripts,
  outputs, tables or the manuscript. Discuss before adopting revised results.
- Three inactive figure bundles (27 files) were moved intact into
  `archive/inactive_figures_20261019/` with a verified hash/path manifest.
  Active dependencies were retained; older stage-named scripts can still be
  required helpers. Do not move them merely based on their names.
- The user subsequently supplied `/Users/sbilmis/developer/physics/Bc_mixing`.
  Inspection found the Bc heavy-heavy calculation, not a verified copy of
  the original Ds/Bs implementation. Its perturbative densities were checked
  algebraically against the present code and agree up to a common factor 3
  after matching flavor roles; this does not settle light-quark condensates.
- The earlier paper's document project was located at
  `/Users/sbilmis/developer/physics/LateX/papers/mixing_of_D_and_B_mesons`.
  Its `.nb` and `.py` files are plotting routines; its available raw curves
  correspond to an older draft and must not be used as the final published
  baseline. The user recalls Mathematica and Python but does not remember
  the calculation filename. Do not repeatedly ask the same filename question.
  See `review/correctness_audit_20261019/bc_folder_followup/README.md` for
  search coverage, exact paths and the added perturbative cross-check.
