(* Read-only use of the author-designated Bc package. All output is local.
   This checks only its two-mass perturbative matrix, not its heavy-heavy
   condensate expansion or the original heavy-strange angle calculation. *)
Get["/Users/sbilmis/developer/physics/Bc_mixing/BcMixingMomentum.wl"];

auditDirectory = "/Users/sbilmis/Dsgamma/papers/Ds_Bs_gamma/review/correctness_audit_20261019/bc_folder_followup";
Clear[mh, ml, z];
auditLambda = (z - (mh + ml)^2) (z - (mh - ml)^2);
auditCommon = Sqrt[auditLambda]/(8 Pi^2);
auditReference = <|
  "AA" -> auditCommon (2 - (mh^2 + ml^2 + 6 mh ml)/z - (mh^2 - ml^2)^2/z^2),
  "AB" -> auditCommon 3 (ml - mh) ((mh + ml)^2 - z)/((mh + ml) z),
  "BB" -> -auditCommon ((mh + ml)^2 - z) (2 ml^2 - 4 ml mh + 2 mh^2 + z)/((mh + ml)^2 z)
|>;

(* bbar c -> sbar Q: mb becomes the antiquark mass ms; mc becomes mQ. *)
auditRatios = AssociationMap[
  FullSimplify[
    (PerturbativeSpectralDensity[#, z] /. {mb -> ml, mc -> mh})/auditReference[#],
    Assumptions -> mh > ml > 0 && z > (mh + ml)^2
  ] &, {"AA", "AB", "BB"}];
auditPass = And @@ (# === 3 & /@ Values[auditRatios]);
auditResult = <|
  "source" -> "/Users/sbilmis/developer/physics/Bc_mixing/BcMixingMomentum.wl",
  "mapping" -> "mb -> ms (antiquark), mc -> mQ (quark)",
  "Bc_density_divided_by_current_spin_averaged_density" -> auditRatios,
  "pass" -> auditPass,
  "interpretation" -> "All three perturbative entries agree up to one common factor 3, which cancels in tan(2 theta). No mass-dependent tensor-current rescaling is applied.",
  "limitation" -> "This does not supply strange-quark or mixed-condensate terms, or reproduce the original Ds/Bs numerical implementation."
|>;
Export[FileNameJoin[{auditDirectory, "bc_perturbative_mapping.json"}], auditResult, "RawJSON"];
Print[InputForm[auditResult]];
If[!auditPass, Exit[1]];
