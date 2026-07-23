(* ::Package:: *)

(* Momentum-routing audit for the perturbative hard-photon diagrams in
   Ds1/Bs1 -> Ds/Bs gamma.

   The reference correlator is the ordering used in Colangelo-De Fazio-
   Ozpineci, hep-ph/0505195, Eq. (4.2):

     T_mu(p,q) = i int d^4x exp(i p.x)
       <gamma(q)| T{ J_P^dagger(x) J_X,mu(0) } |0>.

   Thus the pseudoscalar current at x carries the final channel p, while the
   axial current at 0 carries pPrime=p+q.  This script derives both photon-emission
   routings from those two vertex constraints, verifies their fermion-line
   Ward identities, computes the E1-projected triangle cores, and compares
   them with the legacy repository routing.
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, p, q, pPrime, k, l, ms, mQ, p2, pq,
  dH1, dH2, dH3, dS1, dS2, dS3,
  gammaA, gammaB, gammaP, num, reduceTrace,
  standardHeavyTrace, standardStrangeTrace,
  legacyHeavyTrace, legacyStrangeTrace,
  e1Tensor, e1Projector, projectE1, commonMomentumRules
];

scriptDir = DirectoryName[$InputFileName];
paperDir = DirectoryName[scriptDir];
outDir = FileNameJoin[{paperDir, "outputs"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "mathematica_correlator_routing_audit.txt"}];

gammaA[index_] := GA[index] . GA[5];
gammaB[index_] :=
  I FV[pPrime, alpha] DiracSigma[GA[index], GA[alpha]] . GA[5]/(mQ + ms);
gammaP[] := I GA[5];
num[momentum_, mass_] := GS[momentum] + mass;

reduceTrace[expression_] :=
  Collect2[
    DiracSimplify[expression, DiracSubstitute67 -> True],
    {mQ, ms}
  ];

(* Standard correlator routing.

   Heavy emission:
     final vertex x:  k - (k-p) = p,
     initial vertex 0: (k+q) - (k-p) = p+q.

   Strange emission:
     final vertex x:  k - (k-p) = p,
     initial vertex 0: k - (k-p-q) = p+q.
*)

standardHeavyTrace[current_] :=
  reduceTrace[
    DiracTrace[
      gammaP[] . num[k, mQ] . GA[nu] . num[k + q, mQ] .
      current[mu] . num[k - p, ms]
    ]
  ];

standardStrangeTrace[current_] :=
  reduceTrace[
    DiracTrace[
      gammaP[] . num[k, mQ] . current[mu] . num[k - p - q, ms] .
      GA[nu] . num[k - p, ms]
    ]
  ];

(* Legacy routing used by step2_hard_photon_numerators.wl. *)
legacyHeavyTrace[current_] :=
  reduceTrace[
    DiracTrace[
      current[mu] . num[k + q, mQ] . GA[nu] . num[k, mQ] .
      gammaP[] . num[k - p, ms]
    ]
  ];

legacyStrangeTrace[current_] :=
  reduceTrace[
    DiracTrace[
      current[mu] . num[k, mQ] . gammaP[] . num[k - p, ms] .
      GA[nu] . num[k - p + q, ms]
    ]
  ];

standardAHeavy = standardHeavyTrace[gammaA];
standardAStrange = standardStrangeTrace[gammaA];
standardBHeavy = standardHeavyTrace[gammaB];
standardBStrange = standardStrangeTrace[gammaB];

legacyAHeavy = legacyHeavyTrace[gammaA];
legacyAStrange = legacyStrangeTrace[gammaA];
legacyBHeavy = legacyHeavyTrace[gammaB];
legacyBStrange = legacyStrangeTrace[gammaB];

(* Fermion-line Ward identities at numerator level.  For B=A+q,

     (slash A+m) slash q (slash B+m)
       = D_B (slash A+m) - D_A (slash B+m).

   The identities below check this independently for the heavy line and for
   the strange line with A=l-q, B=l. *)

wardHeavyLHS =
  DiracGammaExpand[num[k, mQ] . GS[q] . num[k + q, mQ]];
wardHeavyRHS =
  SP[k + q, k + q] - mQ^2;
wardHeavyRHS =
  wardHeavyRHS num[k, mQ] - (SP[k, k] - mQ^2) num[k + q, mQ];
wardHeavyResidual =
  FCE[
    DiracSimplify[
      DiracGammaExpand[FCI[wardHeavyLHS - wardHeavyRHS]],
      DiracSubstitute67 -> True,
      DiracOrder -> True
    ]
  ] // ScalarProductExpand // Simplify;

wardStrangeLHS =
  DiracGammaExpand[num[l - q, ms] . GS[q] . num[l, ms]];
wardStrangeRHS =
  (SP[l, l] - ms^2) num[l - q, ms] -
  (SP[l - q, l - q] - ms^2) num[l, ms];
wardStrangeResidual =
  FCE[
    DiracSimplify[
      DiracGammaExpand[FCI[wardStrangeLHS - wardStrangeRHS]],
      DiracSubstitute67 -> True,
      DiracOrder -> True
    ]
  ] // ScalarProductExpand // Simplify;

e1Tensor = FV[p, nu] FV[q, mu] - SP[p, q] MT[mu, nu];
e1Projector = e1Tensor/(2 SP[p, q]^2);

commonMomentumRules = {
  SP[q, q] -> 0,
  Pair[Momentum[q], Momentum[q]] -> 0,
  Pair[Momentum[k], Momentum[pPrime]] ->
    Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
  Pair[Momentum[p], Momentum[pPrime]] ->
    Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
  Pair[Momentum[pPrime], Momentum[q]] -> Pair[Momentum[p], Momentum[q]]
};

projectE1[expression_] :=
  Collect2[
    ScalarProductExpand[Contract[expression e1Projector] /. commonMomentumRules],
    {mQ, ms}
  ];

stdProjAHeavy = projectE1[standardAHeavy];
stdProjAStrange = projectE1[standardAStrange];
stdProjBHeavy = projectE1[standardBHeavy];
stdProjBStrange = projectE1[standardBStrange];

legacyProjAHeavy = projectE1[legacyAHeavy];
legacyProjAStrange = projectE1[legacyAStrange];
legacyProjBHeavy = projectE1[legacyBHeavy];
legacyProjBStrange = projectE1[legacyBStrange];

heavyRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dH2 + mQ^2,
  Pair[Momentum[k], Momentum[q]] -> (dH1 - dH2)/2,
  Pair[Momentum[k], Momentum[p]] ->
    (dH2 + mQ^2 + p2 - ms^2 - dH3)/2
};

(* Standard strange routing uses dS3=(k-p-q)^2-ms^2. *)
standardStrangeRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dS1 + mQ^2,
  Pair[Momentum[k], Momentum[p]] ->
    (dS1 + mQ^2 + p2 - ms^2 - dS2)/2,
  Pair[Momentum[k], Momentum[q]] -> pq + (dS2 - dS3)/2
};

(* Legacy strange routing uses dS3=(k-p+q)^2-ms^2. *)
legacyStrangeRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dS1 + mQ^2,
  Pair[Momentum[k], Momentum[p]] ->
    (dS1 + mQ^2 + p2 - ms^2 - dS2)/2,
  Pair[Momentum[k], Momentum[q]] -> pq + (dS3 - dS2)/2
};

triangleCore[projected_, rules_, denominators_] :=
  Collect2[
    Expand[projected /. rules] /. Thread[denominators -> 0],
    {mQ, ms, p2, pq}
  ];

stdCoreAHeavy = triangleCore[stdProjAHeavy, heavyRules, {dH1, dH2, dH3}];
stdCoreAStrange = triangleCore[stdProjAStrange, standardStrangeRules, {dS1, dS2, dS3}];
stdCoreBHeavy = triangleCore[stdProjBHeavy, heavyRules, {dH1, dH2, dH3}];
stdCoreBStrange = triangleCore[stdProjBStrange, standardStrangeRules, {dS1, dS2, dS3}];

legacyCoreAHeavy = triangleCore[legacyProjAHeavy, heavyRules, {dH1, dH2, dH3}];
legacyCoreAStrange = triangleCore[legacyProjAStrange, legacyStrangeRules, {dS1, dS2, dS3}];
legacyCoreBHeavy = triangleCore[legacyProjBHeavy, heavyRules, {dH1, dH2, dH3}];
legacyCoreBStrange = triangleCore[legacyProjBStrange, legacyStrangeRules, {dS1, dS2, dS3}];

coreDifferences = {
  "A heavy" -> Together[stdCoreAHeavy - legacyCoreAHeavy] // Simplify,
  "A strange" -> Together[stdCoreAStrange - legacyCoreAStrange] // Simplify,
  "B heavy" -> Together[stdCoreBHeavy - legacyCoreBHeavy] // Simplify,
  "B strange" -> Together[stdCoreBStrange - legacyCoreBStrange] // Simplify
};

routingChecks = {
  {"heavy final vertex momentum", (k) - (k - p), p},
  {"heavy initial vertex momentum", (k + q) - (k - p), p + q},
  {"strange final vertex momentum", k - (k - p), p},
  {"strange initial vertex momentum", k - (k - p - q), p + q}
};

legacyRoutingChecks = {
  {"legacy heavy final/initial", {k - (k - p), (k + q) - (k - p)}, {p, p + q}},
  {"legacy strange final/initial", {k - (k - p), k - (k - p + q)}, {p, p + q}}
};

standardRoutingPassed = And @@ (TrueQ[Simplify[#[[2]] - #[[3]]] === 0] & /@ routingChecks);
wardPassed = TrueQ[wardHeavyResidual === 0] && TrueQ[wardStrangeResidual === 0];

stream = OpenWrite[outFile];
writeLine[text_] := WriteString[stream, text, "\n"];
writeExpr[label_, expression_] :=
  WriteString[stream, label, "\n", ToString[expression, InputForm], "\n\n"];

writeLine["Correlator and momentum-routing audit"];
writeLine["====================================="];
writeLine[""];
writeLine["Reference ordering: Rohrwild three-point/background-field form with pPrime=p+q."];
writeLine["This is the ordering of hep-ph/0505195, Eq. (4.2)."];
writeLine[""];
writeLine["STANDARD VERTEX MOMENTA"];
Do[
  writeExpr[check[[1]] <> ":", Simplify[check[[2]]]],
  {check, routingChecks}
];
writeLine["Standard routing status: " <> If[standardRoutingPassed, "PASS", "FAIL"]];
writeLine[""];
writeLine["LEGACY VERTEX MOMENTA"];
Do[
  writeExpr[check[[1]] <> ": actual then required", {Simplify[check[[2]]], check[[3]]}],
  {check, legacyRoutingChecks}
];
writeLine[""];
writeExpr["Heavy-line Ward residual:", wardHeavyResidual];
writeExpr["Strange-line Ward residual:", wardStrangeResidual];
writeLine["Fermion-line Ward status: " <> If[wardPassed, "PASS", "FAIL"]];
writeLine[""];

writeLine["STANDARD E1 TRIANGLE CORES"];
writeExpr["A heavy:", stdCoreAHeavy];
writeExpr["A strange:", stdCoreAStrange];
writeExpr["B heavy:", stdCoreBHeavy];
writeExpr["B strange:", stdCoreBStrange];

writeLine["LEGACY E1 TRIANGLE CORES"];
writeExpr["A heavy:", legacyCoreAHeavy];
writeExpr["A strange:", legacyCoreAStrange];
writeExpr["B heavy:", legacyCoreBHeavy];
writeExpr["B strange:", legacyCoreBStrange];

writeLine["STANDARD MINUS LEGACY"];
Do[writeExpr[ToString[item[[1]]] <> ":", item[[2]]], {item, coreDifferences}];

writeLine["CONCLUSION"];
writeLine["The standard correlator routing is momentum-consistent at both vertices."];
writeLine["The legacy strange-line routing k-p+q carries p-q, not p+q, at the axial vertex."];
writeLine["Therefore its strange-line triangle core must not be used as an independent validation of the standard correlator."];

Close[stream];

Print["Wrote routing audit to: ", outFile];
Print["Standard routing: ", If[standardRoutingPassed, "PASS", "FAIL"]];
Print["Fermion-line Ward identities: ", If[wardPassed, "PASS", "FAIL"]];
Print["Core differences (standard - legacy): ", coreDifferences];

If[! standardRoutingPassed || ! wardPassed, Exit[1]];
Quit[];
