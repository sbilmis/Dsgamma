(* ::Package:: *)

(* Independent central Mathematica check of the direct normalized-current
   AA/AB/BB two-point QCD sum rule for the Bs1 system.

   OPE: exact-mass LO perturbative + local d=3 through O(ms^2)
   + local d=5 mixed condensate.  The local terms belong to this ordinary
   two-point SVZ sum rule, not to the external-photon transition LCSR. *)

ClearAll["Global`*"];

scriptDir = DirectoryName[$InputFileName];
paperDir = DirectoryName[scriptDir];
outputDir = FileNameJoin[{paperDir, "outputs"}];
If[! DirectoryQ[outputDir], CreateDirectory[outputDir]];

mb = 4.18;
ms = 0.093;
qq = -(0.24)^3;
kappaS = 0.8;
ss = kappaS qq;
m0sq = 0.8;
mixedSS = m0sq ss;
m1 = 5.750;
m2 = 5.82870;
thetaExternal = 38.5 Degree;
sThreshold = (mb + ms)^2;

kallen[s_] := (s - (mb + ms)^2) (s - (mb - ms)^2);

rhoAA[s_] := Sqrt[kallen[s]]/(8 Pi^2) (
  2 - (mb^2 + ms^2 + 6 mb ms)/s - (mb^2 - ms^2)^2/s^2
);

rhoAB[s_] := Sqrt[kallen[s]]/(8 Pi^2) (
  3 (ms - mb) ((mb + ms)^2 - s)/((mb + ms) s)
);

rhoBB[s_] := Sqrt[kallen[s]]/(8 Pi^2) (
  -((mb + ms)^2 - s) (2 ms^2 - 4 mb ms + 2 mb^2 + s)/
    ((mb + ms)^2 s)
);

rhoMatrix[s_] := {{rhoAA[s], rhoAB[s]}, {rhoAB[s], rhoBB[s]}};

localPieces[M2_] := Module[{den = mb + ms, e, q0, q1, q2, d5},
  e = Exp[-mb^2/M2];
  q0 = ss {{mb, mb^2/den}, {mb^2/den, mb^3/den^2}};
  q1 = ss {{
      -ms mb^2/(2 M2),
      ms mb (1 - mb^2/M2)/(2 den)
    }, {
      ms mb (1 - mb^2/M2)/(2 den),
      ms mb^2 (1 - mb^2/(2 M2))/den^2
    }};
  q2 = ss {{
      ms^2 mb^3/(2 M2^2),
      ms^2 (mb^4/M2^2 - mb^2/M2 - 1)/(2 den)
    }, {
      ms^2 (mb^4/M2^2 - mb^2/M2 - 1)/(2 den),
      ms^2 (mb^5/(2 M2^2) - mb^3/M2)/den^2
    }};
  d5 = mixedSS {{
      -mb^3/(4 M2^2),
      -(mb^4/M2^2 - mb^2/M2 - 1)/(4 den)
    }, {
      -(mb^4/M2^2 - mb^2/M2 - 1)/(4 den),
      -(mb^5/(4 M2^2) - mb^3/(2 M2))/den^2
    }};
  <|"d3_ms0" -> e q0, "d3_ms1" -> e q1,
    "d3_ms2" -> e q2, "d5_mixed" -> e d5|>
];

localMatrix[M2_] := Total[Values[localPieces[M2]]];

buildPerturbativeInterpolations[M2_?NumericQ] := Module[
  {sGrid, integrand, increments, cumulative},
  sGrid = N@Subdivide[sThreshold, 80.0, 5000];
  integrand = (rhoMatrix[#] Exp[-#/M2]) & /@ sGrid;
  increments = MapThread[0.5 (#1 + #2) #3 &,
    {Most[integrand], Rest[integrand], Differences[sGrid]}];
  cumulative = Prepend[Accumulate[increments], ConstantArray[0., {2, 2}]];
  Table[
    Interpolation[Transpose[{sGrid, cumulative[[All, i, j]]}],
      InterpolationOrder -> 1],
    {i, 1, 2}, {j, 1, 2}
  ]
];

perturbativeInterpolations[M2_?NumericQ] :=
  perturbativeInterpolations[M2] = buildPerturbativeInterpolations[M2];

perturbativeMatrix[M2_?NumericQ, s0_?NumericQ] :=
  Map[#[s0] &, perturbativeInterpolations[M2], {2}];

opeMatrix[M2_?NumericQ, s0_?NumericQ] :=
  perturbativeMatrix[M2, s0] + localMatrix[M2];

rotationMatrix[theta_?NumericQ] :=
  {{Sin[theta], Cos[theta]}, {Cos[theta], -Sin[theta]}};

mixingAngle[matrix_?MatrixQ] := Module[{aa, ab, bb, theta},
  {aa, ab, bb} = {matrix[[1, 1]], matrix[[1, 2]], matrix[[2, 2]]};
  theta = ArcTan[aa - bb, -2 ab]/2;
  Mod[theta, Pi/2]
];

projectedValue[matrix_?MatrixQ, theta_?NumericQ, channel_Integer] :=
  Module[{v = rotationMatrix[theta][[channel]]}, v . matrix . v];

projectedOPE[M2_?NumericQ, s0_?NumericQ, theta_?NumericQ,
    channel_Integer] :=
  projectedValue[opeMatrix[M2, s0], theta, channel];

projectedEffectiveMass[M2_?NumericQ, s0_?NumericQ, theta_?NumericQ,
    channel_Integer, deltaTau_: 0.0002] := Module[
  {tau, plus, minus, massSq},
  tau = 1/M2;
  plus = projectedOPE[1/(tau + deltaTau), s0, theta, channel];
  minus = projectedOPE[1/(tau - deltaTau), s0, theta, channel];
  If[plus <= 0 || minus <= 0, Return[Indeterminate]];
  massSq = -(Log[plus] - Log[minus])/(2 deltaTau);
  If[massSq > 0, Sqrt[massSq], Indeterminate]
];

fittedThreshold[M2_?NumericQ, theta_?NumericQ, channel_Integer,
    targetMass_?NumericQ] := Module[{candidates, pairs},
  candidates = N[Range[40.0, 50.0, 0.025]];
  pairs = Select[
    {#, projectedEffectiveMass[M2, #, theta, channel]} & /@ candidates,
    NumericQ[#[[2]]] &
  ];
  First@MinimalBy[pairs, Abs[#[[2]] - targetMass] &]
];

centralPoint[M2_?NumericQ] := Module[
  {fit1, fit2, pi1, pi2, f1, f2, pert1, pert2, full1, full2,
   pieces1, pieces2, pole1, pole2, d51, d52, angleSamples,
   thetaDiagnostic, rotated, offdiag},
  fit1 = fittedThreshold[M2, thetaExternal, 1, m1];
  fit2 = fittedThreshold[M2, thetaExternal, 2, m2];
  pi1 = projectedOPE[M2, fit1[[1]], thetaExternal, 1];
  pi2 = projectedOPE[M2, fit2[[1]], thetaExternal, 2];
  f1 = Sqrt[pi1 Exp[m1^2/M2]/m1^2];
  f2 = Sqrt[pi2 Exp[m2^2/M2]/m2^2];
  pert1 = projectedValue[
    perturbativeMatrix[M2, fit1[[1]]], thetaExternal, 1];
  pert2 = projectedValue[
    perturbativeMatrix[M2, fit2[[1]]], thetaExternal, 2];
  full1 = projectedValue[
    perturbativeMatrix[M2, 80.0], thetaExternal, 1];
  full2 = projectedValue[
    perturbativeMatrix[M2, 80.0], thetaExternal, 2];
  pieces1 = localPieces[M2];
  pieces2 = localPieces[M2];
  pole1 = pert1/full1;
  pole2 = pert2/full2;
  d51 = Abs[projectedValue[pieces1["d5_mixed"], thetaExternal, 1]/pi1];
  d52 = Abs[projectedValue[pieces2["d5_mixed"], thetaExternal, 2]/pi2];
  angleSamples = mixingAngle[opeMatrix[M2, #]] & /@
    N[Range[42.0, 48.0, 0.25]];
  thetaDiagnostic = Median[angleSamples];
  rotated = rotationMatrix[thetaExternal] . opeMatrix[M2, 45.0] .
    Transpose[rotationMatrix[thetaExternal]];
  offdiag = Abs[rotated[[1, 2]]]/
    Sqrt[Abs[rotated[[1, 1]] rotated[[2, 2]]]];
  <|
    "M2_GeV2" -> M2,
    "theta_deg" -> thetaExternal 180/Pi,
    "theta_matrix_median_deg" -> thetaDiagnostic 180/Pi,
    "Pi12_normalized_at_s0mix45" -> offdiag,
    "s01_GeV2" -> fit1[[1]],
    "s02_GeV2" -> fit2[[1]],
    "mEff1_GeV" -> fit1[[2]],
    "mEff2_GeV" -> fit2[[2]],
    "f1_GeV" -> f1,
    "f2_GeV" -> f2,
    "pole_fraction1" -> pole1,
    "pole_fraction2" -> pole2,
    "d5_fraction1" -> d51,
    "d5_fraction2" -> d52
  |>
];

Print["Starting direct Bs central two-point extraction"];
result = centralPoint[7.0];
Print["Direct Bs central extraction complete"];

pythonReference = <|
  "theta_deg" -> 38.5,
  "theta_matrix_median_deg" -> 39.4966226,
  "Pi12_normalized_at_s0mix45" -> 0.103551881,
  "s01_GeV2" -> 44.35,
  "s02_GeV2" -> 43.025,
  "mEff1_GeV" -> 5.75019813,
  "mEff2_GeV" -> 5.82910879,
  "f1_GeV" -> 0.538110737,
  "f2_GeV" -> 0.0888786545,
  "pole_fraction1" -> 0.804771326,
  "pole_fraction2" -> 0.556171437,
  "d5_fraction1" -> 0.00294528619,
  "d5_fraction2" -> 0.0347572524
|>;

crossCheck = AssociationMap[
  result[#] - pythonReference[#] &,
  Keys[pythonReference]
];
regressionTolerance = 1.*^-6;
regressionStatus =
  If[Max[Abs[Values[crossCheck]]] <= regressionTolerance, "PASS", "FAIL"];

summaryLines = {
  "Mathematica direct Bs1 AA/AB/BB two-point cross-check",
  "======================================================",
  "STATUS=" <> regressionStatus,
  "Absolute regression tolerance: 1e-6",
  "OPE: exact-mass LO perturbative + local d=3 through ms^2 + local d=5.",
  "theta_Bs=38.5 deg is external; no external f1, f2, fA, fB, or overlap is used.",
  "",
  "Central point M2=7.0 GeV^2:",
  ToString[result, InputForm],
  "",
  "Mathematica minus Python:",
  ToString[crossCheck, InputForm],
  ""
};

summaryPath = FileNameJoin[
  {outputDir, "mathematica_twopoint_bs1_matrix_check.txt"}];
stream = OpenWrite[summaryPath];
WriteString[stream, StringRiffle[summaryLines, "\n"], "\n"];
Close[stream];

Print[StringRiffle[summaryLines, "\n"]];
If[regressionStatus =!= "PASS", Quit[1]];
