(* ::Package:: *)

(* Numerical Mathematica cross-check of the normalized-current AA/AB/BB
   two-point QCD sum rule implemented in twopoint_ds1_matrix_sumrule.py.

   OPE truncation: exact-mass LO perturbative + local d=3 through O(ms^2)
   + local d=5 mixed condensate.  These local terms belong to the ordinary
   two-point SVZ sum rule, not to the external-photon transition LCSR. *)

ClearAll["Global`*"];

scriptDir = DirectoryName[$InputFileName];
paperDir = DirectoryName[scriptDir];
outputDir = FileNameJoin[{paperDir, "outputs"}];
If[! DirectoryQ[outputDir], CreateDirectory[outputDir]];

mc = 1.27;
ms = 0.093;
qq = -(0.24)^3;
kappaS = 0.8;
ss = kappaS qq;
m0sq = 0.8;
mixedSS = m0sq ss;
m1 = 2.4595;
m2 = 2.53511;
sThreshold = (mc + ms)^2;

kallen[s_] := (s - (mc + ms)^2) (s - (mc - ms)^2);

rhoAA[s_] := Sqrt[kallen[s]]/(8 Pi^2) (
  2 - (mc^2 + ms^2 + 6 mc ms)/s - (mc^2 - ms^2)^2/s^2
);

rhoAB[s_] := Sqrt[kallen[s]]/(8 Pi^2) (
  3 (ms - mc) ((mc + ms)^2 - s)/((mc + ms) s)
);

rhoBB[s_] := Sqrt[kallen[s]]/(8 Pi^2) (
  -((mc + ms)^2 - s) (2 ms^2 - 4 mc ms + 2 mc^2 + s)/
    ((mc + ms)^2 s)
);

rhoMatrix[s_] := {{rhoAA[s], rhoAB[s]}, {rhoAB[s], rhoBB[s]}};

rhoAAColangelo[s_] := Sqrt[1 - 2 (mc^2 + ms^2)/s +
  (mc^2 - ms^2)^2/s^2]/(8 Pi^2) (
  2 s - mc^2 - ms^2 - 6 mc ms - (mc^2 - ms^2)^2/s
);

rhoAASymbolic[z_, mh_, ml_] :=
  Sqrt[(z - (mh + ml)^2) (z - (mh - ml)^2)]/(8 Pi^2) (
    2 - (mh^2 + ml^2 + 6 mh ml)/z - (mh^2 - ml^2)^2/z^2);
rhoAAColangeloSymbolic[z_, mh_, ml_] :=
  Sqrt[1 - 2 (mh^2 + ml^2)/z + (mh^2 - ml^2)^2/z^2]/(8 Pi^2) (
    2 z - mh^2 - ml^2 - 6 mh ml - (mh^2 - ml^2)^2/z);

analyticChecks = <|
  "rhoAA divided by Colangelo" ->
    FullSimplify[
      rhoAASymbolic[z, mh, ml]/rhoAAColangeloSymbolic[z, mh, ml],
      Assumptions -> z > (mh + ml)^2 && mh > ml > 0],
  "rhoAB symmetry" -> FullSimplify[rhoMatrix[s][[1, 2]] - rhoMatrix[s][[2, 1]]]
|>;
Print["[1/4] perturbative analytic checks complete"];

localPieces[M2_] := Module[{den = mc + ms, e, q0, q1, q2, d5},
  e = Exp[-mc^2/M2];
  q0 = ss {{mc, mc^2/den}, {mc^2/den, mc^3/den^2}};
  q1 = ss {{
      -ms mc^2/(2 M2),
      ms mc (1 - mc^2/M2)/(2 den)
    }, {
      ms mc (1 - mc^2/M2)/(2 den),
      ms mc^2 (1 - mc^2/(2 M2))/den^2
    }};
  q2 = ss {{
      ms^2 mc^3/(2 M2^2),
      ms^2 (mc^4/M2^2 - mc^2/M2 - 1)/(2 den)
    }, {
      ms^2 (mc^4/M2^2 - mc^2/M2 - 1)/(2 den),
      ms^2 (mc^5/(2 M2^2) - mc^3/M2)/den^2
    }};
  d5 = mixedSS {{
      -mc^3/(4 M2^2),
      -(mc^4/M2^2 - mc^2/M2 - 1)/(4 den)
    }, {
      -(mc^4/M2^2 - mc^2/M2 - 1)/(4 den),
      -(mc^5/(4 M2^2) - mc^3/(2 M2))/den^2
    }};
  <|"d3_ms0" -> e q0, "d3_ms1" -> e q1,
    "d3_ms2" -> e q2, "d5_mixed" -> e d5|>
];

localMatrix[M2_] := Total[Values[localPieces[M2]]];

aaLocalColangelo[M2_] := Exp[-mc^2/M2] (
  ss (mc - mc^2 ms/(2 M2) + mc^3 ms^2/(2 M2^2))
  - mixedSS mc^3/(4 M2^2)
);

analyticChecks["AA local minus Colangelo"] =
  Chop@FullSimplify[localMatrix[M2][[1, 1]] - aaLocalColangelo[M2]];
Print["[2/4] local OPE analytic check complete"];

buildPerturbativeInterpolations[M2_?NumericQ] := Module[
  {sGrid, integrand, increments, cumulative},
  sGrid = N@Subdivide[sThreshold, 11.0, 2000];
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

projectedValue[matrix_?MatrixQ, theta_?NumericQ, channel_Integer] := Module[{v},
  v = rotationMatrix[theta][[channel]];
  v . matrix . v
];

projectedOPE[M2_?NumericQ, s0_?NumericQ, theta_?NumericQ, channel_Integer] :=
  projectedValue[opeMatrix[M2, s0], theta, channel];

projectedEffectiveMass[M2_?NumericQ, s0_?NumericQ, theta_?NumericQ,
    channel_Integer, deltaTau_: 0.0002] := Module[{tau, plus, minus, massSq},
  tau = 1/M2;
  plus = projectedOPE[1/(tau + deltaTau), s0, theta, channel];
  minus = projectedOPE[1/(tau - deltaTau), s0, theta, channel];
  If[plus <= 0 || minus <= 0, Return[Indeterminate]];
  massSq = -(Log[plus] - Log[minus])/(2 deltaTau);
  If[massSq > 0, Sqrt[massSq], Indeterminate]
];

fittedThreshold[M2_?NumericQ, theta_?NumericQ, channel_Integer,
    targetMass_?NumericQ] := Module[{candidates, pairs, best},
  candidates = N[Range[6.55, 10.50, 0.025]];
  pairs = Select[
    {#, projectedEffectiveMass[M2, #, theta, channel]} & /@ candidates,
    NumericQ[#[[2]]] &
  ];
  best = First@MinimalBy[pairs, Abs[#[[2]] - targetMass] &];
  best
];

physicalPoint[M2_?NumericQ] := Module[
  {angleSamples, thetaDiagnostic, theta, fit1, fit2, pi1, pi2, f1, f2,
   rotatedCommon, offdiagNorm},
  angleSamples = mixingAngle[opeMatrix[M2, #]] & /@ N[Range[9.0, 11.0, 0.1]];
  thetaDiagnostic = Median[angleSamples];
  theta = 26.6 Degree;
  fit1 = fittedThreshold[M2, theta, 1, m1];
  fit2 = fittedThreshold[M2, theta, 2, m2];
  pi1 = projectedOPE[M2, fit1[[1]], theta, 1];
  pi2 = projectedOPE[M2, fit2[[1]], theta, 2];
  f1 = Sqrt[pi1 Exp[m1^2/M2]/m1^2];
  f2 = Sqrt[pi2 Exp[m2^2/M2]/m2^2];
  rotatedCommon = rotationMatrix[theta] . opeMatrix[M2, 10.0] .
    Transpose[rotationMatrix[theta]];
  offdiagNorm = Abs[rotatedCommon[[1, 2]]]/
    Sqrt[Abs[rotatedCommon[[1, 1]] rotatedCommon[[2, 2]]]];
  <|"M2_GeV2" -> M2, "theta_deg" -> theta 180/Pi,
    "theta_matrix_diagnostic_deg" -> thetaDiagnostic 180/Pi,
    "Pi12_normalized" -> offdiagNorm,
    "s01_GeV2" -> fit1[[1]], "s02_GeV2" -> fit2[[1]],
    "mEff1_GeV" -> fit1[[2]], "mEff2_GeV" -> fit2[[2]],
    "f1_GeV" -> f1, "f2_GeV" -> f2|>
];

(* One central point is sufficient for an independent numerical regression.
   The dense Python scan supplies the Borel and input uncertainty distribution;
   the interactive notebook below exposes the full Mathematica scan cell. *)
scanM2 = {2.35};
Print["[3/4] starting central physical-pole extraction"];
scanRows = physicalPoint /@ scanM2;
centralPoint = First[scanRows];
Print["[4/4] central physical-pole extraction complete"];

(* Regression reference: the M2=2.35 row written by
   outputs/twopoint_ds1_physical_residue_grid.csv.  Values are copied here
   explicitly because CSV Import hangs in the available headless kernel. *)
pythonCentral = <|
  "theta_deg" -> 26.6,
  "theta_matrix_diagnostic_deg" -> 30.58596158596751,
  "Pi12_normalized" -> 0.18100724350506212,
  "s01_fitted_GeV2" -> 7.925,
  "s02_fitted_GeV2" -> 9.975,
  "f1_GeV" -> 0.40058930523172204,
  "f2_GeV" -> 0.1666703674111434
|>;
Print["[5/7] Python regression row loaded"];

crossCheck = <|
  "theta_deg_difference_Mathematica_minus_Python" ->
    centralPoint["theta_deg"] - pythonCentral["theta_deg"],
  "theta_matrix_diagnostic_difference" ->
    centralPoint["theta_matrix_diagnostic_deg"] -
      pythonCentral["theta_matrix_diagnostic_deg"],
  "Pi12_normalized_difference" ->
    centralPoint["Pi12_normalized"] - pythonCentral["Pi12_normalized"],
  "s01_difference_GeV2" ->
    centralPoint["s01_GeV2"] - pythonCentral["s01_fitted_GeV2"],
  "s02_difference_GeV2" ->
    centralPoint["s02_GeV2"] - pythonCentral["s02_fitted_GeV2"],
  "f1_difference_GeV" ->
    centralPoint["f1_GeV"] - pythonCentral["f1_GeV"],
  "f2_difference_GeV" ->
    centralPoint["f2_GeV"] - pythonCentral["f2_GeV"]
|>;
Print["[6/7] Python/Mathematica differences formed"];

summaryLines = Join[
  {"Mathematica Ds1 AA/AB/BB two-point cross-check",
   "================================================",
   "OPE: exact-mass LO perturbative + local d=3 through ms^2 + local d=5.",
   "theta_Ds=26.6 deg is external; no external f1, f2, fA, fB, or overlap parameter is used.",
   "The matrix-diagonalization angle and Pi12 residual are diagnostics.", "",
   "Analytic checks:", ToString[analyticChecks, InputForm], "",
   "Central point M2=2.35 GeV^2:", ToString[centralPoint, InputForm], "",
   "Mathematica minus Python:", ToString[crossCheck, InputForm]},
  {""}
];

scanTable = Prepend[(Values[#] & /@ scanRows), Keys[First[scanRows]]];
scanText = StringRiffle[
  (StringRiffle[ToString[#, InputForm] & /@ #, ","] & /@ scanTable), "\n"];
scanStream = OpenWrite[
  FileNameJoin[{outputDir, "mathematica_twopoint_ds1_matrix_scan.csv"}]];
WriteString[scanStream, scanText, "\n"];
Close[scanStream];

summaryStream = OpenWrite[
  FileNameJoin[{outputDir, "mathematica_twopoint_ds1_matrix_check.txt"}]];
WriteString[summaryStream, StringRiffle[summaryLines, "\n"], "\n"];
Close[summaryStream];
Print["[7/7] regression outputs exported"];

Print[StringRiffle[summaryLines, "\n"]];
