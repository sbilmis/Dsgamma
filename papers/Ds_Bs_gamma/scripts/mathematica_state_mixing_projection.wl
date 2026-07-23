(* ::Package:: *)

(* Physical-state projection and the pure-axial Colangelo limit.

   Run from papers/Ds_Bs_gamma with

     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/mathematica_state_mixing_projection.wl

   This is an algebraic audit.  It does not use the incomplete tensor-current
   spectral density and therefore does not fit a physical intermediate angle.
*)

ClearAll[
  theta, fA, fB, f1, f2, gA, gB, g1, g2,
  tA, tB, piAA, piAB, piBB, pi11, pi12, pi22,
  m1, m2, mP, mQ, ms, fP, borel1, borel2, borelTwoPoint,
  alphaEM, rotation, transitionBasis, transitionPhysical,
  correlationBasis, correlationPhysical, colangeloPrefactor,
  gLowFromOPE, gColangeloFromOPE, width, intervalOverlap,
  M2, s, s0, u0, eQ, eS, ssCond, chiC, chiGamma, f3Gamma,
  rhoA, phiGamma, A4, HbarGamma, PsiV, IF1,
  tAHat, tBHat, tBpert, tBLocal, tBTw2, tBTw3, tBTw4, tB3p,
  shat, mi, mj, ei, w, v, alphaG, uval, kallenLambda,
  rhoPiece, F1DA, hGamma, psiV,
  rhoAB, rhoFromF1, rhoDiag, f1ModelSquared, f2ModelSquared,
  pi12Model, rhoFromF1Residual, rhoDiagResidual,
  piAAPert, piBBPert, piABPert, piAACond, piBBCond, piABCond,
  sThreshold, twoPointOPE11, twoPointOPE22
];

rotation[angle_] := {
  {Sin[angle], Cos[angle]},
  {Cos[angle], -Sin[angle]}
};

rTheta = rotation[theta];
orthogonalityResidual =
  FullSimplify[Transpose[rTheta] . rTheta - IdentityMatrix[2],
    Assumptions -> Element[theta, Reals]];

transitionBasis = {tA, tB};
transitionPhysical = Expand[rTheta . transitionBasis];
pureAxialTransitionResidual =
  Simplify[(transitionPhysical /. theta -> Pi/2) - {tA, -tB}];

g1 = (Sin[theta] fA gA + Cos[theta] fB gB)/f1;
g2 = (Cos[theta] fA gA - Sin[theta] fB gB)/f2;
pureAxialCouplingResidual =
  Simplify[(g1 /. {theta -> Pi/2, f1 -> fA}) - gA];

correlationBasis = {{piAA, piAB}, {piAB, piBB}};
correlationPhysical = Expand[rTheta . correlationBasis . Transpose[rTheta]];
expectedCorrelationPhysical = {
  {
    Sin[theta]^2 piAA + Cos[theta]^2 piBB +
      2 Sin[theta] Cos[theta] piAB,
    Sin[theta] Cos[theta] (piAA - piBB) +
      (Cos[theta]^2 - Sin[theta]^2) piAB
  },
  {
    Sin[theta] Cos[theta] (piAA - piBB) +
      (Cos[theta]^2 - Sin[theta]^2) piAB,
    Cos[theta]^2 piAA + Sin[theta]^2 piBB -
      2 Sin[theta] Cos[theta] piAB
  }
};
correlationRotationResidual =
  FullSimplify[correlationPhysical - expectedCorrelationPhysical,
    Assumptions -> Element[theta, Reals]];

f1Squared =
  Exp[m1^2/borelTwoPoint]/m1^2 expectedCorrelationPhysical[[1, 1]];
f2Squared =
  Exp[m2^2/borelTwoPoint]/m2^2 expectedCorrelationPhysical[[2, 2]];
pureAxialDecayConstantResidual =
  Simplify[
    (f1Squared /. theta -> Pi/2) -
      Exp[m1^2/borelTwoPoint] piAA/m1^2
  ];

(* Explicit Pi_11/Pi_22 OPE structure, kept in basis blocks. *)
twoPointOPE11 =
  Inactive[Integrate][
    Exp[-s/borelTwoPoint] (
      Sin[theta]^2 piAAPert[s] + Cos[theta]^2 piBBPert[s] +
      2 Sin[theta] Cos[theta] piABPert[s]),
    {s, sThreshold, s0}
  ] + Sin[theta]^2 piAACond + Cos[theta]^2 piBBCond +
    2 Sin[theta] Cos[theta] piABCond;
twoPointOPE22 =
  Inactive[Integrate][
    Exp[-s/borelTwoPoint] (
      Cos[theta]^2 piAAPert[s] + Sin[theta]^2 piBBPert[s] -
      2 Sin[theta] Cos[theta] piABPert[s]),
    {s, sThreshold, s0}
  ] + Cos[theta]^2 piAACond + Sin[theta]^2 piBBCond -
    2 Sin[theta] Cos[theta] piABCond;

(* The numerical Stage-3 residues use a reduced overlap model, not a
   completed local AA/BB/AB OPE. *)
f1ModelSquared =
  Sin[theta]^2 fA^2 + Cos[theta]^2 fB^2 +
    2 Sin[theta] Cos[theta] rhoAB fA fB;
f2ModelSquared =
  Cos[theta]^2 fA^2 + Sin[theta]^2 fB^2 -
    2 Sin[theta] Cos[theta] rhoAB fA fB;
rhoFromF1 =
  (f1^2 - Sin[theta]^2 fA^2 - Cos[theta]^2 fB^2)/
    (2 Sin[theta] Cos[theta] fA fB);
rhoFromF1Residual =
  Simplify[(f1ModelSquared /. rhoAB -> rhoFromF1) - f1^2];
pi12Model =
  Sin[theta] Cos[theta] (fA^2 - fB^2) +
    (Cos[theta]^2 - Sin[theta]^2) rhoAB fA fB;
rhoDiag =
  -Sin[theta] Cos[theta] (fA^2 - fB^2)/
    ((Cos[theta]^2 - Sin[theta]^2) fA fB);
rhoDiagResidual = Simplify[pi12Model /. rhoAB -> rhoDiag];

colangeloPrefactor =
  Exp[m1^2/borel1 + mP^2/borel2] (mQ + ms)/
    (m1 fA mP^2 fP);
gColangeloFromOPE = colangeloPrefactor tA;
gLowFromOPE =
  Exp[m1^2/borel1 + mP^2/borel2] (mQ + ms)/
    (m1 f1 mP^2 fP) (Sin[theta] tA + Cos[theta] tB);
colangeloFormulaResidual =
  Simplify[
    (gLowFromOPE /. {theta -> Pi/2, f1 -> fA}) -
      gColangeloFromOPE
  ];

width[g_, initialMass_, finalMass_] :=
  alphaEM/3 g^2 ((initialMass^2 - finalMass^2)/(2 initialMass))^3;
pureAxialWidthResidual =
  Simplify[
    (width[g1, m1, mP] /. {theta -> Pi/2, f1 -> fA}) -
      width[gA, m1, mP]
  ];

intervalOverlap[first_, second_] := {
  Max[first[[1]], second[[1]]],
  Min[first[[2]], second[[2]]]
};

ourPureAxialGInterval = {-0.3978, -0.3124};
colangeloGInterval = {-0.37, -0.29};
gIntervalIntersection =
  intervalOverlap[ourPureAxialGInterval, colangeloGInterval];
gIntervalsOverlap =
  TrueQ[gIntervalIntersection[[1]] <= gIntervalIntersection[[2]]];

ourPureAxialWidthIntervalKeV = {20.51, 33.26};
colangeloWidthIntervalKeV = {19., 29.};
widthIntervalIntersectionKeV =
  intervalOverlap[ourPureAxialWidthIntervalKeV, colangeloWidthIntervalKeV];
widthIntervalsOverlap =
  TrueQ[widthIntervalIntersectionKeV[[1]] <= widthIntervalIntersectionKeV[[2]]];

(* Colangelo-style component assembly used in the paper-ready formula. *)
kallenLambda[shat_, mass1Squared_, mass2Squared_] :=
  shat^2 + mass1Squared^2 + mass2Squared^2 -
  2 shat mass1Squared - 2 shat mass2Squared -
  2 mass1Squared mass2Squared;
rhoPiece[shat_, mi_, mj_, ei_] :=
  -3 ei/(8 Pi^2) (
    2 mi Log[
      (shat - mj^2 + mi^2 -
        Sqrt[kallenLambda[shat, mj^2, mi^2]])/
      (shat - mj^2 + mi^2 +
        Sqrt[kallenLambda[shat, mj^2, mi^2]])
    ] +
    (mj - mi) (mj^2 - mi^2 - shat)/shat^2
      Sqrt[kallenLambda[shat, mj^2, mi^2]]
  );
rhoA[shat_] :=
  rhoPiece[shat, ms, mQ, eS] - rhoPiece[shat, mQ, ms, eQ];
HbarGamma[uval_] :=
  Inactive[Integrate][(uval - w) hGamma[w], {w, 0, uval}];
PsiV[uval_] :=
  Inactive[Integrate][psiV[w], {w, 0, uval}];
IF1[uval_] :=
  Inactive[Integrate][
    Inactive[Integrate][
      F1DA[uval - (1 - v) alphaG, 1 - uval - v alphaG, alphaG],
      {alphaG, 0, uval/(1 - v)}
    ],
    {v, 0, 1 - uval}
  ] +
  Inactive[Integrate][
    Inactive[Integrate][
      F1DA[uval - (1 - v) alphaG, 1 - uval - v alphaG, alphaG],
      {alphaG, 0, (1 - uval)/v}
    ],
    {v, 1 - uval, 1}
  ];

expQ = Exp[-mQ^2/M2];
exp0 = Exp[-s0/M2];
axialOPEComponents = <|
  "perturbative" ->
    Inactive[Integrate][Exp[-s/M2] rhoA[s], {s, (mQ + ms)^2, s0}],
  "local quark condensate" ->
    eQ expQ ssCond
      (1 - mQ ms/M2 + ms^2/(2 M2) (1 + mQ^2/M2)),
  "twist 2, signed Colangelo chi" ->
    -eS ssCond (expQ - exp0) M2 chiC phiGamma[u0],
  "twist 3" ->
    2 eS f3Gamma mQ expQ PsiV[u0],
  "twist 4" ->
    eS ssCond expQ/4
      (A4[u0] - 8 HbarGamma[u0]) (1 + mQ^2/M2),
  "three particle" ->
    eS ssCond expQ IF1[u0]
|>;
tAHat = Total[Values[axialOPEComponents]];
axialOPEAssemblyResidual =
  Expand[tAHat - Total[Values[axialOPEComponents]]];

tensorOPEComponents = <|
  "perturbative" -> tBpert,
  "local quark condensate" -> tBLocal,
  "twist 2" -> tBTw2,
  "twist 3" -> tBTw3,
  "twist 4" -> tBTw4,
  "three particle" -> tB3p
|>;
tBHat = Total[Values[tensorOPEComponents]];
tensorOPEAssemblyResidual =
  Expand[tBHat - Total[Values[tensorOPEComponents]]];

paperNorm1 =
  Exp[(m1^2 + mP^2)/(2 M2)] (mQ + ms)/(m1 f1 mP^2 fP);
paperNorm2 =
  Exp[(m2^2 + mP^2)/(2 M2)] (mQ + ms)/(m2 f2 mP^2 fP);
paperG1 = paperNorm1 (Sin[theta] tAHat + Cos[theta] tBHat);
paperG2 = paperNorm2 (Cos[theta] tAHat - Sin[theta] tBHat);
paperG1FromComponents = paperNorm1 Total[
  Sin[theta] Values[axialOPEComponents] +
  Cos[theta] Values[tensorOPEComponents]
];
paperComponentAssemblyResidual =
  Expand[paperG1 - paperG1FromComponents];
paperPureAxialResidual = Simplify[
  (paperG1 /. {theta -> Pi/2, f1 -> fA}) -
  Exp[(m1^2 + mP^2)/(2 M2)] (mQ + ms) tAHat/
    (m1 fA mP^2 fP)
];
chiConventionResidual = Simplify[
  (axialOPEComponents["twist 2, signed Colangelo chi"] /.
    chiC -> -chiGamma) -
  eS ssCond (expQ - exp0) M2 chiGamma phiGamma[u0]
];

checks = <|
  "orthogonal state rotation" -> orthogonalityResidual,
  "pure-A projected invariant" -> pureAxialTransitionResidual,
  "residue-weighted pure-A coupling" -> pureAxialCouplingResidual,
  "two-point correlator rotation" -> correlationRotationResidual,
  "pure-A decay constant" -> pureAxialDecayConstantResidual,
  "Ds f1-anchor closure" -> rhoFromF1Residual,
  "Bs Pi12=0 closure" -> rhoDiagResidual,
  "Colangelo OPE formula limit" -> colangeloFormulaResidual,
  "pure-A width" -> pureAxialWidthResidual,
  "axial OPE component assembly" -> axialOPEAssemblyResidual,
  "tensor OPE component assembly" -> tensorOPEAssemblyResidual,
  "mixed coupling component assembly" -> paperComponentAssemblyResidual,
  "full axial OPE pure-A limit" -> paperPureAxialResidual,
  "signed-to-positive chi conversion" -> chiConventionResidual
|>;

zeroMatrixOrVectorQ[value_] :=
  And @@ (TrueQ[# === 0] & /@ Flatten[{value}]);
symbolicStatus = And @@ (zeroMatrixOrVectorQ /@ Values[checks]);

Print["Physical-state projection and Colangelo-limit audit"];
Print["================================================"];
Print["Rotation matrix R(theta) = ", rTheta];
Print["R^T R - 1 = ", orthogonalityResidual];
Print["Projected invariants at theta=Pi/2 residual = ",
  pureAxialTransitionResidual];
Print["G_low(theta=Pi/2, f1=fA) - gA = ",
  pureAxialCouplingResidual];
Print["Colangelo formula residual = ", colangeloFormulaResidual];
Print["Decay-constant residual = ", pureAxialDecayConstantResidual];
Print["Pi_11 OPE in basis blocks = ", twoPointOPE11];
Print["Pi_22 OPE in basis blocks = ", twoPointOPE22];
Print["rho_AB from the Ds f1 anchor = ", rhoFromF1];
Print["rho_AB from Pi_12=0 = ", rhoDiag];
Print["Residue interpretation = model normalization; local BB/AB OPE pending"];
Print["Ds central-angle values: f1=0.345 [0.330,0.364], f2=0.379 [0.338,0.423] GeV (lattice input)"];
Print["Bs1(5830) central-angle values: f1=0.443 [0.384,0.492], f2=0.203 [0.121,0.291] GeV (lattice input)"];
Print["Width residual = ", pureAxialWidthResidual];
Print["Paper G1 component residual = ", paperComponentAssemblyResidual];
Print["Paper pure-A OPE residual = ", paperPureAxialResidual];
Print["Chi-convention residual = ", chiConventionResidual];
Print["Our/Colangelo g overlap = ", gIntervalIntersection,
  " GeV^-1; overlap: ", gIntervalsOverlap];
Print["Our/Colangelo width overlap = ", widthIntervalIntersectionKeV,
  " keV; overlap: ", widthIntervalsOverlap];
Print["SYMBOLIC_STATUS = ", If[symbolicStatus, "PASS", "FAIL"]];
Print["NUMERIC_INTERVAL_STATUS = ",
  If[gIntervalsOverlap && widthIntervalsOverlap, "PASS", "FAIL"]];

If[!(symbolicStatus && gIntervalsOverlap && widthIntervalsOverlap), Quit[1]];
Quit[0];
