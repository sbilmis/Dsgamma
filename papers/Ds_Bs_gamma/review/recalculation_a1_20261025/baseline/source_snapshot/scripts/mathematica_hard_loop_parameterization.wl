(* ::Package:: *)

(* Feynman-parameter audit for the perturbative hard-photon triangles in
   Ds1/Bs1 -> Ds/Bs gamma.

   This starts from the standard correlator routing verified in
   mathematica_correlator_routing_audit.wl and derives the shifted loop
   denominators before loop integration.  It intentionally stops before
   discontinuities, continuum subtraction, and the double Borel transform.
*)

<< FeynCalc`

ClearAll[
  k, p, q, pPrime, x, y, z, mQ, ms, p2, pPrime2,
  dH1, dH2, dH3, dS1, dS2, dS3,
  deltaH, deltaS, shiftH, shiftS, expandSP, simplexReduce
];

(* Tell FeynCalc that the Feynman parameters are scalar coefficients rather
   than momentum labels, so scalar products of shifted momenta expand. *)
(DataType[#, FCVariable] = True) & /@ {x, y, z};

scriptDir = DirectoryName[$InputFileName];
paperDir = DirectoryName[scriptDir];
outDir = FileNameJoin[{paperDir, "outputs"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "mathematica_hard_loop_parameterization.txt"}];

(* q^2=0, p^2=p2, and (p+q)^2=pPrime2. *)
kinematicRules = {
  Pair[Momentum[q], Momentum[q]] -> 0,
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> (pPrime2 - p2)/2,
  Pair[Momentum[q], Momentum[p]] -> (pPrime2 - p2)/2
};

expandSP[expression_] :=
  ExpandScalarProduct[FCI[expression]] /. kinematicRules // FCE // Expand;

simplexReduce[expression_] :=
  Together[expandSP[expression] /. z -> 1 - x - y] // Simplify;

(* Heavy-line emission.  Parameters follow the denominator labels:

     x dH1 + y dH2 + z dH3,

   dH1=(k+q)^2-mQ^2, dH2=k^2-mQ^2,
   dH3=(k-p)^2-ms^2. *)
dH1 = SP[k + q, k + q] - mQ^2;
dH2 = SP[k, k] - mQ^2;
dH3 = SP[k - p, k - p] - ms^2;

shiftH = x q - z p;
deltaH = (x + y) mQ^2 + z ms^2 - x z pPrime2 - y z p2;

weightedH = x dH1 + y dH2 + z dH3;
shiftedH = SP[k + shiftH, k + shiftH] - deltaH;
heavyResidual = simplexReduce[weightedH - shiftedH];

(* Strange-line emission.  Again parameters follow the denominator labels:

     x dS1 + y dS2 + z dS3,

   dS1=k^2-mQ^2, dS2=(k-p)^2-ms^2,
   dS3=(k-p-q)^2-ms^2=(k-pPrime)^2-ms^2. *)
dS1 = SP[k, k] - mQ^2;
dS2 = SP[k - p, k - p] - ms^2;
dS3 = SP[k - p - q, k - p - q] - ms^2;

shiftS = -y p - z (p + q);
deltaS = x mQ^2 + (y + z) ms^2 - x y p2 - x z pPrime2;

weightedS = x dS1 + y dS2 + z dS3;
shiftedS = SP[k + shiftS, k + shiftS] - deltaS;
strangeResidual = simplexReduce[weightedS - shiftedS];

(* A general pairwise-invariant reconstruction gives an independent check:

   Delta = Sum_i x_i m_i^2
           - Sum_{i<j} x_i x_j (r_i-r_j)^2. *)
deltaHPairwise =
  (x + y) mQ^2 + z ms^2 - x y 0 - x z pPrime2 - y z p2;
deltaSPairwise =
  x mQ^2 + (y + z) ms^2 - x y p2 - x z pPrime2 - y z 0;

checks = {
  {"heavy shifted-denominator identity", heavyResidual},
  {"strange shifted-denominator identity", strangeResidual},
  {"heavy pairwise-invariant reconstruction", Simplify[deltaH - deltaHPairwise]},
  {"strange pairwise-invariant reconstruction", Simplify[deltaS - deltaSPairwise]}
};

allChecksPassed = And @@ (TrueQ[#[[2]] === 0] & /@ checks);

stream = OpenWrite[outFile];
writeLine[text_] := WriteString[stream, text, "\n"];
writeExpr[label_, expression_] :=
  WriteString[stream, label, "\n", ToString[expression, InputForm], "\n\n"];

writeLine["Hard-triangle Feynman-parameter audit"];
writeLine["========================================"];
writeLine[""];
writeLine["Kinematics: q^2=0, p^2=p2, (p+q)^2=pPrime2, x+y+z=1."];
writeLine[""];

writeLine["HEAVY-LINE EMISSION"];
writeLine["x multiplies dH1=(k+q)^2-mQ^2;"];
writeLine["y multiplies dH2=k^2-mQ^2;"];
writeLine["z multiplies dH3=(k-p)^2-ms^2."];
writeExpr["Loop shift R_H in ell=k+R_H:", shiftH];
writeExpr["Delta_H:", deltaH];
writeLine["x dH1+y dH2+z dH3 = ell^2-Delta_H."];
writeExpr["Identity residual:", heavyResidual];

writeLine["STRANGE-LINE EMISSION"];
writeLine["x multiplies dS1=k^2-mQ^2;"];
writeLine["y multiplies dS2=(k-p)^2-ms^2;"];
writeLine["z multiplies dS3=(k-p-q)^2-ms^2."];
writeExpr["Loop shift R_S in ell=k+R_S:", shiftS];
writeExpr["Delta_S:", deltaS];
writeLine["x dS1+y dS2+z dS3 = ell^2-Delta_S."];
writeExpr["Identity residual:", strangeResidual];

writeLine["FEYNMAN-PARAMETER MEASURE"];
writeLine["1/(D1 D2 D3) = 2 Integral_0^1 dx Integral_0^(1-x) dy"];
writeLine["                  1/[x D1+y D2+(1-x-y)D3]^3."];
writeLine[""];
writeLine["These Delta polynomials are the starting point for the loop integral,"];
writeLine["double discontinuity, continuum subtraction, and double Borel transform."];
writeLine[""];

Do[
  writeLine[
    If[TrueQ[check[[2]] === 0], "PASS  ", "FAIL  "] <> check[[1]] <>
      If[TrueQ[check[[2]] === 0], "", "; residual = " <> ToString[check[[2]], InputForm]]
  ],
  {check, checks}
];
writeLine[""];
writeLine["OVERALL STATUS: " <> If[allChecksPassed, "PASS", "FAIL"]];

Close[stream];

Print["Wrote audit to: ", outFile];
Print["Checks passed: ", Count[checks, {_, 0}], "/", Length[checks]];
Print["Overall status: ", If[allChecksPassed, "PASS", "FAIL"]];

If[! allChecksPassed, Exit[1]];
Quit[];
