(* ::Package:: *)

(* Hard-photon epsilon projection and triangle-core reduction for
   Ds1 -> Ds* gamma.

   Run:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script Dsstar_gamma/scripts/hard_tensor_triangle_reduction.wl

   The output gives the genuine three-denominator triangle cores for photon
   emission from the heavy and strange lines, for both J_A and J_B.
*)

<< FeynCalc`

ClearAll[
  mu, nu, rho, alpha, p, q, P, k, ms, mQ, p2, pq,
  dH0, dH1, dH2, dS0, dS1, dS2
];

gammaV[nu_] := GA[nu];
gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] := I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
numS[mom_, mass_] := GS[mom] + mass;

replaceP = {
  Pair[Momentum[k], Momentum[P]] -> Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
  Pair[Momentum[p], Momentum[P]] -> Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
  Pair[Momentum[P], Momentum[q]] -> Pair[Momentum[p], Momentum[q]],
  Pair[Momentum[q], Momentum[q]] -> 0,
  Pair[LorentzIndex[rho], Momentum[P]] -> Pair[LorentzIndex[rho], Momentum[p]] + Pair[LorentzIndex[rho], Momentum[q]]
};

target = Eps[LorentzIndex[mu], LorentzIndex[nu], LorentzIndex[rho], Momentum[q]];
dualTarget = Eps[LorentzIndex[mu], LorentzIndex[nu], LorentzIndex[rho], Momentum[p]];
targetNorm = Contract[target dualTarget] /. {
  Pair[Momentum[q], Momentum[q]] -> 0,
  SP[q, q] -> 0
};

projectEps[expr_] :=
  Collect2[
    ScalarProductExpand[
      Contract[expr dualTarget/targetNorm] /. {
        Pair[Momentum[q], Momentum[q]] -> 0,
        SP[q, q] -> 0
      } /. replaceP
    ],
    {mQ, ms, Pair[Momentum[k], Momentum[k]], Pair[Momentum[k], Momentum[p]], Pair[Momentum[k], Momentum[q]]}
  ];

heavyTrace[current_] :=
  DiracSimplify[
    DiracTrace[
      gammaV[nu] . numS[k + p, mQ] . GA[rho] .
      numS[k + p + q, mQ] . current[mu] . numS[k, ms]
    ],
    DiracTraceEvaluate -> True,
    DiracSubstitute67 -> True
  ];

strangeTrace[current_] :=
  DiracSimplify[
    DiracTrace[
      gammaV[nu] . numS[k + p, mQ] . current[mu] .
      numS[k - q, ms] . GA[rho] . numS[k, ms]
    ],
    DiracTraceEvaluate -> True,
    DiracSubstitute67 -> True
  ];

projHA = projectEps[heavyTrace[gammaA]];
projHB = projectEps[heavyTrace[gammaB]];
projSA = projectEps[strangeTrace[gammaA]];
projSB = projectEps[strangeTrace[gammaB]];

(* Heavy emission denominators:
   dH0 = k^2 - ms^2
   dH1 = (k+p)^2 - mQ^2
   dH2 = (k+p+q)^2 - mQ^2
*)
heavyRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dH0 + ms^2,
  Pair[Momentum[k], Momentum[p]] -> (dH1 - dH0 - ms^2 - p2 + mQ^2)/2,
  Pair[Momentum[k], Momentum[q]] -> (dH2 - dH1 - 2 pq)/2
};

(* Strange emission denominators:
   dS0 = (k+p)^2 - mQ^2
   dS1 = (k-q)^2 - ms^2
   dS2 = k^2 - ms^2
*)
strangeRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dS2 + ms^2,
  Pair[Momentum[k], Momentum[q]] -> (dS2 - dS1)/2,
  Pair[Momentum[k], Momentum[p]] -> (dS0 - dS2 - ms^2 - p2 + mQ^2)/2
};

reduce[expr_, rules_, ds_] := Collect[Expand[expr /. rules], ds];
triangle[expr_, ds_] := Collect2[expr /. Thread[ds -> ConstantArray[0, Length[ds]]], {mQ, ms, p2, pq}];
cancelPart[red_, tri_] := Collect2[Expand[red - tri], {dH0, dH1, dH2, dS0, dS1, dS2, mQ, ms, p2, pq}];

redHA = reduce[projHA, heavyRules, {dH0, dH1, dH2}];
redHB = reduce[projHB, heavyRules, {dH0, dH1, dH2}];
redSA = reduce[projSA, strangeRules, {dS0, dS1, dS2}];
redSB = reduce[projSB, strangeRules, {dS0, dS1, dS2}];

triHA = triangle[redHA, {dH0, dH1, dH2}];
triHB = triangle[redHB, {dH0, dH1, dH2}];
triSA = triangle[redSA, {dS0, dS1, dS2}];
triSB = triangle[redSB, {dS0, dS1, dS2}];

canHA = cancelPart[redHA, triHA];
canHB = cancelPart[redHB, triHB];
canSA = cancelPart[redSA, triSA];
canSB = cancelPart[redSB, triSB];

outDir = FileNameJoin[{Directory[], "Dsstar_gamma", "outputs"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "hard_tensor_triangle_reduction.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Hard-photon epsilon projection: triangle reduction\n"];
WriteString[stream, "==================================================\n\n"];
WriteString[stream, "Target norm:\n", ToString[targetNorm, InputForm], "\n\n"];
WriteString[stream, "A heavy triangle core:\n", ToString[triHA, InputForm], "\n\n"];
WriteString[stream, "B heavy triangle core:\n", ToString[triHB, InputForm], "\n\n"];
WriteString[stream, "A strange triangle core:\n", ToString[triSA, InputForm], "\n\n"];
WriteString[stream, "B strange triangle core:\n", ToString[triSB, InputForm], "\n\n"];
WriteString[stream, "A heavy cancellation:\n", ToString[canHA, InputForm], "\n\n"];
WriteString[stream, "B heavy cancellation:\n", ToString[canHB, InputForm], "\n\n"];
WriteString[stream, "A strange cancellation:\n", ToString[canSA, InputForm], "\n\n"];
WriteString[stream, "B strange cancellation:\n", ToString[canSB, InputForm], "\n\n"];
Close[stream];

Print["Wrote hard tensor triangle reduction to: ", outFile];
Print["A heavy core = ", triHA];
Print["B heavy core = ", triHB];
Print["A strange core = ", triSA];
Print["B strange core = ", triSB];

Quit[];
