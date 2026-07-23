(* ::Package:: *)

(* Bc1 -> Bcstar gamma: Levi-Civita structure projection.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step2_vec_epsilon_triangle_reduction.wl

   For the vector final state the traces contain Levi-Civita tensors.  A useful
   gauge-invariant structure is

     V_mnr = p_n eps(m,r,p,q) - (p.q) eps(m,n,r,p).

   It is transverse in the photon index n for q^2=0 and transverse in the
   vector-meson index r.  This script projects onto that structure and reduces
   loop scalar products to denominator variables, as in the pseudoscalar case.
*)

<< FeynCalc`

ClearAll[
  mu, nu, rho, alpha, p, q, P, k, mc, mb, ec, ebbar,
  dC1, dC2, dC3, dB1, dB2, dB3, p2, pq,
  gammaA, gammaB, gammaV, numS,
  charmInsertion, antibInsertion,
  vTensor, vNorm, vProjector, vProject, reduceWithRules
];

gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] :=
  I/(mb + mc) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaV[rho_] := GA[rho];
numS[mom_, mass_] := GS[mom] + mass;

charmInsertion[current_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      numS[k + q, mc] .
      GA[nu] .
      numS[k, mc] .
      gammaV[rho] .
      numS[k - p, mb]
    ],
    DiracSubstitute67 -> True
  ];

antibInsertion[current_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      numS[k, mc] .
      gammaV[rho] .
      numS[k - p, mb] .
      GA[nu] .
      numS[k - p + q, mb]
    ],
    DiracSubstitute67 -> True
  ];

vTensor =
  Pair[LorentzIndex[nu], Momentum[p]] *
    Eps[LorentzIndex[mu], LorentzIndex[rho], Momentum[p], Momentum[q]] -
  Pair[Momentum[p], Momentum[q]] *
    Eps[LorentzIndex[mu], LorentzIndex[nu], LorentzIndex[rho], Momentum[p]];

vNorm =
  ScalarProductExpand[
    Contract[vTensor vTensor] /. {
      SP[q, q] -> 0,
      Pair[Momentum[q], Momentum[q]] -> 0,
      Pair[Momentum[p], Momentum[p]] -> p2,
      Pair[Momentum[p], Momentum[q]] -> pq
    }
  ];

vProjector = vTensor/vNorm;

vProject[expr_] :=
  Collect2[
    ScalarProductExpand[
      Contract[expr vProjector] /. {
        SP[q, q] -> 0,
        Pair[Momentum[q], Momentum[q]] -> 0,
        Pair[Momentum[k], Momentum[P]] ->
          Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
        Pair[Momentum[p], Momentum[P]] ->
          Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
        Pair[Momentum[P], Momentum[p]] ->
          Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
        Pair[Momentum[P], Momentum[q]] -> Pair[Momentum[p], Momentum[q]],
        Pair[Momentum[P], Momentum[P]] ->
          Pair[Momentum[p], Momentum[p]] + 2 Pair[Momentum[p], Momentum[q]]
      }
    ],
    {mc, mb, ec, ebbar}
  ];

projectorCheck =
  ScalarProductExpand[
    Contract[vTensor vProjector] /. {
      SP[q, q] -> 0,
      Pair[Momentum[q], Momentum[q]] -> 0,
      Pair[Momentum[p], Momentum[p]] -> p2,
      Pair[Momentum[p], Momentum[q]] -> pq
    }
  ];

charmRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dC2 + mc^2,
  Pair[Momentum[k], Momentum[q]] -> (dC1 - dC2)/2,
  Pair[Momentum[k], Momentum[p]] -> (dC2 + mc^2 + p2 - mb^2 - dC3)/2
};

antibRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dB1 + mc^2,
  Pair[Momentum[k], Momentum[p]] -> (dB1 + mc^2 + p2 - mb^2 - dB2)/2,
  Pair[Momentum[k], Momentum[q]] -> pq + (dB3 - dB2)/2
};

reduceWithRules[expr_, rules_, denoms_] :=
  Collect[Expand[expr /. rules], denoms];

projAC = vProject[charmInsertion[gammaA]];
projAB = vProject[antibInsertion[gammaA]];
projBC = vProject[charmInsertion[gammaB]];
projBB = vProject[antibInsertion[gammaB]];

redAC = reduceWithRules[projAC, charmRules, {dC1, dC2, dC3}];
redAB = reduceWithRules[projAB, antibRules, {dB1, dB2, dB3}];
redBC = reduceWithRules[projBC, charmRules, {dC1, dC2, dC3}];
redBB = reduceWithRules[projBB, antibRules, {dB1, dB2, dB3}];

triAC = Collect2[redAC /. {dC1 -> 0, dC2 -> 0, dC3 -> 0}, {mc, mb, p2, pq}];
triAB = Collect2[redAB /. {dB1 -> 0, dB2 -> 0, dB3 -> 0}, {mc, mb, p2, pq}];
triBC = Collect2[redBC /. {dC1 -> 0, dC2 -> 0, dC3 -> 0}, {mc, mb, p2, pq}];
triBB = Collect2[redBB /. {dB1 -> 0, dB2 -> 0, dB3 -> 0}, {mc, mb, p2, pq}];

canAC = Collect2[Expand[redAC - triAC], {dC1, dC2, dC3, mc, mb}];
canAB = Collect2[Expand[redAB - triAB], {dB1, dB2, dB3, mc, mb}];
canBC = Collect2[Expand[redBC - triBC], {dC1, dC2, dC3, mc, mb}];
canBB = Collect2[Expand[redBB - triBB], {dB1, dB2, dB3, mc, mb}];

totalTriangleA = Collect2[ec triAC + ebbar triAB, {ec, ebbar, mc, mb, p2, pq}];
totalTriangleB = Collect2[ec triBC + ebbar triBB, {ec, ebbar, mc, mb, p2, pq}];

outDir = FileNameJoin[{Directory[], "outputs"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "step2_vec_epsilon_triangle_reduction.txt"}];

stream = OpenWrite[outFile];
WriteString[stream, "Bc1 -> Bcstar gamma epsilon projection and triangle reduction\n"];
WriteString[stream, "=============================================================\n\n"];
WriteString[stream, "Tensor norm = ", ToString[vNorm, InputForm], "\n"];
WriteString[stream, "Projector check V.P = ", ToString[projectorCheck, InputForm], "\n\n"];

WriteString[stream, "J_A c-line projected numerator:\n", ToString[projAC, InputForm], "\n\n"];
WriteString[stream, "J_A anti-b-line projected numerator:\n", ToString[projAB, InputForm], "\n\n"];
WriteString[stream, "J_B c-line projected numerator:\n", ToString[projBC, InputForm], "\n\n"];
WriteString[stream, "J_B anti-b-line projected numerator:\n", ToString[projBB, InputForm], "\n\n"];

WriteString[stream, "J_A c-line triangle core:\n", ToString[triAC, InputForm], "\n\n"];
WriteString[stream, "J_A anti-b-line triangle core:\n", ToString[triAB, InputForm], "\n\n"];
WriteString[stream, "J_A total triangle core:\n", ToString[totalTriangleA, InputForm], "\n\n"];

WriteString[stream, "J_B c-line triangle core:\n", ToString[triBC, InputForm], "\n\n"];
WriteString[stream, "J_B anti-b-line triangle core:\n", ToString[triBB, InputForm], "\n\n"];
WriteString[stream, "J_B total triangle core:\n", ToString[totalTriangleB, InputForm], "\n\n"];

WriteString[stream, "J_A c-line denominator-cancellation part:\n", ToString[canAC, InputForm], "\n\n"];
WriteString[stream, "J_A anti-b-line denominator-cancellation part:\n", ToString[canAB, InputForm], "\n\n"];
WriteString[stream, "J_B c-line denominator-cancellation part:\n", ToString[canBC, InputForm], "\n\n"];
WriteString[stream, "J_B anti-b-line denominator-cancellation part:\n", ToString[canBB, InputForm], "\n\n"];
Close[stream];

Print["Tensor norm = ", vNorm];
Print["Projector check = ", projectorCheck];
Print["Wrote vector epsilon triangle reduction to: ", outFile];
Print["J_A total triangle core = ", totalTriangleA];
Print["J_B total triangle core = ", totalTriangleB];

Quit[];
