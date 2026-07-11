(* ::Package:: *)

(* Bc1 -> Bc gamma: E1 projection and triangle-core reduction.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script Bc_gamma/scripts/step2_ps_e1_triangle_reduction.wl

   This script projects the pseudoscalar-final-state traces onto
     S_mn = p_n q_m - (p.q) g_mn
   and rewrites loop scalar products in terms of inverse propagator
   denominators.  Setting those denominators to zero isolates the genuine
   three-denominator triangle core.  The residual terms are contact or
   denominator-cancellation pieces and must be handled separately in the
   dispersion calculation.
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, p, q, P, k, mc, mb, ec, ebbar,
  dC1, dC2, dC3, dB1, dB2, dB3, p2, pq,
  gammaA, gammaB, gammaP, numS,
  charmInsertion, antibInsertion,
  e1Tensor, e1Projector, e1Project, reduceWithRules
];

gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] :=
  I/(mb + mc) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaP[] := I GA[5];
numS[mom_, mass_] := GS[mom] + mass;

charmInsertion[current_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      numS[k + q, mc] .
      GA[nu] .
      numS[k, mc] .
      gammaP[] .
      numS[k - p, mb]
    ],
    DiracSubstitute67 -> True
  ];

antibInsertion[current_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      numS[k, mc] .
      gammaP[] .
      numS[k - p, mb] .
      GA[nu] .
      numS[k - p + q, mb]
    ],
    DiracSubstitute67 -> True
  ];

e1Tensor = FV[p, nu] FV[q, mu] - SP[p, q] MT[mu, nu];
e1Projector = e1Tensor/(2 SP[p, q]^2);

e1Project[expr_] :=
  Collect2[
    ScalarProductExpand[
      Contract[expr e1Projector] /. {
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
    Contract[e1Tensor e1Projector] /. {
      SP[q, q] -> 0,
      Pair[Momentum[q], Momentum[q]] -> 0
    }
  ];

(* c-line photon denominators:
     dC1 = (k+q)^2 - mc^2,
     dC2 = k^2 - mc^2,
     dC3 = (k-p)^2 - mb^2.
*)
charmRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dC2 + mc^2,
  Pair[Momentum[k], Momentum[q]] -> (dC1 - dC2)/2,
  Pair[Momentum[k], Momentum[p]] -> (dC2 + mc^2 + p2 - mb^2 - dC3)/2
};

(* anti-b-line photon denominators:
     dB1 = k^2 - mc^2,
     dB2 = (k-p)^2 - mb^2,
     dB3 = (k-p+q)^2 - mb^2.
*)
antibRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dB1 + mc^2,
  Pair[Momentum[k], Momentum[p]] -> (dB1 + mc^2 + p2 - mb^2 - dB2)/2,
  Pair[Momentum[k], Momentum[q]] -> pq + (dB3 - dB2)/2
};

reduceWithRules[expr_, rules_, denoms_] :=
  Collect[Expand[expr /. rules], denoms];

projAC = e1Project[charmInsertion[gammaA]];
projAB = e1Project[antibInsertion[gammaA]];
projBC = e1Project[charmInsertion[gammaB]];
projBB = e1Project[antibInsertion[gammaB]];

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

outDir = FileNameJoin[{Directory[], "Bc_gamma", "outputs"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "step2_ps_e1_triangle_reduction.txt"}];

stream = OpenWrite[outFile];
WriteString[stream, "Bc1 -> Bc gamma E1 projection and triangle reduction\n"];
WriteString[stream, "====================================================\n\n"];
WriteString[stream, "Projector check S.P = ", ToString[projectorCheck, InputForm], "\n\n"];

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

Print["Projector check = ", projectorCheck];
Print["Wrote pseudoscalar E1 triangle reduction to: ", outFile];
Print["J_A total triangle core = ", totalTriangleA];
Print["J_B total triangle core = ", totalTriangleB];

Quit[];

