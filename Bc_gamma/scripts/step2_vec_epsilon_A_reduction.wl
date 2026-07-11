(* ::Package:: *)

(* Bc1 -> Bcstar gamma: J_A epsilon-structure reduction only.

   This split file avoids expanding the much larger tensor-current J_B traces.
*)

<< FeynCalc`

ClearAll[
  mu, nu, rho, p, q, k, mc, mb, ec, ebbar,
  dC1, dC2, dC3, dB1, dB2, dB3, p2, pq,
  gammaA, gammaV, numS, charmInsertion, antibInsertion,
  vTensor, vNorm, vProjector, vProject, reduceWithRules
];

gammaA[mu_] := GA[mu] . GA[5];
gammaV[rho_] := GA[rho];
numS[mom_, mass_] := GS[mom] + mass;

charmInsertion[] :=
  DiracSimplify[
    DiracTrace[
      gammaA[mu] .
      numS[k + q, mc] .
      GA[nu] .
      numS[k, mc] .
      gammaV[rho] .
      numS[k - p, mb]
    ],
    DiracSubstitute67 -> True
  ];

antibInsertion[] :=
  DiracSimplify[
    DiracTrace[
      gammaA[mu] .
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
        Pair[Momentum[q], Momentum[q]] -> 0
      }
    ],
    {mc, mb}
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

projC = vProject[charmInsertion[]];
projB = vProject[antibInsertion[]];

redC = reduceWithRules[projC, charmRules, {dC1, dC2, dC3}];
redB = reduceWithRules[projB, antibRules, {dB1, dB2, dB3}];

triC = Collect2[redC /. {dC1 -> 0, dC2 -> 0, dC3 -> 0}, {mc, mb, p2, pq}];
triB = Collect2[redB /. {dB1 -> 0, dB2 -> 0, dB3 -> 0}, {mc, mb, p2, pq}];
totalTriangleA = Collect2[ec triC + ebbar triB, {ec, ebbar, mc, mb, p2, pq}];

outDir = FileNameJoin[{Directory[], "Bc_gamma", "outputs"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "step2_vec_epsilon_A_reduction.txt"}];

stream = OpenWrite[outFile];
WriteString[stream, "Bc1 -> Bcstar gamma J_A epsilon-structure reduction\n"];
WriteString[stream, "===================================================\n\n"];
WriteString[stream, "Tensor norm = ", ToString[vNorm, InputForm], "\n"];
WriteString[stream, "Projector check V.P = ", ToString[projectorCheck, InputForm], "\n\n"];
WriteString[stream, "J_A c-line projected numerator:\n", ToString[projC, InputForm], "\n\n"];
WriteString[stream, "J_A anti-b-line projected numerator:\n", ToString[projB, InputForm], "\n\n"];
WriteString[stream, "J_A c-line triangle core:\n", ToString[triC, InputForm], "\n\n"];
WriteString[stream, "J_A anti-b-line triangle core:\n", ToString[triB, InputForm], "\n\n"];
WriteString[stream, "J_A total triangle core:\n", ToString[totalTriangleA, InputForm], "\n\n"];
Close[stream];

Print["Tensor norm = ", vNorm];
Print["Projector check = ", projectorCheck];
Print["J_A total triangle core = ", totalTriangleA];
Print["Wrote vector J_A reduction to: ", outFile];

Quit[];
