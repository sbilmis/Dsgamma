(* ::Package:: *)

(* Bc1 -> Bc gamma: open-gluon-pair dimension-4 G2 insertions for the
   photon-on-c topology.

   This chunk uses one S_Q^(G) insertion on each of two propagator segments and
   contracts their open Lorentz indices with the vacuum tensor.  It writes the
   E1-projected numerator inventory only; double-Borel reduction comes later.
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, al, be, rh, si, p, q, P, k, mc, mb,
  gammaA, gammaB, gammaP, S0Num, SGNum, GGVacuumTensor,
  e1Tensor, e1Projector, e1Project, segmentNum, traceCPhotonPair,
  currentList, pairList, rows, outDir, outFile, stream
];

gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] := I/(mb + mc) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaP[] := I GA[5];

S0Num[mom_, mass_] := GS[mom] + mass;

SGNum[mom_, mass_, a_, b_] :=
  -1/4 (
    DiracSigma[GA[a], GA[b]] . (GS[mom] + mass) +
    (GS[mom] + mass) . DiracSigma[GA[a], GA[b]]
  );

GGVacuumTensor[a_, b_, r_, t_] := MT[a, r] MT[b, t] - MT[a, t] MT[b, r];

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
    {mc, mb}
  ];

segmentNum[name_, pairA_, pairB_] := Module[
  {mom, mass, indices},
  {mom, mass} = Switch[
    name,
    "c_after_photon", {k + q, mc},
    "c_before_photon", {k, mc},
    "b_spectator", {k - p, mb}
  ];
  indices = If[name === pairA, {al, be}, If[name === pairB, {rh, si}, None]];
  If[indices === None,
    S0Num[mom, mass],
    SGNum[mom, mass, indices[[1]], indices[[2]]]
  ]
];

traceCPhotonPair[current_, pairA_, pairB_] := Module[
  {cAfter, cBefore, bSpec},
  cAfter = segmentNum["c_after_photon", pairA, pairB];
  cBefore = segmentNum["c_before_photon", pairA, pairB];
  bSpec = segmentNum["b_spectator", pairA, pairB];
  Contract[
    GGVacuumTensor[al, be, rh, si]
    DiracSimplify[
      DiracTrace[
        current[mu] . cAfter . GA[nu] . cBefore . gammaP[] . bSpec
      ],
      DiracSubstitute67 -> True
    ]
  ] // Contract // FCE // Simplify
];

currentList = {{"A", gammaA}, {"B", gammaB}};
pairList = {
  {"c_after_photon", "c_before_photon"},
  {"c_after_photon", "b_spectator"},
  {"c_before_photon", "b_spectator"}
};

rows = Flatten[
  Table[
    Module[{expr, projected, expanded},
      expr = traceCPhotonPair[current[[2]], pair[[1]], pair[[2]]];
      projected = e1Project[expr] // Contract // FCE // Simplify;
      expanded = Expand[projected];
      {
        current[[1]],
        pair[[1]],
        pair[[2]],
        Length[List @@ expanded],
        ToString[projected, InputForm]
      }
    ],
    {current, currentList},
    {pair, pairList}
  ],
  1
];

outDir = FileNameJoin[{Directory[], "Bc_gamma", "outputs"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "step3_ps_g2_open_pair_cphoton.txt"}];

stream = OpenWrite[outFile];
WriteString[stream, "Bc1 -> Bc gamma open-gluon-pair G2 chunk: photon on c line\n"];
WriteString[stream, "==============================================================\n\n"];
WriteString[stream, "Overall prefactor: the vacuum-average/color normalization must be attached in the Borel stage.\n"];
WriteString[stream, "Each listed pair has one S_Q^(G) denominator squared on each chosen segment.\n\n"];
Do[
  WriteString[
    stream,
    "current=", row[[1]], ", pair=(", row[[2]], ",", row[[3]], ")",
    ", term_count=", ToString[row[[4]]], "\n",
    row[[5]], "\n\n"
  ],
  {row, rows}
];
Close[stream];

Print["Wrote ", outFile];
Do[
  Print["current=", row[[1]], " pair=(", row[[2]], ",", row[[3]], ") terms=", row[[4]]],
  {row, rows}
];

Quit[];
