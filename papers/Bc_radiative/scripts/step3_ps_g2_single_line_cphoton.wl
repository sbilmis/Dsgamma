(* ::Package:: *)

(* Bc1 -> Bc gamma: single-line dimension-4 G2 insertions for the photon-on-c
   topology.

   This is a deliberately small FeynCalc chunk.  It projects only the
   single-propagator S_Q^(G2) pieces for the topology

     J_i -- S_c(k+q) -- gamma_em -- S_c(k) -- j_5 -- S_b(k-p).

   It does not yet perform the double-Borel transform.  The output is the E1
   projected numerator inventory that must be reduced topology by topology.
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, p, q, P, k, mc, mb, G2,
  gammaA, gammaB, gammaP, S0Num, SG2Num,
  e1Tensor, e1Projector, e1Project,
  traceCPhoton, currentList, segmentList, rows, outDir, outFile, stream
];

gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] := I/(mb + mc) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaP[] := I GA[5];

S0Num[mom_, mass_] := GS[mom] + mass;

(* Same numerator convention as the Bc-mixing workbench:
   S_Q^(G2) numerator is m_Q (k^2 + m_Q slash{k}); the overall G2/12 and
   denominator power are attached later. *)
SG2Num[mom_, mass_] := mass (SP[mom, mom] + mass GS[mom]);

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

traceCPhoton[current_, segment_] := Module[
  {cAfter, cBefore, bSpec},
  cAfter = If[segment === "c_after_photon", SG2Num[k + q, mc], S0Num[k + q, mc]];
  cBefore = If[segment === "c_before_photon", SG2Num[k, mc], S0Num[k, mc]];
  bSpec = If[segment === "b_spectator", SG2Num[k - p, mb], S0Num[k - p, mb]];
  DiracSimplify[
    DiracTrace[
      current[mu] . cAfter . GA[nu] . cBefore . gammaP[] . bSpec
    ],
    DiracSubstitute67 -> True
  ]
];

currentList = {{"A", gammaA}, {"B", gammaB}};
segmentList = {
  {"c_after_photon", "(k+q)^2-m_c^2 raised from 1 to 4"},
  {"c_before_photon", "k^2-m_c^2 raised from 1 to 4"},
  {"b_spectator", "(k-p)^2-m_b^2 raised from 1 to 4"}
};

rows = Flatten[
  Table[
    Module[{expr, projected, expanded},
      expr = traceCPhoton[current[[2]], segment[[1]]];
      projected = e1Project[expr];
      expanded = Expand[projected];
      {
        current[[1]],
        segment[[1]],
        segment[[2]],
        Length[List @@ expanded],
        ToString[projected, InputForm]
      }
    ],
    {current, currentList},
    {segment, segmentList}
  ],
  1
];

outDir = FileNameJoin[{Directory[], "outputs"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "step3_ps_g2_single_line_cphoton.txt"}];

stream = OpenWrite[outFile];
WriteString[stream, "Bc1 -> Bc gamma single-line G2 chunk: photon on c line\n"];
WriteString[stream, "========================================================\n\n"];
WriteString[stream, "Overall prefactor per row: G2/12 times the usual color and loop factors.\n"];
WriteString[stream, "The listed segment denominator is raised from power 1 to power 4.\n\n"];
Do[
  WriteString[
    stream,
    "current=", row[[1]], ", segment=", row[[2]], ", denominator note=", row[[3]],
    ", term_count=", ToString[row[[4]]], "\n",
    row[[5]], "\n\n"
  ],
  {row, rows}
];
Close[stream];

Print["Wrote ", outFile];
Do[
  Print["current=", row[[1]], " segment=", row[[2]], " terms=", row[[4]]],
  {row, rows}
];

Quit[];
