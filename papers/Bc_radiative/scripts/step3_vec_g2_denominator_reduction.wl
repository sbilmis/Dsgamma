(* ::Package:: *)

(* Bc1 -> Bcstar gamma: denominator reduction for projected G2 numerators.

   This is the vector-final-state analogue of
   step3_ps_g2_denominator_reduction.wl.  It projects every dimension-4
   background-field numerator onto the gauge-invariant

     V_{mu nu rho}
       = p_nu eps_{mu rho p q} - (p.q) eps_{mu nu rho p}

   structure and rewrites scalar products in terms of propagator denominators.
   Rows with positive effective denominator powers are ordinary three-propagator
   terms.  Rows with nonpositive powers are contact/derivative-sector terms.
*)

<< FeynCalc`

ClearAll[
  mu, nu, rho, alpha, al, be, rh, si, p, q, P, k, mc, mb,
  d1, d2, d3, p2, pq,
  gammaA, gammaB, gammaV, S0Num, SG2Num, SGNum, GGVacuumTensor,
  vTensor, vNorm, vProjector, vProject,
  cSegmentNumSingle, bSegmentNumSingle, cSegmentNumPair, bSegmentNumPair,
  traceCPhotonSingle, traceBbarPhotonSingle,
  traceCPhotonPair, traceBbarPhotonPair,
  reduceExpr, termRows, appendRows, currentList,
  cSegments, bSegments, cPairs, bPairs, allRows, outDir, outFile, stream,
  csvField
];

gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] := I/(mb + mc) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaV[rho_] := GA[rho];

S0Num[mom_, mass_] := GS[mom] + mass;
SG2Num[mom_, mass_] := mass (SP[mom, mom] + mass GS[mom]);
SGNum[mom_, mass_, a_, b_] :=
  -1/4 (
    DiracSigma[GA[a], GA[b]] . (GS[mom] + mass) +
    (GS[mom] + mass) . DiracSigma[GA[a], GA[b]]
  );
GGVacuumTensor[a_, b_, r_, t_] := MT[a, r] MT[b, t] - MT[a, t] MT[b, r];

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
    {mc, mb}
  ];

cSegmentNumSingle[name_, selected_] := Switch[
  name,
  "c_after_photon", If[selected === name, SG2Num[k + q, mc], S0Num[k + q, mc]],
  "c_before_photon", If[selected === name, SG2Num[k, mc], S0Num[k, mc]],
  "b_spectator", If[selected === name, SG2Num[k - p, mb], S0Num[k - p, mb]]
];

bSegmentNumSingle[name_, selected_] := Switch[
  name,
  "c_spectator", If[selected === name, SG2Num[k, mc], S0Num[k, mc]],
  "b_before_photon", If[selected === name, SG2Num[k - p, mb], S0Num[k - p, mb]],
  "b_after_photon", If[selected === name, SG2Num[k - p + q, mb], S0Num[k - p + q, mb]]
];

cSegmentNumPair[name_, pairA_, pairB_] := Module[{mom, mass, inds},
  {mom, mass} = Switch[
    name,
    "c_after_photon", {k + q, mc},
    "c_before_photon", {k, mc},
    "b_spectator", {k - p, mb}
  ];
  inds = If[name === pairA, {al, be}, If[name === pairB, {rh, si}, None]];
  If[inds === None, S0Num[mom, mass], SGNum[mom, mass, inds[[1]], inds[[2]]]]
];

bSegmentNumPair[name_, pairA_, pairB_] := Module[{mom, mass, inds},
  {mom, mass} = Switch[
    name,
    "c_spectator", {k, mc},
    "b_before_photon", {k - p, mb},
    "b_after_photon", {k - p + q, mb}
  ];
  inds = If[name === pairA, {al, be}, If[name === pairB, {rh, si}, None]];
  If[inds === None, S0Num[mom, mass], SGNum[mom, mass, inds[[1]], inds[[2]]]]
];

traceCPhotonSingle[current_, selected_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      cSegmentNumSingle["c_after_photon", selected] .
      GA[nu] .
      cSegmentNumSingle["c_before_photon", selected] .
      gammaV[rho] .
      cSegmentNumSingle["b_spectator", selected]
    ],
    DiracSubstitute67 -> True
  ];

traceBbarPhotonSingle[current_, selected_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      bSegmentNumSingle["c_spectator", selected] .
      gammaV[rho] .
      bSegmentNumSingle["b_before_photon", selected] .
      GA[nu] .
      bSegmentNumSingle["b_after_photon", selected]
    ],
    DiracSubstitute67 -> True
  ];

traceCPhotonPair[current_, pairA_, pairB_] :=
  Contract[
    GGVacuumTensor[al, be, rh, si]
    DiracSimplify[
      DiracTrace[
        current[mu] .
        cSegmentNumPair["c_after_photon", pairA, pairB] .
        GA[nu] .
        cSegmentNumPair["c_before_photon", pairA, pairB] .
        gammaV[rho] .
        cSegmentNumPair["b_spectator", pairA, pairB]
      ],
      DiracSubstitute67 -> True
    ]
  ] // Contract // FCE // Simplify;

traceBbarPhotonPair[current_, pairA_, pairB_] :=
  Contract[
    GGVacuumTensor[al, be, rh, si]
    DiracSimplify[
      DiracTrace[
        current[mu] .
        bSegmentNumPair["c_spectator", pairA, pairB] .
        gammaV[rho] .
        bSegmentNumPair["b_before_photon", pairA, pairB] .
        GA[nu] .
        bSegmentNumPair["b_after_photon", pairA, pairB]
      ],
      DiracSubstitute67 -> True
    ]
  ] // Contract // FCE // Simplify;

reduceExpr[expr_, topology_] := Module[{rules},
  rules = If[topology === "c_photon",
    {
      Pair[Momentum[p], Momentum[p]] -> p2,
      Pair[Momentum[p], Momentum[q]] -> pq,
      Pair[Momentum[k], Momentum[k]] -> d2 + mc^2,
      Pair[Momentum[k], Momentum[q]] -> (d1 - d2)/2,
      Pair[Momentum[k], Momentum[p]] -> (d2 + mc^2 + p2 - mb^2 - d3)/2
    },
    {
      Pair[Momentum[p], Momentum[p]] -> p2,
      Pair[Momentum[p], Momentum[q]] -> pq,
      Pair[Momentum[k], Momentum[k]] -> d1 + mc^2,
      Pair[Momentum[k], Momentum[p]] -> (d1 + mc^2 + p2 - mb^2 - d2)/2,
      Pair[Momentum[k], Momentum[q]] -> pq + (d3 - d2)/2
    }
  ];
  Collect[Expand[vProject[expr] /. rules], {d1, d2, d3}]
];

termRows[expr_, basePowers_] := Module[{expanded, terms, powers, coeff},
  expanded = List @@ Expand[expr];
  If[Head[Expand[expr]] =!= Plus, expanded = {Expand[expr]}];
  Table[
    powers = Exponent[term, {d1, d2, d3}];
    coeff = term/(d1^powers[[1]] d2^powers[[2]] d3^powers[[3]]) // Together // Simplify;
    {
      powers[[1]], powers[[2]], powers[[3]],
      basePowers[[1]] - powers[[1]],
      basePowers[[2]] - powers[[2]],
      basePowers[[3]] - powers[[3]],
      ToString[coeff, InputForm]
    },
    {term, expanded}
  ]
];

appendRows[rows_, topology_, class_, currentName_, object_, expr_, basePowers_] := Module[
  {trs = termRows[reduceExpr[expr, topology], basePowers]},
  Join[
    rows,
    ({topology, class, currentName, object} ~Join~ #) & /@ trs
  ]
];

currentList = {{"A", gammaA}, {"B", gammaB}};
cSegments = {"c_after_photon", "c_before_photon", "b_spectator"};
bSegments = {"c_spectator", "b_before_photon", "b_after_photon"};
cPairs = {{"c_after_photon", "c_before_photon"}, {"c_after_photon", "b_spectator"}, {"c_before_photon", "b_spectator"}};
bPairs = {{"c_spectator", "b_before_photon"}, {"c_spectator", "b_after_photon"}, {"b_before_photon", "b_after_photon"}};

allRows = {};
Do[
  allRows = appendRows[
    allRows, "c_photon", "single_line_G2", cur[[1]], seg,
    traceCPhotonSingle[cur[[2]], seg],
    If[seg === "c_after_photon", {4, 1, 1}, If[seg === "c_before_photon", {1, 4, 1}, {1, 1, 4}]]
  ],
  {cur, currentList}, {seg, cSegments}
];
Do[
  allRows = appendRows[
    allRows, "bbar_photon", "single_line_G2", cur[[1]], seg,
    traceBbarPhotonSingle[cur[[2]], seg],
    If[seg === "c_spectator", {4, 1, 1}, If[seg === "b_before_photon", {1, 4, 1}, {1, 1, 4}]]
  ],
  {cur, currentList}, {seg, bSegments}
];
Do[
  allRows = appendRows[
    allRows, "c_photon", "open_pair_GG", cur[[1]], pair[[1]] <> "," <> pair[[2]],
    traceCPhotonPair[cur[[2]], pair[[1]], pair[[2]]],
    If[pair === {"c_after_photon", "c_before_photon"}, {2, 2, 1},
      If[pair === {"c_after_photon", "b_spectator"}, {2, 1, 2}, {1, 2, 2}]]
  ],
  {cur, currentList}, {pair, cPairs}
];
Do[
  allRows = appendRows[
    allRows, "bbar_photon", "open_pair_GG", cur[[1]], pair[[1]] <> "," <> pair[[2]],
    traceBbarPhotonPair[cur[[2]], pair[[1]], pair[[2]]],
    If[pair === {"c_spectator", "b_before_photon"}, {2, 2, 1},
      If[pair === {"c_spectator", "b_after_photon"}, {2, 1, 2}, {1, 2, 2}]]
  ],
  {cur, currentList}, {pair, bPairs}
];

outDir = FileNameJoin[{Directory[], "outputs"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "step3_vec_g2_denominator_reduction.csv"}];
csvField[x_] := "\"" <> StringReplace[ToString[x], "\"" -> "\"\""] <> "\"";

stream = OpenWrite[outFile];
WriteString[stream, "topology,class,current,object,num_d1,num_d2,num_d3,eff_n1,eff_n2,eff_n3,coefficient\n"];
Do[
  WriteString[
    stream,
    StringRiffle[
      Join[csvField /@ row[[1 ;; 4]], ToString /@ row[[5 ;; 10]], {csvField[row[[11]]]}],
      ","
    ],
    "\n"
  ],
  {row, allRows}
];
Close[stream];

Print["Vector tensor norm = ", vNorm];
Print["Wrote vector denominator-reduced G2 rows: ", Length[allRows]];

Quit[];
