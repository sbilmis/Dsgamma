(* ::Package:: *)

(* Axial-current hard contribution: triangle-core reduction.

   This mirrors step3_tensor_hard_triangle_reduction.wl for J_A.  Since the
   axial hard contribution has a published spectral density, this reduction is
   a calibration tool: if the denominator-cancellation pieces are nontrivial
   for J_A, then triangle-core-only estimates are not sufficient for J_B.
*)

<< FeynCalc`

ClearAll[
  mu, nu, p, q, k, ms, mQ, dH1, dH2, dH3, dS1, dS2, dS3, p2, pq,
  gammaA, gammaP, numS, heavyInsertion, strangeInsertion, e1Project
];

gammaA[mu_] := GA[mu] . GA[5];
gammaP[] := I GA[5];
numS[mom_, mass_] := GS[mom] + mass;

heavyInsertion[] :=
  DiracSimplify[
    DiracTrace[
      gammaA[mu] .
      numS[k + q, mQ] .
      GA[nu] .
      numS[k, mQ] .
      gammaP[] .
      numS[k - p, ms]
    ],
    DiracSubstitute67 -> True
  ];

strangeInsertion[] :=
  DiracSimplify[
    DiracTrace[
      gammaA[mu] .
      numS[k, mQ] .
      gammaP[] .
      numS[k - p, ms] .
      GA[nu] .
      numS[k - p + q, ms]
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
        Pair[Momentum[q], Momentum[q]] -> 0
      }
    ],
    {mQ, ms}
  ];

projAHeavy = e1Project[heavyInsertion[]];
projAStrange = e1Project[strangeInsertion[]];

heavyRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dH2 + mQ^2,
  Pair[Momentum[k], Momentum[q]] -> (dH1 - dH2)/2,
  Pair[Momentum[k], Momentum[p]] -> (dH2 + mQ^2 + p2 - ms^2 - dH3)/2
};

strangeRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dS1 + mQ^2,
  Pair[Momentum[k], Momentum[p]] -> (dS1 + mQ^2 + p2 - ms^2 - dS2)/2,
  Pair[Momentum[k], Momentum[q]] -> pq + (dS3 - dS2)/2
};

reduceHeavy = Collect[Expand[projAHeavy /. heavyRules], {dH1, dH2, dH3}];
reduceStrange = Collect[Expand[projAStrange /. strangeRules], {dS1, dS2, dS3}];

triangleHeavy = Collect2[reduceHeavy /. {dH1 -> 0, dH2 -> 0, dH3 -> 0}, {mQ, ms, p2, pq}];
triangleStrange = Collect2[reduceStrange /. {dS1 -> 0, dS2 -> 0, dS3 -> 0}, {mQ, ms, p2, pq}];

cancellationHeavy = Collect2[Expand[reduceHeavy - triangleHeavy], {dH1, dH2, dH3, mQ, ms}];
cancellationStrange = Collect2[Expand[reduceStrange - triangleStrange], {dS1, dS2, dS3, mQ, ms}];

outFile = FileNameJoin[{Directory[], "outputs", "step3_axial_hard_triangle_reduction.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Axial-current hard E1 projection: triangle reduction\n"];
WriteString[stream, "====================================================\n\n"];
WriteString[stream, "Heavy-line triangle core:\n", ToString[triangleHeavy, InputForm], "\n\n"];
WriteString[stream, "Strange-line triangle core:\n", ToString[triangleStrange, InputForm], "\n\n"];
WriteString[stream, "Heavy-line denominator-cancellation part:\n", ToString[cancellationHeavy, InputForm], "\n\n"];
WriteString[stream, "Strange-line denominator-cancellation part:\n", ToString[cancellationStrange, InputForm], "\n\n"];
Close[stream];

Print["Wrote axial hard triangle reduction to: ", outFile];
Print["Heavy triangle core = ", triangleHeavy];
Print["Strange triangle core = ", triangleStrange];

Quit[];
