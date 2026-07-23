(* ::Package:: *)

(* Tensor-current hard contribution: triangle-core reduction.

   Purpose:
     The E1-projected hard numerators for J_B contain k.k, k.p and k.q.
     Before deriving a spectral density we rewrite these scalar products in
     terms of inverse propagator denominators.  Setting those denominators to
     zero isolates the genuine three-denominator triangle core.  The remaining
     denominator-cancellation pieces must be handled separately.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step3_tensor_hard_triangle_reduction.wl
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, p, q, P, k, ms, mQ, eS, eQ,
  dH1, dH2, dH3, dS1, dS2, dS3, p2, pq,
  gammaB, gammaP, numS, heavyInsertion, strangeInsertion,
  e1Tensor, e1Projector, e1Project
];

gammaB[mu_] :=
  I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaP[] := I GA[5];
numS[mom_, mass_] := GS[mom] + mass;

heavyInsertion[] :=
  DiracSimplify[
    DiracTrace[
      gammaB[mu] .
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
      gammaB[mu] .
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
        Pair[Momentum[q], Momentum[q]] -> 0,
        Pair[Momentum[k], Momentum[P]] -> Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
        Pair[Momentum[p], Momentum[P]] -> Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
        Pair[Momentum[P], Momentum[q]] -> Pair[Momentum[p], Momentum[q]]
      }
    ],
    {mQ, ms, eQ, eS}
  ];

projBHeavy = e1Project[heavyInsertion[]];
projBStrange = e1Project[strangeInsertion[]];

(* Heavy-line photon denominators:
     dH1 = (k+q)^2 - mQ^2,
     dH2 = k^2 - mQ^2,
     dH3 = (k-p)^2 - ms^2.
*)
heavyRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dH2 + mQ^2,
  Pair[Momentum[k], Momentum[q]] -> (dH1 - dH2)/2,
  Pair[Momentum[k], Momentum[p]] -> (dH2 + mQ^2 + p2 - ms^2 - dH3)/2
};

(* Strange-line photon denominators:
     dS1 = k^2 - mQ^2,
     dS2 = (k-p)^2 - ms^2,
     dS3 = (k-p+q)^2 - ms^2.
*)
strangeRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dS1 + mQ^2,
  Pair[Momentum[k], Momentum[p]] -> (dS1 + mQ^2 + p2 - ms^2 - dS2)/2,
  Pair[Momentum[k], Momentum[q]] -> pq + (dS3 - dS2)/2
};

reduceHeavy = Collect[Expand[projBHeavy /. heavyRules], {dH1, dH2, dH3}];
reduceStrange = Collect[Expand[projBStrange /. strangeRules], {dS1, dS2, dS3}];

triangleHeavy = Collect2[reduceHeavy /. {dH1 -> 0, dH2 -> 0, dH3 -> 0}, {mQ, ms, p2, pq}];
triangleStrange = Collect2[reduceStrange /. {dS1 -> 0, dS2 -> 0, dS3 -> 0}, {mQ, ms, p2, pq}];

cancellationHeavy = Collect2[Expand[reduceHeavy - triangleHeavy], {dH1, dH2, dH3, mQ, ms}];
cancellationStrange = Collect2[Expand[reduceStrange - triangleStrange], {dS1, dS2, dS3, mQ, ms}];

outFile = FileNameJoin[{Directory[], "outputs", "step3_tensor_hard_triangle_reduction.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Tensor-current hard E1 projection: triangle reduction\n"];
WriteString[stream, "======================================================\n\n"];
WriteString[stream, "Heavy-line triangle core:\n", ToString[triangleHeavy, InputForm], "\n\n"];
WriteString[stream, "Strange-line triangle core:\n", ToString[triangleStrange, InputForm], "\n\n"];
WriteString[stream, "Heavy-line denominator-cancellation part:\n", ToString[cancellationHeavy, InputForm], "\n\n"];
WriteString[stream, "Strange-line denominator-cancellation part:\n", ToString[cancellationStrange, InputForm], "\n\n"];
Close[stream];

Print["Wrote tensor hard triangle reduction to: ", outFile];
Print["Heavy triangle core = ", triangleHeavy];
Print["Strange triangle core = ", triangleStrange];

Quit[];
