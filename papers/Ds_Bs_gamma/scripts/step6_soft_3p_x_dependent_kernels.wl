(* ::Package:: *)

(* E1 projection of x-dependent three-particle photon-DA kernels.

   This script keeps x as an external vector and contracts the T3/T4 photon DA
   Lorentz structures with the sigma_{ab} part of the heavy-quark background
   propagator.  It does not yet perform the Fourier/Borel replacement of
   x_alpha/(q.x); the output is the raw projected kernel needed for that step.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step6_soft_3p_x_dependent_kernels.wl
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, beta, rho, lambda, a, b, p, q, P, k, x, mQ, ms,
  gammaA, gammaB, gammaP, numQ, traceSigmaG, e1Tensor, e1Projector,
  projectE1, T3Struct, T4Struct, trAT, trBT, raw
];

gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] :=
  I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaP[] := I GA[5];
numQ[mom_] := GS[mom] + mQ;

traceSigmaG[current_, basis_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      numQ[k] .
      DiracSigma[GA[a], GA[b]] .
      gammaP[] .
      basis
    ],
    DiracSubstitute67 -> True
  ] /. {
    Pair[Momentum[k], Momentum[P]] -> Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
    Pair[Momentum[p], Momentum[P]] -> Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
    Pair[Momentum[P], Momentum[q]] -> Pair[Momentum[p], Momentum[q]],
    Pair[Momentum[q], Momentum[q]] -> 0,
    Pair[LorentzIndex[a], Momentum[P]] -> Pair[LorentzIndex[a], Momentum[p]] + Pair[LorentzIndex[a], Momentum[q]],
    Pair[LorentzIndex[b], Momentum[P]] -> Pair[LorentzIndex[b], Momentum[p]] + Pair[LorentzIndex[b], Momentum[q]],
    Pair[LorentzIndex[rho], Momentum[P]] -> Pair[LorentzIndex[rho], Momentum[p]] + Pair[LorentzIndex[rho], Momentum[q]],
    Pair[LorentzIndex[lambda], Momentum[P]] -> Pair[LorentzIndex[lambda], Momentum[p]] + Pair[LorentzIndex[lambda], Momentum[q]]
  };

e1Tensor = FV[p, nu] FV[q, mu] - SP[p, q] MT[mu, nu];
e1Projector = e1Tensor/(2 SP[p, q]^2);

projectE1[expr_] :=
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
    {mQ, ms}
  ];

(* Appendix structures:
   T3 ~ (q_a x_b - q_b x_a)(eps_rho q_lambda - eps_lambda q_rho)/(q.x)
   T4 ~ (q_rho x_lambda - q_lambda x_rho)(eps_a q_b - eps_b q_a)/(q.x)
   where a,b are G indices and rho,lambda are sigma-bilinear indices. *)
T3Struct =
  ((FV[q, a] FV[x, b] - FV[q, b] FV[x, a]) *
   (MT[rho, nu] FV[q, lambda] - MT[lambda, nu] FV[q, rho])) /
  SP[q, x];

T4Struct =
  ((FV[q, rho] FV[x, lambda] - FV[q, lambda] FV[x, rho]) *
   (MT[a, nu] FV[q, b] - MT[b, nu] FV[q, a])) /
  SP[q, x];

trAT = traceSigmaG[gammaA, DiracSigma[GA[rho], GA[lambda]]];
trBT = traceSigmaG[gammaB, DiracSigma[GA[rho], GA[lambda]]];

kernels = {
  {"A_T3_x_raw", projectE1[trAT T3Struct]},
  {"A_T4_x_raw", projectE1[trAT T4Struct]},
  {"B_T3_x_raw", projectE1[trBT T3Struct]},
  {"B_T4_x_raw", projectE1[trBT T4Struct]}
};

ratioRules = {
  Pair[Momentum[k], Momentum[q]] -> Pair[Momentum[p], Momentum[q]]
};

ratioT3 =
  FullSimplify[
    (kernels[[3, 2]] /. ratioRules) /
    ((mQ/(mQ + ms)) (kernels[[1, 2]] /. ratioRules))
  ];

ratioT4 =
  FullSimplify[
    (kernels[[4, 2]] /. ratioRules) /
    ((mQ/(mQ + ms)) (kernels[[2, 2]] /. ratioRules))
  ];

outFile = FileNameJoin[{Directory[], "outputs", "step6_soft_3p_x_dependent_kernels.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "x-dependent three-particle soft photon-DA E1 kernels\n"];
WriteString[stream, "====================================================\n\n"];
WriteString[stream, "Scope: sigma_{ab} part of S_Q^G; x kept explicit.\n\n"];
Do[
  WriteString[stream, kernels[[i, 1]], ":\n", ToString[kernels[[i, 2]], InputForm], "\n\n"],
  {i, Length[kernels]}
];
WriteString[stream, "Ratio check after k.q -> p.q:\n"];
WriteString[stream, "B_T3 / [(mQ/(mQ+ms)) A_T3] = ", ToString[ratioT3, InputForm], "\n"];
WriteString[stream, "B_T4 / [(mQ/(mQ+ms)) A_T4] = ", ToString[ratioT4, InputForm], "\n"];
Close[stream];

Print["Wrote x-dependent three-particle kernels to: ", outFile];
Do[Print[kernels[[i, 1]], " = ", kernels[[i, 2]]], {i, Length[kernels]}];
Print["T3 ratio check = ", ratioT3];
Print["T4 ratio check = ", ratioT4];

Quit[];
