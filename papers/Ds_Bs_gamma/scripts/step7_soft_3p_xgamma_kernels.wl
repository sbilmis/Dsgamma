(* ::Package:: *)

(* E1 projection of the x_alpha gamma_beta part of S_Q^G.

   The second term in the heavy-quark background-gluon propagator is

     u x_alpha G^{alpha beta}(ux) gamma_beta /(m_Q^2-k^2).

   This script computes its E1-projected trace kernels with x kept explicit.
   As in step6, no Fourier/Borel replacement of x_alpha is performed here.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step7_soft_3p_xgamma_kernels.wl
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, beta, rho, lambda, a, b, p, q, P, k, x, mQ, ms,
  gammaA, gammaB, gammaP, numQ, traceXGammaG, e1Tensor, e1Projector,
  projectE1, photonG, SStruct, StStruct, T1Struct, T2Struct, T3Struct,
  T4Struct, trAS, trASt, trAT, trBS, trBSt, trBT
];

gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] :=
  I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaP[] := I GA[5];
numQ[mom_] := GS[mom] + mQ;

traceXGammaG[current_, basis_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      numQ[k] .
      GA[b] .
      gammaP[] .
      basis
    ],
    DiracSubstitute67 -> True
  ] FV[x, a] /. {
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

photonG = FV[q, b] MT[a, nu] - FV[q, a] MT[b, nu];
SStruct = photonG;
StStruct = photonG;

T1Struct =
  (MT[a, lambda] MT[rho, nu] FV[q, b] -
   MT[b, lambda] MT[rho, nu] FV[q, a] -
   MT[a, rho] MT[lambda, nu] FV[q, b] +
   MT[b, rho] MT[lambda, nu] FV[q, a]);

T2Struct =
  (MT[rho, b] MT[a, nu] FV[q, lambda] -
   MT[lambda, b] MT[a, nu] FV[q, rho] -
   MT[rho, a] MT[b, nu] FV[q, lambda] +
   MT[lambda, a] MT[b, nu] FV[q, rho]);

T3Struct =
  ((FV[q, a] FV[x, b] - FV[q, b] FV[x, a]) *
   (MT[rho, nu] FV[q, lambda] - MT[lambda, nu] FV[q, rho])) /
  SP[q, x];

T4Struct =
  ((FV[q, rho] FV[x, lambda] - FV[q, lambda] FV[x, rho]) *
   (MT[a, nu] FV[q, b] - MT[b, nu] FV[q, a])) /
  SP[q, x];

trAS = traceXGammaG[gammaA, 1];
trASt = traceXGammaG[gammaA, I GA[5]];
trAT = traceXGammaG[gammaA, DiracSigma[GA[rho], GA[lambda]]];

trBS = traceXGammaG[gammaB, 1];
trBSt = traceXGammaG[gammaB, I GA[5]];
trBT = traceXGammaG[gammaB, DiracSigma[GA[rho], GA[lambda]]];

kernels = {
  {"A_S_xgamma_raw", projectE1[trAS SStruct]},
  {"A_Stilde_xgamma_raw", projectE1[trASt StStruct]},
  {"A_T1_xgamma_raw", projectE1[trAT T1Struct]},
  {"A_T2_xgamma_raw", projectE1[trAT T2Struct]},
  {"A_T3_xgamma_raw", projectE1[trAT T3Struct]},
  {"A_T4_xgamma_raw", projectE1[trAT T4Struct]},
  {"B_S_xgamma_raw", projectE1[trBS SStruct]},
  {"B_Stilde_xgamma_raw", projectE1[trBSt StStruct]},
  {"B_T1_xgamma_raw", projectE1[trBT T1Struct]},
  {"B_T2_xgamma_raw", projectE1[trBT T2Struct]},
  {"B_T3_xgamma_raw", projectE1[trBT T3Struct]},
  {"B_T4_xgamma_raw", projectE1[trBT T4Struct]}
};

ratioRules = {
  Pair[Momentum[k], Momentum[q]] -> Pair[Momentum[p], Momentum[q]]
};

ratio[labelA_, labelB_] := Module[{aexpr, bexpr},
  aexpr = SelectFirst[kernels, #[[1]] == labelA &][[2]] /. ratioRules;
  bexpr = SelectFirst[kernels, #[[1]] == labelB &][[2]] /. ratioRules;
  If[aexpr === 0,
    If[bexpr === 0, "both zero", "A zero"],
    FullSimplify[bexpr/((mQ/(mQ + ms)) aexpr)]
  ]
];

ratioChecks = {
  {"S", ratio["A_S_xgamma_raw", "B_S_xgamma_raw"]},
  {"Stilde", ratio["A_Stilde_xgamma_raw", "B_Stilde_xgamma_raw"]},
  {"T1", ratio["A_T1_xgamma_raw", "B_T1_xgamma_raw"]},
  {"T2", ratio["A_T2_xgamma_raw", "B_T2_xgamma_raw"]},
  {"T3", ratio["A_T3_xgamma_raw", "B_T3_xgamma_raw"]},
  {"T4", ratio["A_T4_xgamma_raw", "B_T4_xgamma_raw"]}
};

outFile = FileNameJoin[{Directory[], "outputs", "step7_soft_3p_xgamma_kernels.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "x_alpha gamma_beta three-particle soft photon-DA E1 kernels\n"];
WriteString[stream, "============================================================\n\n"];
WriteString[stream, "Scope: x_alpha G^{alpha beta} gamma_beta part of S_Q^G; x kept explicit.\n\n"];
Do[
  WriteString[stream, kernels[[i, 1]], ":\n", ToString[kernels[[i, 2]], InputForm], "\n\n"],
  {i, Length[kernels]}
];
WriteString[stream, "Ratio checks after k.q -> p.q:\n"];
Do[
  WriteString[stream, ratioChecks[[i, 1]], ": ", ToString[ratioChecks[[i, 2]], InputForm], "\n"],
  {i, Length[ratioChecks]}
];
Close[stream];

Print["Wrote xgamma three-particle kernels to: ", outFile];
Do[Print[kernels[[i, 1]], " = ", kernels[[i, 2]]], {i, Length[kernels]}];
Do[Print["ratio ", ratioChecks[[i, 1]], " = ", ratioChecks[[i, 2]]], {i, Length[ratioChecks]}];

Quit[];
