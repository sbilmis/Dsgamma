(* ::Package:: *)

(* Electromagnetic S_gamma and T4^gamma current kernels.

   This script compares the E1 kernels for the axial/tensor 1^+ currents with
   the magnetic D* -> D reference current used by Rohrwild.  The reference
   projection is included only as a normalization check; the axial/tensor
   kernels are the ingredients required by the present paper.

   Run from papers/Ds_Bs_gamma with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step12_em_da_current_kernels.wl
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, a, b, rho, lambda, p, q, P, k, x, mQ, ms,
  gammaV, gammaA, gammaB, gammaP, numQ, traceSigma, traceXGamma,
  e1Tensor, e1Projector, m1Tensor, m1Projector, projectE1, projectM1,
  photonF, T4Struct, tr, kernels
];

gammaV[mu_] := GA[mu];
gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] :=
  I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaP[] := I GA[5];
numQ[mom_] := GS[mom] + mQ;

routingRules = {
  Pair[Momentum[k], Momentum[P]] ->
    Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
  Pair[Momentum[p], Momentum[P]] ->
    Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
  Pair[Momentum[P], Momentum[q]] -> Pair[Momentum[p], Momentum[q]],
  Pair[Momentum[q], Momentum[q]] -> 0,
  Pair[LorentzIndex[a], Momentum[P]] ->
    Pair[LorentzIndex[a], Momentum[p]] + Pair[LorentzIndex[a], Momentum[q]],
  Pair[LorentzIndex[b], Momentum[P]] ->
    Pair[LorentzIndex[b], Momentum[p]] + Pair[LorentzIndex[b], Momentum[q]],
  Pair[LorentzIndex[rho], Momentum[P]] ->
    Pair[LorentzIndex[rho], Momentum[p]] + Pair[LorentzIndex[rho], Momentum[q]],
  Pair[LorentzIndex[lambda], Momentum[P]] ->
    Pair[LorentzIndex[lambda], Momentum[p]] + Pair[LorentzIndex[lambda], Momentum[q]]
};

traceSigma[current_, basis_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] . numQ[k] . DiracSigma[GA[a], GA[b]] .
      gammaP[] . basis
    ],
    DiracSubstitute67 -> True
  ] /. routingRules;

traceXGamma[current_, basis_] :=
  (DiracSimplify[
     DiracTrace[current[mu] . numQ[k] . GA[b] . gammaP[] . basis],
     DiracSubstitute67 -> True
   ] FV[x, a]) /. routingRules;

e1Tensor = FV[p, nu] FV[q, mu] - SP[p, q] MT[mu, nu];
e1Projector = e1Tensor/(2 SP[p, q]^2);

(* Eps[p,q,mu,nu] is the FeynCalc internal form produced by the trace. *)
m1Tensor = Eps[Momentum[p], Momentum[q], LorentzIndex[mu], LorentzIndex[nu]];
m1Projector = m1Tensor/(2 SP[p, q]^2);

projectE1[expr_] :=
  Collect2[ScalarProductExpand[Contract[expr e1Projector] /. routingRules], {mQ, ms}];

projectM1[expr_] :=
  Collect2[ScalarProductExpand[Contract[expr m1Projector] /. routingRules], {mQ, ms}];

(* Open photon-polarization index nu and electromagnetic indices a,b. *)
photonF = FV[q, b] MT[a, nu] - FV[q, a] MT[b, nu];

(* Lorentz tensor multiplying T4^gamma in Rohrwild Eq. (42), with the
   transverse-polarization terms expressed in the gauge-invariant projector. *)
T4Struct =
  ((FV[q, rho] FV[x, lambda] - FV[q, lambda] FV[x, rho]) *
   (MT[a, nu] FV[q, b] - MT[b, nu] FV[q, a])) / SP[q, x];

kernels = {
  {"V_Sg_sigma", projectM1[traceSigma[gammaV, 1] photonF]},
  {"V_T4g_sigma", projectM1[
    traceSigma[gammaV, DiracSigma[GA[rho], GA[lambda]]] T4Struct]},
  {"V_Sg_xgamma", projectM1[traceXGamma[gammaV, 1] photonF]},
  {"V_T4g_xgamma", projectM1[
    traceXGamma[gammaV, DiracSigma[GA[rho], GA[lambda]]] T4Struct]},
  {"A_Sg_sigma", projectE1[traceSigma[gammaA, 1] photonF]},
  {"A_T4g_sigma", projectE1[
    traceSigma[gammaA, DiracSigma[GA[rho], GA[lambda]]] T4Struct]},
  {"A_Sg_xgamma", projectE1[traceXGamma[gammaA, 1] photonF]},
  {"A_T4g_xgamma", projectE1[
    traceXGamma[gammaA, DiracSigma[GA[rho], GA[lambda]]] T4Struct]},
  {"B_Sg_sigma", projectE1[traceSigma[gammaB, 1] photonF]},
  {"B_T4g_sigma", projectE1[
    traceSigma[gammaB, DiracSigma[GA[rho], GA[lambda]]] T4Struct]},
  {"B_Sg_xgamma", projectE1[traceXGamma[gammaB, 1] photonF]},
  {"B_T4g_xgamma", projectE1[
    traceXGamma[gammaB, DiracSigma[GA[rho], GA[lambda]]] T4Struct]}
};

outFile = FileNameJoin[{Directory[], "outputs", "step12_em_da_current_kernels.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Electromagnetic S_gamma and T4^gamma current kernels\n"];
WriteString[stream, "===================================================\n\n"];
WriteString[stream, "V is Rohrwild's D* reference current; A and B are the present E1 basis currents.\n\n"];
Do[
  WriteString[stream, kernels[[i, 1]], ":\n", ToString[kernels[[i, 2]], InputForm], "\n\n"],
  {i, Length[kernels]}
];
Close[stream];

Print["Wrote electromagnetic-DA current kernels to: ", outFile];
Do[Print[kernels[[i, 1]], " = ", kernels[[i, 2]]], {i, Length[kernels]}];

Quit[];
