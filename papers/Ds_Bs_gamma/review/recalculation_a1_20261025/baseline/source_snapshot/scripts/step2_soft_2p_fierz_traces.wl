(* ::Package:: *)

(* Two-particle photon-DA Fierz traces for Ds1 -> Ds gamma.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step2_soft_2p_fierz_traces.wl

   The soft two-particle DA contribution after contracting the heavy quark is

     -i/4 Sum_i Tr[ Gamma_mu S_Q(x) i gamma5 Gamma_i ]
       <gamma| sbar(x) Gamma_i s(0) |0>.

   This script computes the Dirac traces for the DA-relevant bilinears:
     scalar, pseudoscalar, vector, axial vector, tensor.

   It is a trace-level artifact; x-space Fourier/Borel transforms are handled
   in later scripts.
*)

<< FeynCalc`

ClearAll[mu, alpha, beta, rho, lambda, p, q, P, k, mQ, ms];

gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] :=
  I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaP[] := I GA[5];
numQ[mom_] := GS[mom] + mQ;

traceFor[current_, basis_] :=
  DiracSimplify[
    DiracTrace[current[mu] . numQ[k] . gammaP[] . basis],
    DiracSubstitute67 -> True
  ] /. {
    Pair[Momentum[k], Momentum[P]] -> Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
    Pair[Momentum[p], Momentum[P]] -> Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
    Pair[Momentum[P], Momentum[q]] -> Pair[Momentum[p], Momentum[q]],
    Pair[Momentum[q], Momentum[q]] -> 0,
    Pair[LorentzIndex[rho], Momentum[P]] -> Pair[LorentzIndex[rho], Momentum[p]] + Pair[LorentzIndex[rho], Momentum[q]],
    Pair[LorentzIndex[lambda], Momentum[P]] -> Pair[LorentzIndex[lambda], Momentum[p]] + Pair[LorentzIndex[lambda], Momentum[q]]
  };

basisList = {
  {"scalar_1", 1},
  {"pseudoscalar_i_gamma5", I GA[5]},
  {"vector_gamma_rho", GA[rho]},
  {"axial_gamma_rho_gamma5", GA[rho] . GA[5]},
  {"tensor_sigma_rho_lambda", DiracSigma[GA[rho], GA[lambda]]}
};

outFile = FileNameJoin[{Directory[], "outputs", "step2_soft_2p_fierz_traces.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Two-particle photon-DA Fierz traces\n"];
WriteString[stream, "===================================\n\n"];
WriteString[stream, "Overall soft prefactor not included: -i/4 and heavy denominator.\n\n"];

Do[
  name = item[[1]];
  basis = item[[2]];
  trA = traceFor[gammaA, basis];
  trB = traceFor[gammaB, basis];
  WriteString[stream, "Basis: ", name, "\n"];
  WriteString[stream, "  A-current trace:\n    ", ToString[trA, InputForm], "\n"];
  WriteString[stream, "  B-current trace:\n    ", ToString[trB, InputForm], "\n\n"];
  ,
  {item, basisList}
];

Close[stream];
Print["Wrote soft two-particle Fierz traces to: ", outFile];
Print["Step 2 soft two-particle trace generation completed."];

Quit[];
