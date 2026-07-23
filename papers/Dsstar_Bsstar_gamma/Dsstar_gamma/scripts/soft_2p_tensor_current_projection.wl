(* ::Package:: *)

(* Soft two-particle photon-DA projection for Ds1 -> Ds* gamma.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script Dsstar_gamma/scripts/soft_2p_tensor_current_projection.wl

   The target invariant amplitude is the epsilon-sector structure

     eps_{mu nu rho sigma} e^rho q^sigma,

   where mu is the initial axial-current index, nu is the final vector-current
   index, and rho is the photon polarization index.
*)

<< FeynCalc`

ClearAll[mu, nu, rho, lambda, alpha, beta, p, q, P, k, mQ, ms];

gammaV[nu_] := GA[nu];
gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] := I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
numQ[mom_] := GS[mom] + mQ;

replaceP = {
  Pair[Momentum[k], Momentum[P]] -> Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
  Pair[Momentum[p], Momentum[P]] -> Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
  Pair[Momentum[P], Momentum[q]] -> Pair[Momentum[p], Momentum[q]],
  Pair[Momentum[q], Momentum[q]] -> 0,
  Pair[LorentzIndex[rho], Momentum[P]] -> Pair[LorentzIndex[rho], Momentum[p]] + Pair[LorentzIndex[rho], Momentum[q]],
  Pair[LorentzIndex[lambda], Momentum[P]] -> Pair[LorentzIndex[lambda], Momentum[p]] + Pair[LorentzIndex[lambda], Momentum[q]],
  Pair[LorentzIndex[beta], Momentum[P]] -> Pair[LorentzIndex[beta], Momentum[p]] + Pair[LorentzIndex[beta], Momentum[q]]
};

traceFor[current_, basis_] :=
  DiracSimplify[
    DiracTrace[
      gammaV[nu] . numQ[k] . current[mu] . basis
    ],
    DiracTraceEvaluate -> True,
    DiracSubstitute67 -> True
  ] /. replaceP;

(* Photon-DA Lorentz structures with open photon-polarization index beta. *)
tensorPhoton = MT[rho, beta] FV[q, lambda] - MT[lambda, beta] FV[q, rho];
vectorPhoton = MT[rho, beta];
axialPhoton = Eps[LorentzIndex[rho], LorentzIndex[beta], Momentum[q], Momentum[k]];

target = Eps[LorentzIndex[mu], LorentzIndex[nu], LorentzIndex[beta], Momentum[q]];
dualTarget = Eps[LorentzIndex[mu], LorentzIndex[nu], LorentzIndex[beta], Momentum[p]];
targetNorm =
  Contract[target dualTarget] /. {
    Pair[Momentum[q], Momentum[q]] -> 0,
    SP[q, q] -> 0
  };

projectEps[expr_] :=
  Collect2[
    ScalarProductExpand[
      Contract[expr dualTarget/targetNorm] /. {
        Pair[Momentum[q], Momentum[q]] -> 0,
        SP[q, q] -> 0
      }
    ],
    {mQ, ms, Pair[Momentum[k], Momentum[q]], Pair[Momentum[k], Momentum[p]], Pair[Momentum[p], Momentum[q]]}
  ];

trATensor = traceFor[gammaA, DiracSigma[GA[rho], GA[lambda]]];
trBTensor = traceFor[gammaB, DiracSigma[GA[rho], GA[lambda]]];

trAVector = traceFor[gammaA, GA[rho]];
trBVector = traceFor[gammaB, GA[rho]];

trAAxial = traceFor[gammaA, GA[rho] . GA[5]];
trBAxial = traceFor[gammaB, GA[rho] . GA[5]];

projATensor = projectEps[trATensor tensorPhoton];
projBTensor = projectEps[trBTensor tensorPhoton];
projAVector = projectEps[trAVector vectorPhoton];
projBVector = projectEps[trBVector vectorPhoton];
projAAxial = projectEps[trAAxial axialPhoton];
projBAxial = projectEps[trBAxial axialPhoton];

ratioTensor = FullSimplify[projBTensor/projATensor];
ratioVector = FullSimplify[projBVector/projAVector];
ratioAxial = FullSimplify[projBAxial/projAAxial];

outDir = FileNameJoin[{Directory[], "Dsstar_gamma", "outputs"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "soft_2p_tensor_current_projection.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Soft two-particle DA epsilon projection for Ds1 -> Ds* gamma\n"];
WriteString[stream, "============================================================\n\n"];
WriteString[stream, "Target tensor: Eps[mu,nu,beta,q]\n"];
WriteString[stream, "Target norm:\n", ToString[targetNorm, InputForm], "\n\n"];
WriteString[stream, "Tensor DA A projection:\n", ToString[projATensor, InputForm], "\n\n"];
WriteString[stream, "Tensor DA B projection:\n", ToString[projBTensor, InputForm], "\n\n"];
WriteString[stream, "Tensor DA B/A ratio:\n", ToString[ratioTensor, InputForm], "\n\n"];
WriteString[stream, "Vector DA A projection:\n", ToString[projAVector, InputForm], "\n\n"];
WriteString[stream, "Vector DA B projection:\n", ToString[projBVector, InputForm], "\n\n"];
WriteString[stream, "Vector DA B/A ratio:\n", ToString[ratioVector, InputForm], "\n\n"];
WriteString[stream, "Axial DA A projection:\n", ToString[projAAxial, InputForm], "\n\n"];
WriteString[stream, "Axial DA B projection:\n", ToString[projBAxial, InputForm], "\n\n"];
WriteString[stream, "Axial DA B/A ratio:\n", ToString[ratioAxial, InputForm], "\n"];
Close[stream];

Print["Wrote soft two-particle tensor-current projection to: ", outFile];
Print["Tensor DA B/A ratio = ", ratioTensor];
Print["Vector DA B/A ratio = ", ratioVector];
Print["Axial DA B/A ratio = ", ratioAxial];

Quit[];
