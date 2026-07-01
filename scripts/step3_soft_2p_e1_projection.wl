(* ::Package:: *)

(* E1 projection of two-particle photon-DA trace structures.

   This script contracts the Fierz traces with the Lorentz structures appearing
   in the two-particle photon DA matrix elements and projects onto the E1 tensor.

   It is most useful for checking the tensor-current leading-twist contribution
   from <gamma|sbar sigma_{rho lambda} s|0>.
*)

<< FeynCalc`

ClearAll[mu, nu, rho, lambda, alpha, p, q, P, k, mQ, ms];

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

e1Tensor = FV[p, nu] FV[q, mu] - SP[p, q] MT[mu, nu];
e1Projector = e1Tensor/(2 SP[p, q]^2);

projectE1[expr_] :=
  Collect2[
    ScalarProductExpand[
      Contract[expr e1Projector] /. {
        SP[q, q] -> 0,
        Pair[Momentum[q], Momentum[q]] -> 0
      }
    ],
    {mQ, ms}
  ];

(* Photon DA Lorentz structures with open photon-polarization index nu. *)
tensorPhoton = MT[rho, nu] FV[q, lambda] - MT[lambda, nu] FV[q, rho];
vectorPhoton = MT[rho, nu];

trATensor = traceFor[gammaA, DiracSigma[GA[rho], GA[lambda]]];
trBTensor = traceFor[gammaB, DiracSigma[GA[rho], GA[lambda]]];

trAVector = traceFor[gammaA, GA[rho]];
trBVector = traceFor[gammaB, GA[rho]];

projATensorDA = projectE1[trATensor tensorPhoton];
projBTensorDA = projectE1[trBTensor tensorPhoton];
projAVectorDA = projectE1[trAVector vectorPhoton];
projBVectorDA = projectE1[trBVector vectorPhoton];

outFile = FileNameJoin[{Directory[], "outputs", "step3_soft_2p_e1_projection.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "E1 projection of two-particle soft photon-DA traces\n"];
WriteString[stream, "===================================================\n\n"];
WriteString[stream, "Tensor DA contraction, A current:\n", ToString[projATensorDA, InputForm], "\n\n"];
WriteString[stream, "Tensor DA contraction, B current:\n", ToString[projBTensorDA, InputForm], "\n\n"];
WriteString[stream, "Vector DA contraction, A current:\n", ToString[projAVectorDA, InputForm], "\n\n"];
WriteString[stream, "Vector DA contraction, B current:\n", ToString[projBVectorDA, InputForm], "\n\n"];
Close[stream];

Print["Wrote soft two-particle E1 projections to: ", outFile];
Print["B tensor-DA E1 coefficient = ", projBTensorDA];
Print["B vector-DA E1 coefficient = ", projBVectorDA];

Quit[];
