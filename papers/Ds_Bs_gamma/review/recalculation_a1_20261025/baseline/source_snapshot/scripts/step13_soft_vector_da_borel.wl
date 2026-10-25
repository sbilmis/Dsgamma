(* ::Package:: *)

(* Complete x-space/Fourier reduction of the two-particle vector photon DA.

   The photon matrix element contains

     x^alpha F_{alpha rho} -> (x.q g_{rho nu} - x_nu q_rho) epsilon^nu.

   We first project the A- and B-current traces onto the E1 tensor, retaining
   x explicitly.  With phase Exp[I (p + a q - k).x], integration by parts
   maps r.x to -I r.d/dk.  The result is then evaluated at k=p+a q without
   replacing p^2 or p.q by hadron masses.

   Run from papers/Ds_Bs_gamma with

     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step13_soft_vector_da_borel.wl
*)

<< FeynCalc`

ClearAll[
  mu, nu, rho, alpha, p, q, pPrime, k, xCoord, mQ, ms, a,
  p2, pq, k2, kp, kq, kx, px, qx, Dden, Dim,
  gammaA, gammaB, gammaP, numQ, traceFor, e1Tensor, e1Projector,
  vectorPhotonX, projectE1, rawA, rawB, scalarRules,
  coeff, dirK, dirP, dirQ, reduceKernel, reducedA, reducedB,
  saddleRules, saddleA, saddleB, ratioBA
];

gammaA[index_] := GA[index] . GA[5];
gammaB[index_] :=
  I FV[pPrime, alpha] DiracSigma[GA[index], GA[alpha]] . GA[5]/(mQ + ms);
gammaP[] := I GA[5];
numQ[mom_] := GS[mom] + mQ;

traceFor[current_] :=
  DiracSimplify[
    DiracTrace[current[mu] . numQ[k] . gammaP[] . GA[rho]],
    DiracSubstitute67 -> True
  ] /. {
    Pair[Momentum[k], Momentum[pPrime]] ->
      Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
    Pair[Momentum[p], Momentum[pPrime]] ->
      Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
    Pair[Momentum[pPrime], Momentum[q]] -> Pair[Momentum[p], Momentum[q]],
    Pair[LorentzIndex[rho], Momentum[pPrime]] ->
      Pair[LorentzIndex[rho], Momentum[p]] +
        Pair[LorentzIndex[rho], Momentum[q]],
    Pair[Momentum[q], Momentum[q]] -> 0
  };

e1Tensor = FV[p, nu] FV[q, mu] - SP[p, q] MT[mu, nu];
e1Projector = e1Tensor/(2 SP[p, q]^2);
vectorPhotonX =
  SP[xCoord, q] MT[rho, nu] - FV[xCoord, nu] FV[q, rho];

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

rawA = projectE1[traceFor[gammaA] vectorPhotonX];
rawB = projectE1[traceFor[gammaB] vectorPhotonX];

scalarRules = {
  Pair[Momentum[k], Momentum[k]] -> k2,
  Pair[Momentum[k], Momentum[p]] -> kp,
  Pair[Momentum[k], Momentum[q]] -> kq,
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[xCoord]] -> kx,
  Pair[Momentum[xCoord], Momentum[k]] -> kx,
  Pair[Momentum[p], Momentum[xCoord]] -> px,
  Pair[Momentum[xCoord], Momentum[p]] -> px,
  Pair[Momentum[q], Momentum[xCoord]] -> qx,
  Pair[Momentum[xCoord], Momentum[q]] -> qx
};
rawA = Expand[rawA /. scalarRules];
rawB = Expand[rawB /. scalarRules];

Dim = 4;
Dden = mQ^2 - k2;
coeff[expr_, var_] := Coefficient[Expand[expr], var];
dirK[f_] := 2 k2 D[f, k2] + kp D[f, kp] + kq D[f, kq];
dirP[f_] := 2 kp D[f, k2] + p2 D[f, kp] + pq D[f, kq];
dirQ[f_] := 2 kq D[f, k2] + pq D[f, kp];

reduceKernel[expr_] := Module[{ck, cp, cq, fk, fp, fq},
  ck = coeff[expr, kx];
  cp = coeff[expr, px];
  cq = coeff[expr, qx];
  fk = ck/Dden;
  fp = cp/Dden;
  fq = cq/Dden;
  FullSimplify[-I (Dim fk + dirK[fk] + dirP[fp] + dirQ[fq])]
];

reducedA = reduceKernel[rawA];
reducedB = reduceKernel[rawB];
saddleRules = {
  kq -> pq,
  kp -> p2 + a pq,
  k2 -> p2 + 2 a pq
};
saddleA = FullSimplify[reducedA /. saddleRules];
saddleB = FullSimplify[reducedB /. saddleRules];
ratioBA = FullSimplify[saddleB/saddleA];

outFile = FileNameJoin[
  {Directory[], "outputs", "step13_soft_vector_da_borel.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Two-particle vector photon-DA Fourier reduction\n"];
WriteString[stream, "================================================\n\n"];
WriteString[stream, "Raw A E1 x-kernel:\n",
  ToString[rawA, InputForm], "\n\n"];
WriteString[stream, "Raw B E1 x-kernel:\n",
  ToString[rawB, InputForm], "\n\n"];
WriteString[stream, "Reduced A before k=p+a q:\n",
  ToString[reducedA, InputForm], "\n\n"];
WriteString[stream, "Reduced B before k=p+a q:\n",
  ToString[reducedB, InputForm], "\n\n"];
WriteString[stream, "Reduced A after k=p+a q:\n",
  ToString[saddleA, InputForm], "\n\n"];
WriteString[stream, "Reduced B after k=p+a q:\n",
  ToString[saddleB, InputForm], "\n\n"];
WriteString[stream, "B/A ratio after k=p+a q:\n",
  ToString[ratioBA, InputForm], "\n"];
Close[stream];

Print["A vector-DA reduced kernel = ", saddleA];
Print["B vector-DA reduced kernel = ", saddleB];
Print["B/A vector-DA ratio = ", ratioBA];
Print["Wrote: ", outFile];
Quit[];
