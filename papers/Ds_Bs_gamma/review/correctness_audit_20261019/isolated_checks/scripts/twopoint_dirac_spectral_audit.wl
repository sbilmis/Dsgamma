(* ::Package:: *)

(* Dirac-trace audit for the perturbative AA/AB/BB two-point matrix.

   The normalized basis currents are

     J_A^mu = qbar gamma^mu gamma5 Q,
     J_B^mu = i/(mQ+mq) qbar sigma^(mu alpha) p_alpha gamma5 Q.

   The script projects the cut one-loop numerator onto the transverse
   structure g_(mu nu)-p_mu p_nu/p^2.  It deliberately leaves the universal
   two-body phase-space factor outside the trace so that the result can be
   checked against the published Colangelo axial spectral density.
*)

<< FeynCalc`

ClearAll[mu, nu, alpha, p, k, s, mQ, mq, lam, gammaA, gammaB,
  numerator, transverseNumerator, localNumerator, cutRules];

scriptDir = DirectoryName[$InputFileName];
paperDir = DirectoryName[scriptDir];
outDir = FileNameJoin[{paperDir, "outputs"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "twopoint_dirac_spectral_audit.txt"}];

gammaA[index_] := GA[index] . GA[5];
gammaB[index_] := Module[{a = Unique["a"]},
  I/(mQ + mq) FV[p, a] DiracSigma[GA[index], GA[a]] . GA[5]
];

(* The same momentum-space matrix is used for the conjugate flavor current.
   This convention is the one already used in the radiative-current audit. *)
numerator[left_, right_] :=
  DiracSimplify[
    DiracTrace[
      left[mu] . (GS[k] + mQ) . right[nu] . (GS[k - p] + mq)
    ],
    DiracSubstitute67 -> True
  ];

cutRules = {
  Pair[Momentum[p], Momentum[p]] -> s,
  Pair[Momentum[k], Momentum[k]] -> mQ^2,
  Pair[Momentum[k], Momentum[p]] -> (s + mQ^2 - mq^2)/2
};

transverseNumerator[left_, right_] :=
  Factor[
    ScalarProductExpand[
      Contract[
        (MT[mu, nu] - FV[p, mu] FV[p, nu]/SP[p, p])
          numerator[left, right]/3
      ]
    ] /. cutRules
  ];

nAA = transverseNumerator[gammaA, gammaA];
nAB = transverseNumerator[gammaA, gammaB];
nBA = transverseNumerator[gammaB, gammaA];
nBB = transverseNumerator[gammaB, gammaB];

localNumerator[left_, right_] :=
  Factor[
    ScalarProductExpand[
      Contract[
        (MT[mu, nu] - FV[p, mu] FV[p, nu]/SP[p, p])
          DiracSimplify[
            DiracTrace[left[mu] . (GS[p] + mQ) . right[nu]],
            DiracSubstitute67 -> True
          ]/3
      ]
    ] /. Pair[Momentum[p], Momentum[p]] -> s
  ];

qAA = localNumerator[gammaA, gammaA];
qAB = localNumerator[gammaA, gammaB];
qBA = localNumerator[gammaB, gammaA];
qBB = localNumerator[gammaB, gammaB];

lam = (s - (mQ + mq)^2) (s - (mQ - mq)^2);

(* The common cut factor is fixed below by matching rho_AA to the independent
   Colangelo expression.  The candidate normalization is Nc sqrt(lam)/(16
   Pi^2 s), up to the overall Wick/metric sign printed in the audit. *)
rhoCandidate[n_] := Factor[3 Sqrt[lam] n/(16 Pi^2 s)];

rhoAAColangelo =
  Sqrt[lam]/(8 Pi^2) *
    (2 - (mQ^2 + mq^2 + 6 mQ mq)/s - (mQ^2 - mq^2)^2/s^2);

lines = {
  "Normalized-current perturbative two-point Dirac audit",
  "======================================================",
  "",
  "nAA = " <> ToString[nAA, InputForm],
  "nAB = " <> ToString[nAB, InputForm],
  "nBA = " <> ToString[nBA, InputForm],
  "nBB = " <> ToString[nBB, InputForm],
  "",
  "AB symmetry check nAB-nBA = " <> ToString[Factor[nAB - nBA], InputForm],
  "",
  "candidate rhoAA = " <> ToString[rhoCandidate[nAA], InputForm],
  "Colangelo rhoAA = " <> ToString[rhoAAColangelo, InputForm],
  "candidate/Colangelo = " <>
    ToString[FullSimplify[rhoCandidate[nAA]/rhoAAColangelo], InputForm],
  "",
  "candidate rhoAB = " <> ToString[rhoCandidate[nAB], InputForm],
  "candidate rhoBB = " <> ToString[rhoCandidate[nBB], InputForm],
  "",
  "local heavy-pole transverse traces",
  "qAA = " <> ToString[qAA, InputForm],
  "qAB = " <> ToString[qAB, InputForm],
  "qBA = " <> ToString[qBA, InputForm],
  "qBB = " <> ToString[qBB, InputForm]
};

Export[outFile, StringRiffle[lines, "\n"] <> "\n", "String"];
Print[StringRiffle[lines, "\n"]];
