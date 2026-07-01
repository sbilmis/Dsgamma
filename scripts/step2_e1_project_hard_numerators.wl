(* ::Package:: *)

(* E1 projection of hard-photon numerator traces.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step2_e1_project_hard_numerators.wl

   This script does not perform the loop integration.  It contracts the open
   axial index mu and photon index nu with the E1 projector

     S_{mu nu}/(2 (p.q)^2),  S_{mu nu}=p_nu q_mu - (p.q) g_{mu nu},

   so that a pure E1 tensor T_{mu nu}=C S_{mu nu} returns C in four dimensions
   with q^2=0.  Terms containing loop momentum k remain as scalar products and
   are reduced in the next step.
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, p, q, P, k, ms, mQ, eS, eQ,
  gammaA, gammaB, gammaP, numS, heavyInsertion, strangeInsertion,
  e1Tensor, e1Projector, e1Project
];

gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] :=
  I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];
gammaP[] := I GA[5];
numS[mom_, mass_] := GS[mom] + mass;

heavyInsertion[current_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      numS[k + q, mQ] .
      GA[nu] .
      numS[k, mQ] .
      gammaP[] .
      numS[k - p, ms]
    ],
    DiracSubstitute67 -> True
  ];

strangeInsertion[current_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      numS[k, mQ] .
      gammaP[] .
      numS[k - p, ms] .
      GA[nu] .
      numS[k - p + q, ms]
    ],
    DiracSubstitute67 -> True
  ];

e1Tensor =
  FV[p, nu] FV[q, mu] - SP[p, q] MT[mu, nu];

e1Projector =
  e1Tensor/(2 SP[p, q]^2);

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

projectorCheck =
  ScalarProductExpand[
    Contract[e1Tensor e1Projector] /. {
      SP[q, q] -> 0,
      Pair[Momentum[q], Momentum[q]] -> 0
    }
  ];

NAHeavy = heavyInsertion[gammaA];
NAStrange = strangeInsertion[gammaA];
NBHeavy = heavyInsertion[gammaB];
NBStrange = strangeInsertion[gammaB];

projAHeavy = e1Project[NAHeavy];
projAStrange = e1Project[NAStrange];
projBHeavy = e1Project[NBHeavy];
projBStrange = e1Project[NBStrange];

projATotal = e1Project[eQ NAHeavy + eS NAStrange];
projBTotal = e1Project[eQ NBHeavy + eS NBStrange];

outFile = FileNameJoin[{Directory[], "outputs", "step2_e1_projected_hard_numerators.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "E1-projected hard-photon numerator traces\n"];
WriteString[stream, "========================================\n\n"];
WriteString[stream, "Projector check S.P = ", ToString[projectorCheck, InputForm], "\n\n"];
WriteString[stream, "A heavy projected numerator:\n", ToString[projAHeavy, InputForm], "\n\n"];
WriteString[stream, "A strange projected numerator:\n", ToString[projAStrange, InputForm], "\n\n"];
WriteString[stream, "A total projected numerator:\n", ToString[projATotal, InputForm], "\n\n"];
WriteString[stream, "B heavy projected numerator:\n", ToString[projBHeavy, InputForm], "\n\n"];
WriteString[stream, "B strange projected numerator:\n", ToString[projBStrange, InputForm], "\n\n"];
WriteString[stream, "B total projected numerator:\n", ToString[projBTotal, InputForm], "\n\n"];
Close[stream];

Print["Projector check = ", projectorCheck];
Print["Wrote E1-projected hard numerators to: ", outFile];
Print["A total projected term count: ", Length[List @@ Expand[projATotal]]];
Print["B total projected term count: ", Length[List @@ Expand[projBTotal]]];

Quit[];
