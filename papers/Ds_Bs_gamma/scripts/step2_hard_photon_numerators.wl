(* ::Package:: *)

(* Hard-photon numerator traces for Ds1 -> Ds gamma.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step2_hard_photon_numerators.wl

   Scope of this file:
     - free massive propagator numerators only;
     - one photon vertex on either the heavy line or the strange line;
     - separate axial-current and tensor-current numerators;
     - no condensate-expanded propagators.

   The output is written to outputs/step2_hard_photon_numerators.txt for review.
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, p, q, P, k, ms, mQ, eS, eQ,
  gammaA, gammaB, gammaP, numS, heavyInsertion, strangeInsertion
];

(* Momentum labels:
     p = final Ds momentum,
     q = photon momentum,
     P = p + q = initial Ds1 momentum.

   nu is the photon Lorentz index.  It will later be contracted with the photon
   polarization epsilon_nu.  The hard amplitude is proportional to
     eQ * N_heavy^nu + eS * N_strange^nu.
*)

gammaA[mu_] := GA[mu] . GA[5];

gammaB[mu_] :=
  I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];

gammaP[] := I GA[5];

numS[mom_, mass_] := GS[mom] + mass;

(* Momentum routing for this numerator-level skeleton:

   No-photon trace:
     Tr[ Gamma_Ds1,mu S_Q(k) Gamma_Ds S_s(k-p) ].

   Heavy-line photon insertion:
     Tr[ Gamma_Ds1,mu S_Q(k+q) gamma_nu S_Q(k) Gamma_Ds S_s(k-p) ].

   Strange-line photon insertion:
     Tr[ Gamma_Ds1,mu S_Q(k) Gamma_Ds S_s(k-p) gamma_nu S_s(k-p+q) ].

   Overall factors of i, denominators, color, charges, and loop integration are
   attached in the next step.  Here we keep only the Dirac numerator.
*)

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

NAHeavy = Collect2[heavyInsertion[gammaA], {mQ, ms}];
NAStrange = Collect2[strangeInsertion[gammaA], {mQ, ms}];

NBHeavy = Collect2[heavyInsertion[gammaB], {mQ, ms}];
NBStrange = Collect2[strangeInsertion[gammaB], {mQ, ms}];

totalA = Collect2[eQ NAHeavy + eS NAStrange, {eQ, eS, mQ, ms}];
totalB = Collect2[eQ NBHeavy + eS NBStrange, {eQ, eS, mQ, ms}];

outFile = FileNameJoin[{Directory[], "outputs", "step2_hard_photon_numerators.txt"}];

stream = OpenWrite[outFile];
WriteString[stream, "Step 2 hard-photon numerator traces\n"];
WriteString[stream, "===================================\n\n"];
WriteString[stream, "Momentum routing: P = p + q; tensor current uses P^alpha.\n"];
WriteString[stream, "No condensate-expanded propagators are used.\n\n"];

WriteString[stream, "N_A_heavy(mu,nu) =\n", ToString[NAHeavy, InputForm], "\n\n"];
WriteString[stream, "N_A_strange(mu,nu) =\n", ToString[NAStrange, InputForm], "\n\n"];
WriteString[stream, "N_A_total = eQ N_A_heavy + eS N_A_strange =\n",
  ToString[totalA, InputForm], "\n\n"];

WriteString[stream, "N_B_heavy(mu,nu) =\n", ToString[NBHeavy, InputForm], "\n\n"];
WriteString[stream, "N_B_strange(mu,nu) =\n", ToString[NBStrange, InputForm], "\n\n"];
WriteString[stream, "N_B_total = eQ N_B_heavy + eS N_B_strange =\n",
  ToString[totalB, InputForm], "\n\n"];
Close[stream];

Print["Wrote hard-photon numerator traces to: ", outFile];
Print["A-heavy term count: ", Length[List @@ Expand[NAHeavy]]];
Print["A-strange term count: ", Length[List @@ Expand[NAStrange]]];
Print["B-heavy term count: ", Length[List @@ Expand[NBHeavy]]];
Print["B-strange term count: ", Length[List @@ Expand[NBStrange]]];
Print["Step 2 hard-photon numerator generation completed."];

Quit[];
