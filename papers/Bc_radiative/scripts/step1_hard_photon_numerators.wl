(* ::Package:: *)

(* Heavy-heavy hard-photon numerator traces for Bc1 -> Bc and Bcstar gamma.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step1_hard_photon_numerators.wl

   Scope:
     - free massive c and b propagator numerators only;
     - one photon vertex on either the c line or the anti-b line;
     - axial basis currents J_A and J_B;
     - pseudoscalar Bc and vector Bc* final currents;
     - no photon DAs and no condensate-expanded propagators.
*)

<< FeynCalc`

ClearAll[
  mu, nu, rho, alpha, p, q, P, k, mc, mb, ec, ebbar,
  gammaA, gammaB, gammaP, gammaV, numS, charmInsertion, antibInsertion
];

(* Momentum labels:
     p = final Bc or Bc* momentum,
     q = photon momentum,
     P = p + q = initial Bc1 momentum.

   nu is the photon Lorentz index.  rho is the Bc* current index.
*)

gammaA[mu_] := GA[mu] . GA[5];

gammaB[mu_] :=
  I/(mb + mc) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];

gammaP[] := I GA[5];
gammaV[rho_] := GA[rho];

numS[mom_, mass_] := GS[mom] + mass;

(* Current convention:
     j_final = bbar Gamma_f c,  J_initial^\dagger closes the same heavy-heavy
     trace.  The two photon insertions are written as c-line and anti-b-line
     insertions.  Overall signs, charges, denominators, color factor, and loop
     integrations are attached in the dispersion step.
*)

charmInsertion[current_, final_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      numS[k + q, mc] .
      GA[nu] .
      numS[k, mc] .
      final .
      numS[k - p, mb]
    ],
    DiracSubstitute67 -> True
  ];

antibInsertion[current_, final_] :=
  DiracSimplify[
    DiracTrace[
      current[mu] .
      numS[k, mc] .
      final .
      numS[k - p, mb] .
      GA[nu] .
      numS[k - p + q, mb]
    ],
    DiracSubstitute67 -> True
  ];

(* Bc final state *)
NAPSCharm = Collect2[charmInsertion[gammaA, gammaP[]], {mc, mb}];
NAPSAntib = Collect2[antibInsertion[gammaA, gammaP[]], {mc, mb}];
NBPSCharm = Collect2[charmInsertion[gammaB, gammaP[]], {mc, mb}];
NBPSAntib = Collect2[antibInsertion[gammaB, gammaP[]], {mc, mb}];

(* Bc* final state *)
NAVCharm = Collect2[charmInsertion[gammaA, gammaV[rho]], {mc, mb}];
NAVAntib = Collect2[antibInsertion[gammaA, gammaV[rho]], {mc, mb}];
NBVCharm = Collect2[charmInsertion[gammaB, gammaV[rho]], {mc, mb}];
NBVAntib = Collect2[antibInsertion[gammaB, gammaV[rho]], {mc, mb}];

totalAPS = Collect2[ec NAPSCharm + ebbar NAPSAntib, {ec, ebbar, mc, mb}];
totalBPS = Collect2[ec NBPSCharm + ebbar NBPSAntib, {ec, ebbar, mc, mb}];
totalAV = Collect2[ec NAVCharm + ebbar NAVAntib, {ec, ebbar, mc, mb}];
totalBV = Collect2[ec NBVCharm + ebbar NBVAntib, {ec, ebbar, mc, mb}];

outDir = FileNameJoin[{Directory[], "outputs"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "step1_hard_photon_numerators.txt"}];

stream = OpenWrite[outFile];
WriteString[stream, "Bc1 -> Bc and Bcstar gamma hard-photon numerator traces\n"];
WriteString[stream, "=================================================\n\n"];
WriteString[stream, "Momentum routing: P = p + q; tensor current uses P^alpha.\n"];
WriteString[stream, "No photon DAs and no condensate-expanded propagators are included.\n\n"];

WriteString[stream, "Bc final state, J_A, c-line:\n", ToString[NAPSCharm, InputForm], "\n\n"];
WriteString[stream, "Bc final state, J_A, anti-b-line:\n", ToString[NAPSAntib, InputForm], "\n\n"];
WriteString[stream, "Bc final state, J_A total:\n", ToString[totalAPS, InputForm], "\n\n"];

WriteString[stream, "Bc final state, J_B, c-line:\n", ToString[NBPSCharm, InputForm], "\n\n"];
WriteString[stream, "Bc final state, J_B, anti-b-line:\n", ToString[NBPSAntib, InputForm], "\n\n"];
WriteString[stream, "Bc final state, J_B total:\n", ToString[totalBPS, InputForm], "\n\n"];

WriteString[stream, "Bc* final state, J_A, c-line:\n", ToString[NAVCharm, InputForm], "\n\n"];
WriteString[stream, "Bc* final state, J_A, anti-b-line:\n", ToString[NAVAntib, InputForm], "\n\n"];
WriteString[stream, "Bc* final state, J_A total:\n", ToString[totalAV, InputForm], "\n\n"];

WriteString[stream, "Bc* final state, J_B, c-line:\n", ToString[NBVCharm, InputForm], "\n\n"];
WriteString[stream, "Bc* final state, J_B, anti-b-line:\n", ToString[NBVAntib, InputForm], "\n\n"];
WriteString[stream, "Bc* final state, J_B total:\n", ToString[totalBV, InputForm], "\n\n"];
Close[stream];

Print["Wrote heavy-heavy numerator traces to: ", outFile];
Print["Term counts: "];
Print["  A -> Bc c-line: ", Length[List @@ Expand[NAPSCharm]]];
Print["  A -> Bc anti-b-line: ", Length[List @@ Expand[NAPSAntib]]];
Print["  B -> Bc c-line: ", Length[List @@ Expand[NBPSCharm]]];
Print["  B -> Bc anti-b-line: ", Length[List @@ Expand[NBPSAntib]]];
Print["  A -> Bc* c-line: ", Length[List @@ Expand[NAVCharm]]];
Print["  A -> Bc* anti-b-line: ", Length[List @@ Expand[NAVAntib]]];
Print["  B -> Bc* c-line: ", Length[List @@ Expand[NBVCharm]]];
Print["  B -> Bc* anti-b-line: ", Length[List @@ Expand[NBVAntib]]];

Quit[];
