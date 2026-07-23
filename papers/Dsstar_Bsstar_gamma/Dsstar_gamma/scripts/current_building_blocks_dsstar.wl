(* ::Package:: *)

(* Current building blocks for Ds1 -> Ds* gamma.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script Dsstar_gamma/scripts/current_building_blocks_dsstar.wl

   Purpose:
     - define the final vector current and the two initial axial basis currents;
     - compute toy two-point trace numerators;
     - compute hard-photon numerator skeletons for photon emission from the
       heavy and strange quark lines.

   Denominators, loop integration, continuum subtraction, and photon-DA terms
   are not introduced here.
*)

<< FeynCalc`

ClearAll[mu, nu, rho, alpha, beta, p, q, P, k, ms, mQ];

gammaV[nu_] := GA[nu];
gammaA[mu_] := GA[mu] . GA[5];
gammaB[mu_] := I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];

numS[mom_, mass_] := GS[mom] + mass;

(* A toy non-radiative trace:
      Tr[ final-vector vertex * heavy numerator * initial-current vertex
          * strange numerator ].
   This is only a syntax and parity/epsilon-sector check.
*)

toyA =
  DiracSimplify[
    DiracTrace[
      gammaV[nu] . numS[k + p, mQ] . gammaA[mu] . numS[k, ms]
    ],
    DiracTraceEvaluate -> True
  ];

toyB =
  DiracSimplify[
    DiracTrace[
      gammaV[nu] . numS[k + p, mQ] . gammaB[mu] . numS[k, ms]
    ],
    DiracTraceEvaluate -> True
  ];

(* Hard photon emission skeletons.

   rho is the external photon vertex index.  These traces are not yet
   projected onto S1/S2.  The full hard contribution has denominators from the
   two adjacent propagators and the spectator propagator.
*)

heavyA =
  DiracSimplify[
    DiracTrace[
      gammaV[nu] . numS[k + p, mQ] . GA[rho] .
      numS[k + p + q, mQ] . gammaA[mu] . numS[k, ms]
    ],
    DiracTraceEvaluate -> True
  ];

heavyB =
  DiracSimplify[
    DiracTrace[
      gammaV[nu] . numS[k + p, mQ] . GA[rho] .
      numS[k + p + q, mQ] . gammaB[mu] . numS[k, ms]
    ],
    DiracTraceEvaluate -> True
  ];

strangeA =
  DiracSimplify[
    DiracTrace[
      gammaV[nu] . numS[k + p, mQ] . gammaA[mu] .
      numS[k - q, ms] . GA[rho] . numS[k, ms]
    ],
    DiracTraceEvaluate -> True
  ];

strangeB =
  DiracSimplify[
    DiracTrace[
      gammaV[nu] . numS[k + p, mQ] . gammaB[mu] .
      numS[k - q, ms] . GA[rho] . numS[k, ms]
    ],
    DiracTraceEvaluate -> True
  ];

outDir = FileNameJoin[{Directory[], "Dsstar_gamma", "outputs"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "current_building_blocks_dsstar.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Ds1 -> Ds* gamma current-building-block traces\n"];
WriteString[stream, "================================================\n\n"];
WriteString[stream, "Final vector current: gamma_nu\n"];
WriteString[stream, "Initial A current: gamma_mu gamma5\n"];
WriteString[stream, "Initial B current: i sigma_{mu alpha} P^alpha gamma5/(mQ+ms)\n\n"];
WriteString[stream, "TOY_A = ", ToString[toyA, InputForm], "\n\n"];
WriteString[stream, "TOY_B = ", ToString[toyB, InputForm], "\n\n"];
WriteString[stream, "HEAVY_A = ", ToString[heavyA, InputForm], "\n\n"];
WriteString[stream, "HEAVY_B = ", ToString[heavyB, InputForm], "\n\n"];
WriteString[stream, "STRANGE_A = ", ToString[strangeA, InputForm], "\n\n"];
WriteString[stream, "STRANGE_B = ", ToString[strangeB, InputForm], "\n"];
Close[stream];

Print["Wrote current-building-block traces to: ", outFile];
Print["Toy A trace length: ", StringLength[ToString[toyA, InputForm]]];
Print["Toy B trace length: ", StringLength[ToString[toyB, InputForm]]];
Print["Heavy A trace length: ", StringLength[ToString[heavyA, InputForm]]];
Print["Heavy B trace length: ", StringLength[ToString[heavyB, InputForm]]];
Print["Strange A trace length: ", StringLength[ToString[strangeA, InputForm]]];
Print["Strange B trace length: ", StringLength[ToString[strangeB, InputForm]]];

Quit[];

