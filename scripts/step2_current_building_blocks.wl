(* ::Package:: *)

(* Step 2 current building blocks for Ds1 -> Ds gamma.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step2_current_building_blocks.wl

   This file intentionally does not compute the full sum rule yet.  Its job is
   to lock our FeynCalc conventions, define the two Ds1 currents, and perform
   small algebra checks before we build the hard-photon traces.
*)

<< FeynCalc`

Print["Mathematica version: ", $Version];
Print["FeynCalc loaded from: ", $FeynCalcDirectory];

(* We use the mostly-minus metric and the gamma5/epsilon convention stated in
   conventions.md:
     Tr[gamma5 gamma_mu gamma_nu gamma_alpha gamma_beta]
       = -4 I epsilon_{mu nu alpha beta}.
   FeynCalc's standard 4-dimensional Dirac algebra with GA[5] matches this
   convention when the Levi-Civita sign is kept consistently in later scripts.
*)

ClearAll[mu, nu, alpha, beta, rho, p, q, P, k, ms, mQ];

(* Momentum labels:
     p = final pseudoscalar Ds momentum,
     q = real photon momentum,
     P = p + q = initial Ds1 momentum.
   The tensor-axial current contracts sigma_{mu alpha} with P^alpha.
*)

gammaA[mu_] := GA[mu] . GA[5];

gammaB[mu_] :=
  I/(mQ + ms) FV[P, alpha] DiracSigma[GA[mu], GA[alpha]] . GA[5];

gammaP[] := I GA[5];

(* Massive free-propagator numerators for the base hard-emission OPE.
   Denominators and loop integrations are introduced in the next script.
*)

numS[mom_, mass_] := GS[mom] + mass;

(* Basic convention checks.  These are deliberately simple: if these fail, the
   later hard-emission trace would be meaningless.
*)

traceTwoGamma =
  DiracSimplify[DiracTrace[GA[mu] . GA[nu]]];

traceGamma5 =
  DiracSimplify[DiracTrace[GA[5] . GA[mu] . GA[nu] . GA[alpha] . GA[beta]]];

sigmaAntisymCheck =
  DiracSimplify[
    DiracSigma[GA[mu], GA[nu]] + DiracSigma[GA[nu], GA[mu]]
  ];

(* Toy two-point numerators.  These are not yet the radiative correlator; they
   only verify that the A, B, and pseudoscalar vertices can be combined by
   FeynCalc without syntax or convention problems.
*)

toyA =
  DiracSimplify[
    DiracTrace[gammaA[mu] . numS[k, mQ] . gammaP[] . numS[k - p, ms]]
  ];

toyB =
  DiracSimplify[
    DiracTrace[gammaB[mu] . numS[k, mQ] . gammaP[] . numS[k - p, ms]]
  ];

Print["TRACE_TWO_GAMMA = ", traceTwoGamma];
Print["TRACE_GAMMA5_4GAMMA = ", traceGamma5];
Print["SIGMA_ANTISYM_CHECK = ", sigmaAntisymCheck];
Print["TOY_A_NUMERATOR = ", toyA];
Print["TOY_B_NUMERATOR = ", toyB];

Print["Step 2 current-building-block checks completed."];

Quit[];
