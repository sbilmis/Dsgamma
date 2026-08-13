(* ::Package:: *)

(* Explicit double-Borel images of the off-shell E1 kernels.

   The heavy-line momentum after the Fourier transform is k = p + a q and

     D_a = mQ^2 - a pPrime^2 - (1-a) p^2 .

   The two Borel transforms localize a at a0.  For the equal-Borel choice,
   a0 = 1/2.  This file performs no pole-mass substitution on the QCD side.
   These are raw double-Borel identities and therefore contain EQ.  In the
   channel-specific final sum rule only the leading twist-2 term is assigned
   DeltaEQ = EQ-E0.  Twist-3, twist-4, and the gluonic three-particle term
   retain EQ.  The separate gauge-completion electromagnetic three-particle
   term follows Rohrwild's IF-type DeltaEQ subtraction.

   Run from papers/Ds_Bs_gamma with

     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step14_explicit_double_borel_forms.wl
*)

ClearAll[
  a, a0, mQ, ms, M2, EQ, F, d,
  b0, bPq, bP2, bLinear,
  vectorA, vectorB,
  xAS, xAT1, xAT2, xAT3, xAT4,
  xBS, xBT1, xBT2, xBT3, xBT4,
  expectedVectorB, expectedXAS, expectedXAT1, expectedXAT2,
  expectedXAT3, checks, outFile, stream
];

d = mQ + ms;
EQ = Exp[-mQ^2/M2];

(* Exact n=2 double-Borel maps.  F and F' are evaluated at a0. *)
b0[g_] := EQ (g /. a -> a0);
bPq[g_] := -M2 EQ/2 (D[g, a] /. a -> a0);
bP2[g_] := EQ (
  mQ^2 (g /. a -> a0) + M2 a0 (D[g, a] /. a -> a0));
bLinear[c0_, cP2_, cPq_, g_] :=
  c0 b0[g] + cP2 bP2[g] + bPq[cPq g];

(* Two-particle vector photon DA. *)
vectorA = -8 mQ b0[F[a]];
vectorB = -8/d bLinear[mQ^2, 0, 1 - a, F[a]];
expectedVectorB =
  -8 EQ/d (
    mQ^2 F[a0] -
    M2/2 (D[(1 - a) F[a], a] /. a -> a0));

(* x_alpha gamma_beta parts of the gluonic/electromagnetic 3-particle
   kernels.  The propagator term has no (slash k + m_Q) numerator.
   Overall Fierz, charge and DA normalizations are common factors and are
   intentionally not included here. *)
xAS = -8 b0[F[a]];
xAT1 = 48 I b0[F[a]];
xAT2 = 16 I b0[F[a]];
xAT3 = -16 I b0[F[a]];
xAT4 = 0;

xBS = 0;
xBT1 = 0;
xBT2 = 0;
xBT3 = 0;
xBT4 = 0;

expectedXAS = -8 EQ F[a0];
expectedXAT1 = 48 I EQ F[a0];
expectedXAT2 = 16 I EQ F[a0];
expectedXAT3 = -16 I EQ F[a0];

checks = FullSimplify /@ {
  vectorB - expectedVectorB,
  xAS - expectedXAS,
  xAT1 - expectedXAT1,
  xAT2 - expectedXAT2,
  xAT3 - expectedXAT3,
  xBS, xBT1, xBT2, xBT3, xBT4
};

outFile = FileNameJoin[
  {Directory[], "outputs", "step14_explicit_double_borel_forms.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Explicit post-double-Borel E1 kernels\n"];
WriteString[stream, "=====================================\n\n"];
WriteString[stream,
  "RAW BOREL MAPS: EQ appears below. Final channel-specific factors are: ",
  "tw2 -> DeltaEQ=EQ-E0; tw3, tw4 and 3p,g -> EQ; ",
  "gauge-completion 3p,em -> DeltaEQ.\n\n"];
WriteString[stream, "B[Integral F(a)/D_a^2] = ",
  ToString[b0[F[a]], InputForm], "\n"];
WriteString[stream, "B[Integral (p.q) F(a)/D_a^2] = ",
  ToString[bPq[F[a]], InputForm], "\n"];
WriteString[stream, "B[Integral p^2 F(a)/D_a^2] = ",
  ToString[bP2[F[a]], InputForm], "\n\n"];
WriteString[stream, "Two-particle vector A:\n",
  ToString[vectorA, InputForm], "\n\n"];
WriteString[stream, "Two-particle vector B:\n",
  ToString[FullSimplify[vectorB], InputForm], "\n\n"];
WriteString[stream, "xgamma A kernels:\n"];
MapThread[
  WriteString[stream, #1, ": ", ToString[#2, InputForm], "\n"] &,
  {{"S", "T1", "T2", "T3", "T4"}, {xAS, xAT1, xAT2, xAT3, xAT4}}
];
WriteString[stream, "\nxgamma B kernels:\n"];
MapThread[
  WriteString[stream, #1, ": ",
    ToString[FullSimplify[#2], InputForm], "\n"] &,
  {{"S", "T1", "T2", "T3", "T4"}, {xBS, xBT1, xBT2, xBT3, xBT4}}
];
WriteString[stream, "\nClosed-form residuals (all must vanish):\n",
  ToString[checks, InputForm], "\n"];
Close[stream];

Print["closed-form residuals = ", checks];
Print["Wrote: ", outFile];
Quit[];
