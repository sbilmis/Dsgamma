(* ::Package:: *)

(* Polarization-summed width matrix for A(1+) -> V(1-) gamma.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script Dsstar_gamma/scripts/width_matrix_1p_to_1m_gamma.wl

   We use explicit rest-frame polarizations:
     P = (mA,0,0,0),
     q = (k,0,0,k),
     p = (EV,0,0,-k),
     k = (mA^2-mV^2)/(2 mA).

   This gives the helicity-summed matrix for an amplitude
     M = c1 S1 + c2 S2 + cEps Seps.
*)

ClearAll[mA, mV, k, ev, dotM, s1, s2, seps, eps4, polA, polV, polG, ampVec, h];

k = (mA^2 - mV^2)/(2 mA);
ev = (mA^2 + mV^2)/(2 mA);

dotM[a_, b_] := a[[1]] b[[1]] - a[[2]] b[[2]] - a[[3]] b[[3]] - a[[4]] b[[4]];

P = {mA, 0, 0, 0};
q = {k, 0, 0, k};
p = {ev, 0, 0, -k};

polA = {
  {0, 1, 0, 0},
  {0, 0, 1, 0},
  {0, 0, 0, 1}
};

polG = {
  {0, 1, 0, 0},
  {0, 0, 1, 0}
};

polV = {
  {0, 1, 0, 0},
  {0, 0, 1, 0},
  {k/mV, 0, 0, -ev/mV}
};

s1[eta_, xi_, eps_] := dotM[eta, eps] dotM[xi, q] - dotM[eta, q] dotM[xi, eps];

s2[eta_, xi_, eps_] :=
  dotM[eta, q] (dotM[xi, eps] dotM[P, q] - dotM[xi, q] dotM[P, eps]);

eps4[a_, b_, c_, d_] := Det[{a, b, c, d}];
seps[eta_, xi_, eps_] := eps4[eta, xi, eps, q];

ampVec[eta_, xi_, eps_] := {s1[eta, xi, eps], s2[eta, xi, eps], seps[eta, xi, eps]};

h = Simplify[
  1/3 Sum[
    Outer[Times, ampVec[eta, xi, eps], ampVec[eta, xi, eps]],
    {eta, polA}, {xi, polV}, {eps, polG}
  ],
  Assumptions -> {mA > mV > 0}
];

widthPrefactor = Simplify[k/(8 Pi mA^2), Assumptions -> {mA > mV > 0}];

outDir = FileNameJoin[{Directory[], "Dsstar_gamma", "outputs"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "width_matrix_1p_to_1m_gamma.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "A(1+) -> V(1-) gamma polarization-summed width matrix\n"];
WriteString[stream, "======================================================\n\n"];
WriteString[stream, "Amplitude basis: M = c1 S1 + c2 S2 + cEps Seps.\n"];
WriteString[stream, "Rows/columns are {S1,S2,Seps}.\n\n"];
WriteString[stream, "(1/3) sum_pol S_i S_j =\n", ToString[h, InputForm], "\n\n"];
WriteString[stream, "Phase-space prefactor |q|/(8 pi mA^2) = ",
  ToString[widthPrefactor, InputForm], "\n\n"];
WriteString[stream, "Therefore Gamma = prefactor * {c1,c2,cEps}.H.{c1,c2,cEps}.\n"];
WriteString[stream, "The off-diagonal even/epsilon entries vanish in this basis.\n"];
Close[stream];

Print["Width matrix H = ", h];
Print["Width prefactor = ", widthPrefactor];
Print["Wrote width-matrix checkpoint to: ", outFile];

Quit[];

