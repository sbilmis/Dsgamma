(* ::Package:: *)

(* Fourier-derivative reduction of the x_alpha gamma_beta kernels.

   The E1-projected xgamma kernels from step7 are linear in
     k.x, p.x, q.x
   and multiply the scalar heavy-quark denominator 1/D, D = mQ^2-k^2.

   With phase Exp[I (p + a q - k).x], integration over x turns
     r.x  ->  -I r.d/dk
   after integrating by parts.  For k.x the derivative also acts on k_mu:
     k.x F(k) -> -I d/dk_mu [k_mu F(k)]
              = -I [Dim F + k.d/dk F].

   This script applies those derivative rules to the scalar kernels and then
   substitutes the light-cone routing k = p + a q, q^2=0, while keeping the
   denominator D symbolic.  It is the first actual Fourier/Borel-side reduction
   of the x_alpha G^{alpha beta} gamma_beta term.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/step9_xgamma_derivative_reduction.wl
*)

ClearAll[
  mQ, ms, pq, p2, k2, kp, kq, a, Dim, Dden,
  kx, px, qx, Iu, coeff, dirK, dirP, dirQ, reduceKernel,
  AS, AT1, AT2, AT3, AT4, BS, BSt, BT1, BT2, BT3, BT4,
  saddleRules, onShellRules, ratio
];

Dim = 4;
Dden = mQ^2 - k2;

dirK[f_] := 2 k2 D[f, k2] + kp D[f, kp] + kq D[f, kq];
dirP[f_] := 2 kp D[f, k2] + p2 D[f, kp] + pq D[f, kq];
dirQ[f_] := 2 kq D[f, k2] + pq D[f, kp];

coeff[expr_, var_] := Coefficient[Expand[expr], var];

reduceKernel[expr_] := Module[{ck, cp, cq, fk, fp, fq},
  ck = coeff[expr, kx];
  cp = coeff[expr, px];
  cq = coeff[expr, qx];
  fk = ck/Dden;
  fp = cp/Dden;
  fq = cq/Dden;
  FullSimplify[
    -I (Dim fk + dirK[fk] + dirP[fp] + dirQ[fq])
  ]
];

(* Raw step7 kernels with Pair[Momentum[p],Momentum[q]] -> pq, etc. *)
AS = (-4 I mQ qx)/pq;
AT1 = (-24 mQ qx)/pq;
AT2 = (-8 mQ qx)/pq;
AT3 = (8 mQ qx)/pq;
AT4 = 0;

BS = (-2 I) (kx pq^2 - kq pq px + kq p2 qx + kp pq qx + 2 kq pq qx)/((mQ + ms) pq^2);
BSt = 0; (* epsilon structure handled separately; not part of E1 scalar channel. *)
BT1 = -4 (kx pq^2 - kq pq px + kq p2 qx + 5 kp pq qx + 4 kq pq qx)/((mQ + ms) pq^2);
BT2 = 4 (3 kx pq^2 - 3 kq pq px - kq p2 qx - kp pq qx - 4 kq pq qx)/((mQ + ms) pq^2);
BT3 = 8 (kx pq - kq px + kp qx)/((mQ + ms) pq);
BT4 = -4 (kx pq^2 - kq pq px - kq p2 qx + kp pq qx - 2 kq pq qx)/((mQ + ms) pq^2);

saddleRules = {
  kq -> pq,
  kp -> p2 + a pq,
  k2 -> p2 + 2 a pq
};

onShellRules = {
  mQ^2 -> p2 + 2 a pq
};

raw = {
  {"A_S", AS}, {"A_T1", AT1}, {"A_T2", AT2}, {"A_T3", AT3}, {"A_T4", AT4},
  {"B_S", BS}, {"B_T1", BT1}, {"B_T2", BT2}, {"B_T3", BT3}, {"B_T4", BT4}
};

reduced = raw /. {name_, expr_} :> {name, FullSimplify[reduceKernel[expr]]};
saddleReduced = reduced /. {name_, expr_} :> {name, FullSimplify[expr /. saddleRules]};

get[list_, name_] := SelectFirst[list, #[[1]] == name &][[2]];
ratio[name_] := Module[{aexpr, bexpr},
  aexpr = get[saddleReduced, "A_" <> name];
  bexpr = get[saddleReduced, "B_" <> name];
  If[aexpr === 0,
    If[bexpr === 0, "both zero", "A zero"],
    FullSimplify[bexpr/((mQ/(mQ + ms)) aexpr)]
  ]
];

ratios = {
  {"S", ratio["S"]},
  {"T1", ratio["T1"]},
  {"T2", ratio["T2"]},
  {"T3", ratio["T3"]},
  {"T4", ratio["T4"]}
};

ratioOnShell = ratios /. {name_, expr_} :> {
  name,
  If[StringQ[expr], expr, FullSimplify[expr /. onShellRules]]
};

outFile = FileNameJoin[{Directory[], "outputs", "step9_xgamma_derivative_reduction.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Fourier-derivative reduction of xgamma kernels\n"];
WriteString[stream, "==============================================\n\n"];
WriteString[stream, "Derivative rules use phase Exp[I(p+a q-k).x]. Overall common phases are not physical here.\n\n"];
WriteString[stream, "Reduced kernels before saddle substitution:\n"];
Do[
  WriteString[stream, reduced[[i, 1]], ":\n", ToString[reduced[[i, 2]], InputForm], "\n\n"],
  {i, Length[reduced]}
];
WriteString[stream, "Reduced kernels after k=p+a q:\n"];
Do[
  WriteString[stream, saddleReduced[[i, 1]], ":\n", ToString[saddleReduced[[i, 2]], InputForm], "\n\n"],
  {i, Length[saddleReduced]}
];
WriteString[stream, "Ratio checks B / [(mQ/(mQ+ms)) A] after k=p+a q:\n"];
Do[
  WriteString[stream, ratios[[i, 1]], ": ", ToString[ratios[[i, 2]], InputForm], "\n"],
  {i, Length[ratios]}
];
WriteString[stream, "\nRatio checks with mQ^2=p^2+2 a p.q:\n"];
Do[
  WriteString[stream, ratioOnShell[[i, 1]], ": ", ToString[ratioOnShell[[i, 2]], InputForm], "\n"],
  {i, Length[ratioOnShell]}
];
Close[stream];

Print["Wrote xgamma derivative reduction to: ", outFile];
Do[Print[ratios[[i, 1]], " ratio = ", ratios[[i, 2]]], {i, Length[ratios]}];
Do[Print[ratioOnShell[[i, 1]], " ratio on shell = ", ratioOnShell[[i, 2]]], {i, Length[ratioOnShell]}];

Quit[];
