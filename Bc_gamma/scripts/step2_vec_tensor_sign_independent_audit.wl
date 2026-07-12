(* ::Package:: *)

(* Independent tensor-current sign audit for Bc1 -> Bcstar gamma.

   Purpose:
   - compare the working J_B implementation,
       i DiracSigma[gamma_mu,gamma_alpha] P^alpha gamma_5/(mb+mc),
     with an explicit commutator representation of the same physics convention;
   - identify which simple convention changes flip the radiative J_B sign;
   - keep the vector projector and charge-flow bookkeeping identical to the
     step2 vector-reduction scripts.

   FeynCalc's DiracSigma contains the usual i/2 commutator.  Therefore
       i sigma_{mu alpha}^{FC} = -1/2 [gamma_mu,gamma_alpha].
*)

<< FeynCalc`

ClearAll[
  mu, nu, rho, alpha, p, q, k, mc, mb, ec, ebbar,
  dC1, dC2, dC3, dB1, dB2, dB3, p2, pq,
  gammaV, numS, pAlpha, sigmaComm, gammaB,
  charmInsertion, antibInsertion, vTensor, vNorm, vProjector, vProject,
  reduceWithRules, triangleC, triangleB, totalTriangle, stripCommonI,
  coefficientVector, signRatio, variants, rows, outDir, outTxt, outCsv,
  stream, csvStream
];

gammaV[rho_] := GA[rho];
numS[mom_, mass_] := GS[mom] + mass;
pAlpha := Pair[LorentzIndex[alpha], Momentum[p]] +
  Pair[LorentzIndex[alpha], Momentum[q]];
sigmaComm[mu_, alpha_] := GA[mu] . GA[alpha] - GA[alpha] . GA[mu];

gammaB["working_dirac_sigma", mu_] :=
  I/(mb + mc) pAlpha DiracSigma[GA[mu], GA[alpha]] . GA[5];

gammaB["explicit_commutator", mu_] :=
  -1/(2 (mb + mc)) pAlpha sigmaComm[mu, alpha] . GA[5];

gammaB["opposite_current_phase", mu_] :=
  -I/(mb + mc) pAlpha DiracSigma[GA[mu], GA[alpha]] . GA[5];

gammaB["reversed_sigma_indices", mu_] :=
  I/(mb + mc) pAlpha DiracSigma[GA[alpha], GA[mu]] . GA[5];

gammaB["no_raw_i_dirac_sigma", mu_] :=
  1/(mb + mc) pAlpha DiracSigma[GA[mu], GA[alpha]] . GA[5];

charmInsertion[tag_] :=
  DiracSimplify[
    DiracTrace[
      gammaB[tag, mu] .
      numS[k + q, mc] .
      GA[nu] .
      numS[k, mc] .
      gammaV[rho] .
      numS[k - p, mb]
    ],
    DiracSubstitute67 -> True
  ];

antibInsertion[tag_] :=
  DiracSimplify[
    DiracTrace[
      gammaB[tag, mu] .
      numS[k, mc] .
      gammaV[rho] .
      numS[k - p, mb] .
      GA[nu] .
      numS[k - p + q, mb]
    ],
    DiracSubstitute67 -> True
  ];

vTensor =
  Pair[LorentzIndex[nu], Momentum[p]] *
    Eps[LorentzIndex[mu], LorentzIndex[rho], Momentum[p], Momentum[q]] -
  Pair[Momentum[p], Momentum[q]] *
    Eps[LorentzIndex[mu], LorentzIndex[nu], LorentzIndex[rho], Momentum[p]];

vNorm =
  ScalarProductExpand[
    Contract[vTensor vTensor] /. {
      SP[q, q] -> 0,
      Pair[Momentum[q], Momentum[q]] -> 0,
      Pair[Momentum[p], Momentum[p]] -> p2,
      Pair[Momentum[p], Momentum[q]] -> pq
    }
  ];

vProjector = vTensor/vNorm;

vProject[expr_] :=
  Collect2[
    ScalarProductExpand[
      Contract[expr vProjector] /. {
        SP[q, q] -> 0,
        Pair[Momentum[q], Momentum[q]] -> 0
      }
    ],
    {mc, mb}
  ];

charmRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dC2 + mc^2,
  Pair[Momentum[k], Momentum[q]] -> (dC1 - dC2)/2,
  Pair[Momentum[k], Momentum[p]] -> (dC2 + mc^2 + p2 - mb^2 - dC3)/2
};

antibRules = {
  Pair[Momentum[p], Momentum[p]] -> p2,
  Pair[Momentum[p], Momentum[q]] -> pq,
  Pair[Momentum[k], Momentum[k]] -> dB1 + mc^2,
  Pair[Momentum[k], Momentum[p]] -> (dB1 + mc^2 + p2 - mb^2 - dB2)/2,
  Pair[Momentum[k], Momentum[q]] -> pq + (dB3 - dB2)/2
};

reduceWithRules[expr_, rules_, denoms_] :=
  Collect[Expand[expr /. rules], denoms];

triangleC[tag_] := Collect2[
  reduceWithRules[vProject[charmInsertion[tag]], charmRules, {dC1, dC2, dC3}] /.
    {dC1 -> 0, dC2 -> 0, dC3 -> 0},
  {mc, mb, p2, pq}
];

triangleB[tag_] := Collect2[
  reduceWithRules[vProject[antibInsertion[tag]], antibRules, {dB1, dB2, dB3}] /.
    {dB1 -> 0, dB2 -> 0, dB3 -> 0},
  {mc, mb, p2, pq}
];

totalTriangle[tag_] :=
  Collect2[ec triangleC[tag] + ebbar triangleB[tag], {ec, ebbar, mc, mb, p2, pq}];

stripCommonI[expr_] := Collect2[FullSimplify[expr/I], {ec, ebbar, mc, mb, p2, pq}];

coefficientVector[expr_] := Coefficient[
  Expand[stripCommonI[expr]],
  {ec, ebbar},
  {1, 0}
];

signRatio[expr_, ref_] := FullSimplify[stripCommonI[expr]/stripCommonI[ref]];

variants = {
  "working_dirac_sigma",
  "explicit_commutator",
  "opposite_current_phase",
  "reversed_sigma_indices",
  "no_raw_i_dirac_sigma"
};

ref = totalTriangle["working_dirac_sigma"];
rows = Table[
  With[{tot = totalTriangle[tag]},
    {
      tag,
      ToString[FullSimplify[tot - ref], InputForm],
      ToString[FullSimplify[tot + ref], InputForm],
      ToString[signRatio[tot, ref], InputForm],
      ToString[stripCommonI[tot], InputForm]
    }
  ],
  {tag, variants}
];

outDir = FileNameJoin[{Directory[], "Bc_gamma", "outputs"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outTxt = FileNameJoin[{outDir, "step2_vec_tensor_sign_independent_audit.txt"}];
outCsv = FileNameJoin[{outDir, "step2_vec_tensor_sign_independent_audit.csv"}];

stream = OpenWrite[outTxt];
WriteString[stream, "Bc1 -> Bcstar gamma independent J_B tensor-current sign audit\n"];
WriteString[stream, "================================================================\n\n"];
WriteString[stream, "Projector norm V.V = ", ToString[vNorm, InputForm], "\n"];
WriteString[stream, "Convention reminder: FeynCalc DiracSigma = i/2 [gamma_mu,gamma_alpha].\n"];
WriteString[stream, "Therefore i DiracSigma equals -1/2 times the explicit commutator.\n\n"];
WriteString[stream, "Reference total triangle, working_dirac_sigma:\n"];
WriteString[stream, ToString[ref, InputForm], "\n\n"];
WriteString[stream, "Reference after stripping the common gamma5/epsilon phase I:\n"];
WriteString[stream, ToString[stripCommonI[ref], InputForm], "\n\n"];
Do[
  WriteString[stream, "Variant: ", row[[1]], "\n"];
  WriteString[stream, "  total - reference = ", row[[2]], "\n"];
  WriteString[stream, "  total + reference = ", row[[3]], "\n"];
  WriteString[stream, "  stripped total/reference = ", row[[4]], "\n"];
  WriteString[stream, "  stripped total = ", row[[5]], "\n\n"],
  {row, rows}
];
WriteString[stream,
  "Interpretation:\n",
  "- explicit_commutator agrees with working_dirac_sigma if the tensor current is i sigma_phys P gamma5.\n",
  "- opposite_current_phase and reversed_sigma_indices flip the J_B radiative amplitude.\n",
  "- no_raw_i_dirac_sigma leaves a different overall phase and is not compatible with the real AB/BB convention in the Bc mixing code.\n"
];
Close[stream];

csvStream = OpenWrite[outCsv];
WriteString[csvStream, "variant,total_minus_reference,total_plus_reference,stripped_ratio,stripped_total\n"];
Do[
  WriteString[csvStream, StringRiffle[StringReplace[#, {"," -> ";", "\n" -> " "}] & /@ row, ","], "\n"],
  {row, rows}
];
Close[csvStream];

Print["Wrote independent tensor sign audit to: ", outTxt];
Print["Reference stripped total = ", stripCommonI[ref]];
Do[Print[row[[1]], " ratio = ", row[[4]]], {row, rows}];

Quit[];
