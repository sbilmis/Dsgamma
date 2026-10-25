(* Independent tree-level, dimension-five Wilson-coefficient calculation.
   Light-mass-zero Wilson coefficients; d=mQ+ms is retained as the specified
   current normalization. Two topologies: light bilocal x^2 and heavy-line G.
   Everything is written in this review directory, never production outputs. *)
ClearAll["Global`*"];
auditDirectory = DirectoryName[$InputFileName];
zero2 = ConstantArray[0, {2, 2}]; eye2 = IdentityMatrix[2]; eye4 = IdentityMatrix[4];
pauli = {{{0, 1}, {1, 0}}, {{0, -I}, {I, 0}}, {{1, 0}, {0, -1}}};
g[0] = ArrayFlatten[{{eye2, zero2}, {zero2, -eye2}}];
Do[g[j] = ArrayFlatten[{{zero2, pauli[[j]]}, {-pauli[[j]], zero2}}], {j, 1, 3}];
g5 = I g[0].g[1].g[2].g[3];
eta = {1, -1, -1, -1};
sig[a_, b_] := I (g[a].g[b] - g[b].g[a])/2;
slash = g[0] k0 - g[1] k1 - g[2] k2 - g[3] k3;
heavy = (slash + m eye4)/(k0^2-k1^2-k2^2-k3^2-m^2);
lapHeavy = (D[heavy, {k0,2}] - D[heavy, {k1,2}]
  - D[heavy, {k2,2}] - D[heavy, {k3,2}]);
heavyRest = (p0 g[0]+m eye4)/(p0^2-m^2);
numRest = p0 g[0]+m eye4;
currents = <|"A" -> g[1].g5, "B" -> I p0 sig[1,0].g5/d,
             "P" -> I g5, "V" -> g[1]|>;
channels = <|"AA" -> {"A","A"}, "AB" -> {"A","B"},
  "BA" -> {"B","A"}, "BB" -> {"B","B"},
  "P" -> {"P","P"}, "V" -> {"V","V"}|>;
toS[expr_] := FullSimplify[expr /. p0 -> Sqrt[s], Assumptions -> s>0 && m>0 && d>0];
borel[expr_] := FullSimplify[Exp[-m^2/K] Sum[
  SeriesCoefficient[expr /. s -> m^2-z, {z,0,-n}]/(Factorial[n-1] K^(n-1)),
  {n,1,4}], Assumptions -> m>0 && d>0 && K>0];

(* <g q_i bar q_j G^a_ab> = -C5/192 (t^a)_ij sigma_ab.
   Tr_color sum_a t^a t^a = 4 and Tr_spin sigma^ab sigma_ab = 48.
   The heavy background propagator coefficient is -g t.G {sigma,N}/4D^2.
   Therefore the contracted heavy-line topology has coefficient +C5/192.
   Color factor in the light-bilocal topology is 3: Fourier x^2 -> -box_k
   converts -x^2 C5/192 into +C5 box_k/64. *)
sigmaNorm = Sum[eta[[a+1]] eta[[b+1]] Tr[sig[a,b].sig[a,b]], {a,0,3},{b,0,3}];
If[sigmaNorm =!= 48, Print["Bad mixed-condensate projector normalization"]; Exit[1]];

results = Association@Table[
  Module[{left=currents[channels[ch][[1]]], right=currents[channels[ch][[2]]],
    q3, light5, background5, full5},
    q3 = toS[-Tr[left.heavyRest.right]/4];
    light5 = toS[(Tr[left.lapHeavy.right]/64) /.
      {k0->p0,k1->0,k2->0,k3->0}];
    background5 = toS[Sum[eta[[a+1]] eta[[b+1]]
      Tr[left.(sig[a,b].numRest+numRest.sig[a,b]).right.sig[a,b]],
      {a,0,3},{b,0,3}]/(192 (p0^2-m^2)^2)];
    full5 = Factor[light5+background5];
    ch -> <|"q3_rational"->q3, "light_d5_rational"->light5,
      "heavy_G_d5_rational"->background5, "full_d5_rational"->full5,
      "q3_borel"->borel[q3], "light_d5_borel"->borel[light5],
      "heavy_G_d5_borel"->borel[background5], "full_d5_borel"->borel[full5]|>
  ], {ch, Keys[channels]}];

checks = <|
  "AB_BA_symmetry" -> FullSimplify[results["AB"]["full_d5_rational"]-results["BA"]["full_d5_rational"]]==0,
  "AA_matches_published_axial_control" -> FullSimplify[
    results["AA"]["full_d5_borel"] + Exp[-m^2/K] m^3/(4 K^2)]==0,
  "V_matches_published_vector_control" -> FullSimplify[
    results["V"]["full_d5_borel"] - Exp[-m^2/K] m^3/(4 K^2)]==0,
  "P_matches_published_pseudoscalar_control" -> FullSimplify[
    results["P"]["full_d5_borel"] - Exp[-m^2/K] (-m/(2 K)+m^3/(4 K^2))]==0,
  "AA_light_equals_current_code" -> FullSimplify[
    results["AA"]["light_d5_borel"] + Exp[-m^2/K] m^3/(4 K^2)]==0,
  "AB_light_equals_current_code" -> FullSimplify[
    results["AB"]["light_d5_borel"] + Exp[-m^2/K] (m^4/K^2-m^2/K-1)/(4 d)]==0,
  "BB_light_equals_current_code" -> FullSimplify[
    results["BB"]["light_d5_borel"] + Exp[-m^2/K] (m^5/(4 K^2)-m^3/(2 K))/d^2]==0
|>;
fullDifference = FullSimplify[(results["AA"]["full_d5_borel"]-
  results["BB"]["full_d5_borel"]) /. d->m];
lightDifference = FullSimplify[(results["AA"]["light_d5_borel"]-
  results["BB"]["light_d5_borel"]) /. d->m];
printedDifference = -Exp[-m^2/K] m/K;
comparison = <|"current_light_only_AA_minus_BB" -> lightDifference,
  "independent_full_AA_minus_BB" -> fullDifference,
  "printed_assignment_AA_minus_BB" -> printedDifference,
  "full_minus_printed" -> FullSimplify[fullDifference-printedDifference]|>;

checks = Map[TrueQ, checks];
If[!And@@Values[checks], Print[InputForm[checks]]; Exit[1]];
jsonResults = Map[Map[ToString[#,InputForm]&,#]&,results];
Export[FileNameJoin[{auditDirectory,"dimension5_symbolic.json"}],
  <|"coefficient_convention"->"C3 or C5 factored out; positive transverse invariant",
    "scope"->"tree-level dimension five, ms=0 Wilson coefficient, d normalization retained",
    "sigma_trace_normalization"->sigmaNorm,"checks"->checks,
    "channels"->jsonResults,
    "massless_light_comparison"->Map[ToString[#,InputForm]&,comparison]|>,"RawJSON"];
Put[results,FileNameJoin[{auditDirectory,"dimension5_symbolic.wdx"}]];
Print[InputForm[checks]];
Print[InputForm[comparison]];
Print[InputForm[results]];
If[!And@@Values[checks],Exit[1]];
