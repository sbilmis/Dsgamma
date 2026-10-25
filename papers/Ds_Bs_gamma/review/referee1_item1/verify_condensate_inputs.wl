(* Independent transcription and high-precision quadrature of printed RHSs.
   This checks the numerical reconstruction, not the physical assignment
   of the appendix labels or completeness of its Wilson coefficients. *)
auditDirectory = DirectoryName[$InputFileName];
Clear[coefficients, reconstructed];

coefficients[k_, m_, l_, z_, g2_, g3_] := Module[{g4 = g2^2, c1, c2, b1, b2},
  (* Group all terms linear in <sbar s> separately from q-independent terms. *)
  c1 = m^5 (4 g4 m^2 + 24 z l m g3)/(2^9 3^4 k^7)
    - m^3 ((52 g4 - 72 g3 z) m^2 + 48 z l m (4 g3 + g2 m^2))/(2^9 3^4 k^6)
    + 2 m^2 (m (14 g4 - 51 g3 z - 12 g2 z m^2)
        + 6 l g3 (z + 2 m^2) - 48 g2 z l m^2)/(2^8 3^3 k^5)
    - (5 g2 - 6 m (7 z m - 2 m^3 - 2 l (6 z - m^2))) m g2/(2^6 3^3 k^4)
    + 6 g3 (23 z m - 8 m^3 + l (30 z + 8 m^2))/(2^8 3^3 k^4)
    - (g2 (4 z m - 6 m^3 - l (4 z - 7 m^2)) - 12 z l m^4)/(144 k^3)
    + g3 (2 z + 12 l m + 7 m^2)/(288 m k^3)
    - (2 m (g2 - 6 z m^2) - 3 l (g2 - 4 z m^2))/(48 k^2)
    + g3/(96 m k^2) - m (z + 2 l m)/(4 k)
    - (z - 2 l m + 4 m^2)/(4 m);
  b1 = g3 l m/(192 k^2 Pi^2) - m g2 l/(96 k Pi^2)
    - g3 (l + 2 m)/(384 m k Pi^2) + g2 (l + m)/(96 m Pi^2)
    - g3 (3 l + 8 m)/(384 m^3 Pi^2);
  c2 = -m^3 (14 g4 m^2 + 96 g3 z l m)/(2^8 3^4 k^6)
    - m (-(22 g4 - 36 g3 z) m^2 - 7 z l m (17 g3 + 4 g2 m^2))/(2^7 3^3 k^5)
    - g2 m (25 g2 - 6 m (5 z m - l (14 z - 4 m^2)))/(2^6 3^3 k^4)
    + g3 (5 (g2 + z m^2) - l m (23 z - 10 m^2))/(2^6 3^2 m k^4)
    + (3 g2 (g2 - 6 z m^2) + l m g2 (20 z - 39 m^2)
        + 72 z l m^5)/(432 m k^3)
    + g3 (108 z m - l (20 z + 81 m^2))/(2^6 3^3 m^2 k^3)
    - l (3 g3 - 4 m^2 (3 g2 - 14 z m^2))/(96 m^2 k^2)
    + (2 g2 - 3 m (4 z m - l (z + 4 m^2)))/(12 m k)
    + l (z + 9 m^2)/(6 m^2);
  b2 = -g3 (2 l - m)/(96 m k Pi^2) - 35 g2/(96 Pi^2)
    + l g2/(48 m Pi^2) - g3 (5 l - 2 m)/(96 m^3 Pi^2);
  Exp[-m^2/k] {b1, b2, c1, c2}
];

reconstructed[k_, s0_, m_] := Module[
  {l = 63/500, z = 4/5, g2 = 4 Pi^2 (3/250), g3 = 57/100,
   terms, integrals, qOld = -(33/2000) (4/5), qNew = -(6/25)^3 (4/5),
   old, new, angle, oldAngle, newAngle, x},
  terms = coefficients[k, m, l, z, g2, g3];
  integrals = {
    NIntegrate[-Exp[-x/k] (g2/(32 x Pi^2)
      + 3 (x-m^2) (x m-m^3-l (x+m^2))/(8 m x Pi^2)),
      {x, m^2, s0}, WorkingPrecision -> 50, AccuracyGoal -> 35, PrecisionGoal -> 35],
    NIntegrate[Exp[-x/k] (g2 (4 x-m^2)/(48 x^2 Pi^2)
      - g3 (2 x+m^2)/(384 m^2 x^2 Pi^2)
      - (x-m^2)^2 (x^2-6 x l m-m^4)/(8 Pi^2 m^2 x^2)),
      {x, m^2, s0}, WorkingPrecision -> 50, AccuracyGoal -> 35, PrecisionGoal -> 35]
  };
  old = N[terms[[1;;2]] + qOld terms[[3;;4]] + integrals, 40];
  new = N[terms[[1;;2]] + qNew terms[[3;;4]] + integrals, 40];
  angle[v_] := Mod[ArcTan[v[[2]], 2 v[[1]]] 90/Pi, 90];
  oldAngle = angle[old]; newAngle = angle[new];
  <|"M2_GeV2" -> N[k], "s0_GeV2" -> N[s0],
    "A1_old" -> N[old[[1]]], "A2_old" -> N[old[[2]]],
    "A1_new" -> N[new[[1]]], "A2_new" -> N[new[[2]]],
    "dA1_dss" -> N[terms[[3]]], "dA2_dss" -> N[terms[[4]]],
    "theta_old_deg" -> N[oldAngle], "theta_new_deg" -> N[newAngle],
    "delta_theta_deg" -> N[newAngle-oldAngle]|>
];

auditRows = Flatten[Table[
  Table[Join[<|"sector" -> sector[[1]]|>, reconstructed[k, s0, sector[[2]]]],
    {k, sector[[3]]}, {s0, sector[[4]]}],
  {sector, {{"Ds", 127/100, {3,9/2,6}, {9,10,11}},
            {"Bs", 209/50, {9,11,13}, {42,43,44}}}}], 2];
Export[FileNameJoin[{auditDirectory, "wolfram_high_precision.json"}], auditRows, "RawJSON"];
Print["Exported ", Length[auditRows], " independent high-precision points."];
