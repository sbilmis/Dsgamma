(* ::Package:: *)

(* Build a self-contained, step-by-step Mathematica notebook for the complete
   normalized-current AA/AB/BB two-point Ds1 sum-rule calculation. *)

ClearAll[inputCell, textCell, sectionCell, subsectionCell];

inputCell[held_HoldComplete] :=
  Cell[BoxData[ToBoxes[held /. HoldComplete -> Defer, StandardForm]], "Input"];
textCell[text_String] := Cell[text, "Text"];
sectionCell[text_String] := Cell[text, "Section"];
subsectionCell[text_String] := Cell[text, "Subsection"];

scriptDir = DirectoryName[$InputFileName];
paperDir = DirectoryName[scriptDir];
notebookDir = FileNameJoin[{paperDir, "notebooks"}];
If[! DirectoryQ[notebookDir],
  CreateDirectory[notebookDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{notebookDir, "Ds1_AA_AB_BB_two_point_sumrule.nb"}];

notebook = Notebook[
  {
    Cell["Ds1 AA/AB/BB two-point QCD sum rule", "Title"],
    Cell["Fixed-angle projection and decay constants without an external f1 anchor", "Subtitle"],
    textCell[
      "This notebook replaces the former normalization shortcut by a single " <>
      "2 x 2 two-point QCD sum-rule calculation. Evaluate the input cells in " <>
      "order. The OPE truncation is exact-mass leading-order perturbation theory " <>
      "+ local dimension-3 strange condensate through O(ms^2) + local " <>
      "dimension-5 mixed condensate. Gluon-condensate and alpha_s corrections " <>
      "are not included and are not hidden in fitted constants."
    ],
    textCell[
      "Scope warning: local condensates are necessary here because this is an " <>
      "ordinary two-point SVZ sum rule. They are not inserted into the " <>
      "external-photon Rohrwild transition LCSR, where photon matrix elements " <>
      "are organized through non-local distribution amplitudes."
    ],

    sectionCell["1. Basis currents and physical-state rotation"],
    textCell[
      "The normalized basis is J_A^mu = sbar gamma^mu gamma_5 c and " <>
      "J_B^mu = i sbar sigma^(mu nu) p'_nu gamma_5 c/(mc+ms). " <>
      "The physical currents are J_1 = sin(theta) J_A + cos(theta) J_B and " <>
      "J_2 = cos(theta) J_A - sin(theta) J_B."
    ],
    inputCell[HoldComplete[
      ClearAll[theta, piAA, piAB, piBB];
      piBasis = {{piAA, piAB}, {piAB, piBB}};
      rotation[theta_] := {{Sin[theta], Cos[theta]},
        {Cos[theta], -Sin[theta]}};
      piPhysical = FullSimplify[rotation[theta] . piBasis . Transpose[rotation[theta]]];
      MatrixForm[piPhysical]
    ]],
    inputCell[HoldComplete[
      <|
        "Pi11" -> Expand[piPhysical[[1, 1]]],
        "Pi22" -> Expand[piPhysical[[2, 2]]],
        "Pi12" -> TrigReduce[Expand[piPhysical[[1, 2]]]]
      |>
    ]],
    textCell[
      "With the currents above, Pi11 = sin^2(theta) PiAA + cos^2(theta) PiBB " <>
      "+ 2 sin(theta) cos(theta) PiAB. Any formula that interchanges the AA " <>
      "and BB coefficients for state 1 is inconsistent with the stated current definition."
    ],

    sectionCell["2. Exact-mass LO perturbative spectral matrix"],
    inputCell[HoldComplete[
      ClearAll[s, mc, ms, kallen, rhoAA, rhoAB, rhoBB, rhoMatrix];
      kallen[s_] := (s - (mc + ms)^2) (s - (mc - ms)^2);
      rhoAA[s_] := Sqrt[kallen[s]]/(8 Pi^2) (
        2 - (mc^2 + ms^2 + 6 mc ms)/s - (mc^2 - ms^2)^2/s^2);
      rhoAB[s_] := Sqrt[kallen[s]]/(8 Pi^2) (
        3 (ms - mc) ((mc + ms)^2 - s)/((mc + ms) s));
      rhoBB[s_] := Sqrt[kallen[s]]/(8 Pi^2) (
        -((mc + ms)^2 - s) (2 ms^2 - 4 mc ms + 2 mc^2 + s)/
          ((mc + ms)^2 s));
      rhoMatrix[s_] := {{rhoAA[s], rhoAB[s]}, {rhoAB[s], rhoBB[s]}};
      MatrixForm[rhoMatrix[s]]
    ]],
    textCell[
      "All three entries follow from the same cut one-loop correlator and the " <>
      "same normalized tensor current. The off-diagonal entry is symmetric."
    ],

    subsectionCell["2.1 Analytic checks"],
    inputCell[HoldComplete[
      ClearAll[z, mh, ml, rhoAASymbolic, rhoAAColangeloSymbolic];
      rhoAASymbolic[z_, mh_, ml_] :=
        Sqrt[(z - (mh + ml)^2) (z - (mh - ml)^2)]/(8 Pi^2) (
          2 - (mh^2 + ml^2 + 6 mh ml)/z - (mh^2 - ml^2)^2/z^2);
      rhoAAColangeloSymbolic[z_, mh_, ml_] :=
        Sqrt[1 - 2 (mh^2 + ml^2)/z + (mh^2 - ml^2)^2/z^2]/(8 Pi^2) (
          2 z - mh^2 - ml^2 - 6 mh ml - (mh^2 - ml^2)^2/z);
      FullSimplify[
        rhoAASymbolic[z, mh, ml]/rhoAAColangeloSymbolic[z, mh, ml],
        Assumptions -> z > (mh + ml)^2 && mh > ml > 0]
    ]],
    inputCell[HoldComplete[
      FullSimplify[rhoMatrix[s][[1, 2]] - rhoMatrix[s][[2, 1]]]
    ]],
    textCell[
      "The first cell must return 1: the AA spectral density is algebraically " <>
      "identical to the Colangelo axial-current expression. The second must return 0."
    ],

    sectionCell["3. Local two-point OPE matrix"],
    textCell[
      "The coordinate-space local light-quark terms used to generate this matrix are " <>
      "-<sbar s>/12 [1 - i ms xslash/4 - ms^2 x^2/8] " <>
      "- x^2 <sbar g sigma.G s>/192. They are displayed here through their " <>
      "Borel-transformed matrix contributions, which are what enter the sum rule."
    ],
    inputCell[HoldComplete[
      ClearAll[M2, ss, mixedSS, localPieces, localMatrix];
      localPieces[M2_] := Module[{den = mc + ms, e, q0, q1, q2, d5},
        e = Exp[-mc^2/M2];
        q0 = ss {{mc, mc^2/den}, {mc^2/den, mc^3/den^2}};
        q1 = ss {{
            -ms mc^2/(2 M2),
            ms mc (1 - mc^2/M2)/(2 den)}, {
            ms mc (1 - mc^2/M2)/(2 den),
            ms mc^2 (1 - mc^2/(2 M2))/den^2}};
        q2 = ss {{
            ms^2 mc^3/(2 M2^2),
            ms^2 (mc^4/M2^2 - mc^2/M2 - 1)/(2 den)}, {
            ms^2 (mc^4/M2^2 - mc^2/M2 - 1)/(2 den),
            ms^2 (mc^5/(2 M2^2) - mc^3/M2)/den^2}};
        d5 = mixedSS {{
            -mc^3/(4 M2^2),
            -(mc^4/M2^2 - mc^2/M2 - 1)/(4 den)}, {
            -(mc^4/M2^2 - mc^2/M2 - 1)/(4 den),
            -(mc^5/(4 M2^2) - mc^3/(2 M2))/den^2}};
        <|"d3_ms0" -> e q0, "d3_ms1" -> e q1,
          "d3_ms2" -> e q2, "d5_mixed" -> e d5|>
      ];
      localMatrix[M2_] := Total[Values[localPieces[M2]]];
      localPieces[M2] // Map[MatrixForm]
    ]],

    subsectionCell["3.1 AA local-term check against Colangelo"],
    inputCell[HoldComplete[
      aaLocalColangelo[M2_] := Exp[-mc^2/M2] (
        ss (mc - mc^2 ms/(2 M2) + mc^3 ms^2/(2 M2^2))
        - mixedSS mc^3/(4 M2^2));
      FullSimplify[localMatrix[M2][[1, 1]] - aaLocalColangelo[M2]]
    ]],
    textCell["This cell must return 0 before numerical inputs are assigned."],

    sectionCell["4. Numerical inputs"],
    inputCell[HoldComplete[
      mc = 1.27;
      ms = 0.093;
      qq = -(0.24)^3;
      kappaS = 0.8;
      ss = kappaS qq;
      m0sq = 0.8;
      mixedSS = m0sq ss;
      m1 = 2.4595;
      m2 = 2.53511;
      sThreshold = (mc + ms)^2;
      <|"mc_GeV" -> mc, "ms_GeV" -> ms,
        "ss_GeV3" -> ss, "mixedSS_GeV5" -> mixedSS,
        "m1_GeV" -> m1, "m2_GeV" -> m2|>
    ]],

    sectionCell["5. Borel matrix and continuum subtraction"],
    textCell[
      "The interpolated cumulative trapezoid is a fast numerical realization " <>
      "of Integrate[rho_ij(s) Exp[-s/M2], {s, threshold, s0}]. Clear the " <>
      "memoized interpolations after changing a numerical input."
    ],
    inputCell[HoldComplete[
      ClearAll[buildPerturbativeInterpolations, perturbativeInterpolations,
        perturbativeMatrix, opeMatrix];
      buildPerturbativeInterpolations[M2_?NumericQ] := Module[
        {sGrid, integrand, increments, cumulative},
        sGrid = N@Subdivide[sThreshold, 11.0, 2000];
        integrand = (rhoMatrix[#] Exp[-#/M2]) & /@ sGrid;
        increments = MapThread[0.5 (#1 + #2) #3 &,
          {Most[integrand], Rest[integrand], Differences[sGrid]}];
        cumulative = Prepend[Accumulate[increments], ConstantArray[0., {2, 2}]];
        Table[Interpolation[
          Transpose[{sGrid, cumulative[[All, i, j]]}], InterpolationOrder -> 1],
          {i, 1, 2}, {j, 1, 2}]
      ];
      perturbativeInterpolations[M2_?NumericQ] :=
        perturbativeInterpolations[M2] = buildPerturbativeInterpolations[M2];
      perturbativeMatrix[M2_?NumericQ, s0_?NumericQ] :=
        Map[#[s0] &, perturbativeInterpolations[M2], {2}];
      opeMatrix[M2_?NumericQ, s0_?NumericQ] :=
        perturbativeMatrix[M2, s0] + localMatrix[M2];
      MatrixForm[opeMatrix[2.35, 10.0]]
    ]],

    sectionCell["6. External mixing angle, mass matching, and residues"],
    textCell[
      "The nominal projection uses theta_Ds = 26.6 +/- 0.6 degrees from the " <>
      "dedicated previous QCD-sum-rule study (arXiv:2408.08014v2). The angle " <>
      "obtained by diagonalizing the present LO+d3+d5 matrix is retained as a " <>
      "truncation diagnostic; it is not used as the nominal input."
    ],
    inputCell[HoldComplete[
      ClearAll[mixingAngle, projectedValue, projectedOPE,
        projectedEffectiveMass, fittedThreshold, physicalPoint];
      mixingAngle[matrix_?MatrixQ] := Module[{aa, ab, bb, angle},
        {aa, ab, bb} = {matrix[[1, 1]], matrix[[1, 2]], matrix[[2, 2]]};
        angle = ArcTan[aa - bb, -2 ab]/2;
        Mod[angle, Pi/2]
      ];
      projectedValue[matrix_?MatrixQ, angle_?NumericQ, channel_Integer] :=
        Module[{v = rotation[angle][[channel]]}, v . matrix . v];
      projectedOPE[M2_?NumericQ, s0_?NumericQ, angle_?NumericQ,
          channel_Integer] := projectedValue[opeMatrix[M2, s0], angle, channel];
      projectedEffectiveMass[M2_?NumericQ, s0_?NumericQ, angle_?NumericQ,
          channel_Integer, deltaTau_: 0.0002] := Module[
        {tau = 1/M2, plus, minus, massSq},
        plus = projectedOPE[1/(tau + deltaTau), s0, angle, channel];
        minus = projectedOPE[1/(tau - deltaTau), s0, angle, channel];
        If[plus <= 0 || minus <= 0, Return[Indeterminate]];
        massSq = -(Log[plus] - Log[minus])/(2 deltaTau);
        If[massSq > 0, Sqrt[massSq], Indeterminate]
      ];
      fittedThreshold[M2_?NumericQ, angle_?NumericQ, channel_Integer,
          targetMass_?NumericQ] := Module[{candidates, pairs},
        candidates = N[Range[6.55, 10.50, 0.025]];
        pairs = Select[
          {#, projectedEffectiveMass[M2, #, angle, channel]} & /@ candidates,
          NumericQ[#[[2]]] &];
        First@MinimalBy[pairs, Abs[#[[2]] - targetMass] &]
      ];
    ]],
    inputCell[HoldComplete[
      thetaDsInput = 26.6 Degree;
      thetaDsSigma = 0.6 Degree;
      physicalPoint[M2_?NumericQ] := Module[
        {angleSamples, angleDiagnostic, angle, fit1, fit2, pi1, pi2,
         rotatedCommon, offdiagNorm},
        angleSamples = mixingAngle[opeMatrix[M2, #]] & /@
          N[Range[9.0, 11.0, 0.1]];
        angleDiagnostic = Median[angleSamples];
        angle = thetaDsInput;
        fit1 = fittedThreshold[M2, angle, 1, m1];
        fit2 = fittedThreshold[M2, angle, 2, m2];
        pi1 = projectedOPE[M2, fit1[[1]], angle, 1];
        pi2 = projectedOPE[M2, fit2[[1]], angle, 2];
        rotatedCommon = rotation[angle] . opeMatrix[M2, 10.0] .
          Transpose[rotation[angle]];
        offdiagNorm = Abs[rotatedCommon[[1, 2]]]/
          Sqrt[Abs[rotatedCommon[[1, 1]] rotatedCommon[[2, 2]]]];
        <|"M2_GeV2" -> M2, "theta_deg" -> angle 180/Pi,
          "theta_matrix_diagnostic_deg" -> angleDiagnostic 180/Pi,
          "Pi12_normalized" -> offdiagNorm,
          "s01_GeV2" -> fit1[[1]], "s02_GeV2" -> fit2[[1]],
          "mEff1_GeV" -> fit1[[2]], "mEff2_GeV" -> fit2[[2]],
          "f1_GeV" -> Sqrt[pi1 Exp[m1^2/M2]/m1^2],
          "f2_GeV" -> Sqrt[pi2 Exp[m2^2/M2]/m2^2]|>
      ];
      centralPoint = physicalPoint[2.35]
    ]],
    textCell[
      "Expected central output: theta input = 26.6 degrees, matrix diagnostic " <>
      "theta = 30.58596 degrees, s01 = 7.925 GeV^2, s02 = 9.975 GeV^2, " <>
      "f1 = 0.4005893 GeV, and f2 = 0.1666704 GeV. Neither residue is externally anchored."
    ],

    sectionCell["7. Python/Mathematica central-point regression"],
    inputCell[HoldComplete[
      pythonReference = <|
        "theta_deg" -> 26.6,
        "theta_matrix_diagnostic_deg" -> 30.58596158596751,
        "Pi12_normalized" -> 0.18100724350506212,
        "s01_GeV2" -> 7.925,
        "s02_GeV2" -> 9.975,
        "f1_GeV" -> 0.40058930523172204,
        "f2_GeV" -> 0.1666703674111434|>;
      AssociationMap[centralPoint[#] - pythonReference[#] &,
        Keys[pythonReference]]
    ]],
    textCell[
      "With the default 2000-point cumulative integration grid, all displayed " <>
      "differences are below 10^-7. Increase the grid in Section 5 for an " <>
      "arbitrarily tighter integration comparison."
    ],

    sectionCell["8. Optional full deterministic Borel scan"],
    textCell[
      "This cell is deliberately separate because it builds three interpolation " <>
      "tables for every M2 value. It is the Mathematica analogue of the dense " <>
      "Python physical-pole grid."
    ],
    inputCell[HoldComplete[
      fullScan = physicalPoint /@ N[Range[2.0, 2.7, 0.05]];
      Dataset[fullScan]
    ]],
    inputCell[HoldComplete[
      scanSummary = AssociationMap[
        Quantile[Lookup[fullScan, #], {0.16, 0.50, 0.84}] &,
        {"theta_deg", "theta_matrix_diagnostic_deg", "Pi12_normalized",
         "s01_GeV2", "s02_GeV2", "f1_GeV", "f2_GeV"}]
    ]],

    sectionCell["9. Interpretation"],
    textCell[
      "The Python Monte Carlo samples theta_Ds = 26.6 +/- 0.6 degrees as an " <>
      "external input and gives f1 = 0.4046 [0.3799, 0.4291] GeV and " <>
      "f2 = 0.1677 [0.1594, 0.1756] GeV. Diagonalizing the same truncated " <>
      "matrix gives the diagnostic theta = 30.63 [29.94, 31.29] degrees. " <>
      "The normalized off-diagonal residual after projection at the external " <>
      "angle is about 0.18, and is retained as an OPE-truncation diagnostic."
    ],
    textCell[
      "For the derivation-level Dirac checks, run scripts/" <>
      "twopoint_dirac_spectral_audit.wl and scripts/" <>
      "twopoint_local_condensate_audit.wl. For the independent numerical " <>
      "regression, run scripts/mathematica_twopoint_ds1_matrix_check.wl."
    ]
  },
  WindowTitle -> "Ds1 AA/AB/BB two-point sum rule",
  Saveable -> True,
  StyleDefinitions -> "Default.nb"
];

Put[notebook, outFile];
Print["Wrote notebook to: ", outFile];
Print["Input cells: ", Count[notebook, Cell[_, "Input", ___], Infinity]];
Quit[];
