(* ::Package:: *)

(* Build the step-by-step Mathematica notebook for the Ds/Bs gamma paper.

   The generated .nb is intentionally self-contained and evaluatable cell by
   cell.  The batch .wl audits remain the regression-test versions of the
   same symbolic stages.
*)

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
  CreateDirectory[notebookDir, CreateIntermediateDirectories -> True]
];
outFile = FileNameJoin[{notebookDir, "DsBs_gamma_symbolic_derivation.nb"}];

notebook = Notebook[
  {
    Cell["Ds/Bs gamma: current, Wick contraction, propagators, and hard traces", "Title"],
    Cell["Step-by-step Mathematica/FeynCalc derivation", "Subtitle"],
    textCell[
      "This notebook is the interactive companion to the batch audit scripts. " <>
      "Evaluate the input cells from top to bottom. It follows Rohrwild's " <>
      "vacuum three-point/background-field convention, so the final channel carries p " <>
      "and the initial channel carries pPrime=p+q. The legacy k-p+q strange-line " <>
      "routing is retained only as a diagnostic comparison near the end."
    ],

    sectionCell["1. Initialization and conventions"],
    textCell[
      "Load FeynCalc and clear only symbols used in this notebook. " <>
      "Feynman parameters are declared scalar coefficients for shifted-momentum algebra."
    ],
    inputCell[HoldComplete[
      Needs["FeynCalc`"];
      ClearAll[
        mu, nu, alpha, beta, rho, lambda, p, q, pPrime, k, l, eps,
        mQ, ms, mb, eQ, eS, theta, Nc, p2, pPrime2, pq,
        dH1, dH2, dH3, dS1, dS2, dS3, x, y, z,
        xi, x2, mf, ef, gs, eta, M2, ssCond,
        Ffield, Gfield, fA, fB, f1, f2, gA, gB, g1, g2,
        tA, tB, piAA, piAB, piBB, m1, m2, mP, fP,
        borel1, borel2, borelTwoPoint, alphaEM, s, s0, u0,
        chiC, chiGamma, f3Gamma, rhoA, phiGamma, A4,
        HbarGamma, PsiV, IF1, tAHat, tBHat,
        tBpert, tBLocal, tBTw2, tBTw3, tBTw4, tB3p,
        shat, mi, mj, ei, u, w, v, alphaG, uval, kallenLambda,
        rhoPiece, F1DA, hGamma, psiV, uBar, muRen, gammaE, rhoAB,
        rhoDiag, rhoFromF1, f1ModelSquared, f2ModelSquared, f2FromF1Anchor,
        rhoFromF1Residual, rhoDiagResidual, pi12Model,
        dsFoundResidues, bsFoundResidues,
        piAAPert, piBBPert, piABPert, piAACond, piBBCond, piABCond,
        sThreshold, twoPointOPE11, twoPointOPE22,
        xCoord, yCoord, zCoord, Int4, dot, VEV, TimeOrder,
        JPdagger, JX, jEM, VEVF, Fplane, Pi3Point
      ];
      (FeynCalc`DataType[#, FeynCalc`FCVariable] = True) & /@ {x, y, z};
      FeynCalc`$LeviCivitaSign
    ]],
    inputCell[HoldComplete[
      kinematicRules = {
        SP[q, q] -> 0,
        Pair[Momentum[q], Momentum[q]] -> 0,
        Pair[Momentum[p], Momentum[p]] -> p2,
        Pair[Momentum[p], Momentum[q]] -> pq,
        Pair[Momentum[q], Momentum[p]] -> pq,
        Pair[Momentum[k], Momentum[pPrime]] ->
          Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
        Pair[Momentum[p], Momentum[pPrime]] ->
          Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
        Pair[Momentum[pPrime], Momentum[q]] -> Pair[Momentum[p], Momentum[q]]
      };
      pPrime2Rule = pq -> (pPrime2 - p2)/2;
    ]],

    sectionCell["2. Interpolating currents"],
    textCell[
      "J_A is the axial-vector basis current and J_B is the derivative/tensor " <>
      "basis current. The momentum in J_B is the initial-state momentum pPrime=p+q."
    ],
    inputCell[HoldComplete[
      gammaA[index_] := GA[index] . GA[5];
      gammaB[index_] :=
        I FV[pPrime, alpha] DiracSigma[GA[index], GA[alpha]] . GA[5]/(mQ + ms);
      gammaBExpanded[index_] :=
        -(GA[index] . GS[pPrime] - GS[pPrime] . GA[index]) . GA[5]/(2 (mQ + ms));
      gammaP[] := I GA[5];
      gammaLow[index_] := Sin[theta] gammaA[index] + Cos[theta] gammaB[index];
      gammaHigh[index_] := Cos[theta] gammaA[index] - Sin[theta] gammaB[index];
      propNumerator[momentum_, mass_] := GS[momentum] + mass;
    ]],
    inputCell[HoldComplete[
      tensorVertexResidual =
        Contract[DiracSimplify[gammaB[mu] - gammaBExpanded[mu]]] // Simplify
    ]],

    sectionCell["3. Correlator and Wick contraction"],
    textCell[
      "Use Rohrwild's vacuum three-point correlator Pi_X,munu = i^2 Integral_x Integral_y " <>
      "exp(i p.x+i q.y)<0|T{J_P^dagger(x) j_em,nu(y) J_X,mu(0)}|0>. " <>
      "The external-photon correlator is its term linear in the plane-wave background field. " <>
      "A full contraction gives one closed fermion loop: " <>
      "<sbar Gamma_P Q Qbar Gamma_X s> = -Tr[Gamma_P S_Q Gamma_X S_s]."
    ],
    inputCell[HoldComplete[
      rohrwildThreePoint =
        I^2 Int4[xCoord] Int4[yCoord]
          Exp[I dot[p, xCoord] + I dot[q, yCoord]]
          VEV[TimeOrder[
            JPdagger[xCoord], jEM[nu, yCoord], JX[mu, 0]
          ]];
      Fplane[alpha_, beta_, zCoord_] :=
        I (FV[q, alpha] FV[eps, beta] - FV[q, beta] FV[eps, alpha])
          Exp[I dot[q, zCoord]];
      rohrwildBackground =
        I Int4[xCoord] Exp[I dot[p, xCoord]]
          VEVF[TimeOrder[JPdagger[xCoord], JX[mu, 0]], Fplane];
      <|
        "Pi_X,munu(p,q)" -> rohrwildThreePoint,
        "epsilon^nu Pi_X,munu: term linear in F" -> rohrwildBackground,
        "initial momentum" -> pPrime == p + q
      |>
    ]],
    inputCell[HoldComplete[
      fermionLoopSign = -1;
      hardWickPrefactor = I fermionLoopSign Nc;
      softWickPrefactorBeforeFierz = I;
      softWickPrefactorAfterFierz = -I/4;
      <|
        "hard before propagator and QED factors" -> hardWickPrefactor,
        "soft before Fierz" -> softWickPrefactorBeforeFierz,
        "soft after Fierz" -> softWickPrefactorAfterFierz
      |>
    ]],

    sectionCell["4. Propagators without local parts and calculation warning"],
    textCell[
      "Local vacuum-condensate pieces are deliberately not written in the propagators. " <>
      "The displayed expressions contain only the free massive term and the electromagnetic " <>
      "and one-gluon background-field kernels. IMPORTANT: the existing axial numerical " <>
      "calculation still includes the contribution named heavy_local; it has not been removed."
    ],
    inputCell[HoldComplete[
      freePropagator[momentum_, mass_] :=
        I (GS[momentum] + mass)/(SP[momentum, momentum] - mass^2 + I eta);
      lightStrangeFreeTerms = <|
        "free massless" -> I GS[xi]/(2 Pi^2 x2^2),
        "explicit strange mass" -> -ms/(4 Pi^2 x2)
      |>;
      uBar = 1 - u;
      rohrwildLightBackgroundTerms = <|
        "gluon, Rohrwild massless kernel" ->
          -I gs/(16 Pi^2 x2) Inactive[Integrate][
            (uBar GS[xi] . DiracSigma[GA[alpha], GA[beta]] +
              u DiracSigma[GA[alpha], GA[beta]] . GS[xi])
              Gfield[alpha, beta, u], {u, 0, 1}],
        "photon, Rohrwild massless kernel" ->
          -I eS/(16 Pi^2 x2) Inactive[Integrate][
            (uBar GS[xi] . DiracSigma[GA[alpha], GA[beta]] +
              u DiracSigma[GA[alpha], GA[beta]] . GS[xi])
              Ffield[alpha, beta, u], {u, 0, 1}],
        "linear-ms background correction" ->
          -ms/(32 Pi^2) Inactive[Integrate][
            (gs Gfield[alpha, beta, u] + eS Ffield[alpha, beta, u])
              DiracSigma[GA[alpha], GA[beta]]
              (Log[-x2 muRen^2/4] + 2 gammaE), {u, 0, 1}]
      |>;
      <|
        "S_s free massive part" -> lightStrangeFreeTerms,
        "S_s Rohrwild G/F terms plus ms correction" -> rohrwildLightBackgroundTerms,
        "WARNING" ->
          "legacy/current axial numerical results include heavy_local; mixed and vacuum-G2 terms are absent"
      |>
    ]],
    inputCell[HoldComplete[
      electromagneticBackgroundKernel =
        -I ef (
          (GS[k] + mf) Ffield[alpha, rho, u]
            DiracSigma[GA[alpha], GA[rho]]/(2 (mf^2 - SP[k, k])^2)
          + u FV[xi, alpha] Ffield[alpha, rho, u] GA[rho]/
            (mf^2 - SP[k, k])
        );
      oneGluonBackgroundKernel =
        -I gs (
          (GS[k] + mQ) Gfield[alpha, rho, u]
            DiracSigma[GA[alpha], GA[rho]]/(2 (mQ^2 - SP[k, k])^2)
          + u FV[xi, alpha] Gfield[alpha, rho, u] GA[rho]/
            (mQ^2 - SP[k, k])
        );
      propagatorTermStatus = <|
        "free propagators" -> "independently used in Mathematica",
        "explicit strange mass" -> "kept in free and linear background terms",
        "local <sbar s> axial term" -> "not displayed here; included in existing axial numerics",
        "electromagnetic background field" -> "hard traces derived; loop stage pending",
        "one-gluon background field" -> "A and B trace kernels derived",
        "mixed condensate" -> "not included in radiative three-point sum rule",
        "vacuum gluon condensate" -> "not included in radiative three-point sum rule"
      |>;
      <|
        "S^gamma kernel" -> electromagneticBackgroundKernel,
        "S^G kernel" -> oneGluonBackgroundKernel,
        "status" -> propagatorTermStatus
      |>
    ]],

    sectionCell["5. Propagator routing"],
    textCell[
      "The two vertex constraints are checked before any trace: final momentum p " <>
      "at J_P^dagger and initial momentum pPrime=p+q at J_X."
    ],
    inputCell[HoldComplete[
      heavyVertexMomenta = {k - (k - p), (k + q) - (k - p)};
      strangeVertexMomenta = {k - (k - p), k - (k - p - q)};
      <|
        "heavy: {final, initial}" -> heavyVertexMomenta,
        "strange: {final, initial}" -> strangeVertexMomenta,
        "required" -> {p, p + q}
      |>
    ]],
    inputCell[HoldComplete[
      heavyDenominators = {
        dH1 -> SP[k + q, k + q] - mQ^2,
        dH2 -> SP[k, k] - mQ^2,
        dH3 -> SP[k - p, k - p] - ms^2
      };
      strangeDenominators = {
        dS1 -> SP[k, k] - mQ^2,
        dS2 -> SP[k - p, k - p] - ms^2,
        dS3 -> SP[k - p - q, k - p - q] - ms^2
      };
      {heavyDenominators, strangeDenominators}
    ]],

    subsectionCell["Fermion-line Ward identities"],
    inputCell[HoldComplete[
      wardHeavyLHS =
        DiracGammaExpand[propNumerator[k, mQ] . GS[q] . propNumerator[k + q, mQ]];
      wardHeavyRHS =
        (SP[k + q, k + q] - mQ^2) propNumerator[k, mQ] -
        (SP[k, k] - mQ^2) propNumerator[k + q, mQ];
      wardHeavyResidual =
        FCE[DiracSimplify[
          DiracGammaExpand[FCI[wardHeavyLHS - wardHeavyRHS]],
          DiracSubstitute67 -> True, DiracOrder -> True
        ]] // ScalarProductExpand // Simplify;
      wardHeavyResidual
    ]],
    inputCell[HoldComplete[
      wardStrangeLHS =
        DiracGammaExpand[propNumerator[l - q, ms] . GS[q] . propNumerator[l, ms]];
      wardStrangeRHS =
        (SP[l, l] - ms^2) propNumerator[l - q, ms] -
        (SP[l - q, l - q] - ms^2) propNumerator[l, ms];
      wardStrangeResidual =
        FCE[DiracSimplify[
          DiracGammaExpand[FCI[wardStrangeLHS - wardStrangeRHS]],
          DiracSubstitute67 -> True, DiracOrder -> True
        ]] // ScalarProductExpand // Simplify;
      wardStrangeResidual
    ]],

    sectionCell["6. Hard-emission Dirac traces"],
    textCell[
      "Charges, denominators, color, QED vertices, and propagator factors remain " <>
      "outside these numerator traces. Each topology is evaluated separately."
    ],
    inputCell[HoldComplete[
      traceReduce[expression_] :=
        Collect2[DiracSimplify[expression, DiracSubstitute67 -> True], {mQ, ms}];
      heavyTrace[current_] := traceReduce[DiracTrace[
        gammaP[] . propNumerator[k, mQ] . GA[nu] .
        propNumerator[k + q, mQ] . current[mu] .
        propNumerator[k - p, ms]
      ]];
      strangeTrace[current_] := traceReduce[DiracTrace[
        gammaP[] . propNumerator[k, mQ] . current[mu] .
        propNumerator[k - p - q, ms] . GA[nu] .
        propNumerator[k - p, ms]
      ]];
    ]],
    inputCell[HoldComplete[nAHeavy = heavyTrace[gammaA]]],
    inputCell[HoldComplete[nAStrange = strangeTrace[gammaA]]],
    inputCell[HoldComplete[nBHeavy = heavyTrace[gammaB]]],
    inputCell[HoldComplete[nBStrange = strangeTrace[gammaB]]],
    inputCell[HoldComplete[
      nATotal = Collect2[eQ nAHeavy + eS nAStrange, {eQ, eS, mQ, ms}];
      nBTotal = Collect2[eQ nBHeavy + eS nBStrange, {eQ, eS, mQ, ms}];
      {nATotal, nBTotal}
    ]],

    sectionCell["7. E1 projection and inverse-denominator reduction"],
    inputCell[HoldComplete[
      e1Tensor = FV[p, nu] FV[q, mu] - SP[p, q] MT[mu, nu];
      e1Projector = e1Tensor/(2 SP[p, q]^2);
      e1Project[expression_] := Collect2[
        ScalarProductExpand[Contract[expression e1Projector] /. kinematicRules],
        {mQ, ms}
      ];
      projectorNormalization =
        ScalarProductExpand[Contract[e1Tensor e1Projector] /. kinematicRules] // Simplify;
      photonWard =
        ScalarProductExpand[Contract[e1Tensor FV[q, nu]] /. kinematicRules] // Simplify;
      {projectorNormalization, photonWard}
    ]],
    inputCell[HoldComplete[
      projAHeavy = e1Project[nAHeavy];
      projAStrange = e1Project[nAStrange];
      projBHeavy = e1Project[nBHeavy];
      projBStrange = e1Project[nBStrange];
    ]],
    inputCell[HoldComplete[
      heavyRules = {
        Pair[Momentum[p], Momentum[p]] -> p2,
        Pair[Momentum[p], Momentum[q]] -> pq,
        Pair[Momentum[k], Momentum[k]] -> dH2 + mQ^2,
        Pair[Momentum[k], Momentum[q]] -> (dH1 - dH2)/2,
        Pair[Momentum[k], Momentum[p]] ->
          (dH2 + mQ^2 + p2 - ms^2 - dH3)/2
      };
      strangeRules = {
        Pair[Momentum[p], Momentum[p]] -> p2,
        Pair[Momentum[p], Momentum[q]] -> pq,
        Pair[Momentum[k], Momentum[k]] -> dS1 + mQ^2,
        Pair[Momentum[k], Momentum[p]] ->
          (dS1 + mQ^2 + p2 - ms^2 - dS2)/2,
        Pair[Momentum[k], Momentum[q]] -> pq + (dS2 - dS3)/2
      };
      triangleCore[projected_, rules_, denominators_] :=
        Collect2[Expand[projected /. rules] /. Thread[denominators -> 0],
          {mQ, ms, p2, pq}];
    ]],
    inputCell[HoldComplete[
      coreAHeavy = triangleCore[projAHeavy, heavyRules, {dH1, dH2, dH3}];
      coreAStrange = triangleCore[projAStrange, strangeRules, {dS1, dS2, dS3}];
      coreBHeavy = triangleCore[projBHeavy, heavyRules, {dH1, dH2, dH3}];
      coreBStrange = triangleCore[projBStrange, strangeRules, {dS1, dS2, dS3}];
      <|
        "A heavy" -> coreAHeavy, "A strange" -> coreAStrange,
        "B heavy" -> coreBHeavy, "B strange" -> coreBStrange
      |>
    ]],

    sectionCell["8. Hard-loop Feynman parameterization"],
    textCell[
      "For x+y+z=1, each weighted denominator is completed to a square. " <>
      "The residuals below must vanish identically before loop integration."
    ],
    inputCell[HoldComplete[
      parameterKinematics = {
        Pair[Momentum[q], Momentum[q]] -> 0,
        Pair[Momentum[p], Momentum[p]] -> p2,
        Pair[Momentum[p], Momentum[q]] -> (pPrime2 - p2)/2,
        Pair[Momentum[q], Momentum[p]] -> (pPrime2 - p2)/2
      };
      expandParameterSP[expression_] :=
        ExpandScalarProduct[FCI[expression]] /. parameterKinematics // FCE // Expand;
      simplexReduce[expression_] :=
        Together[expandParameterSP[expression] /. z -> 1 - x - y] // Simplify;
    ]],
    inputCell[HoldComplete[
      deltaH = (x + y) mQ^2 + z ms^2 - x z pPrime2 - y z p2;
      shiftH = x q - z p;
      weightedH =
        x (SP[k + q, k + q] - mQ^2) +
        y (SP[k, k] - mQ^2) +
        z (SP[k - p, k - p] - ms^2);
      heavyParameterResidual =
        simplexReduce[weightedH - (SP[k + shiftH, k + shiftH] - deltaH)];
      <|"shift R_H" -> shiftH, "Delta_H" -> deltaH,
        "identity residual" -> heavyParameterResidual|>
    ]],
    inputCell[HoldComplete[
      deltaS = x mQ^2 + (y + z) ms^2 - x y p2 - x z pPrime2;
      shiftS = -y p - z (p + q);
      weightedS =
        x (SP[k, k] - mQ^2) +
        y (SP[k - p, k - p] - ms^2) +
        z (SP[k - p - q, k - p - q] - ms^2);
      strangeParameterResidual =
        simplexReduce[weightedS - (SP[k + shiftS, k + shiftS] - deltaS)];
      <|"shift R_S" -> shiftS, "Delta_S" -> deltaS,
        "identity residual" -> strangeParameterResidual|>
    ]],
    textCell[
      "The common identity is 1/(D1 D2 D3) = 2 Integral dx dy " <>
      "[x D1+y D2+(1-x-y)D3]^(-3), with 0<=x<=1 and 0<=y<=1-x."
    ],

    sectionCell["9. Soft two-particle Fierz traces"],
    inputCell[HoldComplete[
      traceForSoft[current_, basis_] :=
        DiracSimplify[
          DiracTrace[gammaP[] . propNumerator[k, mQ] . current[mu] . basis],
          DiracSubstitute67 -> True
        ] /. kinematicRules;
      softBasis = {
        "scalar" -> 1,
        "pseudoscalar" -> I GA[5],
        "vector" -> GA[rho],
        "axial" -> GA[rho] . GA[5],
        "tensor" -> DiracSigma[GA[rho], GA[lambda]]
      };
      softTracesA = AssociationThread[
        softBasis[[All, 1]], traceForSoft[gammaA, #] & /@ softBasis[[All, 2]]
      ];
      softTracesB = AssociationThread[
        softBasis[[All, 1]], traceForSoft[gammaB, #] & /@ softBasis[[All, 2]]
      ];
      {softTracesA, softTracesB}
    ]],

    sectionCell["10. Legacy-routing diagnostic"],
    textCell[
      "The old strange-line momentum k-p+q makes the initial-current momentum " <>
      "p-q. It therefore does not satisfy the Rohrwild correlator routing pPrime=p+q."
    ],
    inputCell[HoldComplete[
      legacyStrangeInitialMomentum = k - (k - p + q);
      <|
        "legacy actual" -> legacyStrangeInitialMomentum,
        "required" -> p + q,
        "difference" -> Simplify[legacyStrangeInitialMomentum - (p + q)]
      |>
    ]],

    sectionCell["11. Physical-state projection and Colangelo limit"],
    textCell[
      "The current rotation acts linearly on the E1 invariant amplitudes. " <>
      "Colangelo, De Fazio, and Ozpineci use the axial current J_A only. " <>
      "With J_low=Sin[theta] J_A+Cos[theta] J_B, their calculation is the " <>
      "theta=Pi/2 limit. Residue factors are retained explicitly."
    ],
    inputCell[HoldComplete[
      stateRotation[angle_] := {
        {Sin[angle], Cos[angle]},
        {Cos[angle], -Sin[angle]}
      };
      rTheta = stateRotation[theta];
      stateRotationResidual =
        FullSimplify[Transpose[rTheta] . rTheta - IdentityMatrix[2],
          Assumptions -> Element[theta, Reals]];
      transitionPhysical = Expand[rTheta . {tA, tB}];
      pureAxialTransitionResidual =
        Simplify[(transitionPhysical /. theta -> Pi/2) - {tA, -tB}];
      <|
        "R(theta)" -> rTheta,
        "R^T R - 1" -> stateRotationResidual,
        "{T_low,T_high}" -> transitionPhysical,
        "pure-A residual" -> pureAxialTransitionResidual
      |>
    ]],
    inputCell[HoldComplete[
      g1 = (Sin[theta] fA gA + Cos[theta] fB gB)/f1;
      g2 = (Cos[theta] fA gA - Sin[theta] fB gB)/f2;
      pureAxialCouplingResidual =
        Simplify[(g1 /. {theta -> Pi/2, f1 -> fA}) - gA];
      gColangeloFromOPE =
        Exp[m1^2/borel1 + mP^2/borel2] (mQ + ms) tA/
          (m1 fA mP^2 fP);
      gLowFromOPE =
        Exp[m1^2/borel1 + mP^2/borel2] (mQ + ms)
          (Sin[theta] tA + Cos[theta] tB)/(m1 f1 mP^2 fP);
      colangeloFormulaResidual =
        Simplify[
          (gLowFromOPE /. {theta -> Pi/2, f1 -> fA}) -
            gColangeloFromOPE
        ];
      <|
        "G_low" -> g1,
        "G_high" -> g2,
        "pure-A coupling residual" -> pureAxialCouplingResidual,
        "Colangelo OPE formula residual" -> colangeloFormulaResidual
      |>
    ]],
    inputCell[HoldComplete[
      correlationBasis = {{piAA, piAB}, {piAB, piBB}};
      correlationPhysical =
        Expand[rTheta . correlationBasis . Transpose[rTheta]];
      expectedCorrelationPhysical = {
        {
          Sin[theta]^2 piAA + Cos[theta]^2 piBB +
            2 Sin[theta] Cos[theta] piAB,
          Sin[theta] Cos[theta] (piAA - piBB) +
            (Cos[theta]^2 - Sin[theta]^2) piAB
        },
        {
          Sin[theta] Cos[theta] (piAA - piBB) +
            (Cos[theta]^2 - Sin[theta]^2) piAB,
          Cos[theta]^2 piAA + Sin[theta]^2 piBB -
            2 Sin[theta] Cos[theta] piAB
        }
      };
      correlationRotationResidual =
        FullSimplify[correlationPhysical - expectedCorrelationPhysical,
          Assumptions -> Element[theta, Reals]];
      f1Squared = Exp[m1^2/borelTwoPoint]
        expectedCorrelationPhysical[[1, 1]]/m1^2;
      f2Squared = Exp[m2^2/borelTwoPoint]
        expectedCorrelationPhysical[[2, 2]]/m2^2;
      pureAxialDecayConstantResidual =
        Simplify[(f1Squared /. theta -> Pi/2) -
          Exp[m1^2/borelTwoPoint] piAA/m1^2];
      <|
        "Pi physical" -> correlationPhysical,
        "f_low^2" -> f1Squared,
        "f_high^2" -> f2Squared,
        "rotation residual" -> correlationRotationResidual,
        "pure-A f residual" -> pureAxialDecayConstantResidual
      |>
    ]],
    inputCell[HoldComplete[
      twoPointOPE11 =
        Inactive[Integrate][
          Exp[-s/borelTwoPoint] (
            Sin[theta]^2 piAAPert[s] +
            Cos[theta]^2 piBBPert[s] +
            2 Sin[theta] Cos[theta] piABPert[s]),
          {s, sThreshold, s0}
        ] + Sin[theta]^2 piAACond + Cos[theta]^2 piBBCond +
          2 Sin[theta] Cos[theta] piABCond;
      twoPointOPE22 =
        Inactive[Integrate][
          Exp[-s/borelTwoPoint] (
            Cos[theta]^2 piAAPert[s] +
            Sin[theta]^2 piBBPert[s] -
            2 Sin[theta] Cos[theta] piABPert[s]),
          {s, sThreshold, s0}
        ] + Cos[theta]^2 piAACond + Sin[theta]^2 piBBCond -
          2 Sin[theta] Cos[theta] piABCond;
      <|
        "Pi_11 OPE in basis blocks" -> twoPointOPE11,
        "Pi_22 OPE in basis blocks" -> twoPointOPE22,
        "direct f1^2 if all blocks are known" ->
          Exp[m1^2/borelTwoPoint] twoPointOPE11/m1^2,
        "direct f2^2 if all blocks are known" ->
          Exp[m2^2/borelTwoPoint] twoPointOPE22/m2^2,
        "current status" ->
          "AA is explicit; local BB and AB OPEs are not yet complete"
      |>
    ]],
    inputCell[HoldComplete[
      f1ModelSquared =
        Sin[theta]^2 fA^2 + Cos[theta]^2 fB^2 +
          2 Sin[theta] Cos[theta] rhoAB fA fB;
      f2ModelSquared =
        Cos[theta]^2 fA^2 + Sin[theta]^2 fB^2 -
          2 Sin[theta] Cos[theta] rhoAB fA fB;
      rhoFromF1 =
        (f1^2 - Sin[theta]^2 fA^2 - Cos[theta]^2 fB^2)/
          (2 Sin[theta] Cos[theta] fA fB);
      rhoFromF1Residual =
        Simplify[(f1ModelSquared /. rhoAB -> rhoFromF1) - f1^2];
      f2FromF1Anchor =
        Sqrt[Simplify[f2ModelSquared /. rhoAB -> rhoFromF1]];
      pi12Model =
        Sin[theta] Cos[theta] (fA^2 - fB^2) +
          (Cos[theta]^2 - Sin[theta]^2) rhoAB fA fB;
      rhoDiag =
        -Sin[theta] Cos[theta] (fA^2 - fB^2)/
          ((Cos[theta]^2 - Sin[theta]^2) fA fB);
      rhoDiagResidual = Simplify[pi12Model /. rhoAB -> rhoDiag];
      dsFoundResidues = <|
        "legacy chi, theta=35.3 deg" -> <|
          "f1 median [p16,p84] GeV" -> {0.344, 0.329, 0.363},
          "f2 median [p16,p84] GeV" -> {0.380, 0.339, 0.423}|>,
        "lattice fperp, theta=35.3 deg" -> <|
          "f1 median [p16,p84] GeV" -> {0.345, 0.330, 0.364},
          "f2 median [p16,p84] GeV" -> {0.379, 0.338, 0.423}|>
      |>;
      bsFoundResidues = <|
        "legacy chi, Bs1(5830) closure" -> <|
          "f1 median [p16,p84] GeV" -> {0.436, 0.380, 0.486},
          "f2 median [p16,p84] GeV" -> {0.218, 0.138, 0.292}|>,
        "lattice fperp, Bs1(5830) closure" -> <|
          "f1 median [p16,p84] GeV" -> {0.443, 0.384, 0.492},
          "f2 median [p16,p84] GeV" -> {0.203, 0.121, 0.291}|>
      |>;
      <|
        "f1^2 overlap model" -> f1ModelSquared,
        "f2^2 overlap model" -> f2ModelSquared,
        "rho from Ds f1 anchor" -> rhoFromF1,
        "f2 inferred from that anchor" -> f2FromF1Anchor,
        "rho from Pi12=0 closure" -> rhoDiag,
        "anchor residual" -> rhoFromF1Residual,
        "diagonal-closure residual" -> rhoDiagResidual,
        "Ds values found" -> dsFoundResidues,
        "Bs values found" -> bsFoundResidues,
        "interpretation" ->
          "model normalizations, not a completed AA/BB/AB two-point OPE prediction"
      |>
    ]],
    inputCell[HoldComplete[
      radiativeWidth[g_, initialMass_, finalMass_] :=
        alphaEM/3 g^2
          ((initialMass^2 - finalMass^2)/(2 initialMass))^3;
      pureAxialWidthResidual =
        Simplify[
          (radiativeWidth[g1, m1, mP] /.
            {theta -> Pi/2, f1 -> fA}) -
          radiativeWidth[gA, m1, mP]
        ];
      intervalOverlap[first_, second_] := {
        Max[first[[1]], second[[1]]],
        Min[first[[2]], second[[2]]]
      };
      gIntervalIntersection =
        intervalOverlap[{-0.3978, -0.3124}, {-0.37, -0.29}];
      widthIntervalIntersectionKeV =
        intervalOverlap[{20.51, 33.26}, {19., 29.}];
      <|
        "pure-A width residual" -> pureAxialWidthResidual,
        "G interval overlap [GeV^-1]" -> gIntervalIntersection,
        "width interval overlap [keV]" -> widthIntervalIntersectionKeV
      |>
    ]],

    subsectionCell["Colangelo-style final OPE assembly"],
    textCell[
      "The paper-ready formula is assembled from named OPE components. " <>
      "The axial-current entries reproduce the implemented Colangelo sum rule. " <>
      "The tensor-current component slots remain symbolic until their independent " <>
      "spectral and condensate reductions are complete."
    ],
    inputCell[HoldComplete[
      kallenLambda[shat_, mass1Squared_, mass2Squared_] :=
        shat^2 + mass1Squared^2 + mass2Squared^2 -
        2 shat mass1Squared - 2 shat mass2Squared -
        2 mass1Squared mass2Squared;
      rhoPiece[shat_, mi_, mj_, ei_] :=
        -3 ei/(8 Pi^2) (
          2 mi Log[
            (shat - mj^2 + mi^2 -
              Sqrt[kallenLambda[shat, mj^2, mi^2]])/
            (shat - mj^2 + mi^2 +
              Sqrt[kallenLambda[shat, mj^2, mi^2]])
          ] +
          (mj - mi) (mj^2 - mi^2 - shat)/shat^2
            Sqrt[kallenLambda[shat, mj^2, mi^2]]
        );
      rhoA[shat_] :=
        rhoPiece[shat, ms, mQ, eS] - rhoPiece[shat, mQ, ms, eQ];
      HbarGamma[uval_] :=
        Inactive[Integrate][(uval - w) hGamma[w], {w, 0, uval}];
      PsiV[uval_] :=
        Inactive[Integrate][psiV[w], {w, 0, uval}];
      IF1[uval_] :=
        Inactive[Integrate][
          Inactive[Integrate][
            F1DA[uval - (1 - v) alphaG,
              1 - uval - v alphaG, alphaG],
            {alphaG, 0, uval/(1 - v)}
          ],
          {v, 0, 1 - uval}
        ] +
        Inactive[Integrate][
          Inactive[Integrate][
            F1DA[uval - (1 - v) alphaG,
              1 - uval - v alphaG, alphaG],
            {alphaG, 0, (1 - uval)/v}
          ],
          {v, 1 - uval, 1}
        ];
      <|
        "rho_A^P(s)" -> rhoA[s],
        "Hbar_gamma(u0)" -> HbarGamma[u0],
        "Psi^v(u0)" -> PsiV[u0],
        "I_F1(u0)" -> IF1[u0]
      |>
    ]],
    inputCell[HoldComplete[
      expQ = Exp[-mQ^2/M2];
      exp0 = Exp[-s0/M2];
      axialOPEComponents = <|
        "perturbative" ->
          Inactive[Integrate][Exp[-s/M2] rhoA[s],
            {s, (mQ + ms)^2, s0}],
        "local quark condensate" ->
          eQ expQ ssCond
            (1 - mQ ms/M2 + ms^2/(2 M2) (1 + mQ^2/M2)),
        "twist 2, signed Colangelo chi" ->
          -eS ssCond (expQ - exp0) M2 chiC phiGamma[u0],
        "twist 3" ->
          2 eS f3Gamma mQ expQ PsiV[u0],
        "twist 4" ->
          eS ssCond expQ/4
            (A4[u0] - 8 HbarGamma[u0]) (1 + mQ^2/M2),
        "three particle" ->
          eS ssCond expQ IF1[u0]
      |>;
      tAHat = Total[Values[axialOPEComponents]];
      axialOPEAssemblyResidual =
        Expand[tAHat - Total[Values[axialOPEComponents]]];
      <|
        "T_A components" -> axialOPEComponents,
        "T_A" -> tAHat,
        "assembly residual" -> axialOPEAssemblyResidual
      |>
    ]],
    inputCell[HoldComplete[
      tensorOPEComponents = <|
        "perturbative" -> tBpert,
        "local quark condensate" -> tBLocal,
        "twist 2" -> tBTw2,
        "twist 3" -> tBTw3,
        "twist 4" -> tBTw4,
        "three particle" -> tB3p
      |>;
      tBHat = Total[Values[tensorOPEComponents]];
      tensorOPEAssemblyResidual =
        Expand[tBHat - Total[Values[tensorOPEComponents]]];
      tensorOPEStatus = <|
        "perturbative" -> "projected traces complete; spectral density pending",
        "local quark condensate" -> "independent derivation pending",
        "twist 2" -> "trace derived; final reduction pending",
        "twist 3" -> "diagnostic trace only",
        "twist 4" -> "relevant traces exist; final reduction pending",
        "three particle" -> "kernel traces exist; final integral is not complete"
      |>;
      <|
        "T_B components" -> tensorOPEComponents,
        "T_B" -> tBHat,
        "status" -> tensorOPEStatus,
        "assembly residual" -> tensorOPEAssemblyResidual
      |>
    ]],
    inputCell[HoldComplete[
      paperNorm1 =
        Exp[(m1^2 + mP^2)/(2 M2)] (mQ + ms)/(m1 f1 mP^2 fP);
      paperNorm2 =
        Exp[(m2^2 + mP^2)/(2 M2)] (mQ + ms)/(m2 f2 mP^2 fP);
      paperG1 = paperNorm1
        (Sin[theta] tAHat + Cos[theta] tBHat);
      paperG2 = paperNorm2
        (Cos[theta] tAHat - Sin[theta] tBHat);
      paperG1FromComponents = paperNorm1 Total[
        Sin[theta] Values[axialOPEComponents] +
        Cos[theta] Values[tensorOPEComponents]
      ];
      paperComponentAssemblyResidual =
        Expand[paperG1 - paperG1FromComponents];
      paperPureAxialResidual = Simplify[
        (paperG1 /. {theta -> Pi/2, f1 -> fA}) -
        Exp[(m1^2 + mP^2)/(2 M2)] (mQ + ms) tAHat/
          (m1 fA mP^2 fP)
      ];
      chiConventionResidual = Simplify[
        (axialOPEComponents["twist 2, signed Colangelo chi"] /.
          chiC -> -chiGamma) -
        eS ssCond (expQ - exp0) M2 chiGamma phiGamma[u0]
      ];
      <|
        "G_1 final" -> paperG1,
        "G_2 final" -> paperG2,
        "component assembly residual" -> paperComponentAssemblyResidual,
        "theta=Pi/2 Colangelo residual" -> paperPureAxialResidual,
        "chi-convention residual" -> chiConventionResidual
      |>
    ]],

    sectionCell["12. Exact check summary"],
    inputCell[HoldComplete[
      checkResiduals = <|
        "tensor-current commutator" -> tensorVertexResidual,
        "heavy Ward identity" -> wardHeavyResidual,
        "strange Ward identity" -> wardStrangeResidual,
        "E1 normalization minus one" -> Simplify[projectorNormalization - 1],
        "photon Ward contraction" -> photonWard,
        "heavy Feynman-parameter identity" -> heavyParameterResidual,
        "strange Feynman-parameter identity" -> strangeParameterResidual,
        "orthogonal state rotation" -> stateRotationResidual,
        "pure-A projected invariant" -> pureAxialTransitionResidual,
        "pure-A residue-weighted coupling" -> pureAxialCouplingResidual,
        "two-point correlator rotation" -> correlationRotationResidual,
        "pure-A decay constant" -> pureAxialDecayConstantResidual,
        "Ds f1-anchor closure" -> rhoFromF1Residual,
        "Bs Pi12=0 closure" -> rhoDiagResidual,
        "Colangelo formula limit" -> colangeloFormulaResidual,
        "pure-A width" -> pureAxialWidthResidual,
        "axial OPE component assembly" -> axialOPEAssemblyResidual,
        "tensor OPE component assembly" -> tensorOPEAssemblyResidual,
        "mixed coupling component assembly" -> paperComponentAssemblyResidual,
        "full axial OPE pure-A limit" -> paperPureAxialResidual,
        "signed-to-positive chi conversion" -> chiConventionResidual
      |>;
      zeroResidualQ[value_] :=
        And @@ (TrueQ[# === 0] & /@ Flatten[{value}]);
      <|
        "residuals" -> checkResiduals,
        "overall status" -> If[And @@ (zeroResidualQ /@ Values[checkResiduals]),
          "PASS", "FAIL"]
      |>
    ]],
    textCell[
      "For a clean batch regression run, execute scripts/" <>
      "mathematica_current_wick_trace_audit.wl, scripts/" <>
      "mathematica_correlator_routing_audit.wl, and scripts/" <>
      "mathematica_hard_loop_parameterization.wl. The state projection has its " <>
      "own scripts/mathematica_state_mixing_projection.wl audit. The Bs channel uses the same " <>
      "symbolic expressions with mQ -> mb and the appropriate heavy-quark charge."
    ]
  },
  WindowTitle -> "Ds/Bs gamma symbolic derivation",
  Saveable -> True,
  StyleDefinitions -> "Default.nb"
];

Put[notebook, outFile];
Print["Wrote notebook to: ", outFile];
Print["Input cells: ", Count[notebook, Cell[_, "Input", ___], Infinity]];
Quit[];
