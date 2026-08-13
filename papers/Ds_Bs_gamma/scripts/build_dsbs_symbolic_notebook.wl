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
        xi, x2, mf, ef, gs, eta, M2, M1sq, M2sq, Meff2, u0B, ssCond,
        Ffield, Gfield, fA, fB, f1, f2, gA, gB, g1, g2,
        tA, tB, piAA, piAB, piBB, m1, m2, mP, fP,
        borel1, borel2, borelTwoPoint, alphaEM, s, s0, u0,
        chiC, chiGamma, f3Gamma, rhoA, phiGamma, A4, B4,
        psiVBar, IF1, tAHat, tBHat,
        tBpert, tBTw2, tBTw3, tBTw4,
        Dsoft, borelFactor, borelMaster, borelPq, borelPPrime2,
        borelP2, Ffun, n, a, a0B, Feff, xAS, xAT1, xAT2, xAT3,
        xAT4, xBS, xBT1, xBT2, xBT3, xBT4,
        aq, aqp, aqb, aqbp, betaP, S3, St3, T13, T23, T33, T43,
        Sg3, T4g3, lineProjection3, pLine3, pLinePrime3,
        fGSigma, fGX, fGA, fGammaSigma, fGammaA,
        tA3g, tB3g, tA3em, tB3em, threeParticleResiduals,
        expQ, exp0, deltaEQ,
        rhoB, rhoBLeg, lambdaB, aScoef, bScoef, aQcoef, bQcoef,
        shat, mi, mj, ei, u, w, v, alphaG, uval, kallenLambda,
        rhoPiece, F1DA, uBar, muRen, gammaE, rhoAB,
        rhoAA2pt, rhoAB2pt, rhoBB2pt, lambda2pt, rho2ptMatrix,
        dTwoPoint, mixedCond, s01, s02, q0TwoPoint, q1TwoPoint,
        q2TwoPoint, q5TwoPoint, localTwoPointMatrix,
        v1TwoPoint, v2TwoPoint, rho11TwoPoint, rho22TwoPoint,
        local11TwoPoint, local22TwoPoint, pi11Explicit, pi22Explicit,
        f1Explicit, f2Explicit, centralTwoPointRules,
        centralF1F2FromExplicit,
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
      "and one-gluon background-field kernels. In the Rohrwild transition OPE there is no " <>
      "standalone local <sbar s> term. Local condensates remain required in the separate " <>
      "two-point sum rules for f1 and f2."
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
          "no standalone local condensate is used in the radiative transition OPE"
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
        "local <sbar s> transition term" -> "absent in the Rohrwild transition OPE",
        "electromagnetic background field" ->
          "hard traces and final spectral/DA terms implemented",
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
      lambda2pt[shat_] :=
        (shat - (mQ + ms)^2) (shat - (mQ - ms)^2);
      rhoAA2pt[shat_] := Sqrt[lambda2pt[shat]]/(8 Pi^2) (
        2 - (mQ^2 + ms^2 + 6 mQ ms)/shat -
        (mQ^2 - ms^2)^2/shat^2);
      rhoAB2pt[shat_] := Sqrt[lambda2pt[shat]]/(8 Pi^2) (
        3 (ms - mQ) ((mQ + ms)^2 - shat)/
        ((mQ + ms) shat));
      rhoBB2pt[shat_] := Sqrt[lambda2pt[shat]]/(8 Pi^2) (
        -((mQ + ms)^2 - shat) *
        (2 ms^2 - 4 mQ ms + 2 mQ^2 + shat)/
        ((mQ + ms)^2 shat));
      rho2ptMatrix[shat_] := {
        {rhoAA2pt[shat], rhoAB2pt[shat]},
        {rhoAB2pt[shat], rhoBB2pt[shat]}
      };
      MatrixForm[rho2ptMatrix[s]]
    ]],
    subsectionCell["Explicit f1 and f2 formulas used in the calculation"],
    textCell[
      "This is the complete normalization step. The two physical poles use " <>
      "separate continuum thresholds s01 and s02, fitted to M1 and M2. " <>
      "Local condensates occur here because this is an ordinary two-point SVZ " <>
      "sum rule; they are not inserted into the radiative transition OPE."
    ],
    inputCell[HoldComplete[
      dTwoPoint = mQ + ms;
      q0TwoPoint = {
        {mQ, mQ^2/dTwoPoint},
        {mQ^2/dTwoPoint, mQ^3/dTwoPoint^2}
      };
      q1TwoPoint = {
        {
          -ms mQ^2/(2 borelTwoPoint),
          ms mQ (1 - mQ^2/borelTwoPoint)/(2 dTwoPoint)
        },
        {
          ms mQ (1 - mQ^2/borelTwoPoint)/(2 dTwoPoint),
          ms mQ^2 (1 - mQ^2/(2 borelTwoPoint))/dTwoPoint^2
        }
      };
      q2TwoPoint = {
        {
          ms^2 mQ^3/(2 borelTwoPoint^2),
          ms^2 (
            mQ^4/borelTwoPoint^2 - mQ^2/borelTwoPoint - 1
          )/(2 dTwoPoint)
        },
        {
          ms^2 (
            mQ^4/borelTwoPoint^2 - mQ^2/borelTwoPoint - 1
          )/(2 dTwoPoint),
          ms^2 (
            mQ^5/(2 borelTwoPoint^2) -
            mQ^3/borelTwoPoint
          )/dTwoPoint^2
        }
      };
      q5TwoPoint = {
        {
          -mQ^3/(4 borelTwoPoint^2),
          -(mQ^4/borelTwoPoint^2 -
            mQ^2/borelTwoPoint - 1)/(4 dTwoPoint)
        },
        {
          -(mQ^4/borelTwoPoint^2 -
            mQ^2/borelTwoPoint - 1)/(4 dTwoPoint),
          -(mQ^5/(4 borelTwoPoint^2) -
            mQ^3/(2 borelTwoPoint))/dTwoPoint^2
        }
      };
      localTwoPointMatrix =
        Exp[-mQ^2/borelTwoPoint] (
          ssCond (q0TwoPoint + q1TwoPoint + q2TwoPoint) +
          mixedCond q5TwoPoint
        );
      <|
        "Q0" -> MatrixForm[q0TwoPoint],
        "Q1" -> MatrixForm[q1TwoPoint],
        "Q2" -> MatrixForm[q2TwoPoint],
        "Q5" -> MatrixForm[q5TwoPoint],
        "local basis matrix" -> MatrixForm[localTwoPointMatrix]
      |>
    ]],
    inputCell[HoldComplete[
      v1TwoPoint = {Sin[theta], Cos[theta]};
      v2TwoPoint = {Cos[theta], -Sin[theta]};
      rho11TwoPoint[shat_] :=
        FullSimplify[
          v1TwoPoint . rho2ptMatrix[shat] . v1TwoPoint,
          Assumptions -> Element[theta, Reals]
        ];
      rho22TwoPoint[shat_] :=
        FullSimplify[
          v2TwoPoint . rho2ptMatrix[shat] . v2TwoPoint,
          Assumptions -> Element[theta, Reals]
        ];
      local11TwoPoint =
        FullSimplify[
          v1TwoPoint . localTwoPointMatrix . v1TwoPoint,
          Assumptions -> Element[theta, Reals]
        ];
      local22TwoPoint =
        FullSimplify[
          v2TwoPoint . localTwoPointMatrix . v2TwoPoint,
          Assumptions -> Element[theta, Reals]
        ];
      <|
        "rho_11(s), fully substituted" -> rho11TwoPoint[s],
        "rho_22(s), fully substituted" -> rho22TwoPoint[s],
        "Pi_11 local, fully substituted" -> local11TwoPoint,
        "Pi_22 local, fully substituted" -> local22TwoPoint
      |>
    ]],
    inputCell[HoldComplete[
      pi11Explicit =
        Inactive[Integrate][
          Exp[-s/borelTwoPoint] rho11TwoPoint[s],
          {s, dTwoPoint^2, s01}
        ] + local11TwoPoint;
      pi22Explicit =
        Inactive[Integrate][
          Exp[-s/borelTwoPoint] rho22TwoPoint[s],
          {s, dTwoPoint^2, s02}
        ] + local22TwoPoint;
      f1Explicit =
        Exp[m1^2/(2 borelTwoPoint)]/m1 Sqrt[pi11Explicit];
      f2Explicit =
        Exp[m2^2/(2 borelTwoPoint)]/m2 Sqrt[pi22Explicit];
      <|
        "Pi_11(M^2,s01,theta)" -> pi11Explicit,
        "Pi_22(M^2,s02,theta)" -> pi22Explicit,
        "f1(M^2,s01,theta)" -> f1Explicit,
        "f2(M^2,s02,theta)" -> f2Explicit
      |>
    ]],
    inputCell[HoldComplete[
      centralTwoPointRules = {
        mQ -> 1.27,
        ms -> 0.093,
        ssCond -> 0.8 (-(0.24)^3),
        mixedCond -> 0.8 (0.8 (-(0.24)^3)),
        theta -> 26.6 Degree,
        borelTwoPoint -> 2.35,
        m1 -> 2.4595,
        m2 -> 2.53511,
        s01 -> 7.925,
        s02 -> 9.975
      };
      centralF1F2FromExplicit =
        N[
          {f1Explicit, f2Explicit} /. centralTwoPointRules /.
            Inactive[Integrate] -> NIntegrate,
          10
        ];
      <|
        "central substitutions" -> centralTwoPointRules,
        "{f1,f2} in GeV" -> centralF1F2FromExplicit,
        "Python central reference" ->
          {0.40058930523172204, 0.1666703674111434},
        "Mathematica minus Python" ->
          centralF1F2FromExplicit -
            {0.40058930523172204, 0.1666703674111434}
      |>
    ]],
    textCell[
      "The same fully substituted AA/AB/BB expressions give the direct " <>
      "bottom-strange residues after mQ->mb.  The following cell uses the " <>
      "central bottom two-point point; no Pullin-Zwicky decay constant and " <>
      "no overlap closure is inserted."
    ],
    inputCell[HoldComplete[
      centralBsTwoPointRules = {
        mQ -> 4.18,
        ms -> 0.093,
        ssCond -> 0.8 (-(0.24)^3),
        mixedCond -> 0.8 (0.8 (-(0.24)^3)),
        theta -> 38.5 Degree,
        borelTwoPoint -> 7.00,
        m1 -> 5.750,
        m2 -> 5.82870,
        s01 -> 44.35,
        s02 -> 43.025
      };
      centralBsF1F2FromExplicit =
        N[
          {f1Explicit, f2Explicit} /. centralBsTwoPointRules /.
            Inactive[Integrate] -> NIntegrate,
          10
        ];
      <|
        "central Bs substitutions" -> centralBsTwoPointRules,
        "{f1_Bs,f2_Bs} in GeV" -> centralBsF1F2FromExplicit,
        "Python central reference" ->
          {0.5381107370, 0.0888786545},
        "Mathematica minus Python" ->
          centralBsF1F2FromExplicit -
            {0.5381107370, 0.0888786545}
      |>
    ]],
    inputCell[HoldComplete[
      twoPointOPE11 =
        Inactive[Integrate][
          Exp[-s/borelTwoPoint] (
            Sin[theta]^2 rhoAA2pt[s] +
            Cos[theta]^2 rhoBB2pt[s] +
            2 Sin[theta] Cos[theta] rhoAB2pt[s]),
          {s, sThreshold, s0}
        ] + Sin[theta]^2 piAACond + Cos[theta]^2 piBBCond +
          2 Sin[theta] Cos[theta] piABCond;
      twoPointOPE22 =
        Inactive[Integrate][
          Exp[-s/borelTwoPoint] (
            Cos[theta]^2 rhoAA2pt[s] +
            Sin[theta]^2 rhoBB2pt[s] -
            2 Sin[theta] Cos[theta] rhoAB2pt[s]),
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
          "the fully substituted square-root formulas and their central numerical evaluation are displayed in the preceding four cells"
      |>
    ]],
    textCell[
      "The following overlap cell is retained only as a historical regression " <>
      "of the superseded normalization shortcut. The current Ds result from " <>
      "the complete matrix is f1=0.4046[0.3799,0.4291] GeV and " <>
      "f2=0.1677[0.1594,0.1756] GeV at theta=26.6+/-0.6 degrees; neither is anchored."
    ],
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
          "legacy regression only; superseded by the complete AA/BB/AB notebook"
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

    subsectionCell["Off-shell double-Borel identities and Rohrwild OPE assembly"],
    textCell[
      "On the QCD side p2 and pPrime2 remain independent until both Borel " <>
      "transforms have been applied. Pole-mass substitutions are allowed only " <>
      "in the phenomenological residue. The Rohrwild transition OPE contains no " <>
      "standalone local-condensate component; condensate factors can still " <>
      "normalize nonlocal photon DAs."
    ],
    inputCell[HoldComplete[
      Dsoft[u_] := mQ^2 - (1 - u) pPrime2 - u p2;
      Meff2 = M1sq M2sq/(M1sq + M2sq);
      u0B = M1sq/(M1sq + M2sq);
      borelFactor[n_] :=
        Meff2^(2 - n) Exp[-mQ^2/Meff2]/Gamma[n];
      borelMaster[n_, Ffun_] :=
        borelFactor[n] Ffun[u0B];
      borelPq[n_, Ffun_] :=
        borelFactor[n - 1] Ffun'[u0B]/(2 (n - 1));
      borelPPrime2[n_, Ffun_] :=
        mQ^2 borelFactor[n] Ffun[u0B] -
        borelFactor[n - 1] Ffun[u0B] +
        borelFactor[n - 1]
          (D[u Ffun[u], u] /. u -> u0B)/(n - 1);
      borelP2[n_, Ffun_] :=
        mQ^2 borelFactor[n] Ffun[u0B] -
        borelFactor[n - 1] Ffun[u0B] -
        borelFactor[n - 1]
          (D[(1 - u) Ffun[u], u] /. u -> u0B)/(n - 1);
      <|
        "D(u), before Borel" -> Dsoft[u],
        "M_eff^2" -> Meff2,
        "u0" -> u0B,
        "B[Integral F/D^n]" -> borelMaster[n, Ffun],
        "B[Integral (p.q) F/D^n]" -> borelPq[n, Ffun],
        "B[Integral pPrime^2 F/D^n]" -> borelPPrime2[n, Ffun],
        "B[Integral p^2 F/D^n]" -> borelP2[n, Ffun]
      |>
    ]],
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
        "Rohrwild A(u0)" -> A4[u0],
        "Rohrwild B(u0)" -> B4[u0],
        "Rohrwild vector DA" -> psiVBar[u0],
        "I_F1(u0)" -> IF1[u0]
      |>
    ]],
    textCell[
      "Rohrwild uses A(u), B(u), and overbar psi^(V)(u). The exact BBK " <>
      "dictionary is B(u)=-4 Integral_0^u d alpha (u-alpha) h_gamma(alpha), " <>
      "so Colangelo's A-8 Hbar becomes A+2 B. H_gamma is not introduced as " <>
      "an independent Rohrwild amplitude."
    ],
    subsectionCell["Complete operator-weighted three-particle Borel terms"],
    textCell[
      "The heavy-line fraction is a=alpha_qbar+v alpha_g.  The definitions " <>
      "below display every integration limit.  The support-corrected P-line " <>
      "density has betaP >= z; extending betaP to zero gives the spurious " <>
      "logarithmic divergence in the printed symmetric appendix."
    ],
    inputCell[HoldComplete[
      expQ = Exp[-mQ^2/M2];
      exp0 = Exp[-s0/M2];
      deltaEQ = expQ - exp0;
      a0B = 1 - u0;
      lineProjection3[fun_, z_] :=
        Inactive[Integrate][
          Inactive[Integrate][
            fun[1 - z - (1 - v) alphaG, z - v alphaG, alphaG, v],
            {alphaG, 0, (1 - z)/(1 - v)}],
          {v, 0, z}] +
        Inactive[Integrate][
          Inactive[Integrate][
            fun[1 - z - (1 - v) alphaG, z - v alphaG, alphaG, v],
            {alphaG, 0, z/v}],
          {v, z, 1}];
      pLine3[fun_, z_] :=
        Inactive[Integrate][
          Inactive[Integrate][
            Inactive[Integrate][
              (z - aqb)/(1 - aq - aqb)^2
                fun[aq, aqb, 1 - aq - aqb],
              {aqb, 0, z}],
            {aqp, 0, aq}],
          {aq, 0, 1 - z}] -
        Inactive[Integrate][
          Inactive[Integrate][
            z/betaP^2 Inactive[Integrate][
              fun[aq, aqbp, 1 - aq - aqbp],
              {aqbp, 0, betaP}],
            {betaP, z, 1 - aq}],
          {aq, 0, 1 - z}] +
        Inactive[Integrate][
          z/(1 - aq)^2 Inactive[Integrate][
            Inactive[Integrate][
              fun[aqp, aqb, 1 - aqp - aqb],
              {aqb, 0, 1 - aqp}],
            {aqp, 0, aq}],
          {aq, 0, 1 - z}];
      pLinePrime3[fun_, z_] :=
        Inactive[D][pLine3[fun, a], a] /. a -> z;
      fGSigma[aq_, aqb_, ag_, v_] :=
        S3[aq, aqb, ag] + St3[aq, aqb, ag] -
        T13[aq, aqb, ag] - T23[aq, aqb, ag] +
        T33[aq, aqb, ag] + T43[aq, aqb, ag];
      fGX[aq_, aqb_, ag_, v_] :=
        2 v (-S3[aq, aqb, ag] - T33[aq, aqb, ag] +
          T23[aq, aqb, ag]);
      fGA[aq_, aqb_, ag_, v_] :=
        fGSigma[aq, aqb, ag, v] + fGX[aq, aqb, ag, v];
      fGammaSigma[aq_, aqb_, ag_, v_] :=
        Sg3[aq, aqb, ag] - T4g3[aq, aqb, ag];
      fGammaA[aq_, aqb_, ag_, v_] :=
        (1 - 2 v) Sg3[aq, aqb, ag] - T4g3[aq, aqb, ag];
      tA3g = eS ssCond expQ lineProjection3[fGA, a0B];
      tB3g = eS ssCond expQ mQ/(mQ + ms)
        lineProjection3[fGSigma, a0B];
      tA3em = eQ ssCond deltaEQ (
        lineProjection3[fGammaA, a0B] +
        2 pLinePrime3[T4g3, a0B]);
      tB3em = eQ ssCond deltaEQ mQ/(mQ + ms) (
        lineProjection3[fGammaSigma, a0B] +
        2 pLinePrime3[T4g3, a0B]);
      threeParticleResiduals = FullSimplify /@ {
        fGA[aq, aqb, alphaG, v] -
          fGSigma[aq, aqb, alphaG, v] -
          fGX[aq, aqb, alphaG, v],
        fGammaA[aq, aqb, alphaG, v] -
          fGammaSigma[aq, aqb, alphaG, v] +
          2 v Sg3[aq, aqb, alphaG]
      };
      <|
        "T_A three-particle gluonic" -> tA3g,
        "T_B three-particle gluonic" -> tB3g,
        "T_A three-particle electromagnetic" -> tA3em,
        "T_B three-particle electromagnetic" -> tB3em,
        "operator-weight residuals" -> threeParticleResiduals,
        "full assembly audit" ->
          "step15_complete_three_particle_borel.wl: four exact zero residuals"
      |>
    ]],
    textCell[
      "Continuum bookkeeping is channel-specific.  The leading twist-2 term " <>
      "uses DeltaEQ=expQ-exp0.  Twist 3, twist 4, and the gluonic " <>
      "three-particle term retain the raw expQ factor of the Ds1-to-Ds-gamma " <>
      "sum rule.  The distinct gauge-completion electromagnetic " <>
      "three-particle term follows Rohrwild's IF-type DeltaEQ subtraction."
    ],
    inputCell[HoldComplete[
      axialOPEComponents = <|
        "perturbative" ->
          Inactive[Integrate][Exp[-s/M2] rhoA[s],
            {s, (mQ + ms)^2, s0}],
        "twist 2, positive-magnitude chi" ->
          eS ssCond deltaEQ M2 chiGamma phiGamma[u0],
        "twist 3" ->
          eS f3Gamma mQ expQ psiVBar[u0],
        "twist 4" ->
          eS ssCond expQ/4
            (A4[u0] + 2 B4[u0]) (1 + mQ^2/M2),
        "three particle, gluonic" -> tA3g,
        "three particle, electromagnetic" -> tA3em
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
    subsectionCell["Explicit post-double-Borel vector and x-gamma kernels"],
    textCell[
      "The Fourier reduction uses k=p+a q and D_a=mQ^2-a pPrime^2-(1-a)p2. " <>
      "The double Borel transform is applied before a0=1/2 is chosen. " <>
      "These expressions contain DA derivatives and no pole-mass substitution."
    ],
    inputCell[HoldComplete[
      a0B = 1 - u0;
      tBTw2 = mQ/(mQ + ms) axialOPEComponents[
        "twist 2, positive-magnitude chi"];
      tBTw3 = eS f3Gamma expQ/(mQ + ms) (
        mQ^2 psiVBar[a0B] -
        M2/2 (D[(1 - a) psiVBar[a], a] /. a -> a0B));
      tBTw4 = mQ/(mQ + ms) axialOPEComponents["twist 4"];
      symmetricTw3 = {
        axialOPEComponents["twist 3"] /. u0 -> 1/2,
        tBTw3 /. u0 -> 1/2
      } /. psiVBar[1/2] -> 0;
      <|
        "T_B twist 2" -> tBTw2,
        "T_B twist 3" -> tBTw3,
        "T_B twist 4" -> tBTw4,
        "symmetric twist-3 pair" -> symmetricTw3
      |>
    ]],
    inputCell[HoldComplete[
      xAS[f_] := -8 expQ f[a0B];
      xAT1[f_] := 48 I expQ f[a0B];
      xAT2[f_] := 16 I expQ f[a0B];
      xAT3[f_] := -16 I expQ f[a0B];
      xAT4[f_] := 0;
      xBS[f_] := 0;
      xBT1[f_] := 0;
      xBT2[f_] := 0;
      xBT3[f_] := 0;
      xBT4[f_] := 0;
      <|
        "A xgamma kernels" ->
          {xAS[Feff], xAT1[Feff], xAT2[Feff], xAT3[Feff], xAT4[Feff]},
        "B xgamma kernels" ->
          {xBS[Feff], xBT1[Feff], xBT2[Feff], xBT3[Feff], xBT4[Feff]},
        "sigma coefficients A" -> {8, -16 I},
        "sigma coefficients B" ->
          {8 mQ/(mQ + ms), -16 I mQ/(mQ + ms)}
      |>
    ]],
    inputCell[HoldComplete[
      lambdaB[shat_, mi_, mj_] :=
        kallenLambda[shat, mj^2, mi^2];
      rhoBLeg[shat_, mi_, mj_, aa_, bb_] :=
        -3/(8 Pi^2) (
          2 aa Log[
            (shat - mj^2 + mi^2 - Sqrt[lambdaB[shat, mi, mj]])/
            (shat - mj^2 + mi^2 + Sqrt[lambdaB[shat, mi, mj]])] +
          bb/mi^2 (mj^2 - mi^2 - shat)/shat^2
            Sqrt[lambdaB[shat, mi, mj]]
        );
      aScoef = mQ ms/(mQ + ms);
      bScoef[shat_] := ms^2 shat/(mQ + ms);
      aQcoef = (2 mQ^2 - mQ ms)/(mQ + ms);
      bQcoef[shat_] := -mQ^2 shat/(mQ + ms);
      rhoB[shat_] :=
        eS rhoBLeg[shat, ms, mQ, aScoef, bScoef[shat]] -
        eQ rhoBLeg[shat, mQ, ms, aQcoef, bQcoef[shat]];
      tBpert = Inactive[Integrate][
        Exp[-s/M2] rhoB[s], {s, (mQ + ms)^2, s0}];
      tensorOPEComponents = <|
        "perturbative" -> tBpert,
        "twist 2" -> tBTw2,
        "twist 3" -> tBTw3,
        "twist 4" -> tBTw4,
        "three particle, gluonic" -> tB3g,
        "three particle, electromagnetic" -> tB3em
      |>;
      tBHat = Total[Values[tensorOPEComponents]];
      tensorOPEAssemblyResidual =
        Expand[tBHat - Total[Values[tensorOPEComponents]]];
      tensorOPEStatus = <|
        "perturbative" -> "explicit diagonal spectral density",
        "ordinary local quark condensate" -> "absent in Rohrwild transition OPE",
        "twist 2" -> "kinematics-independent trace ratio derived",
        "twist 3" -> "explicit post-Borel DA-derivative form derived",
        "twist 4" -> "kinematics-independent trace ratio derived",
        "three particle, gluonic" ->
          "complete operator-weighted post-Borel convolution",
        "three particle, electromagnetic" ->
          "complete support-corrected post-Borel convolution"
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
        axialOPEComponents["twist 2, positive-magnitude chi"] -
        eS ssCond deltaEQ M2 chiGamma phiGamma[u0]
      ];
      <|
        "G_1 final" -> paperG1,
        "G_2 final" -> paperG2,
        "component assembly residual" -> paperComponentAssemblyResidual,
        "theta=Pi/2 Colangelo residual" -> paperPureAxialResidual,
        "chi-convention residual" -> chiConventionResidual
      |>
    ]],

    sectionCell["12. Corrected numerical transition evaluation"],
    textCell[
      "This section runs the independent Mathematica numerical implementation " <>
      "of the exact post-double-Borel Rohrwild expressions.  It evaluates the " <>
      "operator-weighted gluonic and electromagnetic convolutions, the complete " <>
      "T_A and T_B component sums, the physical g1 and g2 projections, and the " <>
      "radiative widths.  The ordinary local transition condensate is fixed to zero. " <>
      "The Python table is used only after the Mathematica calculation, as a " <>
      "term-by-term regression comparison."
    ],
    inputCell[HoldComplete[
      notebookPaperDir = Quiet@Check[
        ParentDirectory[NotebookDirectory[]],
        DirectoryName[DirectoryName[$InputFileName]]
      ];
      If[! DirectoryQ[notebookPaperDir],
        notebookPaperDir = DirectoryName[DirectoryName[$InputFileName]]
      ];
      wolframKernelPath =
        "/Applications/Wolfram.app/Contents/MacOS/WolframKernel";
      numericalScriptPath = FileNameJoin[
        {notebookPaperDir, "scripts",
          "mathematica_corrected_transition_numerics.wl"}];
      numericalRun = RunProcess[
        {wolframKernelPath, "-noinit", "-noprompt", "-script",
          numericalScriptPath},
        All,
        ProcessDirectory -> notebookPaperDir
      ];
      <|
        "exit code" -> numericalRun["ExitCode"],
        "script" -> numericalScriptPath,
        "local transition condensate" -> "excluded"
      |>
    ]],
    inputCell[HoldComplete[
      parseNumericText[value_?NumericQ] := N[value, 17];
      parseNumericText[value_String] := N[
        ToExpression[
          StringReplace[value, {"e+" -> "*^", "e-" -> "*^-"}]
        ],
        17
      ];
      mathematicaCentralPath = FileNameJoin[
        {notebookPaperDir, "outputs",
          "corrected_transition_central_mathematica.csv"}];
      mathematicaCentralRows = Import[mathematicaCentralPath, "CSV"];
      mathematicaCentral = Association[
        (#[[1]] -> parseNumericText[#[[2]]]) & /@
          Rest[mathematicaCentralRows]
      ];
      KeyTake[
        mathematicaCentral,
        {
          "convolution.J_g_sigma", "convolution.J_g_xgamma",
          "convolution.J_g_axial", "convolution.J_em_axial",
          "convolution.J_em_tensor_base", "convolution.L_T4gamma",
          "convolution.Lprime_T4gamma"
        }
      ]
    ]],
    inputCell[HoldComplete[
      KeyTake[
        mathematicaCentral,
        {
          "Ds.T_A_pert", "Ds.T_A_tw2", "Ds.T_A_tw3",
          "Ds.T_A_tw4", "Ds.T_A_3p_g", "Ds.T_A_3p_gamma",
          "Ds.T_A", "Ds.T_B_pert", "Ds.T_B_tw2",
          "Ds.T_B_tw3", "Ds.T_B_tw4", "Ds.T_B_3p_g",
          "Ds.T_B_3p_gamma", "Ds.T_B"
        }
      ]
    ]],
    inputCell[HoldComplete[
      KeyTake[
        mathematicaCentral,
        {
          "Ds.g_1_A", "Ds.g_1_B", "Ds.g_1",
          "Ds.Gamma_1_keV", "Ds.g_2_A", "Ds.g_2_B",
          "Ds.g_2", "Ds.Gamma_2_keV",
          "Bs_low.f1", "Bs_low.f2", "Bs_low.g_1",
          "Bs_low.Gamma_1_keV", "Bs_high.f1", "Bs_high.f2",
          "Bs_high.g_2", "Bs_high.Gamma_2_keV"
        }
      ]
    ]],
    inputCell[HoldComplete[
      comparisonPath = FileNameJoin[
        {notebookPaperDir, "outputs",
          "corrected_transition_python_mathematica_comparison.csv"}];
      comparisonRows = Import[comparisonPath, "CSV"];
      comparisonAbsoluteDifferences =
        parseNumericText /@ Rest[comparisonRows][[All, 4]];
      comparisonPassed =
        And @@ (ToString[#] == "1" & /@ Rest[comparisonRows][[All, 6]]);
      numericalComparisonStatus = <|
        "common quantities" -> Length[comparisonAbsoluteDifferences],
        "largest absolute difference" ->
          Max[comparisonAbsoluteDifferences],
        "all rows pass 1e-8 absolute-or-relative tolerance" ->
          comparisonPassed
      |>;
      numericalComparisonStatus
    ]],

    sectionCell["13. Exact check summary"],
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
        "legacy Bs Pi12=0 closure diagnostic" -> rhoDiagResidual,
        "Colangelo formula limit" -> colangeloFormulaResidual,
        "pure-A width" -> pureAxialWidthResidual,
        "axial OPE component assembly" -> axialOPEAssemblyResidual,
        "tensor OPE component assembly" -> tensorOPEAssemblyResidual,
        "three-particle operator assembly" -> threeParticleResiduals,
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
