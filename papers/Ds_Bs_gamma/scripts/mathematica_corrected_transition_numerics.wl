(* ::Package:: *)

(* Independent Mathematica numerical implementation of the exact post-double-
   Borel Rohrwild transition invariants.

   Run from papers/Ds_Bs_gamma:

     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noinit \
       -noprompt -script scripts/mathematica_corrected_transition_numerics.wl

   The output is a flat key/value CSV matching
   outputs/corrected_transition_central_python.csv.
*)

ClearAll["Global`*"];

wp = 30;
agOf[aq_, aqb_] := 1 - aq - aqb;

kappa = 3/20;
kappaP = -1/20;
zeta1 = 2/5;
zeta1P = 0;
zeta2 = 3/10;
zeta2P = 0;
omegaA = -21/10;
omegaV = 19/5;

phiGamma[u_] := 6 u (1 - u);

aTw4[u_] := Module[{ub = 1 - u, bracket},
  bracket =
    u ub (2 + 13 u ub) +
    2 u^3 (10 - 15 u + 6 u^2) Log[u] +
    2 ub^3 (10 - 15 ub + 6 ub^2) Log[ub];
  40 u^2 ub^2 (3 kappa - kappaP + 1) +
    8 (zeta2P - 3 zeta2) bracket
];

bTw4[u_] := Module[{inner},
  inner =
    u (2 u^3 - 3 u^2 + u) -
    (3 u^4/2 - 2 u^3 + u^2/2);
  40 (1 + 2 kappaP) inner
];

barPsiV[u_] :=
  -20 u (1 - u) (2 u - 1) +
  15/16 (omegaA - 3 omegaV) u (1 - u) (2 u - 1)
    (7 (2 u - 1)^2 - 3);

barPsiVPrime[u_] := Module[{t, duub, basePrime, shapePrime, coeff},
  t = 2 u - 1;
  duub = 1 - 2 u;
  basePrime = duub t + 2 u (1 - u);
  shapePrime =
    basePrime (7 t^2 - 3) +
    28 u (1 - u) t^2;
  coeff = 15/16 (omegaA - 3 omegaV);
  -20 basePrime + coeff shapePrime
];

s3[aq_, aqb_] := Module[{ag = agOf[aq, aqb]},
  30 ag^2 (
    (kappa + kappaP) (1 - ag) +
    (zeta1 + zeta1P) (1 - ag) (1 - 2 ag) +
    zeta2 (3 (aq - aqb)^2 - ag (1 - ag))
  )
];

st3[aq_, aqb_] := Module[{ag = agOf[aq, aqb]},
  -30 ag^2 (
    (kappa - kappaP) (1 - ag) +
    (zeta1 - zeta1P) (1 - ag) (1 - 2 ag) +
    zeta2 (3 (aq - aqb)^2 - ag (1 - ag))
  )
];

t13[aq_, aqb_] := Module[{ag = agOf[aq, aqb]},
  -120 (3 zeta2 + zeta2P) (aq - aqb) aq aqb ag
];

t23[aq_, aqb_] := Module[{ag = agOf[aq, aqb]},
  30 ag^2 (aq - aqb) (
    (kappa - kappaP) +
    (zeta1 - zeta1P) (1 - 2 ag) +
    zeta2 (3 - 4 ag)
  )
];

t33[aq_, aqb_] := Module[{ag = agOf[aq, aqb]},
  -120 (3 zeta2 - zeta2P) (aq - aqb) aq aqb ag
];

t43[aq_, aqb_] := Module[{ag = agOf[aq, aqb]},
  30 ag^2 (aq - aqb) (
    (kappa + kappaP) +
    (zeta1 + zeta1P) (1 - 2 ag) +
    zeta2 (3 - 4 ag)
  )
];

sg3[aq_, aqb_] := Module[{sum = aq + aqb, ag = agOf[aq, aqb]},
  60 ag^2 sum (4 - 7 sum)
];

t4g3[aq_, aqb_] := Module[{sum = aq + aqb, ag = agOf[aq, aqb]},
  60 ag^2 (aq - aqb) (4 - 7 sum)
];

fGSigma[aq_, aqb_, ag_, v_] :=
  s3[aq, aqb] + st3[aq, aqb] -
  t13[aq, aqb] - t23[aq, aqb] +
  t33[aq, aqb] + t43[aq, aqb];

fGXgamma[aq_, aqb_, ag_, v_] :=
  2 v (-s3[aq, aqb] - t33[aq, aqb] + t23[aq, aqb]);

fGAxial[aq_, aqb_, ag_, v_] :=
  fGSigma[aq, aqb, ag, v] + fGXgamma[aq, aqb, ag, v];

fEMAxial[aq_, aqb_, ag_, v_] :=
  (1 - 2 v) sg3[aq, aqb] - t4g3[aq, aqb];

fEMTensorBase[aq_, aqb_, ag_, v_] :=
  sg3[aq, aqb] - t4g3[aq, aqb];

Needs["NumericalDifferentialEquationAnalysis`"];
lineRule = N[
  NumericalDifferentialEquationAnalysis`GaussianQuadratureWeights[
    20, 0, 1
  ],
  18
];
pRule = N[
  NumericalDifferentialEquationAnalysis`GaussianQuadratureWeights[
    16, 0, 1
  ],
  18
];
Print["RULES_READY"];
fixedQuad[fun_, lo_?NumericQ, hi_?NumericQ, rule_List] := Module[
  {span = hi - lo},
  If[span <= 0, 0,
    span Total[
      (#[[2]] fun[lo + span #[[1]]]) & /@ rule
    ]
  ]
];

lineProjection[fun_, z_?NumericQ] := Module[
  {domain1, domain2, integrateAg},
  integrateAg[v_?NumericQ, upper_?NumericQ] :=
    fixedQuad[
      Function[{ag},
        fun[
          1 - z - (1 - v) ag,
          z - v ag,
          ag,
          v
        ]
      ],
      0,
      upper,
      lineRule
    ];
  domain1 = fixedQuad[
    Function[{v}, integrateAg[v, (1 - z)/(1 - v)]],
    0,
    z,
    lineRule
  ];
  domain2 = fixedQuad[
    Function[{v}, integrateAg[v, z/v]],
    z,
    1,
    lineRule
  ];
  domain1 + domain2
];

(* Direct numerical implementation of the support-corrected L_z density.
   In the first term T4gamma/ag^2 is cancelled analytically. *)
pLineT4[z_?NumericQ] := Module[{term1, term2, term3},
  term1 = fixedQuad[
    Function[{aq},
      fixedQuad[
        Function[{aqb},
          aq (z - aqb) 60 (aq - aqb) (4 - 7 (aq + aqb))
        ],
        0,
        z,
        pRule
      ]
    ],
    0,
    1 - z,
    pRule
  ];
  term2 = -fixedQuad[
    Function[{aq},
      fixedQuad[
        Function[{beta},
          z/beta^2 fixedQuad[
            Function[{aqbp}, t4g3[aq, aqbp]],
            0,
            beta,
            pRule
          ]
        ],
        z,
        1 - aq,
        pRule
      ]
    ],
    0,
    1 - z,
    pRule
  ];
  term3 = fixedQuad[
    Function[{aq},
      z/(1 - aq)^2 fixedQuad[
        Function[{aqp},
          fixedQuad[
            Function[{aqb}, t4g3[aqp, aqb]],
            0,
            1 - aqp,
            pRule
          ]
        ],
        0,
        aq,
        pRule
      ]
    ],
    0,
    1 - z,
    pRule
  ];
  term1 + term2 + term3
];

pLinePrimeT4[z_?NumericQ, h_: 2/10000] := (
  pLineT4[z - 2 h] -
  8 pLineT4[z - h] +
  8 pLineT4[z + h] -
  pLineT4[z + 2 h]
)/(12 h);

a0 = 1/2;
jGSigma = lineProjection[fGSigma, a0]; Print["J_G_SIGMA_READY"];
jGXgamma = lineProjection[fGXgamma, a0]; Print["J_G_X_READY"];
jGAxial = lineProjection[fGAxial, a0]; Print["J_G_A_READY"];
jEMAxial = lineProjection[fEMAxial, a0]; Print["J_EM_A_READY"];
jEMTensor = lineProjection[fEMTensorBase, a0]; Print["J_EM_B_READY"];
lT4 = pLineT4[a0]; Print["L_T4_READY"];
lPrimeT4 = pLinePrimeT4[a0]; Print["LPRIME_T4_READY"];
convolutions = <|
  "z" -> N[a0, 20],
  "J_g_sigma" -> jGSigma,
  "J_g_xgamma" -> jGXgamma,
  "J_g_axial" -> jGAxial,
  "J_em_axial" -> jEMAxial,
  "J_em_tensor_base" -> jEMTensor,
  "L_T4gamma" -> lT4,
  "Lprime_T4gamma" -> lPrimeT4
|>;
convolutions["J_g_closure"] =
  convolutions["J_g_axial"] -
  convolutions["J_g_sigma"] -
  convolutions["J_g_xgamma"];

kallen[s_, mi_, mj_] := (s - (mi + mj)^2) (s - (mi - mj)^2);

rhoPiece[s_?NumericQ, mi_, mj_, charge_] := Module[
  {lam, logTerm, braces},
  lam = Sqrt[Max[kallen[s, mj, mi], 0]];
  logTerm = Log[
    (s - mj^2 + mi^2 - lam)/
    (s - mj^2 + mi^2 + lam)
  ];
  braces =
    2 mi logTerm +
    (mj - mi) (mj^2 - mi^2 - s) lam/s^2;
  -3 charge braces/(8 Pi^2)
];

rhoA[s_?NumericQ, mQ_, ms_, eQ_, es_] :=
  rhoPiece[s, ms, mQ, es] - rhoPiece[s, mQ, ms, eQ];

lineKernel[s_?NumericQ, mi_, mj_, aa_, bb_] := Module[
  {lam, logTerm},
  lam = Sqrt[Max[kallen[s, mj, mi], 0]];
  logTerm = Log[
    (s - mj^2 + mi^2 - lam)/
    (s - mj^2 + mi^2 + lam)
  ];
  -3/(8 Pi^2) (
    2 aa logTerm +
    (bb/mi^2) (mj^2 - mi^2 - s) lam/s^2
  )
];

rhoB[s_?NumericQ, mQ_, ms_, eQ_, es_] := Module[
  {d = mQ + ms, as, bs, aQ, bQ},
  as = mQ ms/d;
  bs = ms^2 s/d;
  aQ = (2 mQ^2 - mQ ms)/d;
  bQ = -mQ^2 s/d;
  es lineKernel[s, ms, mQ, as, bs] -
  eQ lineKernel[s, mQ, ms, aQ, bQ]
];

transitionInvariants[params_Association] := Module[
  {
    M2 = params["M2"], s0 = params["s0"],
    mQ = params["mQ"], ms = params["ms"],
    eQ = params["eQ"], es = params["es"],
    ss = params["ss"], chi = params["chi"],
    f3g = params["f3g"], d, rQ, threshold, EQ, E0, deltaEQ,
    taPert, tbPert, taTw2, taTw3, taTw4,
    tbTw2, tbTw3, tbTw4, ta3g, tb3g, ta3em, tb3em,
    barPsi, barPsiPrime, derivativeProduct
  },
  d = mQ + ms;
  rQ = mQ/d;
  threshold = d^2;
  EQ = Exp[-mQ^2/M2];
  E0 = Exp[-s0/M2];
  deltaEQ = EQ - E0;
  taPert = NIntegrate[
    Exp[-spectralS/M2] rhoA[spectralS, mQ, ms, eQ, es],
    {spectralS, threshold, s0},
    WorkingPrecision -> 25,
    AccuracyGoal -> 16,
    PrecisionGoal -> 14,
    MaxRecursion -> 20,
    Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
  ];
  tbPert = NIntegrate[
    Exp[-spectralS/M2] rhoB[spectralS, mQ, ms, eQ, es],
    {spectralS, threshold, s0},
    WorkingPrecision -> 25,
    AccuracyGoal -> 16,
    PrecisionGoal -> 14,
    MaxRecursion -> 20,
    Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
  ];
  taTw2 =
    es ss deltaEQ M2 chi phiGamma[a0];
  taTw4 =
    es ss EQ/4 (aTw4[a0] + 2 bTw4[a0])
      (1 + mQ^2/M2);
  barPsi = barPsiV[a0];
  barPsiPrime = barPsiVPrime[a0];
  taTw3 = es f3g mQ EQ barPsi;
  tbTw2 = rQ taTw2;
  tbTw4 = rQ taTw4;
  derivativeProduct = -barPsi + (1 - a0) barPsiPrime;
  tbTw3 =
    es f3g EQ/d (mQ^2 barPsi - M2 derivativeProduct/2);
  ta3g =
    es ss EQ convolutions["J_g_axial"];
  tb3g =
    es ss EQ rQ convolutions["J_g_sigma"];
  ta3em =
    eQ ss deltaEQ (
      convolutions["J_em_axial"] +
      2 convolutions["Lprime_T4gamma"]
    );
  tb3em =
    eQ ss deltaEQ rQ (
      convolutions["J_em_tensor_base"] +
      2 convolutions["Lprime_T4gamma"]
    );
  Join[<|
    "ordinary_local_transition_condensate" -> 0,
    "M2" -> M2,
    "s0" -> s0,
    "a0" -> a0,
    "E_Q" -> EQ,
    "E_0" -> E0,
    "Delta_E_Q" -> deltaEQ,
    "T_A_pert" -> taPert,
    "T_A_tw2" -> taTw2,
    "T_A_tw3" -> taTw3,
    "T_A_tw4" -> taTw4,
    "T_A_3p_g" -> ta3g,
    "T_A_3p_gamma" -> ta3em,
    "T_A" -> Total[{taPert, taTw2, taTw3, taTw4, ta3g, ta3em}],
    "T_B_pert" -> tbPert,
    "T_B_tw2" -> tbTw2,
    "T_B_tw3" -> tbTw3,
    "T_B_tw4" -> tbTw4,
    "T_B_3p_g" -> tb3g,
    "T_B_3p_gamma" -> tb3em,
    "T_B" -> Total[{tbPert, tbTw2, tbTw3, tbTw4, tb3g, tb3em}],
    "bar_psi_V" -> barPsi,
    "bar_psi_V_prime" -> barPsiPrime
  |>, convolutions]
];

widthKeV[mInitial_, mFinal_, coupling_] := Module[{alpha, qGamma},
  alpha = 1/137.036;
  qGamma = (mInitial^2 - mFinal^2)/(2 mInitial);
  alpha coupling^2 qGamma^3 10^6/3
];

physicalCouplings[
  inv_Association, thetaDeg_, m1_, m2_, f1_, f2_,
  mP_, fP_, mQ_, ms_
] := Module[
  {theta, nn1, nn2, ta, tb, g1a, g1b, g2a, g2b, g1, g2},
  theta = thetaDeg Degree;
  ta = inv["T_A"];
  tb = inv["T_B"];
  nn1 =
    Exp[(m1^2 + mP^2)/(2 inv["M2"])] (mQ + ms)/
    (m1 f1 mP^2 fP);
  nn2 =
    Exp[(m2^2 + mP^2)/(2 inv["M2"])] (mQ + ms)/
    (m2 f2 mP^2 fP);
  g1a = nn1 Sin[theta] ta;
  g1b = nn1 Cos[theta] tb;
  g2a = nn2 Cos[theta] ta;
  g2b = -nn2 Sin[theta] tb;
  g1 = g1a + g1b;
  g2 = g2a + g2b;
  <|
    "N_1" -> nn1,
    "N_2" -> nn2,
    "g_1" -> g1,
    "g_2" -> g2,
    "g_1_A" -> g1a,
    "g_1_B" -> g1b,
    "g_2_A" -> g2a,
    "g_2_B" -> g2b,
    "Gamma_1_keV" -> widthKeV[m1, mP, g1],
    "Gamma_2_keV" -> widthKeV[m2, mP, g2]
  |>
];

ssCentral = SetPrecision[0.8 (-(0.240)^3), wp];
fPerpS = SetPrecision[-0.0510 1.08, wp];
chiLattice = fPerpS/ssCentral;

dsParams = <|
  "M2" -> SetPrecision[3.75, wp],
  "s0" -> SetPrecision[8.0, wp],
  "mQ" -> SetPrecision[1.27, wp],
  "ms" -> SetPrecision[0.093, wp],
  "eQ" -> 2/3,
  "es" -> -1/3,
  "ss" -> ssCentral,
  "chi" -> chiLattice,
  "f3g" -> SetPrecision[-0.0039, wp]
|>;

dsInv = transitionInvariants[dsParams];
Print["DS_INVARIANTS_READY"];
dsPhys = physicalCouplings[
  dsInv,
  SetPrecision[26.6, wp],
  SetPrecision[2.4595, wp],
  SetPrecision[2.53511, wp],
  SetPrecision[0.40058930493635564, wp],
  SetPrecision[0.16667036765413243, wp],
  SetPrecision[1.96835, wp],
  SetPrecision[0.2499, wp],
  dsParams["mQ"],
  dsParams["ms"]
];

bsTwoPointConstants[] := Module[
  {mb, msb, m0sq, mixed, mLow, mHigh, m2tp, s01, s02, theta,
   den, threshold, kallen, rhoAA, rhoAB, rhoBB, rho, q0, q1, q2,
   d5, local, rotation, v1, v2, pi1, pi2, f1, f2},
  mb = SetPrecision[4.18, wp];
  msb = SetPrecision[0.093, wp];
  m0sq = SetPrecision[0.8, wp];
  mixed = m0sq ssCentral;
  mLow = SetPrecision[5.750, wp];
  mHigh = SetPrecision[5.82870, wp];
  m2tp = SetPrecision[7.0, wp];
  s01 = SetPrecision[44.35, wp];
  s02 = SetPrecision[43.025, wp];
  theta = SetPrecision[38.5, wp] Degree;
  den = mb + msb;
  threshold = den^2;
  kallen[z_] := (z - (mb + msb)^2) (z - (mb - msb)^2);
  rhoAA[z_] := Sqrt[kallen[z]]/(8 Pi^2) (
    2 - (mb^2 + msb^2 + 6 mb msb)/z -
      (mb^2 - msb^2)^2/z^2);
  rhoAB[z_] := Sqrt[kallen[z]]/(8 Pi^2) (
    3 (msb - mb) ((mb + msb)^2 - z)/(den z));
  rhoBB[z_] := Sqrt[kallen[z]]/(8 Pi^2) (
    -((mb + msb)^2 - z)
      (2 msb^2 - 4 mb msb + 2 mb^2 + z)/(den^2 z));
  rho[z_] := {{rhoAA[z], rhoAB[z]}, {rhoAB[z], rhoBB[z]}};
  q0 = ssCentral {{mb, mb^2/den}, {mb^2/den, mb^3/den^2}};
  q1 = ssCentral {{
      -msb mb^2/(2 m2tp),
      msb mb (1 - mb^2/m2tp)/(2 den)
    }, {
      msb mb (1 - mb^2/m2tp)/(2 den),
      msb mb^2 (1 - mb^2/(2 m2tp))/den^2
    }};
  q2 = ssCentral {{
      msb^2 mb^3/(2 m2tp^2),
      msb^2 (mb^4/m2tp^2 - mb^2/m2tp - 1)/(2 den)
    }, {
      msb^2 (mb^4/m2tp^2 - mb^2/m2tp - 1)/(2 den),
      msb^2 (mb^5/(2 m2tp^2) - mb^3/m2tp)/den^2
    }};
  d5 = mixed {{
      -mb^3/(4 m2tp^2),
      -(mb^4/m2tp^2 - mb^2/m2tp - 1)/(4 den)
    }, {
      -(mb^4/m2tp^2 - mb^2/m2tp - 1)/(4 den),
      -(mb^5/(4 m2tp^2) - mb^3/(2 m2tp))/den^2
    }};
  local = Exp[-mb^2/m2tp] (q0 + q1 + q2 + d5);
  rotation = {{Sin[theta], Cos[theta]}, {Cos[theta], -Sin[theta]}};
  v1 = rotation[[1]];
  v2 = rotation[[2]];
  pi1 = NIntegrate[
      v1 . rho[z] . v1 Exp[-z/m2tp], {z, threshold, s01},
      AccuracyGoal -> 12, PrecisionGoal -> 12] +
    v1 . local . v1;
  pi2 = NIntegrate[
      v2 . rho[z] . v2 Exp[-z/m2tp], {z, threshold, s02},
      AccuracyGoal -> 12, PrecisionGoal -> 12] +
    v2 . local . v2;
  f1 = Sqrt[pi1 Exp[mLow^2/m2tp]/mLow^2];
  f2 = Sqrt[pi2 Exp[mHigh^2/m2tp]/mHigh^2];
  <|"f1" -> f1, "f2" -> f2, "two_point_M2" -> m2tp,
    "two_point_s01" -> s01, "two_point_s02" -> s02|>
];

bsConstants = bsTwoPointConstants[];
Print["BS_TWOPOINT_READY"];

computeBs[label_, mass_, s0_] := Module[
  {params, inv, phys},
  params = <|
    "M2" -> SetPrecision[12.0, wp],
    "s0" -> SetPrecision[s0, wp],
    "mQ" -> SetPrecision[4.18, wp],
    "ms" -> SetPrecision[0.093, wp],
    "eQ" -> -1/3,
    "es" -> -1/3,
    "ss" -> ssCentral,
    "chi" -> chiLattice,
    "f3g" -> SetPrecision[-0.0039, wp]
  |>;
  inv = transitionInvariants[params];
  phys = physicalCouplings[
    inv,
    SetPrecision[38.5, wp],
    SetPrecision[5.750, wp],
    SetPrecision[5.82870, wp],
    bsConstants["f1"],
    bsConstants["f2"],
    SetPrecision[5.36692, wp],
    SetPrecision[0.2303, wp],
    params["mQ"],
    params["ms"]
  ];
  <|
    "label" -> label,
    "invariants" -> inv,
    "physical" -> Join[
      phys,
      bsConstants
    ]
  |>
];

bsLow = computeBs["Bs_low", 5.750, 38.5];
Print["BS_LOW_READY"];
bsHigh = computeBs["Bs_high", 5.82870, 39.5];
Print["BS_HIGH_READY"];

associationRows[prefix_, assoc_Association] := Module[{numericKeys},
  numericKeys = Select[Keys[assoc], NumericQ[Lookup[assoc, #]] &];
  Table[
    {prefix <> "." <> ToString[key], N[Lookup[assoc, key], 17]},
    {key, numericKeys}
  ]
];

rows = Join[
  associationRows["convolution", convolutions],
  associationRows["Ds", dsInv],
  associationRows["Ds", dsPhys],
  associationRows["Bs_low", bsLow["invariants"]],
  associationRows["Bs_low", bsLow["physical"]],
  associationRows["Bs_high", bsHigh["invariants"]],
  associationRows["Bs_high", bsHigh["physical"]]
];
Print["ROWS_READY=", Length[rows]];

outDir = FileNameJoin[{Directory[], "outputs"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir]];
csvPath = FileNameJoin[{outDir, "corrected_transition_central_mathematica.csv"}];
txtPath = FileNameJoin[{outDir, "corrected_transition_central_mathematica.txt"}];
csvStream = OpenWrite[csvPath];
WriteString[csvStream, "key,value\n"];
Scan[
  WriteString[
    csvStream,
    #[[1]], ",", ToString[#[[2]], InputForm], "\n"
  ] &,
  rows
];
Close[csvStream];
Print["CSV_READY"];

stream = OpenWrite[txtPath];
WriteString[stream,
  "Corrected Rohrwild central numerical reference (Mathematica)\n",
  "============================================================\n",
  "Transition local condensate: excluded.\n",
  "Tensor/EM numerator prescription: exact off-shell double Borel.\n\n",
  "Convolutions:\n",
  ToString[N[convolutions, 15], InputForm], "\n\n",
  "Ds invariants:\n",
  ToString[N[dsInv, 15], InputForm], "\n\n",
  "Ds physical:\n",
  ToString[N[dsPhys, 15], InputForm], "\n\n",
  "Bs low physical:\n",
  ToString[N[bsLow["physical"], 15], InputForm], "\n\n",
  "Bs high physical:\n",
  ToString[N[bsHigh["physical"], 15], InputForm], "\n"
];
Close[stream];

Print["CONVOLUTIONS=", N[convolutions, 12]];
Print["DS_PHYSICAL=", N[dsPhys, 12]];
Print["BS_LOW_PHYSICAL=", N[bsLow["physical"], 12]];
Print["BS_HIGH_PHYSICAL=", N[bsHigh["physical"], 12]];
Print["WROTE=", csvPath];
Quit[];
