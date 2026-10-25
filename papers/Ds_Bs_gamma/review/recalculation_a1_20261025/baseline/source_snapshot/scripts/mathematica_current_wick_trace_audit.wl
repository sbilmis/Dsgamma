(* ::Package:: *)

(* Independent Mathematica/FeynCalc audit for Ds1/Bs1 -> Ds/Bs gamma.

   This script rebuilds the symbolic beginning of the calculation in one
   place.  It deliberately stops before loop integration, dispersion
   relations, continuum subtraction, and Borel transformation.

   The correlator ordering is the standard one used in hep-ph/0505195,

     T_mu = i int exp(i p.x) <gamma|T{J_P^dagger(x) J_X,mu(0)}|0>,

   so the final pseudoscalar channel carries p and the initial axial channel
   carries pPrime=p+q.  The chain checked here is

     interpolating currents
       -> Wick contraction and fermion signs
       -> free propagators and one-photon insertions
       -> hard-emission Dirac traces
       -> E1 projection and inverse-denominator reduction
       -> soft two-particle Fierz traces.

   No Python source or numerical output is used to construct the algebra.

   Run from any directory with

     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script papers/Ds_Bs_gamma/scripts/mathematica_current_wick_trace_audit.wl
*)

<< FeynCalc`

ClearAll[
  mu, nu, alpha, beta, rho, lambda, p, q, pPrime, k, eps, ms, mQ, mb, eS, eQ,
  theta, Nc, p2, pq, dH1, dH2, dH3, dS1, dS2, dS3,
  gammaA, gammaB, gammaBExpanded, gammaP, gammaLow, gammaHigh,
  propNumerator, heavyTrace, strangeTrace, traceReduce,
  e1Tensor, e1Projector, e1Project, momentumRules,
  traceForSoft, projectSoftE1, tensorPhoton, vectorPhoton
];

scriptDir = DirectoryName[$InputFileName];
paperDir = DirectoryName[scriptDir];
outDir = FileNameJoin[{paperDir, "outputs"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outFile = FileNameJoin[{outDir, "mathematica_current_wick_trace_audit.txt"}];

(* ---------------------------------------------------------------------- *)
(* 1. Conventions and current vertices                                    *)
(* ---------------------------------------------------------------------- *)

(* pPrime is the momentum of the initial 1+ meson.  It is substituted as p+q
   only after the Dirac algebra, which keeps the tensor-current definition
   visible during the trace. *)

gammaA[index_] := GA[index] . GA[5];

gammaB[index_] :=
  I/(mQ + ms) FV[pPrime, alpha] DiracSigma[GA[index], GA[alpha]] . GA[5];

(* Since sigma_(mu alpha) = i/2 [gamma_mu,gamma_alpha], the extra i in J_B
   turns the vertex into -1/2 times the gamma-matrix commutator.  This second
   definition provides an algebraically independent representation. *)
gammaBExpanded[index_] :=
  -(GA[index] . GS[pPrime] - GS[pPrime] . GA[index]) . GA[5]/(2 (mQ + ms));

(* At x we use J_P^dagger = sbar i gamma5 Q, following the correlator ordering
   above.  At 0 the conjugate axial/tensor vertex has the same Dirac matrix. *)
gammaP[] := I GA[5];

gammaLow[index_] := Sin[theta] gammaA[index] + Cos[theta] gammaB[index];
gammaHigh[index_] := Cos[theta] gammaA[index] - Sin[theta] gammaB[index];

propNumerator[momentum_, mass_] := GS[momentum] + mass;

(* ---------------------------------------------------------------------- *)
(* 2. Wick signs and propagator bookkeeping                               *)
(* ---------------------------------------------------------------------- *)

(* The correlator starts with an overall i.  Contracting both flavors makes
   one closed fermion loop:

     <T sbar GammaP Q Qbar Gamma_X s>
       = -Tr[GammaP S_Q Gamma_X S_s].

   Hence the hard Wick prefactor before propagator/vertex factors is -i Nc.
   Contracting only Q leaves the ordered bilocal sbar ... s and has no closed
   loop, so the corresponding soft prefactor is +i.  The standard Fierz
   identity contributes a further -1/4. *)

fermionLoopSign = -1;
hardWickPrefactor = I fermionLoopSign Nc;
softWickPrefactorBeforeFierz = I;
softWickPrefactorAfterFierz = -I/4;

(* Momentum-space denominators.  They are recorded explicitly even though
   the trace calculation below works on the numerators. *)
heavyDenominators = {
  "dH1 = (k+q)^2 - mQ^2",
  "dH2 = k^2 - mQ^2",
  "dH3 = (k-p)^2 - ms^2"
};
strangeDenominators = {
  "dS1 = k^2 - mQ^2",
  "dS2 = (k-p)^2 - ms^2",
  "dS3 = (k-p-q)^2 - ms^2"
};

(* ---------------------------------------------------------------------- *)
(* 3. Hard photon emission traces                                         *)
(* ---------------------------------------------------------------------- *)

(* Photon emission from the heavy line:

   Tr[i gamma5 (k+mQ) gamma_nu (k+q+mQ) Gamma_X (k-p+ms)].

   Photon emission from the strange line:

   Tr[i gamma5 (k+mQ) Gamma_X (k-p-q+ms) gamma_nu
      (k-p+ms)].

   Charges, color, denominators, the QED vertex, and propagator i factors
   remain outside these numerator traces. *)

traceReduce[expression_] :=
  Collect2[
    DiracSimplify[expression, DiracSubstitute67 -> True],
    {mQ, ms, eQ, eS}
  ];

heavyTrace[current_] :=
  traceReduce[
    DiracTrace[
      gammaP[] .
      propNumerator[k, mQ] .
      GA[nu] .
      propNumerator[k + q, mQ] .
      current[mu] .
      propNumerator[k - p, ms]
    ]
  ];

strangeTrace[current_] :=
  traceReduce[
    DiracTrace[
      gammaP[] .
      propNumerator[k, mQ] .
      current[mu] .
      propNumerator[k - p - q, ms] .
      GA[nu] .
      propNumerator[k - p, ms]
    ]
  ];

nAHeavy = heavyTrace[gammaA];
nAStrange = strangeTrace[gammaA];
nBHeavy = heavyTrace[gammaB];
nBStrange = strangeTrace[gammaB];

nATotal = Collect2[eQ nAHeavy + eS nAStrange, {eQ, eS, mQ, ms}];
nBTotal = Collect2[eQ nBHeavy + eS nBStrange, {eQ, eS, mQ, ms}];

(* Cyclic rotations of each trace provide a second construction of the same
   numerator.  This catches vertex-order mistakes without using saved output. *)
nAHeavyCyclic = traceReduce[
  DiracTrace[
    gammaA[mu] . propNumerator[k - p, ms] . gammaP[] .
    propNumerator[k, mQ] . GA[nu] . propNumerator[k + q, mQ]
  ]
];
nAStrangeCyclic = traceReduce[
  DiracTrace[
    gammaA[mu] . propNumerator[k - p - q, ms] . GA[nu] .
    propNumerator[k - p, ms] . gammaP[] . propNumerator[k, mQ]
  ]
];
nBHeavyCyclic = traceReduce[
  DiracTrace[
    gammaB[mu] . propNumerator[k - p, ms] . gammaP[] .
    propNumerator[k, mQ] . GA[nu] . propNumerator[k + q, mQ]
  ]
];
nBStrangeCyclic = traceReduce[
  DiracTrace[
    gammaB[mu] . propNumerator[k - p - q, ms] . GA[nu] .
    propNumerator[k - p, ms] . gammaP[] . propNumerator[k, mQ]
  ]
];

(* ---------------------------------------------------------------------- *)
(* 4. E1 projector and triangle-core reduction                            *)
(* ---------------------------------------------------------------------- *)

e1Tensor = FV[p, nu] FV[q, mu] - SP[p, q] MT[mu, nu];
e1Projector = e1Tensor/(2 SP[p, q]^2);

momentumRules = {
  SP[q, q] -> 0,
  Pair[Momentum[q], Momentum[q]] -> 0,
  Pair[Momentum[k], Momentum[pPrime]] ->
    Pair[Momentum[k], Momentum[p]] + Pair[Momentum[k], Momentum[q]],
  Pair[Momentum[p], Momentum[pPrime]] ->
    Pair[Momentum[p], Momentum[p]] + Pair[Momentum[p], Momentum[q]],
  Pair[Momentum[pPrime], Momentum[q]] -> Pair[Momentum[p], Momentum[q]]
};

e1Project[expression_] :=
  Collect2[
    ScalarProductExpand[Contract[expression e1Projector] /. momentumRules],
    {mQ, ms, eQ, eS}
  ];

projectorNormalization =
  ScalarProductExpand[
    Contract[e1Tensor e1Projector] /. {
      SP[q, q] -> 0,
      Pair[Momentum[q], Momentum[q]] -> 0
    }
  ] // Simplify;

photonWard =
  ScalarProductExpand[
    Contract[e1Tensor FV[q, nu]] /. {
      SP[q, q] -> 0,
      Pair[Momentum[q], Momentum[q]] -> 0
    }
  ] // Simplify;

initialPolarizationWard =
  ScalarProductExpand[
    Contract[e1Tensor (FV[p, mu] + FV[q, mu]) FV[eps, nu]] /. {
      SP[q, q] -> 0,
      SP[q, eps] -> 0,
      SP[eps, q] -> 0,
      Pair[Momentum[q], Momentum[q]] -> 0,
      Pair[Momentum[q], Momentum[eps]] -> 0,
      Pair[Momentum[eps], Momentum[q]] -> 0
    }
  ] // Simplify;

projAHeavy = e1Project[nAHeavy];
projAStrange = e1Project[nAStrange];
projBHeavy = e1Project[nBHeavy];
projBStrange = e1Project[nBStrange];

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
  Collect2[
    Expand[projected /. rules] /. Thread[denominators -> 0],
    {mQ, ms, p2, pq}
  ];

coreAHeavy = triangleCore[projAHeavy, heavyRules, {dH1, dH2, dH3}];
coreAStrange = triangleCore[projAStrange, strangeRules, {dS1, dS2, dS3}];
coreBHeavy = triangleCore[projBHeavy, heavyRules, {dH1, dH2, dH3}];
coreBStrange = triangleCore[projBStrange, strangeRules, {dS1, dS2, dS3}];

referenceCoreAHeavy =
  4 I mQ + 4 I mQ^3/pq - 4 I mQ^2 ms/pq;
referenceCoreAStrange =
  -4 I ms + 4 I mQ ms^2/pq - 4 I ms^3/pq;
referenceCoreBHeavy =
  8 I mQ^2/(mQ + ms) - 4 I mQ ms/(mQ + ms) +
    4 I mQ^2 p2/((mQ + ms) pq);
referenceCoreBStrange =
  -4 I mQ ms/(mQ + ms) + 8 I ms^2/(mQ + ms) +
    4 I ms^2 p2/((mQ + ms) pq);

(* ---------------------------------------------------------------------- *)
(* 5. Soft two-particle Fierz traces                                      *)
(* ---------------------------------------------------------------------- *)

traceForSoft[current_, basis_] :=
  DiracSimplify[
    DiracTrace[gammaP[] . propNumerator[k, mQ] . current[mu] . basis],
    DiracSubstitute67 -> True
  ] /. momentumRules;

traceForSoftCyclic[current_, basis_] :=
  DiracSimplify[
    DiracTrace[basis . gammaP[] . propNumerator[k, mQ] . current[mu]],
    DiracSubstitute67 -> True
  ] /. momentumRules;

softBasis = {
  {"scalar 1", 1},
  {"pseudoscalar i gamma5", I GA[5]},
  {"vector gamma_rho", GA[rho]},
  {"axial gamma_rho gamma5", GA[rho] . GA[5]},
  {"tensor sigma_rho_lambda", DiracSigma[GA[rho], GA[lambda]]}
};

softTracesA = Table[{entry[[1]], traceForSoft[gammaA, entry[[2]]]}, {entry, softBasis}];
softTracesB = Table[{entry[[1]], traceForSoft[gammaB, entry[[2]]]}, {entry, softBasis}];

tensorPhoton = MT[rho, nu] FV[q, lambda] - MT[lambda, nu] FV[q, rho];
vectorPhoton = MT[rho, nu];

projectSoftE1[expression_] :=
  Collect2[
    ScalarProductExpand[
      Contract[expression e1Projector] /. momentumRules
    ],
    {mQ, ms}
  ];

softATensorE1 = projectSoftE1[
  traceForSoft[gammaA, DiracSigma[GA[rho], GA[lambda]]] tensorPhoton
];
softBTensorE1 = projectSoftE1[
  traceForSoft[gammaB, DiracSigma[GA[rho], GA[lambda]]] tensorPhoton
];
softAVectorE1 = projectSoftE1[traceForSoft[gammaA, GA[rho]] vectorPhoton];
softBVectorE1 = projectSoftE1[traceForSoft[gammaB, GA[rho]] vectorPhoton];

(* ---------------------------------------------------------------------- *)
(* 6. Assertions                                                         *)
(* ---------------------------------------------------------------------- *)

checkDifference[description_, lhs_, rhs_] := Module[{difference},
  difference = Together[Contract[DiracSimplify[lhs - rhs, DiracSubstitute67 -> True]]] // Simplify;
  {description, difference, TrueQ[difference === 0]}
];

checks = {
  checkDifference["FeynCalc Levi-Civita convention", FeynCalc`$LeviCivitaSign, -1],
  checkDifference["tensor vertex equals expanded commutator", gammaB[mu], gammaBExpanded[mu]],
  checkDifference["axial heavy trace is cyclic", nAHeavy, nAHeavyCyclic],
  checkDifference["axial strange trace is cyclic", nAStrange, nAStrangeCyclic],
  checkDifference["tensor heavy trace is cyclic", nBHeavy, nBHeavyCyclic],
  checkDifference["tensor strange trace is cyclic", nBStrange, nBStrangeCyclic],
  checkDifference["photon Ward contraction vanishes", photonWard, 0],
  checkDifference["initial physical-polarization contraction vanishes", initialPolarizationWard, 0],
  checkDifference["E1 projector normalization", projectorNormalization, 1],
  checkDifference["axial heavy triangle core matches routing audit", coreAHeavy, referenceCoreAHeavy],
  checkDifference["axial strange triangle core matches routing audit", coreAStrange, referenceCoreAStrange],
  checkDifference["tensor heavy triangle core matches routing audit", coreBHeavy, referenceCoreBHeavy],
  checkDifference["tensor strange triangle core matches routing audit", coreBStrange, referenceCoreBStrange],
  checkDifference["axial tensor-DA soft trace is cyclic", traceForSoft[gammaA, DiracSigma[GA[rho], GA[lambda]]], traceForSoftCyclic[gammaA, DiracSigma[GA[rho], GA[lambda]]]],
  checkDifference["tensor tensor-DA soft trace is cyclic", traceForSoft[gammaB, DiracSigma[GA[rho], GA[lambda]]], traceForSoftCyclic[gammaB, DiracSigma[GA[rho], GA[lambda]]]],
  checkDifference["axial vector-DA soft trace is cyclic", traceForSoft[gammaA, GA[rho]], traceForSoftCyclic[gammaA, GA[rho]]],
  checkDifference["tensor vector-DA soft trace is cyclic", traceForSoft[gammaB, GA[rho]], traceForSoftCyclic[gammaB, GA[rho]]],
  checkDifference[
    "low-current heavy trace is the stated mixing rotation",
    heavyTrace[gammaLow], Sin[theta] nAHeavy + Cos[theta] nBHeavy
  ],
  checkDifference[
    "high-current strange trace is the stated mixing rotation",
    strangeTrace[gammaHigh], Cos[theta] nAStrange - Sin[theta] nBStrange
  ]
};

allChecksPassed = And @@ (#[[3]] & /@ checks);

(* ---------------------------------------------------------------------- *)
(* 7. Human-readable audit output                                         *)
(* ---------------------------------------------------------------------- *)

stream = OpenWrite[outFile];

writeLine[text_] := WriteString[stream, text, "\n"];
writeExpr[label_, expression_] := (
  WriteString[stream, label, "\n", ToString[expression, InputForm], "\n\n"]
);

writeLine["Mathematica/FeynCalc current, Wick-contraction, and trace audit"];
writeLine["================================================================"];
writeLine[""];
writeLine["Scope: currents -> Wick contraction -> propagators -> traces -> E1 projection."];
writeLine["Loop integration, spectral densities, and Borel transforms are not performed here."];
writeLine[""];

writeLine["1. CURRENT DEFINITIONS"];
writeLine["J_A,mu = sbar gamma_mu gamma5 Q"];
writeLine["J_B,mu = i sbar sigma_(mu alpha) pPrime^alpha gamma5 Q/(mQ+ms)"];
writeLine["Correlator: T{J_P^dagger(x) J_X,mu(0)} with J_P^dagger=sbar i gamma5 Q."];
writeLine["pPrime = p+q; q^2=0 is imposed only in the E1 projection."];
writeLine["Physical-current rotations are linear combinations of J_A and J_B."];
writeLine[""];

writeLine["2. WICK CONTRACTION AND GLOBAL SIGN"];
writeLine["Hard double contraction:"];
writeLine["  <T{sbar GammaP Q Qbar Gamma_X s}> = -Tr[GammaP S_Q Gamma_X S_s]."];
writeExpr["Hard correlator prefactor before propagator/vertex factors (-i Nc):", hardWickPrefactor];
writeExpr["Soft correlator prefactor before Fierz (+i):", softWickPrefactorBeforeFierz];
writeExpr["Soft correlator prefactor after the -1/4 Fierz factor:", softWickPrefactorAfterFierz];
writeLine["The old numerator files intentionally omitted these overall factors."];
writeLine[""];

writeLine["3. PROPAGATOR ROUTING"];
writeLine["Free propagator: S_f^0(r)=i (slash(r)+m_f)/(r^2-m_f^2+i0)."];
writeLine["One-photon insertion is oriented separately on each fermion line; see the routing list below."];
Scan[writeLine, heavyDenominators];
Scan[writeLine, strangeDenominators];
writeLine[""];

writeLine["4. HARD NUMERATOR TRACES"];
writeExpr["N_A,heavy(mu,nu) =", nAHeavy];
writeExpr["N_A,strange(mu,nu) =", nAStrange];
writeExpr["eQ N_A,heavy + eS N_A,strange =", nATotal];
writeExpr["N_B,heavy(mu,nu) =", nBHeavy];
writeExpr["N_B,strange(mu,nu) =", nBStrange];
writeExpr["eQ N_B,heavy + eS N_B,strange =", nBTotal];

writeLine["5. E1-PROJECTED TRIANGLE CORES"];
writeExpr["Projector normalization =", projectorNormalization];
writeExpr["Axial heavy core =", coreAHeavy];
writeExpr["Axial strange core =", coreAStrange];
writeExpr["Tensor heavy core =", coreBHeavy];
writeExpr["Tensor strange core =", coreBStrange];

writeLine["6. SOFT TWO-PARTICLE FIERZ TRACES"];
Do[
  writeExpr["A current, " <> softTracesA[[index, 1]] <> ":", softTracesA[[index, 2]]],
  {index, Length[softTracesA]}
];
Do[
  writeExpr["B current, " <> softTracesB[[index, 1]] <> ":", softTracesB[[index, 2]]],
  {index, Length[softTracesB]}
];
writeExpr["A current, tensor photon-DA E1 coefficient =", softATensorE1];
writeExpr["B current, tensor photon-DA E1 coefficient =", softBTensorE1];
writeExpr["A current, vector photon-DA E1 coefficient =", softAVectorE1];
writeExpr["B current, vector photon-DA E1 coefficient =", softBVectorE1];

writeLine["7. ASSERTIONS"];
Do[
  writeLine[
    If[check[[3]], "PASS  ", "FAIL  "] <> check[[1]] <>
      If[check[[3]], "", "; residual = " <> ToString[check[[2]], InputForm]]
  ],
  {check, checks}
];
writeLine[""];
writeLine["OVERALL STATUS: " <> If[allChecksPassed, "PASS", "FAIL"]];
writeLine["The Bs result uses the same symbolic traces with mQ->mb and eQ->-1/3."];

Close[stream];

Print["Wrote audit to: ", outFile];
Print["Checks passed: ", Count[checks, {_, _, True}], "/", Length[checks]];
Print["Overall status: ", If[allChecksPassed, "PASS", "FAIL"]];

If[! allChecksPassed, Exit[1]];
Quit[];
