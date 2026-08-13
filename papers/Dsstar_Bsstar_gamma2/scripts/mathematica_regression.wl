(* Independent Mathematica numerical implementation; no Python result is read. *)
ClearAll["Global`*"];
Print["regression: start"];
scriptDir = DirectoryName[$InputFileName];
projectDir = DirectoryName[scriptDir];
csvDir = FileNameJoin[{projectDir, "outputs", "csv"}];
If[!DirectoryQ[csvDir], CreateDirectory[csvDir, CreateIntermediateDirectories -> True]];

alphaEM = 0.0072973525693;
esCharge = -1/3;
qq = -(0.24)^3;
kappaS = 0.8;
ss = kappaS qq;
m0sq = 0.8;
mixedSS = m0sq ss;
chi = -3.15;
f3 = -0.0039;
omegaV = 3.8;
omegaA = -2.1;
kappaDA = 0.2;
zeta1 = 0.4;
zeta2 = 0.3;
ms0 = 0.093;

(* Native adaptive integration at the declared regression points. *)
gIntegrate[f_, a_?NumericQ, b_?NumericQ] := NIntegrate[
  f[z], {z, a, b}, Method -> "GlobalAdaptive", AccuracyGoal -> 10,
  MaxRecursion -> 12
];

lambdaF[s_, m_, ms_] := (s - (m + ms)^2) (s - (m - ms)^2);
rhoAA[s_, m_, ms_] := Sqrt[lambdaF[s, m, ms]]/(8 Pi^2) *
  (s - (m + ms)^2) ((m - ms)^2 + 2 s)/s^2;
rhoAB[s_, m_, ms_] := Sqrt[lambdaF[s, m, ms]]/(8 Pi^2) *
  3 (m - ms) (s - (m + ms)^2)/((m + ms) s);
rhoBB[s_, m_, ms_] := Sqrt[lambdaF[s, m, ms]]/(8 Pi^2) *
  (s - (m + ms)^2) (s + 2 (m - ms)^2)/((m + ms)^2 s);
rhoV[s_, m_, ms_] := Sqrt[lambdaF[s, m, ms]]/(8 Pi^2) *
  (2 - (m^2 + ms^2 - 6 m ms)/s - (m^2 - ms^2)^2/s^2);

localMatrices[M2_, m_, ms_] := Module[{b = m + ms, e, d3, corr, v, d3ms, d5},
  e = Exp[-m^2/M2];
  d3 = ss e {{m, m^2/b}, {m^2/b, m^3/b^2}};
  corr = ss e (-m^2 ms/(2 M2) + m^3 ms^2/(2 M2^2));
  v = {1, m/b};
  d3ms = corr Outer[Times, v, v];
  d5 = e {
    {-mixedSS m^3/(4 M2^2), -mixedSS/(4 b) (m^4/M2^2 - m^2/M2 - 1)},
    {-mixedSS/(4 b) (m^4/M2^2 - m^2/M2 - 1),
     mixedSS/b^2 (-m^5/(4 M2^2) + m^3/(2 M2))}
  };
  {d3, d3ms, d5}
];

borelBasis[M2_, s0_, m_, ms_] := Module[{sth, p, locals},
  sth = (m + ms)^2;
  p = {
    {gIntegrate[(Exp[-#/M2] rhoAA[#, m, ms])&, sth, s0],
     gIntegrate[(Exp[-#/M2] rhoAB[#, m, ms])&, sth, s0]},
    {gIntegrate[(Exp[-#/M2] rhoAB[#, m, ms])&, sth, s0],
     gIntegrate[(Exp[-#/M2] rhoBB[#, m, ms])&, sth, s0]}
  };
  locals = localMatrices[M2, m, ms];
  p + Total[locals]
];
rot[thetaDeg_] := {{Sin[thetaDeg Degree], Cos[thetaDeg Degree]},
  {Cos[thetaDeg Degree], -Sin[thetaDeg Degree]}};
projected[M2_, s0_, m_, ms_, theta_] := rot[theta].borelBasis[M2, s0, m, ms].Transpose[rot[theta]];

fullProjected[M2_, m_, ms_, theta_] := Module[{sth, hi, p, locals},
  sth = (m + ms)^2; hi = sth + 60 M2;
  p = {
    {gIntegrate[(Exp[-#/M2] rhoAA[#, m, ms])&, sth, hi],
     gIntegrate[(Exp[-#/M2] rhoAB[#, m, ms])&, sth, hi]},
    {gIntegrate[(Exp[-#/M2] rhoAB[#, m, ms])&, sth, hi],
     gIntegrate[(Exp[-#/M2] rhoBB[#, m, ms])&, sth, hi]}
  };
  locals = Total[localMatrices[M2, m, ms]];
  rot[theta].p.Transpose[rot[theta]] + rot[theta].locals.Transpose[rot[theta]]
];
massSR[M2_, s0_, m_, ms_, theta_, state_] := Module[{tau, h, fp, fm, f0, r1},
  tau = 1/M2; h = Max[10^-5, tau 2 10^-4];
  fp = projected[1/(tau + h), s0, m, ms, theta][[state, state]];
  fm = projected[1/(tau - h), s0, m, ms, theta][[state, state]];
  f0 = projected[M2, s0, m, ms, theta][[state, state]];
  r1 = -(fp - fm)/(2 h);
  Sqrt[Max[r1/f0, 0]]
];
fInitial[M2_, s0_, m_, ms_, theta_, state_, mass_] :=
  Exp[mass^2/(2 M2)] Sqrt[projected[M2, s0, m, ms, theta][[state, state]]]/mass;

borelVector[M2_, s0_, m_, ms_] := Module[{sth, pert, locs, local},
  sth = (m + ms)^2;
  pert = gIntegrate[(Exp[-#/M2] rhoV[#, m, ms])&, sth, s0];
  locs = localMatrices[M2, m, ms];
  local = -Total[locs][[1, 1]];
  {pert, local, pert + local}
];
fVector[M2_, s0_, m_, ms_, mass_] :=
  Exp[mass^2/(2 M2)] Sqrt[borelVector[M2, s0, m, ms][[3]]]/mass;

phi[u_] := 6 u (1 - u);
hGamma[u_] := -10 (3 (2 u - 1)^2 - 1)/2;
psiV[u_] := 5 (3 (2 u - 1)^2 - 1) +
  3/64 (15 omegaV - 5 omegaA) (3 - 30 (2 u - 1)^2 + 35 (2 u - 1)^4);
psiA[u_] := (1 - (2 u - 1)^2) (5 (2 u - 1)^2 - 1) 5/2 *
  (1 + 9 omegaV/16 - 3 omegaA/16);
Atw4[u_] := 40 u^2 (1-u)^2 (3 kappaDA + 1) -
  24 zeta2 (u (1-u) (2 + 13 u (1-u)) +
   2 u^3 (10 - 15 u + 6 u^2) Log[u] +
   2 (1-u)^3 (10 - 15 (1-u) + 6 (1-u)^2) Log[1-u]);
HGamma[u_] := Integrate[hGamma[x], {x, 0, u}];
Hbar[u_] := Integrate[(u-x) hGamma[x], {x, 0, u}];
PsiV[u_] := Integrate[psiV[x], {x, 0, u}];
iF2 = 0;
iF3 = 15 (omegaA - 3 omegaV + 4)/32;

rhoT[s_, m_, ms_, eQ_] := Module[{root, aS, aQ},
  root = Sqrt[lambdaF[s, m, ms]];
  aS = s - m^2 + ms^2; aQ = s - ms^2 + m^2;
  3 ms m/(4 Pi^2) (esCharge Log[(aS-root)/(aS+root)] +
    eQ Log[(aQ-root)/(aQ+root)])
];
transitionTerms[M2_, s0_, m_, ms_, eQ_] := Module[
  {u0=1/2, sth, em, es0, diff, hard, local, tw2, tw42, tw32, tw43, tw33, total},
  sth=(m+ms)^2; em=Exp[-m^2/M2]; es0=Exp[-s0/M2]; diff=em-es0;
  hard=gIntegrate[(Exp[-#/M2] rhoT[#,m,ms,eQ])&,sth,s0];
  local=eQ m em ss (1-ms^2/M2 (1-m^2/M2));
  tw2=esCharge m ss diff M2 chi phi[u0];
  tw42=esCharge m ss em (-m^2 Atw4[u0]/(4 M2)-
    HGamma[u0](1-u0)-Hbar[u0](1-2 m^2/M2));
  tw32=esCharge f3 M2 diff ((1-u0) psiA'[u0]/4-psiA[u0]/4-
    PsiV[u0](1+2m^2/M2)+(1-u0)psiV[u0]);
  tw43=m esCharge ss em iF2;
  tw33=-esCharge f3 M2 diff iF3;
  total=hard+tw2+tw42+tw32+tw43+tw33;
  <|"hard"->hard,"tw2"->tw2,"tw4_2p"->tw42,"tw3_2p"->tw32,
    "tw4_3p"->tw43,"tw3_3p"->tw33,"local_heavy_diagnostic"->local,
    "quoted_total_A"->total,"quoted_total_B"->m/(m+ms) total|>
];

(* Reproduce the accepted-grid residue averages used by the publication path. *)
initialAverage[m_, theta_, state_, mass_, window_] := Module[
 {mgrid, sgrid, vals={}, M2, s0, pi0, pc, d5frac, msr, pieces, full},
 mgrid=Subdivide[window[[1]],window[[2]],5];
 sgrid=Subdivide[window[[3]],window[[4]],4];
 Do[
  pi0=projected[M2,s0,m,ms0,theta][[state,state]];
  full=fullProjected[M2,m,ms0,theta][[state,state]];
  pc=pi0/full;
  pieces=localMatrices[M2,m,ms0];
  d5frac=Abs[(rot[theta].pieces[[3]].Transpose[rot[theta]])[[state,state]]]/Abs[pi0];
  msr=massSR[M2,s0,m,ms0,theta,state];
  If[pi0>0 && pc>=0.40 && d5frac<=0.10 && Abs[msr/mass-1]<=0.05,
    AppendTo[vals,fInitial[M2,s0,m,ms0,theta,state,mass]]],
  {M2,mgrid},{s0,sgrid}];
 Mean[Flatten[vals]]
];
vectorAverage[m_, mass_, window_] := Module[{mgrid,sgrid,vals={},M2,s0,bv,sth,fullPert,full,pc,lf},
 mgrid=Subdivide[window[[1]],window[[2]],5];sgrid=Subdivide[window[[3]],window[[4]],4];
 Do[bv=borelVector[M2,s0,m,ms0];sth=(m+ms0)^2;
  fullPert=gIntegrate[(Exp[-#/M2]rhoV[#,m,ms0])&,sth,sth+60M2];
  full=fullPert+bv[[2]];pc=bv[[3]]/full;lf=Abs[bv[[2]]]/Abs[bv[[3]]];
  If[pc>=0.40&&lf<=0.25,AppendTo[vals,fVector[M2,s0,m,ms0,mass]]],
  {M2,mgrid},{s0,sgrid}];Mean[Flatten[vals]]
];

sectors = {
 <|"name"->"D","m"->1.27,"eQ"->2/3,"theta"->38.4,"mf"->2.1122,
   "vwindow"->{2.5,3.5,6.5,7.5},"test"->{3.0,8.0,5.0,6.75},
   "states"->{
    <|"key"->"Ds1_2460","index"->1,"mass"->2.4595,"tp"->{2.,3.,7.5,8.5}|>,
    <|"key"->"Ds1_2536","index"->2,"mass"->2.53512,"tp"->{2.,3.,8.5,10.}|>
   }|>,
 <|"name"->"B","m"->4.18,"eQ"->-1/3,"theta"->35.264,"mf"->5.4154,
   "vwindow"->{7.,9.,36.,38.},"test"->{9.0,42.0,6.0,40.0},
   "states"->{
    <|"key"->"Bs1_lower","index"->1,"mass"->5.75,"tp"->{8.,11.,40.,42.}|>,
    <|"key"->"Bs1_5830","index"->2,"mass"->5.82873,"tp"->{7.,9.,42.,44.}|>
   }|>
};
Print["regression: definitions complete"];

rows={{"sector","quantity","value"}};
Do[
 Print["regression: sector ", sec["name"]];
 name=sec["name"];m=sec["m"];eQ=sec["eQ"];theta=sec["theta"];mf=sec["mf"];
 {m2tp,s0tp,m2tr,s0tr}=sec["test"];
 sTest=(m+ms0)^2+1;
 rows=Join[rows,{{name,"rho_AA",N[rhoAA[sTest,m,ms0],16]},
   {name,"rho_AB",N[rhoAB[sTest,m,ms0],16]},
   {name,"rho_BB",N[rhoBB[sTest,m,ms0],16]}}];
 mat=borelBasis[m2tp,s0tp,m,ms0];
 Print["regression: basis matrix complete"];
 rows=Join[rows,{{name,"Pi_AA",mat[[1,1]]},{name,"Pi_AB",mat[[1,2]]},{name,"Pi_BB",mat[[2,2]]}}];
 terms=transitionTerms[m2tr,s0tr,m,ms0,eQ];
 Print["regression: transition terms complete"];
 Do[AppendTo[rows,{name,"T_"<>key,terms[key]}],
  {key,{"hard","tw2","tw4_2p","tw3_2p","tw4_3p","tw3_3p","quoted_total_A","quoted_total_B"}}];
 fv=fVector[Mean[sec["vwindow"][[1;;2]]],Mean[sec["vwindow"][[3;;4]]],m,ms0,mf];
 Print["regression: vector residue complete"];
 Do[
  fi=fInitial[m2tp,s0tp,m,ms0,theta,st["index"],st["mass"]];
  rr=rot[theta][[st["index"]]];
  tphys=rr[[1]]terms["quoted_total_A"]+rr[[2]]terms["quoted_total_B"];
  pref=Exp[(st["mass"]^2+mf^2)/(2m2tr)]/(fi fv st["mass"] mf);
  coupling=pref tphys;
  kg=(st["mass"]^2-mf^2)/(2st["mass"]);
  width=alphaEM coupling^2 kg^3 (st["mass"]^2+mf^2)/(3st["mass"]^2 mf^2) 10^6;
  rows=Join[rows,{{name,st["key"]<>"_f",fi},{name,st["key"]<>"_g",coupling},
    {name,st["key"]<>"_width_keV",width}}],
  {st,sec["states"]}];
 Print["regression: state residues and observables complete"];
 AppendTo[rows,{name,"f_"<>name<>"sstar",fv}],
 {sec,sectors}];

csvText = StringRiffle[
  (StringRiffle[(ToString[#, InputForm]& /@ #), ","]& /@ rows), "\n"];
Export[FileNameJoin[{csvDir,"mathematica_regression.csv"}],csvText,"Text"];
Print["regression: wrote ", Length[rows]-1, " rows"];
Exit[0];
