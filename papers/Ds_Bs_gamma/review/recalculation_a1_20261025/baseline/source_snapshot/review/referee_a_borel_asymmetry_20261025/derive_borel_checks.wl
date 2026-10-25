(* Isolated scalar checks. No production or manuscript outputs. *)
ClearAll["Global`*"];
here = DirectoryName[$InputFileName];
$Assumptions = 0 < z < 1 && K > 0 && r > 0;
int[expr_, lim_] := Integrate[expr, lim, GenerateConditions -> False];

(* The production implementation's T4gamma orientation is aq-aqb. *)
t4[x_, y_] := 60 (1-x-y)^2 (x-y) (4-7(x+y));
l1 = int[int[60 x (z-y) (x-y) (4-7(x+y)), {y,0,z}], {x,0,1-z}];
l2 = int[int[z/beta^2 int[t4[x,y], {y,0,beta}], {beta,z,1-x}], {x,0,1-z}];
primitive = int[int[t4[xp,y], {y,0,1-xp}], {xp,0,x}];
l3 = int[z primitive/(1-x)^2, {x,0,1-z}];
line = FullSimplify[l1-l2+l3];
prime = FullSimplify[D[line,z]];
Print["L(z) = ", InputForm[line]];
Print["Lprime(z) = ", InputForm[prime]];
Print["Lprime(1/2) = ", N[prime/.z->1/2,30]];

(* Schwinger constraints for D=m^2-a P^2-(1-a)p^2. *)
mAx2 = (1+r) K;
mPs2 = (1+r) K/r;
saddle = Solve[{tau a==1/mAx2, tau (1-a)==1/mPs2},{tau,a}];
checks = <|
  "effective_borel" -> FullSimplify[mAx2 mPs2/(mAx2+mPs2)-K],
  "saddle_fraction" -> FullSimplify[(a/.First[saddle])-1/(1+r)],
  "saddle_tau" -> FullSimplify[(tau/.First[saddle])-1/K],
  "pole_exponents_mass_ratio" -> FullSimplify[(massA^2/mAx2-massP^2/mPs2)/.r->massA^2/massP^2]
|>;
Export[FileNameJoin[{here,"symbolic_checks.json"}], <|
  "checks" -> Map[ToString[#,InputForm]&,checks],
  "L_T4gamma" -> ToString[line,InputForm],
  "Lprime_T4gamma" -> ToString[prime,InputForm],
  "orientation" -> "production aq-aqb; not manuscript aqb-aq",
  "saddle" -> ToString[saddle,InputForm]
|>, "RawJSON"];
Print["Scalar checks: ",checks];
If[!And@@(TrueQ[#==0]&/@Values[checks]),Exit[1]];
(* RawJSON export above is sufficient. Avoid CSV Export's external importer
   initialization; the independent Python quadrature check writes the samples. *)
Exit[0];
