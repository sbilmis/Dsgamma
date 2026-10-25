(* Isolated analytic distribution checks. No transition or width calculation. *)
ClearAll[s, t, v, w, x, y, q, j, n, rhoG, anti, g, boundary, printed, direct];

(* b=1, sigma=2; q=S/(S+T), and t=d(s) defines a linear domain.
   The fourth-order zero at b removes lower-endpoint contributions through
   all the derivatives tested here. x,y are independent inverse Borel masses. *)
rhoG[s_, t_] := (s-1)^4 Exp[-x s-y t];
boundary[s_] := (2-(1-q)s)/q;
anti = Integrate[(s-1)^4 Exp[-(x+y)s], s,
  Assumptions -> x>0 && y>0];
g[v_] := Exp[-y v] ((anti /. s -> 2-q v)-(anti /. s -> 1));
printed[n_] := (
  (-1)^(n-1) Integrate[
    D[rhoG[s,t], {t,n-1}] /. t -> s, {s,1,2},
    Assumptions -> x>0 && y>0]
  + Sum[(-1)^(j-1) q^(n-j)
    (D[(D[rhoG[s,t], {t,j-1}] /. t -> boundary[s]),
      {s,n-j-1}] /. s -> 2), {j,1,n-1}])/Gamma[n];
direct[n_] := (-1)^(n-1) (D[g[v], {v,n-1}] /. v -> 0)/Gamma[n];
linearChecks = Table[
  FullSimplify[direct[n]-printed[n], Assumptions -> x>0 && y>0 && 0<q<1],
  {n,1,4}];
Print["Linear domain C.4, n=1,2,3,4 residuals: ", linearChecks];

(* A curved-domain control: n=3, rho=(s-1)^2, g=1, sigma=2.
   With v=t-s, s^2+t^2<=8 gives upper s=(-v+Sqrt[16-v^2])/2.
   The lower-bound contribution is O(v^3) and does not affect G''(0).
   The printed C.4 expression misses the curvature part in this test. *)
upper[v_] := (-v+Sqrt[16-v^2])/2;
curvedG[v_] := (upper[v]-1)^3/3;
curvedDirect = FullSimplify[(D[curvedG[v],{v,2}] /. v->0)/Gamma[3]];
curvedPrinted = (1/2)^2*2/Gamma[3];
curvedCorrection = FullSimplify[
  ((upper[0]-1)^2*(D[upper[v],{v,2}]/.v->0))/Gamma[3]];
Print["Curved control {direct, printed C.4, missing curvature}: ",
  {curvedDirect,curvedPrinted,curvedCorrection}];

Export[FileNameJoin[{Directory[],"review",
  "referee_a_continuum_planning_20261025","pz_boundary_checks.json"}],
  <|"linear_n1_through_n4_residuals" -> linearChecks,
    "curved_n3_direct_exact" -> ToString[curvedDirect,InputForm],
    "curved_n3_printed_exact" -> ToString[curvedPrinted,InputForm],
    "curved_n3_missing_term_exact" -> ToString[curvedCorrection,InputForm],
    "scope" -> "Analytic test kernels only; no production or physics result changed"|>,
  "RawJSON"];
Quit[];
