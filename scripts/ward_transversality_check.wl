(* ::Package:: *)

(* Ward/transversality check for the E1 invariant tensor.

   Run with:
     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noprompt \
       -script scripts/ward_transversality_check.wl

   The selected E1 tensor is
     S_{mu nu} = p_nu q_mu - (p.q) g_{mu nu}.
   It is the negative of the tensor multiplying the usual scalar structure
     (eta.eps) (p.q) - (eta.q) (p.eps).
*)

<< FeynCalc`

ClearAll[mu, nu, p, q, P, eta, eps, e1Tensor, scalarStructure, rules];

e1Tensor =
  FV[p, nu] FV[q, mu] - SP[p, q] MT[mu, nu];

scalarStructure =
  SP[eta, eps] SP[p, q] - SP[eta, q] SP[p, eps];

rules = {
  SP[q, q] -> 0,
  SP[eps, q] -> 0,
  SP[q, eps] -> 0,
  SP[eta, P] -> 0,
  SP[P, eta] -> 0,
  Pair[Momentum[q], Momentum[q]] -> 0,
  Pair[Momentum[eps], Momentum[q]] -> 0,
  Pair[Momentum[q], Momentum[eps]] -> 0,
  Pair[Momentum[eta], Momentum[P]] -> 0,
  Pair[Momentum[P], Momentum[eta]] -> 0,
  FV[P, mu_] :> FV[p, mu] + FV[q, mu],
  SP[P, q] -> SP[p, q],
  SP[q, P] -> SP[p, q],
  SP[P, eps] -> SP[p, eps],
  SP[eps, P] -> SP[p, eps],
  Pair[Momentum[P], Momentum[q]] -> Pair[Momentum[p], Momentum[q]],
  Pair[Momentum[q], Momentum[P]] -> Pair[Momentum[p], Momentum[q]],
  Pair[Momentum[P], Momentum[eps]] -> Pair[Momentum[p], Momentum[eps]],
  Pair[Momentum[eps], Momentum[P]] -> Pair[Momentum[eps], Momentum[p]]
};

photonWard =
  ScalarProductExpand[
    Contract[e1Tensor FV[q, nu]] //. rules
  ] // Contract // Simplify;

scalarWard =
  ScalarProductExpand[
    scalarStructure /. eps -> q //. rules
  ] // Contract // Simplify;

initialPhysical =
  ScalarProductExpand[
    Contract[e1Tensor FV[P, mu] FV[eps, nu]] //. rules
  ] //. rules // Contract // Simplify;

projectorNormalization =
  ScalarProductExpand[
    Contract[e1Tensor e1Tensor/(2 SP[p, q]^2)] //. rules
  ] //. rules // Contract // Simplify;

outFile = FileNameJoin[{Directory[], "outputs", "ward_transversality_check.txt"}];
stream = OpenWrite[outFile];
WriteString[stream, "Ward/transversality check for the E1 tensor\n"];
WriteString[stream, "==========================================\n\n"];
WriteString[stream, "E1 tensor S_mu_nu = p_nu q_mu - (p.q) g_mu_nu\n"];
WriteString[stream, "Assumptions: q^2=0, eps.q=0, P=p+q, eta.P=0.\n\n"];
WriteString[stream, "Photon-index Ward contraction S_mu_nu q^nu = ",
  ToString[photonWard, InputForm], "\n"];
WriteString[stream, "Scalar structure after eps -> q = ",
  ToString[scalarWard, InputForm], "\n"];
WriteString[stream, "Physical initial-polarization contraction P^mu S_mu_nu eps^nu = ",
  ToString[initialPhysical, InputForm], "\n"];
WriteString[stream, "Projector normalization S_mu_nu S^(mu nu)/(2 (p.q)^2) = ",
  ToString[projectorNormalization, InputForm], "\n\n"];
WriteString[stream, "PASS iff the first three entries are 0 and the normalization is 1.\n"];
Close[stream];

Print["Photon Ward contraction = ", photonWard];
Print["Scalar eps->q check = ", scalarWard];
Print["Physical initial contraction = ", initialPhysical];
Print["Projector normalization = ", projectorNormalization];
Print["Wrote Ward/transversality check to: ", outFile];

Quit[];
