(* ::Package:: *)

(* Evaluate every input cell in the generated Ds/Bs gamma notebook in the
   same top-to-bottom order used by an interactive Mathematica session.

   Run from papers/Ds_Bs_gamma with

     /Applications/Wolfram.app/Contents/MacOS/WolframKernel -noinit \
       -noprompt -script scripts/check_dsbs_symbolic_notebook.wl
*)

ClearAll[
  scriptDir, paperDir, notebookFile, notebookExpression, inputBoxes,
  cellResults, failedCells, finalStatus, centralReference,
  centralDifference, centralStatus
];

scriptDir = DirectoryName[$InputFileName];
paperDir = DirectoryName[scriptDir];
notebookFile = FileNameJoin[
  {paperDir, "notebooks", "DsBs_gamma_symbolic_derivation.nb"}];

notebookExpression = Get[notebookFile];
inputBoxes = Cases[
  notebookExpression,
  Cell[BoxData[boxes_], "Input", ___] :> boxes,
  Infinity
];

(* Parse each box only when its turn is reached.  Parsing all cells before the
   first Needs["FeynCalc`"] would incorrectly create Global` shadows of the
   FeynCalc symbols used in later cells. *)
cellResults = Table[
  Quiet[ToExpression[inputBoxes[[cell]], StandardForm]],
  {cell, Length[inputBoxes]}
];

failedCells = Flatten[Position[cellResults, $Failed]];
finalStatus = If[
  AssociationQ[Last[cellResults]],
  Lookup[Last[cellResults], "overall status", "MISSING"],
  "MISSING"
];
centralReference = {0.40058930523172204, 0.1666703674111434};
centralDifference =
  If[ValueQ[centralF1F2FromExplicit],
    N[centralF1F2FromExplicit - centralReference],
    {Indeterminate, Indeterminate}
  ];
centralStatus =
  VectorQ[centralDifference, NumericQ] &&
  Max[Abs[centralDifference]] < 10^-7;

Print["INPUT_COUNT = ", Length[inputBoxes]];
Print["FAILED_CELLS = ", failedCells];
Print["FINAL_STATUS = ", finalStatus];
Print["THREE_PARTICLE_RESIDUALS = ", threeParticleResiduals];
Print["CENTRAL_F1_F2 = ", centralF1F2FromExplicit];
Print["CENTRAL_F1_F2_MINUS_PYTHON = ", centralDifference];

If[
  failedCells =!= {} || finalStatus =!= "PASS" || ! centralStatus,
  Quit[1]
];
Quit[0];
