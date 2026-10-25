(* Parse proposed full sources without evaluating their main routines, then
   evaluate only localPieces and compare against the independent derivation. *)
ClearAll["Global`*"];
auditDirectory = DirectoryName[$InputFileName];
sources = Import[FileNameJoin[{auditDirectory,"proposed_wolfram_sources.json"}],"RawJSON"];
derived = Get[FileNameJoin[{auditDirectory,"dimension5_symbolic.wdx"}]];
expected = mixedSS {{derived["AA"]["full_d5_borel"], derived["AB"]["full_d5_borel"]},
  {derived["BA"]["full_d5_borel"], derived["BB"]["full_d5_borel"]}} /. K->M2;
checks = Association@KeyValueMap[Function[{name,source},Module[{parsed,start,finish,snippet,difference},
  parsed = Check[ToExpression[source,InputForm,HoldComplete],$Failed];
  If[parsed===$Failed,Print["Syntax failure: ",name];Exit[1]];
  start = StringPosition[source,"localPieces[M2_] :="][[1,1]];
  finish = StringPosition[source,"localMatrix[M2_] :="][[1,1]]-1;
  snippet = StringTake[source,{start,finish}];
  Clear[localPieces];
  ToExpression[snippet,InputForm];
  difference = FullSimplify[(localPieces[M2]["d5_mixed"] /.
    {mc->m,mb->m,ms->d-m})-expected,Assumptions->m>0&&d>0&&M2>0];
  name -> <|"syntax_pass"->True,"d5_matches_independent_derivation"->
    TrueQ[difference==ConstantArray[0,{2,2}]]|>
]],sources];
Export[FileNameJoin[{auditDirectory,"proposed_wolfram_patch_check.json"}],checks,"RawJSON"];
Print[InputForm[checks]];
If[!And@@Flatten[Values/@Values[checks]],Exit[1]];
