(* ::Package:: *)

(* Figure 4: selected literature comparison.
   Run from the repository root or from this directory:

   /Applications/Wolfram.app/Contents/MacOS/WolframKernel \
     -script papers/Ds_Bs_gamma/figures/figure4_literature_comparison/plot_figure4.wl
*)

ClearAll["Global`*"];

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
SetDirectory[scriptDir];

dataFile = FileNameJoin[{scriptDir, "figure4_literature_comparison.csv"}];
outPdf = FileNameJoin[{scriptDir, "figure4_literature_comparison_mathematica.pdf"}];
outPng = FileNameJoin[{scriptDir, "figure4_literature_comparison_mathematica.png"}];

If[! FileExistsQ[dataFile],
  Print["Missing input CSV: ", dataFile];
  Exit[1];
];

matexPath = FileNameJoin[{
  $UserBaseDirectory, "Paclets", "Repository", "MaTeX-1.7.10", "MaTeX.m"
}];
matexAvailable = TimeConstrained[
  Check[
    If[FileExistsQ[matexPath],
      Get[matexPath],
      Needs["MaTeX`"]
    ];
    MaTeX`ConfigureMaTeX[
      "pdfLaTeX" -> "/Library/TeX/texbin/pdflatex",
      "Ghostscript" -> "/opt/homebrew/bin/gs",
      "WorkingDirectory" -> scriptDir,
      "CacheSize" -> 100
    ];
    True,
    False
  ],
  45,
  False
];
Print["MaTeX available: ", matexAvailable];
tex[s_String, size_: 13] := If[
  matexAvailable,
  MaTeX`MaTeX[s, Magnification -> size/12],
  Style[s, FontFamily -> "Times", size]
];

raw = Import[dataFile, "CSV"];
headers = First[raw];
rows = AssociationThread[headers, #] & /@ Rest[raw];
Print["Rows imported: ", Length[rows]];

num[x_?NumericQ] := x;
num[x_String] := ToExpression[x];

log10[x_] := Log[10, x];
xVals = Range[Length[rows]];
this = num /@ rows[[All, "this_work_median_keV"]];
p16 = num /@ rows[[All, "this_work_p16_keV"]];
p84 = num /@ rows[[All, "this_work_p84_keV"]];
bondar = num /@ rows[[All, "bondar_milstein_2025_keV"]];
labels = {
  "D_{s1}(2460)",
  "D_{s1}(2536)",
  "B_{s1}(5750)",
  "B_{s1}(5830)"
};

thisColor = RGBColor[0.000, 0.278, 0.671];
bondarColor = GrayLevel[0.05];
yMin = log10[0.02];
yMax = log10[500.0];

yticks = {
  {log10[0.02], tex["0.02", 8]},
  {log10[0.05], tex["0.05", 8]},
  {log10[0.1], tex["0.1", 8]},
  {log10[0.5], tex["0.5", 8]},
  {log10[1], tex["1", 8]},
  {log10[5], tex["5", 8]},
  {log10[10], tex["10", 8]},
  {log10[50], tex["50", 8]},
  {log10[100], tex["100", 8]},
  {log10[500], tex["500", 8]}
};

legend = {
  {Directive[thisColor, AbsoluteThickness[1.8]], Line[{{3.10, yMax - 0.30}, {3.38, yMax - 0.30}}],
    PointSize[0.012], Point[{3.24, yMax - 0.30}],
    Text[Style[tex["\\mathrm{This\\ work}", 8.5], GrayLevel[0.10]], {3.46, yMax - 0.30}, {-1, 0}]},
  {Directive[bondarColor, AbsoluteThickness[1.8]], PointSize[0.014], Point[{3.24, yMax - 0.55}],
    Text[Style[tex["\\mathrm{Bondar\\!-\!Milstein}", 8.5], GrayLevel[0.10]], {3.46, yMax - 0.55}, {-1, 0}]}
};

primitives = Flatten[
  Table[
    {
      {Directive[thisColor, AbsoluteThickness[1.5]],
        Line[{{xVals[[i]] - 0.08, log10[p16[[i]]]}, {xVals[[i]] - 0.08, log10[p84[[i]]]}}],
        Line[{{xVals[[i]] - 0.15, log10[p16[[i]]]}, {xVals[[i]] - 0.01, log10[p16[[i]]]}}],
        Line[{{xVals[[i]] - 0.15, log10[p84[[i]]]}, {xVals[[i]] - 0.01, log10[p84[[i]]]}}]},
      {Directive[thisColor, PointSize[0.013]], Point[{xVals[[i]] - 0.08, log10[this[[i]]]}]},
      {Directive[bondarColor, PointSize[0.015]], Point[{xVals[[i]] + 0.08, log10[bondar[[i]]]}]}
    },
    {i, Length[xVals]}
  ],
  1
];

fig = Graphics[
  {
    primitives,
    legend
  },
  PlotRange -> {{0.55, 4.45}, {yMin, yMax}},
  Frame -> True,
  Axes -> False,
  Background -> White,
  ImagePadding -> {{60, 16}, {58, 16}},
  FrameStyle -> Directive[GrayLevel[0.15], AbsoluteThickness[0.8]],
  FrameTicksStyle -> Directive[GrayLevel[0.15], 10],
  LabelStyle -> Directive[GrayLevel[0.08], 10],
  FrameTicks -> {
    {yticks, None},
    {Thread[{xVals, tex[#, 9.2] & /@ labels}], None}
  },
  FrameLabel -> {None, tex["\\Gamma\\,[\\mathrm{keV}]", 12]},
  PlotLabel -> tex["\\mathrm{Comparison\\ with\\ selected\\ literature\\ estimates}", 13],
  AspectRatio -> 0.48,
  ImageSize -> 720
];

Print["Exporting PDF"];
Export[outPdf, fig];
Print["Exporting PNG"];
Export[outPng, fig, ImageResolution -> 240];
Print["Wrote ", outPdf];
Print["Wrote ", outPng];
Quit[];
