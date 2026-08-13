(* ::Package:: *)

(* Figure 5: Monte Carlo width distributions.
   Run from the repository root or from this directory:

   /Applications/Wolfram.app/Contents/MacOS/WolframKernel \
     -script papers/Ds_Bs_gamma/figures/figure5_monte_carlo_widths/plot_figure5.wl
*)

ClearAll["Global`*"];

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
SetDirectory[scriptDir];

dataFile = FileNameJoin[{scriptDir, "figure5_monte_carlo_widths.csv"}];
outPdf = FileNameJoin[{scriptDir, "figure5_monte_carlo_widths_mathematica.pdf"}];
outPng = FileNameJoin[{scriptDir, "figure5_monte_carlo_widths_mathematica.png"}];

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

valuesFor[state_] := num /@ Select[rows, #["state_key"] == state &][[All, "Gamma_keV"]];

summaryText[values_] := Module[{med, p16, p84},
  med = Quantile[values, 0.50];
  p16 = Quantile[values, 0.16];
  p84 = Quantile[values, 0.84];
  ToString[NumberForm[med, {4, 2}]] <> "_{-" <>
    ToString[NumberForm[med - p16, {3, 2}]] <> "}^{+" <>
    ToString[NumberForm[p84 - med, {3, 2}]] <> "}\\,\\mathrm{keV}"
];

panel[state_, title_, xrange_] := Module[
  {vals, med, p16, p84, hist},
  vals = valuesFor[state];
  med = Quantile[vals, 0.50];
  p16 = Quantile[vals, 0.16];
  p84 = Quantile[vals, 0.84];
  hist = Histogram[
    vals,
    {xrange[[1]], xrange[[2]], (xrange[[2]] - xrange[[1]])/28.0},
    "PDF",
    ChartStyle -> Directive[RGBColor[0.000, 0.278, 0.671], Opacity[0.55]],
    ChartBaseStyle -> EdgeForm[White],
    Frame -> True,
    Axes -> False,
    PlotRange -> {xrange, All},
    PlotRangePadding -> {Scaled[0.015], Scaled[0.06]},
    GridLines -> None,
    Background -> White,
    ImagePadding -> {{56, 14}, {45, 16}},
    FrameStyle -> Directive[GrayLevel[0.15], AbsoluteThickness[0.8]],
    FrameTicksStyle -> Directive[GrayLevel[0.15], 9],
    LabelStyle -> Directive[GrayLevel[0.08], 9],
    FrameLabel -> {tex["\\Gamma\\,[\\mathrm{keV}]", 11], tex["\\mathrm{density}", 11]},
    PlotLabel -> tex[title, 12],
    ImageSize -> 342
  ];
  Show[
    hist,
    Epilog -> {
      {Directive[GrayLevel[0.05], AbsoluteThickness[1.35]], InfiniteLine[{{med, 0}, {med, 1}}]},
      {Directive[RGBColor[0.800, 0.200, 0.000], AbsoluteThickness[1.1], Dashing[{0.035, 0.018}]],
        InfiniteLine[{{p16, 0}, {p16, 1}}],
        InfiniteLine[{{p84, 0}, {p84, 1}}]},
      Text[Style[tex[summaryText[vals], 7.6], GrayLevel[0.10]], Scaled[{0.72, 0.88}]]
    }
  ]
];

panels = {
  panel["Ds1_2460", "D_{s1}(2460)\\to D_s\\gamma", {0.0, 78.0}],
  panel["Ds1_2536", "D_{s1}(2536)\\to D_s\\gamma", {0.0, 84.0}],
  panel["Bs1_5750", "B_{s1}(5750)\\to B_s\\gamma", {0.0, 27.0}],
  panel["Bs1_5830", "B_{s1}(5830)\\to B_s\\gamma", {0.0, 43.0}]
};

fig = GraphicsGrid[Partition[panels, 2], Spacings -> {0.18, 0.28}, ImageSize -> 720, Background -> White];

Print["Exporting PDF"];
Export[outPdf, fig];
Print["Exporting PNG"];
Export[outPng, fig, ImageResolution -> 240];
Print["Wrote ", outPdf];
Print["Wrote ", outPng];
Quit[];
