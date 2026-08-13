(* ::Package:: *)

(* Figure 1: Borel-parameter dependence of transition form factors.
   Run from the repository root or from this directory:

   /Applications/Wolfram.app/Contents/MacOS/WolframKernel \
     -script papers/Ds_Bs_gamma/figures/figure1_borel_form_factors/plot_figure1.wl
*)

ClearAll["Global`*"];

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
SetDirectory[scriptDir];

dataFile = FileNameJoin[{scriptDir, "figure1_borel_form_factors.csv"}];
outPdf = FileNameJoin[{scriptDir, "figure1_borel_form_factors_mathematica.pdf"}];
outPng = FileNameJoin[{scriptDir, "figure1_borel_form_factors_mathematica.png"}];

If[! FileExistsQ[dataFile],
  Print["Missing input CSV: ", dataFile];
  Exit[1];
];

Print["Reading data from ", dataFile];

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

blankQ[x_] := StringQ[x] && StringLength[x] == 0;
num[x_?NumericQ] := x;
num[x_String] := If[blankQ[x], Missing["Empty"], ToExpression[x]];
num[x_] := x;
nearQ[x_, y_] := Abs[num[x] - y] < 10^-8;

styles = {
  <|"color" -> RGBColor[0.000, 0.278, 0.671], "dash" -> {}, "label" -> "solid"|>,
  <|"color" -> RGBColor[0.800, 0.200, 0.000], "dash" -> {0.035, 0.018}, "label" -> "dashed"|>,
  <|"color" -> GrayLevel[0.05], "dash" -> {0.045, 0.016, 0.010, 0.016}, "label" -> "dot-dashed"|>
};

channelRows[state_, s0_] := SortBy[
  Select[
    rows,
    #["state"] == state && #["scan"] == "M2" && nearQ[#["fixed_s0"], s0] &&
      ! blankQ[#["M2"]] && ! blankQ[#["s0"]] && ! blankQ[#["g_abs"]] &
  ],
  num[#["M2"]] &
];

curveData[state_, s0_] := ({num[#["M2"]], num[#["g_abs"]]} &) /@
  channelRows[state, s0];

yRange[data_, minSpan_] := Module[{ys, ymin, ymax, center, half},
  ys = data[[All, All, 2]] // Flatten;
  ymin = Min[ys]; ymax = Max[ys];
  center = (ymin + ymax)/2;
  half = Max[(ymax - ymin)/2, minSpan/2, 0.10 Abs[center], 0.002];
  {center - 1.45 half, center + 1.45 half}
];

fmtS0[x_] := If[
  Abs[x - Round[x]] < 10^-8,
  ToString[Round[x]],
  ToString[NumberForm[x, {4, 1}, NumberPadding -> {"", "0"}]]
];

compactLegend[s0Values_, spec_] := Module[
  {xr, yr, x0, y0, lineLen, step, ordered},
  xr = spec["xmax"] - spec["xmin"];
  yr = spec["yrange"][[2]] - spec["yrange"][[1]];
  x0 = spec["xmin"] + spec["legendAnchor"][[1]] xr;
  y0 = spec["yrange"][[1]] + spec["legendAnchor"][[2]] yr;
  lineLen = 0.075 xr;
  step = 0.050 yr;
  ordered = Transpose[{s0Values, styles, Range[Length[s0Values]]}];
  Map[
    With[{s0 = #[[1]], style = #[[2]], k = #[[3]] - 1},
      {
        {Directive[style["color"], AbsoluteThickness[2.1], Dashing[style["dash"]]], Line[{{x0, y0 - k step}, {x0 + lineLen, y0 - k step}}]},
        Text[
          Style[tex["s_0=" <> fmtS0[s0], 7.2], GrayLevel[0.12]],
          {x0 + 1.25 lineLen, y0 - k step},
          {-1, 0}
        ]
      }
    ] &,
    ordered
  ]
];

plotPanel[spec_] := Module[
  {data, yr, plotSpec},
  data = curveData[spec["state"], #] & /@ spec["s0"];
  yr = yRange[data, spec["minSpan"]];
  plotSpec = Join[spec, <|"yrange" -> yr|>];
  ListLinePlot[
    data,
    PlotStyle -> (Directive[#["color"], AbsoluteThickness[2.1], Dashing[#["dash"]]] & /@ styles),
    Frame -> True,
    Axes -> False,
    PlotRange -> {{spec["xmin"], spec["xmax"]}, yr},
    PlotRangePadding -> {Scaled[0.015], Scaled[0.03]},
    ImagePadding -> {{58, 14}, {48, 16}},
    GridLines -> None,
    FrameStyle -> Directive[GrayLevel[0.15], AbsoluteThickness[0.8]],
    FrameTicksStyle -> Directive[GrayLevel[0.15], 10],
    LabelStyle -> Directive[GrayLevel[0.08], 10],
    FrameLabel -> {
      tex["M^2\\,[\\mathrm{GeV}^2]", 12],
      tex["|g|\\,[\\mathrm{GeV}^{-1}]", 12]
    },
    PlotLabel -> tex[spec["title"], 13],
    Epilog -> compactLegend[spec["s0"], plotSpec],
    ImageSize -> 330
  ]
];

specs = {
  <|"state" -> "Ds1_2460",
    "s0" -> {8.5, 9.0, 9.5}, "xmin" -> 3.0, "xmax" -> 4.5,
    "minSpan" -> 0.08, "legendAnchor" -> {0.08, 0.88},
    "title" -> "D_{s1}(2460)\\to D_s\\gamma"|>,
  <|"state" -> "Ds1_2536",
    "s0" -> {9.0, 9.5, 10.0}, "xmin" -> 3.0, "xmax" -> 4.5,
    "minSpan" -> 0.04, "legendAnchor" -> {0.70, 0.86},
    "title" -> "D_{s1}(2536)\\to D_s\\gamma"|>,
  <|"state" -> "Bs1_5750",
    "s0" -> {39.0, 40.0, 41.0}, "xmin" -> 10.0, "xmax" -> 14.0,
    "minSpan" -> 0.10, "legendAnchor" -> {0.68, 0.87},
    "title" -> "B_{s1}(5750)\\to B_s\\gamma"|>,
  <|"state" -> "Bs1_5830",
    "s0" -> {40.0, 41.0, 42.0}, "xmin" -> 10.0, "xmax" -> 14.0,
    "minSpan" -> 0.025, "legendAnchor" -> {0.68, 0.87},
    "title" -> "B_{s1}(5830)\\to B_s\\gamma"|>
};

Print["Building panels"];
panels = plotPanel /@ specs;

fig = GraphicsGrid[
  Partition[panels, 2],
  Spacings -> {0.2, 0.35},
  ImageSize -> 720
];

Print["Exporting PDF"];
Export[outPdf, fig];
Print["Exporting PNG"];
Export[outPng, fig, ImageResolution -> 240];

Print["Wrote ", outPdf];
Print["Wrote ", outPng];
Quit[];
