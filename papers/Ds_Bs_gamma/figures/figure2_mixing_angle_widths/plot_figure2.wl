(* ::Package:: *)

(* Figure 2: mixing-angle dependence of radiative widths.
   Run from the repository root or from this directory:

   /Applications/Wolfram.app/Contents/MacOS/WolframKernel \
     -script papers/Ds_Bs_gamma/figures/figure2_mixing_angle_widths/plot_figure2.wl
*)

ClearAll["Global`*"];

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
SetDirectory[scriptDir];

dataFile = FileNameJoin[{scriptDir, "figure2_mixing_angle_widths.csv"}];
outPdf = FileNameJoin[{scriptDir, "figure2_mixing_angle_widths_mathematica.pdf"}];
outPng = FileNameJoin[{scriptDir, "figure2_mixing_angle_widths_mathematica.png"}];

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

blankQ[x_] := StringQ[x] && StringLength[x] == 0;
num[x_?NumericQ] := x;
num[x_String] := If[blankQ[x], Missing["Empty"], ToExpression[x]];
num[x_] := x;

styles = <|
  "low" -> <|"color" -> RGBColor[0.000, 0.278, 0.671], "dash" -> {}|>,
  "high" -> <|"color" -> GrayLevel[0.05], "dash" -> {0.040, 0.018}|>
|>;

curveData[state_] := SortBy[
  ({num[#["theta_deg"]], num[#["Gamma_keV"]]} &) /@
    Select[rows, #["state"] == state && ! blankQ[#["theta_deg"]] &],
  First
];

comboFor[state_] := First[Select[rows, #["state"] == state &]]["combo"];

compactLegend[entries_, spec_] := Module[
  {xr, yr, x0, y0, lineLen, step},
  xr = 45.0 - 25.0;
  yr = spec["yrange"][[2]] - spec["yrange"][[1]];
  x0 = 25.0 + spec["legendAnchor"][[1]] xr;
  y0 = spec["yrange"][[1]] + spec["legendAnchor"][[2]] yr;
  lineLen = 0.080 xr;
  step = 0.085 yr;
  MapIndexed[
    With[{style = styles[comboFor[#[[1]]]], label = #[[2]], k = First[#2] - 1},
      {
        {Directive[style["color"], AbsoluteThickness[2.1], Dashing[style["dash"]]],
          Line[{{x0, y0 - k step}, {x0 + lineLen, y0 - k step}}]},
        Text[
          Style[tex[label, 8.0], GrayLevel[0.10]],
          {x0 + 1.22 lineLen, y0 - k step},
          {-1, 0}
        ]
      }
    ] &,
    entries
  ]
];

plotPanel[spec_] := Module[
  {states, data, plotStyles},
  states = spec["states"][[All, 1]];
  data = curveData /@ states;
  plotStyles = (Directive[
      styles[comboFor[#]]["color"],
      AbsoluteThickness[2.1],
      Dashing[styles[comboFor[#]]["dash"]]
    ] &) /@ states;
  ListLinePlot[
    data,
    PlotStyle -> plotStyles,
    Frame -> True,
    Axes -> False,
    PlotRange -> {{25.0, 45.0}, spec["yrange"]},
    PlotRangePadding -> {Scaled[0.015], Scaled[0.03]},
    ImagePadding -> {{58, 14}, {47, 16}},
    Background -> White,
    GridLines -> None,
    FrameStyle -> Directive[GrayLevel[0.15], AbsoluteThickness[0.8]],
    FrameTicksStyle -> Directive[GrayLevel[0.15], 10],
    LabelStyle -> Directive[GrayLevel[0.08], 10],
    FrameLabel -> {
      tex["\\theta\\,[^\\circ]", 12],
      tex["\\Gamma\\,[\\mathrm{keV}]", 12]
    },
    PlotLabel -> tex[spec["title"], 13],
    Epilog -> Join[
      {
        {Directive[GrayLevel[0.45], AbsoluteThickness[1.25], Dashing[{0.012, 0.018}]],
          Line[{{35.3, spec["yrange"][[1]]}, {35.3, spec["yrange"][[2]]}}]},
        Text[
          Style[tex["35.3^\\circ", 8.0], GrayLevel[0.30]],
          spec["thetaLabel"],
          {0, 0}
        ]
      },
      compactLegend[spec["states"], spec]
    ],
    ImageSize -> 342
  ]
];

specs = {
  <|"states" -> {{"Ds1_2460", "D_{s1}\\,(2460)"}, {"Ds1_2536", "D_{s1}\\,(2536)"}},
    "yrange" -> {0.0, 52.0}, "legendAnchor" -> {0.08, 0.70},
    "thetaLabel" -> {34.45, 26.0},
    "title" -> "D_{s1}\\to D_s\\gamma"|>,
  <|"states" -> {{"Bs1_5750", "B_{s1}\\,(5750)"}, {"Bs1_5830", "B_{s1}\\,(5830)"}},
    "yrange" -> {0.0, 30.0}, "legendAnchor" -> {0.08, 0.72},
    "thetaLabel" -> {34.45, 15.0},
    "title" -> "B_{s1}\\to B_s\\gamma"|>
};

Print["Building panels"];
panels = plotPanel /@ specs;
fig = GraphicsGrid[{panels}, Spacings -> {0.22, 0.0}, ImageSize -> 720, Background -> White];

Print["Exporting PDF"];
Export[outPdf, fig];
Print["Exporting PNG"];
Export[outPng, fig, ImageResolution -> 240];
Print["Wrote ", outPdf];
Print["Wrote ", outPng];
Quit[];
