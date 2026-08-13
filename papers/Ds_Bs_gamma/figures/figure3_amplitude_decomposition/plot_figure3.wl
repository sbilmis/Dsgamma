(* ::Package:: *)

(* Figure 3: amplitude decomposition.
   Run from the repository root or from this directory:

   /Applications/Wolfram.app/Contents/MacOS/WolframKernel \
     -script papers/Ds_Bs_gamma/figures/figure3_amplitude_decomposition/plot_figure3.wl
*)

ClearAll["Global`*"];

scriptDir = DirectoryName[$InputFileName];
If[scriptDir === "", scriptDir = NotebookDirectory[]];
SetDirectory[scriptDir];

dataFile = FileNameJoin[{scriptDir, "figure3_amplitude_decomposition.csv"}];
outPdf = FileNameJoin[{scriptDir, "figure3_amplitude_decomposition_mathematica.pdf"}];
outPng = FileNameJoin[{scriptDir, "figure3_amplitude_decomposition_mathematica.png"}];

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

rowByState[state_] := First[Select[rows, #["state"] == state &]];

components = {
  <|"key" -> "A_component_GeV_inv", "label" -> "A\\;\\mathrm{piece}",
    "color" -> RGBColor[0.000, 0.278, 0.671], "dash" -> {}, "dy" -> 0.22|>,
  <|"key" -> "B_component_GeV_inv", "label" -> "B\\;\\mathrm{piece}",
    "color" -> RGBColor[0.800, 0.200, 0.000], "dash" -> {0.035, 0.018}, "dy" -> 0.0|>,
  <|"key" -> "G_total_GeV_inv", "label" -> "G=A+B",
    "color" -> GrayLevel[0.05], "dash" -> {0.045, 0.016, 0.010, 0.016}, "dy" -> -0.22|>
};

stateLabels = <|
  "Ds1_2460" -> "D_{s1}(2460)",
  "Ds1_2536" -> "D_{s1}(2536)",
  "Bs1_5750" -> "B_{s1}(5750)",
  "Bs1_5830" -> "B_{s1}(5830)"
|>;

legend[] := Module[{x0 = 0.22, y0 = 1.45, lineLen = 0.070, step = 0.17},
  MapIndexed[
    With[{comp = #, k = First[#2] - 1},
      {
        {Directive[comp["color"], AbsoluteThickness[2.2], Dashing[comp["dash"]]],
          Line[{{x0, y0 - k step}, {x0 + lineLen, y0 - k step}}]},
        Text[
          Style[tex[comp["label"], 8.0], GrayLevel[0.10]],
          {x0 + 1.25 lineLen, y0 - k step},
          {-1, 0}
        ]
      }
    ] &,
    components
  ]
];

componentGraphics[states_, title_] := Module[
  {stateRows, yPositions, primitives, xMin = -0.90, xMax = 0.90},
  stateRows = rowByState /@ states;
  yPositions = {1.0, 0.0};
  primitives = Flatten[
    Table[
      With[
        {
          row = stateRows[[i]],
          y = yPositions[[i]],
          comp = components[[j]],
          value = num[stateRows[[i]][components[[j]]["key"]]]
        },
        {
          {Directive[comp["color"], AbsoluteThickness[3.0], Dashing[comp["dash"]]],
            Line[{{0.0, y + comp["dy"]}, {value, y + comp["dy"]}}]},
          {Directive[comp["color"], AbsolutePointSize[5.0]], Point[{value, y + comp["dy"]}]}
        }
      ],
      {i, Length[stateRows]},
      {j, Length[components]}
    ],
    2
  ];
  Show[
    Graphics[
      {
        {Directive[GrayLevel[0.45], AbsoluteThickness[1.0]], Line[{{0.0, -0.55}, {0.0, 1.55}}]},
        primitives,
        legend[]
      },
      PlotRange -> {{xMin, xMax}, {-0.55, 1.55}},
      Frame -> True,
      Axes -> False,
      Background -> White,
      ImagePadding -> {{88, 14}, {48, 18}},
      FrameStyle -> Directive[GrayLevel[0.15], AbsoluteThickness[0.8]],
      FrameTicksStyle -> Directive[GrayLevel[0.15], 10],
      LabelStyle -> Directive[GrayLevel[0.08], 10],
      FrameTicks -> {
        {
          Thread[{yPositions, (tex[#, 10] & /@ (stateLabels /@ states))}],
          None
        },
        {
          Automatic,
          None
        }
      },
      FrameLabel -> {tex["G\\;\\mathrm{contribution}\\,[\\mathrm{GeV}^{-1}]", 12], None},
      PlotLabel -> tex[title, 13],
      ImageSize -> 342
    ]
  ]
];

panels = {
  componentGraphics[{"Ds1_2460", "Ds1_2536"}, "D_{s1}\\to D_s\\gamma"],
  componentGraphics[{"Bs1_5750", "Bs1_5830"}, "B_{s1}\\to B_s\\gamma"]
};

fig = GraphicsGrid[{panels}, Spacings -> {0.22, 0.0}, ImageSize -> 720, Background -> White];

Print["Exporting PDF"];
Export[outPdf, fig];
Print["Exporting PNG"];
Export[outPng, fig, ImageResolution -> 240];
Print["Wrote ", outPdf];
Print["Wrote ", outPng];
Quit[];
