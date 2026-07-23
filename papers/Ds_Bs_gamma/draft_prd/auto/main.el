;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("revtex4-2" "aps" "prd" "preprint" "nofootinbib" "superscriptaddress")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("amsmath" "") ("amssymb" "") ("mathtools" "") ("bm" "") ("booktabs" "") ("graphicx" "") ("xcolor" "") ("hyperref" "colorlinks=true" "linkcolor=blue" "citecolor=blue" "urlcolor=blue")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "setfloatlink")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "homepage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "email")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "revtex4-2"
    "revtex4-210"
    "fontenc"
    "inputenc"
    "amsmath"
    "amssymb"
    "mathtools"
    "bm"
    "booktabs"
    "graphicx"
    "xcolor"
    "hyperref")
   (TeX-add-symbols
    "ii"
    "Ds"
    "Dst"
    "Dsone"
    "Bs"
    "Bst"
    "Bsone"
    "GeV"
    "MeV"
    "keV"
    "qqbar"
    "ssbar")
   (LaTeX-add-labels
    "sec:intro"
    "sec:framework"
    "eq:j5"
    "eq:basis-currents"
    "eq:physical-currents"
    "eq:fp-def"
    "eq:fi-def"
    "eq:e1"
    "eq:width"
    "eq:corr-physical"
    "eq:projector"
    "eq:basis-decomp"
    "eq:hadronic-invariant"
    "eq:heavy-prop"
    "eq:pi-decomp"
    "eq:borel-defs"
    "eq:pert-part"
    "eq:pi1-rotation"
    "eq:pi2-rotation"
    "eq:general-sum-rule"
    "eq:g1-sum-rule"
    "eq:g2-sum-rule"
    "sec:numerics"
    "tab:inputs"
    "eq:numerical-borel-choice"
    "tab:windows"
    "fig:stability"
    "tab:results"
    "eq:bs-low-width"
    "fig:mc"
    "eq:br-ds12536"
    "sec:concl"
    "app:sumrule-expressions"
    "eq:app-phi"
    "eq:app-psiv"
    "eq:app-A"
    "eq:app-B"
    "eq:app-S"
    "eq:app-Stilde"
    "eq:app-T1"
    "eq:app-T2"
    "eq:app-T3"
    "eq:app-T4"
    "eq:app-Psiv-int"
    "eq:app-B-dictionary"
    "eq:app-three-particle-integral"
    "eq:app-F-combination"
    "tab:app-photon-params"
    "app:basis-amplitudes"
    "eq:app-pi-decomp"
    "eq:app-pert")
   (LaTeX-add-bibliographies
    "all"))
 :latex)

