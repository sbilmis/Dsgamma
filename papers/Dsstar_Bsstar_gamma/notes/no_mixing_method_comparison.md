# No-Mixing Same-Method Comparison

This is the intended benchmark check: switch off physical mixing and compare the pure basis-current LCSR result with literature calculations using the same current.

With our convention

\[
J_{1\mu}=\sin\theta\,J^A_\mu+\cos\theta\,J^B_\mu,
\]

the pure axial-current comparison corresponds to the basis quantity `g_A`, equivalently the lower-state current at `theta=90 deg` before adding tensor-current mixing.

## Results

- `D_s1(2460) -> D_s gamma`: Colangelo-like pure axial baseline gives `g_1=[-0.344,-0.245] GeV^-1` and `Gamma=[12.60,24.92] keV`, compared with the published `g_1=[-0.37,-0.29] GeV^-1` and `Gamma=[19,29] keV`.
- `D_s1 -> D_s* gamma`: Colangelo-like pure axial baseline gives `g_A*=[-0.192,-0.141]`, compared with the published `g*=[-0.18,-0.13]`.  This is a direct sign/normalization agreement.

The comparison table is stored in:

- `/Users/sbilmis/Dsgamma/papers/Dsstar_Bsstar_gamma/outputs/no_mixing_method_comparison.csv`
- `/Users/sbilmis/Dsgamma/papers/Dsstar_Bsstar_gamma/outputs/no_mixing_method_comparison_table.tex`
