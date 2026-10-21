# CorPresDefl

Compute three transformed corneal net pressure–volume slopes over
Initial–A1, A1–HC, and HC–A2.

## Install

```r
install.packages("pak")
pak::pak("GulumimiRen/CorPresDefl")
```

Requires R 4.1.0 or later.

## Usage

```r
library(CorPresDefl)

corpd(
  A1t = 7.5, HCt = 17, A2t = 22, bIOP = 14,
  x_a1 = 0.2, x_hc = 4, x_a2 = 0.4
)

# Data frame with the native pipe
df <- data.frame(
  id = c("example_1", "example_2"),
  A1t = c(7.5, 8), HCt = c(17, 17.5), A2t = c(22, 22.5),
  bIOP = c(14, 15),
  A1DeflVol = c(0.2, 0.25), HCDeflVol = c(4, 4.5),
  A2DeflVol = c(0.4, 0.5)
)
df |> corpd()

# Custom volume column names
custom <- data.frame(A1t = 7.5, HCt = 17, A2t = 22, bIOP = 14,
                     volume_a1 = 0.2, volume_hc = 4, volume_a2 = 0.4)
custom |> corpd(
  columns = c(x_a1 = "volume_a1", x_hc = "volume_hc", x_a2 = "volume_a2")
)
```

## Inputs and outputs

| Parameter | Description | Unit |
|---|---|---|
| `A1t`, `HCt`, `A2t` | Device times of first applanation, highest concavity, and second applanation | ms |
| `bIOP` | Biomechanically corrected IOP from the same examination | mmHg |
| `x_a1`, `x_hc`, `x_a2` | Calculated DeflVol at A1, HC, and A2 | mm³ = μL |

Data-frame defaults are `A1t`, `HCt`, `A2t`, `bIOP`, `A1DeflVol`,
`HCDeflVol`, and `A2DeflVol`. Numeric inputs may be vectors of equal
length or scalars. The `measure` argument accepts `"Vol"` (the default).

| Output | Interval | Pressure transformation |
|---|---|---|
| `k_initial_a1` | Initial–A1 | Exponential, scale 4 mmHg |
| `k_a1_hc` | A1–HC | Natural logarithm |
| `k_hc_a2` | HC–A2 | Exponential, scale 20 mmHg |

Numeric calls return these three columns and `status`. Data-frame calls
preserve the input and append `Vol_k_initial_a1`, `Vol_k_a1_hc`,
`Vol_k_hc_a2`, and `Vol_status`. Use `columns` for custom input names
and `prefix` for custom output names. All three slopes have units mm⁻³.

Valid inputs satisfy `0 < A1t < HCt < A2t <= 30`, positive bIOP and event
volumes, HC volume greater than A1 and A2, and positive A1 and HC net
pressures. Invalid inputs produce `NA` slopes and a reason in `status`.
Use `invalid = "warn"`, `"quiet"`, or `"error"` to select validation handling.

## Algorithm

Reference on-eye pressure is interpolated from the bundled pressure–time
table using monotone cubic Hermite interpolation on the rising and falling
branches. Net pressure is
$P(t)=P_{\mathrm{on-eye}}(t)-\mathrm{bIOP}$.
Write $V=\mathrm{DeflVol}$, with $V_0=0$ and $P_0=-\mathrm{bIOP}$.
The three endpoint slopes are:

$$k_{\mathrm{Initial-A1}}=
\frac{\exp[P(A1t)/c_1]-\exp[-\mathrm{bIOP}/c_1]}{V(A1t)},
\qquad c_1=4\,\mathrm{mmHg}.$$

$$k_{\mathrm{A1-HC}}=
\frac{\ln[P(HCt)/(1\,\mathrm{mmHg})]-\ln[P(A1t)/(1\,\mathrm{mmHg})]}
{V(HCt)-V(A1t)}.$$

$$k_{\mathrm{HC-A2}}=
\frac{\exp[P(A2t)/c_2]-\exp[P(HCt)/c_2]}{V(A2t)-V(HCt)},
\qquad c_2=20\,\mathrm{mmHg}.$$

## References

Joda AA, Shervin MMS, Kook D, Elsheikh A. Development and validation of a
correction equation for Corvis tonometry. *Computer Methods in Biomechanics
and Biomedical Engineering*. 2016;19(9):943–953.
[doi:10.1080/10255842.2015.1077515](https://doi.org/10.1080/10255842.2015.1077515).

Ren Q, Deng Z, Chu Z, Pan Z, Wan P. DeflVol: a model-derived volumetric
descriptor of corneal deflection at highest concavity under Corvis ST
air-puff loading. *Frontiers in Bioengineering and Biotechnology*.
2026;14:1933733.
[doi:10.3389/fbioe.2026.1933733](https://doi.org/10.3389/fbioe.2026.1933733).

## License

R code: [MIT](LICENSE). Bundled reference data: CC BY 4.0;
see [NOTICE](inst/NOTICE).
