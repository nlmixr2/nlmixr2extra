# Example single dose Theophylline ODE model

This is a nlmixr2 model that is pre-run so that it can be used in
package testing and development. It is regenerated whenever binaries of
`nlmixr2extra` are created. If there is a binary incompatability between
the fit objects, a simple rerun of the installation will fix this
nlmixr2 fit object.

## Format

A (modified) data frame with 132 rows and 22 columns.

- ID:

  Patient identifier

- TIME:

  Time (hr)

- DV:

  Dependent variable (concentration)

- PRED:

  Predictions without any between subject variability

- RES:

  Population Residual

- WRES:

  Weighted Residuals under the FO assumption

- IPRED:

  Individual Predictions

- IRES:

  Individual Residuals

- IWRES:

  Individual Weighted Residuals

- CPRED:

  Conditional Prediction under the FOCE assumption

- CRES:

  Conditional Residuals under the FOCE assumption

- CWRES:

  Conditional Weighted Residuals under the FOCE assumption

- eta.ka:

  Between subject changes for ka

- eta.cl:

  Between subject changes for v

- depot:

  amount in the depot compartment

- center:

  amount in the central compartment

- ka:

  Individual ka values

- cl:

  Individual cl values

- v:

  Individual volume of distribution

- tad:

  Time after dose

- dosenum:

  Dose number
