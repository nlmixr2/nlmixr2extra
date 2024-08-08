# nlmixr2extra 2.0.9

* `bootstrapFit()` now will be more careful handling `NA` values so
  they do not completely affect results (Issue #59)

* `bootstrapFit()` will now only take the correlation of the non-zero
  diagonals (Issue #59).

* New method for `knit_print()` will generate model equations for LaTeX
  reporting automatically.

* Tests are now skipped if they contain linear compartment models that
  need gradients when the gradients are not compiled (as in the case
  of intel c++).
