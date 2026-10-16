# Description of this portfolio project

This is new AI-assisted portfolio work created on 16 October 2026. It does not
establish that Chris Kisia conducted the original authors' experiments or
performed this work during a previous appointment or degree.

The package contains runnable code, numerical results, tests and a coding
question. Before discussing it in an application or interview, run the
workflow, review the equations and explain the implementation and limitations.
Describe which parts you personally wrote, reviewed, changed or validated.

## Neutral project description

An AI-assisted Python portfolio project implements two-temperature virial
adsorption fitting, heat calculations and pressure-to-uptake inversion. It
includes a digitization and reproducibility audit of a published CO2 adsorption
plot, alongside a separate synthetic coding task with reference answers and
automated validation. Data provenance, original authorship, source licensing
and limits of the graphical reconstruction are documented.

## For the application's research-link field

Use the URL of your own repository containing this implementation, if the
employer accepts portfolio reproductions. The cited paper is the method/data
source, not evidence that you authored the paper or collected its data.

The form requests a calculation from your own research. A new AI-assisted
portfolio reproduction should be identified as such; it does not establish
original doctoral research or satisfy an unstated programming-experience
requirement by itself.

## Interview checks

- Explain why the fit is linear in the coefficients despite exponential pressure predictions.
- Explain why two temperatures are needed to distinguish the a and b terms.
- Explain the pressure and loading reference units inside the logarithm.
- Derive Qst from the derivative of ln(P) with respect to 1/T at fixed loading.
- Explain why inversion needs a bracket and a monotonicity check.
- Explain why a good synthetic fit cannot resolve a published-plot discrepancy.
