# DefectSpeciesAction

!syntax description /DefectDynamics/DefectSpeciesAction

## Description

The `DefectSpeciesAction` is a MOOSE action that describes one defect species (e.g. a vacancy or
an interstitial) under the enclosing `DefectDynamics` block: its diffusivity, capillary
formation-energy law, and (for a loop) dislocation bias. This is invoked using a named sub-block
within [DefectDynamics](syntax/DefectDynamics/index.md), such as `[vacancy]` or `[interstitial]`.
See that page for a description, example usage, and parameters.
