# SamplerFullSolveMultiApp

The [SamplerFullSolveMultiApp](#) simply creates a full-solve type sub application (see [MultiApps])
for each row of each matrix returned from the [Sampler](stochastic_tools/index.md#samplers) object.

This object is capable of running in batch mode by setting the 'mode' parameter. For more
information refer to [batch_mode.md].


## Example Syntax

!listing modules/stochastic_tools/test/tests/multiapps/sampler_full_solve_multiapp/master_full_solve.i block=MultiApps

!syntax parameters /MultiApps/SamplerFullSolveMultiApp

!syntax inputs /MultiApps/SamplerFullSolveMultiApp

!syntax children /MultiApps/SamplerFullSolveMultiApp
