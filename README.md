# JohnsonSrc
Code used for generating precise Johnson Solids

The code found in `src/johnson.jl` was used to generate the files in `Oscar.jl`'s `data/JohnsonMatrices`. These matrices store the precise coordinates for the Johnson solids that require non-quadratic and/or non-simple number fields, which then is used in an active `Oscar` session to construct the corresponding `Polyhedron` object based on the $V$-description.

To replicate the generation on your own device, inject the code into an active `Oscar` session (e.g. by using `include("/path/to/johnson.jl")`) and call the function with the desired index `i`: `_johnson_solid(Val(i))`.
