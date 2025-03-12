function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{ExplicitTaylor2Cache, ExplicitTaylor2ConstantCache
}}
    dense ? "specialized 2th order \"free\" interpolation" : "1st order linear"
end

function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{ExplicitTaylorCache, ExplicitTaylorConstantCache
}}
dense ? "specialized nth order \"free\" interpolation" : "1st order linear"
end
