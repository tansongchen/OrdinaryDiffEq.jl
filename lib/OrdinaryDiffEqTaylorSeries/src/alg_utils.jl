alg_order(::ExplicitTaylor2) = 2
alg_stability_size(alg::ExplicitTaylor2) = 1

alg_order(::ExplicitTaylor{P}) where {P} = P
alg_stability_size(alg::ExplicitTaylor) = 1

alg_order(alg::ExplicitTaylorAdaptiveOrder) = get_value(alg.min_order)
get_current_adaptive_order(::ExplicitTaylorAdaptiveOrder, cache) = cache.current_order[]
get_current_alg_order(::ExplicitTaylorAdaptiveOrder, cache) = cache.current_order[]

TaylorScalar{T, P}(x::TaylorScalar{T, Q}) where {T, P, Q} = TaylorScalar{P}(x)

const JET_CACHE = IdDict()

function make_term(a)
    term(TaylorScalar, Symbolics.unwrap(a.value), map(Symbolics.unwrap, a.partials))
end

function get_value(::Val{P}) where {P}
    return P
end

function build_jet(f::ODEFunction{iip}, p, order, length = nothing) where {iip}
    key = (f, order, p)
    if haskey(JET_CACHE, key)
        return JET_CACHE[key]
    end
    @variables t0::Real
    u0 = isnothing(length) ? Symbolics.variable(:u0) : Symbolics.variables(:u0, 1:length)
    if iip
        @assert length isa Integer
        f0 = similar(u0)
        f(f0, u0, p, t0)
    else
        f0 = f(u0, p, t0)
    end
    u = TaylorDiff.make_seed(u0, f0, Val(1))
    for index in 2:order
        t = TaylorScalar{index - 1}(t0, one(t0))
        if iip
            fu = similar(u)
            f(fu, u, p, t)
        else
            fu = f(u, p, t)
        end
        d = get_coefficient(fu, index - 1) / index
        u = append_coefficient(u, d)
    end
    u_term = make_term.(u)
    jet = build_function(u_term, u0, t0; expression = Val(false), cse = true)
    if !haskey(JET_CACHE, key)
        JET_CACHE[key] = jet
    end
    return jet
end

# evaluate using Qin Jiushao's algorithm
@generated function evaluate_polynomial(t::TaylorScalar{T, P}, z) where {T, P}
    ex = :(v[$(P + 1)])
    for i in P:-1:1
        ex = :(v[$i] + z * $ex)
    end
    return :($(Expr(:meta, :inline)); v = flatten(t); $ex)
end
