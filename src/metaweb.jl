"""
    metawebify(m::T; binary=true) where T <: Metaweb

Build metaweb from scale field of `Metaweb` objects. Optionally convert to
binary using `binary=true.`
"""
function metawebify(m::T; binary::Bool=true) where T <: Metaweb
    m_adj = convert.(Matrix{Int}, adjacency(m.scale))
    m_acc = accumulate(+, vec(m_adj))
    m_end = m_acc[end]
    if binary
        m_end = Matrix{Int}(m_end .>= 1)
    end
    return m_end
end

"""
    extract(f::Function, m::T) where T <: Metaweb

Extract given measure with function `f` across local networks of `Metaweb`
object `m`. Similar to `map` but with different function name to avoid type
piracy.
"""
function extract(f::Function, m::T) where T <: Metaweb
    nets = SIS.networks(m)
    return f.(nets)
end
