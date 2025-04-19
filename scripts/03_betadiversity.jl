using DrWatson
@quickactivate :NetworkMonitoring

includet("main.jl")

d = Dict(
    :ns => 75,
    :nsites => 100,
    :C_exp => 0.2,
    :ra_sigma => 1.2,
    :ra_scaling => 50.0,
    :energy_NFL => 10_000,
    :H_nlm => 0.5,
    :nbon => 100,
    :refmethod => "metawebify",
)

Random.seed!(42)
res = main(d; res=:all)

bon = res[:bon]
m = render(Binary, res[:metaweb].metaweb)
n = res[:pos]
# extract(x -> render(Binary, x), r)

# Should I simplify and when?
# Should I rebuild species list based on ranges?
# networks = monitor(x -> simplify(render(Binary, x)), n, bon)
networks = monitor(x -> render(Binary, x), n, bon)


# betadiversity(βS, m, n)
# betadiversity(βOS, m, n)
# betadiversity(βWN, m, n)
# KGL08(betadiversity(βOS, m, n))

βcomponents = [betadiversity(βOS, m, n) for n in networks]
βosprime = KGL08.(βcomponents)

unique(βcomponents)
unique(βosprime)

p1 = hist(βosprime; axis=(; limits=((0.0, 1.0), (nothing, nothing))))


S, OS, WN = Float64[], Float64[], Float64[]
for i in 1:(length(networks)-1)
  for j in (i+1):length(networks)
    push!(S, KGL08(betadiversity(βS, simplify(networks[i]), simplify(networks[j]))))
    push!(OS, KGL08(betadiversity(βOS, networks[i], networks[j])))
    push!(WN, KGL08(betadiversity(βWN, simplify(networks[i]), simplify(networks[j]))))
  end
end
unique(S)
unique(OS)
unique(WN)

begin
    p2 = scatter(S, OS)
    scatter!(S, WN)
    p2
end
OS == WN # ?
