

#
# Comparing convergence between algorithms
include("unconstrained_methods.jl")

f(x) = exp(-(x[1]-3)/2) + exp((4x[2] + x[1])/10) + exp((-4x[2] + x[1])/10) # standard function
xstart = [10;-4.5] # standard
xopt = [-(5/6)*(-3 + 2*log(2) - 2*log(5)); 0]
fopt = f(xopt)

exact_method = :golden_ratio
ϵ_reference = 1e-7

xg = gradient(f, xstart, exact_method, ϵ = ϵ_reference)
xn = newton_method(f, xstart, exact_method, ϵ = ϵ_reference)
xj = conjugate_gradient(f, xstart, exact_method, ϵ = ϵ_reference)
xb = bfgs(f, xstart, exact_method, ϵ = ϵ_reference)

dist_xg = sqrt.(sum(( xg .- xopt).^2,dims=1)');
dist_xn = sqrt.(sum(( xn .- xopt).^2,dims=1)');
dist_xj = sqrt.(sum(( xj .- xopt).^2,dims=1)');
dist_xb = sqrt.(sum(( xb .- xopt).^2,dims=1)');

plot(dist_xg, yscale=:log10, label = "Gradient")
plot!(dist_xn, yscale=:log10, label = "Newton")
plot!(dist_xj, yscale=:log10, label = "Conjugate")
plot!(dist_xb, yscale=:log10, label = "BFGS",
    xaxis = ("iterations", (1,15)),
    yaxis = (L"$||x_k - \overline{x}||$", ( 10*ϵ_reference, 10)))

savefig("convergence.pdf")


# xg = gradient(f, xstart, exact_method, ϵ = ϵ_reference)
# xn = newton_method(f, xstart, exact_method, ϵ = ϵ_reference)
# xj = conjugate_gradient(f, xstart, exact_method, ϵ = ϵ_reference)
# xb = bfgs(f, xstart, exact_method, ϵ = ϵ_reference)


# fxg = [f(xg[:,k]) for k=1:size(xg)[2]]
# fxn = [f(xn[:,k]) for k=1:size(xn)[2]]
# fxj = [f(xj[:,k]) for k=1:size(xj)[2]]
# fxb = [f(xb[:,k]) for k=1:size(xb)[2]]

# distf_xg = abs.(fxg .- fopt);
# distf_xn = abs.(fxn .- fopt);
# distf_xj = abs.(fxj .- fopt);
# distf_xb = abs.(fxb .- fopt);

# plot(distf_xg, label = "Gradient")
# plot!(distf_xn, label = "Newton")
# plot!(distf_xj,  label = "Conjugate")
# plot!(distf_xb,  label = "BFGS",
#     xaxis = ("iterations", (1,15)),
#     yaxis = (L"$||x_k - \overline{x}||$", (ϵ_reference)))