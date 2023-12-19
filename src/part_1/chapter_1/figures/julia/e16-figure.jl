using LaTeXStrings, Plots
using Random
Random.seed!(2024)
N1 = 50
N2 = 50
x1 = randn(N1,2).+[2.5 1]
x2 = randn(N2,2)
x = [x1; x2]
I1 = 1:N1
I2 = (N1+1):(N1+N2)
scatter(x[I1,1], x[I1,2], label=L"i \in I_1", xlabel="feature 1", ylabel="feature 2", xlim=(-3,6), ylim=(-3,5))
scatter!(x[I2,1], x[I2,2], label=L"i \in I_2")
plot!([-2,5], -[-2,5]./0.45.+3.71, linestyle=:dash, linecolor=:black, label=false)
savefig("figure_e16.pdf")