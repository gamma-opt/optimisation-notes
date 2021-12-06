## Figure 3 - classifier

using JuMP, Ipopt, Distributions, LaTeXStrings, Plots
pyplot()

# Generating random multivariate data.
data11 = rand(Normal(10, 1.5), 100)
data12 = rand(Normal(8, 1.5), 100)
data1 = [data11 data12]

data01 = rand(Normal(2, 1.5), 100)
data02 = rand(Normal(2, 1.5), 100)
data0 = [data01 data02]

scatter(
    data1[:, 1],
    data1[:, 2],
    label = "Pos. obs.",
    xaxis = (L"x_1", (-5,15)),
    yaxis = (L"x_2", (-5,15))
)

scatter!(data0[:, 1], data0[:, 2], label = "Neg. obs.")

savefig("classes_no_classifier.pdf")


function classify(data1, data0, δ)
    m = Model(with_optimizer(Ipopt.Optimizer))

    @variable(m, u[1:100] >= 0)
    @variable(m, v[1:100] >= 0)
    @variable(m, a[1:2])
    @variable(m, b)

    @constraint(m, error1[i = 1:100],
        sum(a[j] * data1[i, j] for j = 1:2) + b - u[i] <= -1
    )
    @constraint(m, error0[i = 1:100],
        sum(a[j] * data0[i, j] for j = 1:2) + b + v[i] >= 1
    )

    # This parameter trades off accuracy and robustness for the classifier
    @objective(m, Min,
        sum(u[i] for i = 1:100) + sum(v[i] for i = 1:100) +
        δ * (sum(a[i]^2 for i = 1:2))
    )

    optimize!(m)

    return value.(a), value(b)
end

# Classify without robustness
a, b = classify(data1, data0, 0.0)

# Calculating the points on plane defined by a and b
x1 = range(-10, stop = 20, length = 100)
x2 = -(b .+ a[1] .* x1) / a[2]

plot!(
    x1,
    x2,
    xaxis = (L"x_1", (-5,15)),
    yaxis = (L"x_2", (-5,15)),
    label = "Classifier",
    color = 3,
)

savefig("classes_with_classifier.pdf")

# Classify with robustness
δ = 0.1
a, b = classify(data1, data0, δ)

# Boundaries of the slabs
x2u = -(b .+ 1 .+ a[1] .* x1) / a[2]
x2d = -(b .- 1 .+ a[1] .* x1) / a[2]
x2 = -(b .+ a[1] .* x1) / a[2]

plot!(
    x1,
    [x2u x2 x2d],
    xaxis = (L"x_1", (-5,15)),
    yaxis = (L"x_2"),
    label = [" " L"Classifier $\delta = 0.1$" " "],
    color = [4 4 4],
    linestyle = [:dot :solid :dot]
)

savefig("classes_with_robust_classifier.pdf")
