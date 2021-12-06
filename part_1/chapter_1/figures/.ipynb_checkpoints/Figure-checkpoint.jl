using Plots, Distributions, LatexStrings

## Figure 2 - robust ellipsoids

a1 = rand(Normal(10,2), 100)
a2 = rand(Normal(5,3), 100)

scatter(a1, a2, label = "observations", title="Distribution of (a1, a2)")
scatter!((10,5), shape =:square,label="nominal")

# create an ellipse in the sky
pts = Plots.partialcircle(0,2π,200)
x1, y1 = Plots.unzip(pts)
x2, y2 = Plots.unzip(pts)
x1 = 2x1 + 10
y1 += 3y1 + 5
x2 = 3x2 + 10
y2 += 5.5y2 + 5

plot!(Shape(x2,y2), label="Gamma2", opacity=.1)
plot!(Shape(x1,y1), label="Gamma1", opacity=.1)

savefig("Figure2.pdf")


## Figure 3 - classifier

using JuMP, Ipopt

# Generating random multivariate data.
data11 = rand(Normal(8, 4), 100)
data12 = rand(Normal(10, 2), 100)
data1 = [data11 data12]

data01 = rand(Normal(2, 6), 100)
data02 = rand(Normal(5, 2), 100)
data0 = [data01 data02]

scatter(data1[:,1], data1[:,2], label="Pos. obs.")
scatter!(data0[:,1], data0[:,2], label="Neg. obs.")

function classify(data1, data0, δ)
    m = Model(solver = IpoptSolver())

    @variable(m, u[1:100] >= 0)
    @variable(m, v[1:100] >= 0)
    @variable(m, a[1:2])
    @variable(m, b)

    @constraint(m, error1[i=1:100], sum(a[j]*data1[i,j] for j = 1:2) + b - u[i] <= -1)
    @constraint(m, error0[i=1:100], sum(a[j]*data0[i,j] for j = 1:2) + b + v[i] >= 1)

# This parameter trades off accuracy and robustness for the classifier
    @objective(m, Min, sum(u[i] for i=1:100) + sum(v[i] for i=1:100) + δ*(sum(a[i]^2 for i=1:2)))

    solve(m);

    return getvalue(a), getvalue(b)
end

a, b = classify(data1, data0, 0.5)

# Calculating the points on plane defined by a and b
x1 = linspace(-10,20,100);
x2 = -(b + a[1]*x1)/ a[2];
# boundaries of the slab
x2u = -(b + 1 + a[1]*x1)/ a[2];
x2d = -(b - 1 + a[1]*x1)/ a[2];
plot!(x1, [x2u x2 x2d],
    xlabel = "x_1",
    ylabel = "x_2",
    label= [ "", "Classifier", ""] ,
    color=[:green :green :green],
    linestyle=[ :dot :solid :dot],
    )

savefig("Figure3.pdf")
