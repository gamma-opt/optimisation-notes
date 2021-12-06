## Figure 2 - robust ellipsoids

using Plots, Distributions, LaTeXStrings
pyplot()

# creates an ellipse
pts = Plots.partialcircle(0, 2π, 200)
x1, y1 = Plots.unzip(pts)
x2, y2 = Plots.unzip(pts)
x1 = 2x1 .+ 10
y1 += 3y1 .+ 5
x2 = 3x2 .+ 10
y2 += 5.5y2 .+ 5

#plot(Shape(x2, y2), label = L"\Gamma_2", opacity = .2)
#plot!(Shape(x1, y1), label = L"\Gamma_1", opacity = .2)

a1 = rand(Normal(10, 2), 100)
a2 = rand(Normal(5, 3), 100)

scatter(
    a1,
    a2,
    label = "Observations",
    title = L"Distribution of ($a_{1}$, $a_2$)",
    xaxis = (L"a_1"),
    yaxis = (L"a_2"),
    color = 1
)
scatter!((10, 5), shape = :square, label = "Nominal value", color = 2)

savefig("Figure2.pdf")

plot!(Shape(x2, y2), label = L"\Gamma_2=1.5", opacity = .2)
plot!(Shape(x1, y1), label = L"\Gamma_1=1", opacity = .2)

savefig("Figure2_with_ellip.pdf")
