using Plots, ForwardDiff, LaTeXStrings, LinearAlgebra
pyplot(size = (600,200))

include("unconstrained_methods.jl");

# Select function
f(x) = (1/2)*(x[1]^2 + 2x[2]^2)
# f(x) = (1/2)*(x[1]^2 + 10x[2]^2)

# Plotting the contours of the function to be optimised
n = 1000
x = range(-20,20,length=n);
y = range(-10,10,length=n);
z = [f([x[i],y[j]]) for j = 1:n, i = 1:n];

# Checking eigenvalues to calculate condition number
H(f,x) = ForwardDiff.hessian(f, x);
xstart = [10; 1] # starting point
exact_method = :golden_ratio
ϵ_reference = 1e-7

# Calculating the conditio number
eigens = eigen(H(f, xstart))
κ = maximum(eigens.values)/ minimum(eigens.values)

contour(x,y,z, # Change title accordingly
        title = L"$f(x) = (1/2)(x_1^2 + 10x_2^2), \kappa=$" * string(κ),
        levels = [0.5, 2, 6, 10, 15, 25, 40, 60],
        xaxis = (L"$x_1$", (-12,12)),
        yaxis = (L"$x_2$", (-4,4)),
        clims = (0,70),
        clabels = true,
        legend = false,
        aspect_ratio = :equal,
        titlefontsize = 9
        )

xg = gradient(f, xstart, exact_method, ϵ = ϵ_reference)

#Plotting progress gradient descent
plot!( xg[1,:], xg[2,:],label = "Gradient", marker=:circle)


savefig("gradient_k2.pdf")
# savefig("gradient_k10.pdf")
