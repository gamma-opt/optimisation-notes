using Plots, LaTeXStrings
pyplot(size=(600,300))

include("unconstrained_methods.jl")

# which_plot = "coordinate_descent"
# which_plot = "gradient"
# which_plot = "gradient_fixed"
# which_plot = "newton"
# which_plot = "conjugate_gradient"
which_plot = "conjugate_gradient_vs_gradient"
# which_plot = "bfgs"
# which_plot = "bfgs_vs_newton"
# Plots all methods
# which_plot = "exact"
# which_plot = "armijo"

exact_method = :golden_ratio
save_pic = true

# function to be optimised
f(x) = exp(-(x[1]-3)/2) + exp((4x[2] + x[1])/10) + exp((-4x[2] + x[1])/10) # standard function
# f(x) = exp(-(0.1x[1]-3)/2) + exp((16x[2] + x[1])/10) + exp((-16x[2] + x[1])/10) # higher conditioning

# Calculating first- and second-order derivatives
∇(f,x) = ForwardDiff.gradient(f, x);
H(f,x) = ForwardDiff.hessian(f, x);

# Plotting the contours of the function to be optimised
n = 1000
x = range(-10, 25, length=n);
y = range(-10, 10, length=n);

contour(x,y,(x,y)-> f([x,y]),
        levels = [3.6, 4, 5, 6 , 8, 10, 15, 20, 30],
        # xaxis = (L"$x_1$", (-5,15)), # standard
        # yaxis = (L"$x_2$", (-5,5)),  
        xaxis = (L"$x_1$", (0, 20)),  # higher conditioning
        yaxis = (L"$x_2$", (-5,5)),
        clims = (0,35),
        contour_labels = true,
        cbar = false,
        aspect_ratio = :equal
        )

# The optimal value of x for standard function:
xopt = [-(5/6)*(-3 + 2*log(2) - 2*log(5)); 0]
fopt = f(xopt)

# xstart = [10;-4.5] # standard
xstart = [18;-1]   # high conditioning

if which_plot == "coordinate_descent"
        # Calling the methods and plotting their progress
        xc = coordinate_descent(f, xstart)
        plot!( xc[1,:], xc[2,:], label = "Coord. desc.", marker=:circle)
        save_pic ? savefig("coordinate_decent_exact.pdf") : 0

elseif which_plot == "gradient"
        xg = gradient(f, xstart, exact_method)
        plot!( xg[1,:], xg[2,:],label = "Gradient (exact)", marker=:circle)

        xg = gradient(f, xstart, :armijo)
        plot!( xg[1,:], xg[2,:],label = "Gradient (Armijo)", marker=:circle, line=:dash)
        save_pic ? savefig("gradient.pdf") : 0

elseif which_plot == "newton"
        xn = newton_method(f, xstart, exact_method)
        plot!(xn[1,:], xn[2,:],label = "Newton's method (exact)", marker=:circle)

        xn = newton_method(f, xstart, :armijo)
        plot!(xn[1,:], xn[2,:],label = "Newton's method (Armijo)", marker=:circle, line=:dash)
        save_pic ? savefig("newton.pdf") : println("Picture not saved.")

elseif which_plot == "conjugate_gradient"
        xj = conjugate_gradient(f, xstart, exact_method)
        plot!( xj[1,:], xj[2,:], label = "Conj. grad. (exact)", marker=:circle)

        xj = conjugate_gradient(f, xstart, :armijo)
        plot!( xj[1,:], xj[2,:], label = "Conj. grad. (Armijo)", marker=:circle, line=:dash)
        save_pic ? savefig("conjugate_gradient.pdf") : 0

elseif which_plot == "conjugate_gradient_vs_gradient"
        xj = conjugate_gradient(f, xstart, exact_method)
        plot!( xj[1,:], xj[2,:], label = "Conj. grad. (exact)", marker=:circle)

        xg = gradient(f, xstart, exact_method)
        plot!( xg[1,:], xg[2,:],label = "Gradient (exact)", marker=:circle, line=:dash)
        save_pic ? savefig("conjugate_vs_gradient.pdf") : 0

elseif which_plot == "bfgs"
        xb = bfgs(f, xstart, exact_method)
        plot!( xb[1,:], xb[2,:], label = "BFGS (exact)", marker=:circle)

        xb = bfgs(f, xstart, :armijo)
        plot!( xb[1,:], xb[2,:], label = "BFGS (Armijo)", marker=:circle, line=:dash)
        save_pic ? savefig("bfgs.pdf") : 0

elseif which_plot == "bfgs_vs_newton"
        xb = bfgs(f, xstart, exact_method)
        plot!( xb[1,:], xb[2,:], label = "BFGS (exact)", marker=:circle)

        xn = newton_method(f, xstart, exact_method)
        plot!(xn[1,:], xn[2,:],label = "Newton's method (exact)", marker=:circle, line=:dash)
        save_pic ? savefig("bfgs_vs_newton.pdf") : 0

elseif which_plot == "exact"
        xg = gradient(f, xstart, exact_method)
        plot!( xg[1,:], xg[2,:],label = "Gradient (exact)", marker=:circle)

        xn = newton_method(f, xstart, exact_method)
        plot!(xn[1,:], xn[2,:],label = "Newton (exact)", marker=:circle)

        xj = conjugate_gradient(f, xstart, exact_method)
        plot!( xj[1,:], xj[2,:], label = "Conj. grad. (exact)", marker=:circle)

        xb = bfgs(f, xstart, exact_method)
        plot!( xb[1,:], xb[2,:], label = "BFGS (exact)", marker=:circle)
        save_pic ? savefig("all_methods_exact.pdf") : 0

elseif which_plot == "armijo"
        xg = gradient(xstart, :armijo)
        plot!( xg[1,:], xg[2,:],label = "Gradient (Armijo)", marker=:circle, line=:dash)

        xn = newton_method(xstart, :armijo)
        plot!(xn[1,:], xn[2,:],label = "Newton (Armijo)", marker=:circle, line=:dash)

        xj = conjugate_gradient(xstart, :armijo)
        plot!( xj[1,:], xj[2,:], label = "Conj. grad. (Armijo)", marker=:circle, line=:dash)

        xb = bfgs(xstart, :armijo)
        plot!( xb[1,:], xb[2,:], label = "BFGS (Armijo)", marker=:circle, line=:dash)
        save_pic ? savefig("all_methods_armijo.pdf") : 0
end