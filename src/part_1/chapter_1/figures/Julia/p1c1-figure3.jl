using Plots, LaTeXStrings
pyplot()

n = 500
x1 = range(0,100,length=n)

# Setting standard colours for annotation
blue = RGBA(0.0,0.605603,0.97868,1.0)
orange = RGBA(0.888873,0.435649,0.278122)
green = RGBA(0.242224,0.643275,0.304449,1.0)
purple = RGBA(0.76444,0.444112,0.824298,1.0)
yellow = RGBA(0.675544,0.555662,0.0942343,1.0)

n = 500
x1 = range(0,2000,length=n)

#6x1 + 4x2 <= 24
plot(x1, (24 .- 6 .* x1)./4, color = :1, fill=(0,0.1),
    xaxis = (L"$x_1$", (0,6)),
        yaxis = (L"$x_2$", (0,6)),
        aspect_ratio = 1,
        size = (300,300),
        legend = false) #fill=(2000,0.1))

annotate!([(0.8, 5.0, text(L"$6x_1 + 4x_2 \leq 24$",8,:left, blue))])

#savefig("figure3_c1.pdf")


plot!(x1, (6 .- x1)./2.0, color = :2, fill=(0, 0.1))
annotate!([(4.5, 0.9, text(L"$x_1 + 2x_2 \leq 6$",8,:left, orange))])

#savefig("figure3_c2.pdf")


#-x1 + x2 <= 1
plot!(x1, 1 .+ x1, color = :3, fill=(0,0.1))
annotate!(3.6, 4.5, text(L"$-x_1 + x_2 \leq 1$", 8, :left, green))

#savefig("figure3_c3.pdf")


#x2 <= 2
plot!(x1, 0 .* x1 .+ 2 , color = :4, fill=(0,0.1))
annotate!(5.1, 1.8, text(L"$x_2 \leq 2$", 8, :left, purple))

savefig("figure3_c4.pdf")


#Clearing the fill of the constraints
#6x1 + 4x2 <= 24
plot(x1, (24 .- 6x1)./4, color = :1,
    xaxis = (L"$x_1$", (0,6)),
    yaxis = (L"$x_2$", (0,6)),
    aspect_ratio = 1,
    size = (300,300),
    legend = false) #fill=(2000,0.1))
annotate!([(0.8, 5.0, text(L"$6x_1 + 4x_2 \leq 24$",8,:left, blue))])

#x1 + 2x2 <= 6
plot!(x1, (6 .- x1)./2.0, color = :2) #, fill=(2000,0.1))
annotate!(4.5, 0.9, text(L"$x_1 + 2x_2 \leq 6$",8,:left, orange))

#-x1 + x2 <= 1
plot!(x1, 1 .+ x1, color = :3) #, fill=(0,0.1))
annotate!(3.6, 4.5, text(L"$-x_1 + x_2 \leq 1$",8,:left, green))

#x2 <= 2
plot!(x1, 0 .* x1 .+ 2 , color = :4) #, fill=(0,0.1))
annotate!(5.1, 2.2, text(L"$x_2 \leq 2$",8,:left, purple))

# Polyhedral set
plot!(Shape([0,0,1,2,3,4], [0,1,2,2,1.5,0]), opacity = 0.1, color = :6)

savefig("figure3_c5.pdf")


# Adding level curves
f(x) = 5x[1] + 4x[2] 

plot!(x1, (5 .- 5x1)/4, color = :5, line = :dash) #, fill=(0,0.1))
annotate!(0, 1.3, text(L"$z=5$",7,:right, yellow))

plot!(x1, (10 .- 5x1)/4, color = :5, line = :dash) #, fill=(0,0.1))
annotate!(-0.1, 2.5, text(L"$z=10$",7,:right, yellow))

plot!(x1, (15 .- 5x1)/4, color = :5, line = :dash) #, fill=(0,0.1))
annotate!(0, 3.7, text(L"$z=15$",7,:right, yellow))

annotate!(3.2, 2.3, text(L"$\nabla z$", 8, :right, yellow))
quiver!([2],[5/4], quiver=([5/5],[4/5]), color = :5)

savefig("figure3_c6.pdf")


# optimal point 
plot!(x1, (21 .- 5x1)/4, color = :5, line = :dash) #, fill=(0,0.1))

scatter!([3],[1.5], color = :orange)
annotate!(3.2, 1.5, text(L"$x^* = (3, 1.5)$", 8,:left))
annotate!(3.2, 1.8, text(L"$z^* = 21$", 8, :left))

savefig("figure3_complete.pdf")


# Animation

#6x1 + 4x2 <= 24
plt = plot(x1, (24 .- 6x1)./4, color = :1,
    xaxis = (L"$x_1$", (0,6)),
    yaxis = (L"$x_2$", (0,6)),
    aspect_ratio = 1,
    size = (300,300),
    legend = false) #fill=(2000,0.1))
annotate!([(0.8, 5.0, text(L"$6x_1 + 4x_2 \leq 24$",8,:left, blue))])

#x1 + 2x2 <= 6
plot!(x1, (6 .- x1)./2.0, color = :2) #, fill=(2000,0.1))
annotate!([(4.5, 0.9, text(L"$x_1 + 2x_2 \leq 6$",8,:left, orange))])

#-x1 + x2 <= 1
plot!(x1, 1 .+ x1, color = :3) #, fill=(0,0.1))
annotate!([(3.6, 4.5, text(L"$-x_1 + x_2 \leq 1$",8,:left, green))])

#x2 <= 2
plot!(x1, 0 .* x1 .+ 2 , color = :4) #, fill=(0,0.1))
annotate!([(5.1, 2.2, text(L"$x_2 \leq 2$",8,:left, purple))])

# f(x) = 5x[1] + 4x[2] 

anim = @animate for i = 1:6
    plot!(x1, (i*3.5 .- 5x1)/4, color = :5, line = :dash) #, fill=(0,0.1))
#     annotate!([(-0.1, 1.3, text(L"$z=5$",7,:right, yellow))])
    if i == 6
        scatter!([3],[1.5], color = :orange)
        annotate!([(3.2, 1.5, text(L"$x^* = (x_1^* = 3, x_2^* = 1.5)$",8,:left))])
        annotate!([(3.2, 1.8, text(L"$z^* = 21$",8,:left))])
    end    
end

#gif(anim, "anim_fps15.gif", fps = 1)