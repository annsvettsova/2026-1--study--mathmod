using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using DataFrames
using CSV

x0 = 0.1
v0 = 1.1

u0 = [x0, v0]

tspan = (0.0, 63.0)
save_step = 0.05

mkpath("data")
mkpath("plots")

function oscillator1!(du, u, p, t)
    x, v = u

    du[1] = v
    du[2] = -3.7x
end

prob1 = ODEProblem(oscillator1!, u0, tspan)

sol1 = solve(
    prob1,
    Tsit5();
    saveat=save_step,
    abstol=1e-9,
    reltol=1e-9,
)

t1 = sol1.t
x1 = getindex.(sol1.u, 1)
v1 = getindex.(sol1.u, 2)

function oscillator2!(du, u, p, t)
    x, v = u

    du[1] = v
    du[2] = -3v - 10x
end

prob2 = ODEProblem(oscillator2!, u0, tspan)

sol2 = solve(
    prob2,
    Tsit5();
    saveat=save_step,
    abstol=1e-9,
    reltol=1e-9,
)

t2 = sol2.t
x2 = getindex.(sol2.u, 1)
v2 = getindex.(sol2.u, 2)

function oscillator3!(du, u, p, t)
    x, v = u

    du[1] = v
    du[2] = 0.9sin(0.9t) - 3v - 11x
end

prob3 = ODEProblem(oscillator3!, u0, tspan)

sol3 = solve(
    prob3,
    Tsit5();
    saveat=save_step,
    abstol=1e-9,
    reltol=1e-9,
)

t3 = sol3.t
x3 = getindex.(sol3.u, 1)
v3 = getindex.(sol3.u, 2)

CSV.write(
    "data/oscillator_case1.csv",
    DataFrame(time=t1, x=x1, velocity=v1),
)

CSV.write(
    "data/oscillator_case2.csv",
    DataFrame(time=t2, x=x2, velocity=v2),
)

CSV.write(
    "data/oscillator_case3.csv",
    DataFrame(time=t3, x=x3, velocity=v3),
)

p1 = plot(
    t1,
    x1;
    label="Без затухания",
    linewidth=2,
    xlabel="t",
    ylabel="x(t)",
    title="Гармонический осциллятор",
)

plot!(p1, t2, x2; label="С затуханием", linewidth=2)
plot!(p1, t3, x3; label="С внешней силой", linewidth=2)

savefig(p1, "plots/oscillator_solutions.png")

p2 = plot(
    x1,
    v1;
    label="Без затухания",
    linewidth=2,
    xlabel="x",
    ylabel="dx/dt",
    title="Фазовые портреты",
)

plot!(p2, x2, v2; label="С затуханием", linewidth=2)
plot!(p2, x3, v3; label="С внешней силой", linewidth=2)

savefig(p2, "plots/oscillator_phase.png")

p3 = plot(
    t1,
    x1;
    label="x(t)",
    linewidth=2,
    xlabel="t",
    ylabel="x",
    title="Без затухания",
)
savefig(p3, "plots/oscillator_case1.png")

p4 = plot(
    t2,
    x2;
    label="x(t)",
    linewidth=2,
    xlabel="t",
    ylabel="x",
    title="С затуханием",
)
savefig(p4, "plots/oscillator_case2.png")

p5 = plot(
    t3,
    x3;
    label="x(t)",
    linewidth=2,
    xlabel="t",
    ylabel="x",
    title="С затуханием и внешней силой",
)
savefig(p5, "plots/oscillator_case3.png")

println("Случай 1: max |x| = ", maximum(abs.(x1)))
println("Случай 2: max |x| = ", maximum(abs.(x2)))
println("Случай 3: max |x| = ", maximum(abs.(x3)))

println()
println("LAB03 BASE OK")
