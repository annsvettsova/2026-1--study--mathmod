using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using DataFrames
using CSV

N = 17_854.0

I0 = 199.0
R0 = 35.0
S0 = N - I0 - R0

u0 = [S0, I0, R0]

alpha = 1.0e-5
beta = 0.02

tspan = (0.0, 120.0)
save_step = 0.1

mkpath("data")
mkpath("plots")

println("Population N = ", N)
println("S(0) = ", S0)
println("I(0) = ", I0)
println("R(0) = ", R0)
println("alpha = ", alpha)
println("beta = ", beta)

function epidemic_isolated!(du, u, p, t)
    S, I, R = u

    du[1] = 0.0
    du[2] = -beta * I
    du[3] = beta * I
end

problem1 = ODEProblem(
    epidemic_isolated!,
    u0,
    tspan,
)

solution1 = solve(
    problem1,
    Tsit5();
    saveat=save_step,
    abstol=1e-9,
    reltol=1e-9,
)

t1 = solution1.t
S1 = getindex.(solution1.u, 1)
I1 = getindex.(solution1.u, 2)
R1 = getindex.(solution1.u, 3)

function epidemic_active!(du, u, p, t)
    S, I, R = u

    infection = alpha * S * I
    recovery = beta * I

    du[1] = -infection
    du[2] = infection - recovery
    du[3] = recovery
end

problem2 = ODEProblem(
    epidemic_active!,
    u0,
    tspan,
)

solution2 = solve(
    problem2,
    Tsit5();
    saveat=save_step,
    abstol=1e-9,
    reltol=1e-9,
)

t2 = solution2.t
S2 = getindex.(solution2.u, 1)
I2 = getindex.(solution2.u, 2)
R2 = getindex.(solution2.u, 3)

CSV.write(
    "data/epidemic_isolated.csv",
    DataFrame(
        time=t1,
        susceptible=S1,
        infected=I1,
        recovered=R1,
    ),
)

CSV.write(
    "data/epidemic_active.csv",
    DataFrame(
        time=t2,
        susceptible=S2,
        infected=I2,
        recovered=R2,
    ),
)

p1 = plot(
    t1,
    S1;
    label="S(t)",
    linewidth=2,
    xlabel="Time",
    ylabel="Population",
    title="Case 1: I(0) <= I*",
)

plot!(p1, t1, I1; label="I(t)", linewidth=2)
plot!(p1, t1, R1; label="R(t)", linewidth=2)

savefig(
    p1,
    "plots/epidemic_isolated.png",
)

p2 = plot(
    t2,
    S2;
    label="S(t)",
    linewidth=2,
    xlabel="Time",
    ylabel="Population",
    title="Case 2: I(0) > I*",
)

plot!(p2, t2, I2; label="I(t)", linewidth=2)
plot!(p2, t2, R2; label="R(t)", linewidth=2)

savefig(
    p2,
    "plots/epidemic_active.png",
)

p3 = plot(
    t1,
    I1;
    label="I(0) <= I*",
    linewidth=2,
    xlabel="Time",
    ylabel="Infected",
    title="Comparison of epidemic regimes",
)

plot!(
    p3,
    t2,
    I2;
    label="I(0) > I*",
    linewidth=2,
)

savefig(
    p3,
    "plots/epidemic_comparison.png",
)

peak_index = argmax(I2)

println()
println("Case 1 final:")
println("S = ", last(S1))
println("I = ", last(I1))
println("R = ", last(R1))

println()
println("Case 2:")
println("Peak infected = ", I2[peak_index])
println("Peak time = ", t2[peak_index])
println("Final S = ", last(S2))
println("Final I = ", last(I2))
println("Final R = ", last(R2))

println()
println("Population check, case 1 = ", last(S1 + I1 + R1))
println("Population check, case 2 = ", last(S2 + I2 + R2))

println()
println("LAB05 BASE OK")
