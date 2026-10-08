using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using DataFrames
using CSV

M10 = 7.2
M20 = 9.1

p_cr = 42.0
N = 45.0
q = 1.0

tau1 = 28.0
tau2 = 22.0

p1 = 8.1
p2 = 10.5

u0 = [M10, M20]

mkpath("data")
mkpath("plots")

a1 =
    p_cr /
    (tau1^2 * p1^2 * N * q)

a2 =
    p_cr /
    (tau2^2 * p2^2 * N * q)

b =
    p_cr /
    (
        tau1^2 *
        tau2^2 *
        p1^2 *
        p2^2 *
        N *
        q
    )

c1 =
    (p_cr - p1) /
    (tau1 * p1)

c2 =
    (p_cr - p2) /
    (tau2 * p2)

println("Model coefficients:")
println("a1 = ", a1)
println("a2 = ", a2)
println("b  = ", b)
println("c1 = ", c1)
println("c2 = ", c2)

function firms_case1!(du, u, p, theta)
    M1, M2 = u

    du[1] =
        M1 -
        (b / c1) * M1 * M2 -
        (a1 / c1) * M1^2

    du[2] =
        (c2 / c1) * M2 -
        (b / c1) * M1 * M2 -
        (a2 / c1) * M2^2
end

problem1 = ODEProblem(
    firms_case1!,
    u0,
    (0.0, 30.0),
)

solution1 = solve(
    problem1,
    Tsit5();
    saveat=0.05,
    abstol=1e-9,
    reltol=1e-9,
)

theta1 = solution1.t
M1_case1 = getindex.(solution1.u, 1)
M2_case1 = getindex.(solution1.u, 2)

determinant =
    a1 * a2 - b^2

M1_stationary =
    (c1 * a2 - b * c2) /
    determinant

M2_stationary =
    (a1 * c2 - b * c1) /
    determinant

println()
println("Case 1 stationary state:")
println("M1* = ", M1_stationary)
println("M2* = ", M2_stationary)

social_factor = 0.00048

function firms_case2!(du, u, p, theta)
    M1, M2 = u

    du[1] =
        M1 -
        (b / c1 + social_factor) *
        M1 * M2 -
        (a1 / c1) * M1^2

    du[2] =
        (c2 / c1) * M2 -
        (b / c1) * M1 * M2 -
        (a2 / c1) * M2^2
end

problem2 = ODEProblem(
    firms_case2!,
    u0,
    (0.0, 30.0),
)

solution2 = solve(
    problem2,
    Tsit5();
    saveat=0.05,
    abstol=1e-9,
    reltol=1e-9,
)

theta2 = solution2.t
M1_case2 = getindex.(solution2.u, 1)
M2_case2 = getindex.(solution2.u, 2)

CSV.write(
    "data/firms_case1.csv",
    DataFrame(
        theta=theta1,
        firm1=M1_case1,
        firm2=M2_case1,
    ),
)

CSV.write(
    "data/firms_case2.csv",
    DataFrame(
        theta=theta2,
        firm1=M1_case2,
        firm2=M2_case2,
    ),
)

p_case1 = plot(
    theta1,
    M1_case1;
    label="Firm 1",
    linewidth=2,
    xlabel="Normalized time",
    ylabel="Working capital",
    title="Competition of firms: case 1",
)

plot!(
    p_case1,
    theta1,
    M2_case1;
    label="Firm 2",
    linewidth=2,
)

savefig(
    p_case1,
    "plots/firms_case1.png",
)

p_case2 = plot(
    theta2,
    M1_case2;
    label="Firm 1",
    linewidth=2,
    xlabel="Normalized time",
    ylabel="Working capital",
    title="Competition of firms: case 2",
)

plot!(
    p_case2,
    theta2,
    M2_case2;
    label="Firm 2",
    linewidth=2,
)

savefig(
    p_case2,
    "plots/firms_case2.png",
)

p_phase = plot(
    M1_case1,
    M2_case1;
    label="Trajectory",
    linewidth=2,
    xlabel="Firm 1",
    ylabel="Firm 2",
    title="Case 1: phase trajectory",
)

scatter!(
    p_phase,
    [M1_stationary],
    [M2_stationary];
    label="Stationary state",
    markersize=6,
)

savefig(
    p_phase,
    "plots/firms_case1_phase.png",
)

p_compare = plot(
    theta1,
    M1_case1;
    label="Firm 1, case 1",
    linewidth=2,
    xlabel="Normalized time",
    ylabel="Working capital",
    title="Comparison of competition models",
)

plot!(
    p_compare,
    theta1,
    M2_case1;
    label="Firm 2, case 1",
    linewidth=2,
)

plot!(
    p_compare,
    theta2,
    M1_case2;
    label="Firm 1, case 2",
    linewidth=2,
    linestyle=:dash,
)

plot!(
    p_compare,
    theta2,
    M2_case2;
    label="Firm 2, case 2",
    linewidth=2,
    linestyle=:dash,
)

savefig(
    p_compare,
    "plots/firms_comparison.png",
)

println()
println("Final values:")
println(
    "Case 1: M1 = ",
    last(M1_case1),
    ", M2 = ",
    last(M2_case1),
)

println(
    "Case 2: M1 = ",
    last(M1_case2),
    ", M2 = ",
    last(M2_case2),
)

println()
println("LAB07 BASE OK")
