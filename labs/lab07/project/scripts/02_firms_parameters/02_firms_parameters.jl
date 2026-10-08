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

social_values = [
    0.0,
    0.00048,
    0.00096,
]

tspan = (0.0, 30.0)

mkpath("data")
mkpath("plots")

results = DataFrame(
    social_factor=Float64[],
    final_firm1=Float64[],
    final_firm2=Float64[],
)

plot_result = plot(
    xlabel="Normalized time",
    ylabel="Working capital",
    title="Influence of social-psychological factor",
)

for social_factor in social_values

    function model!(du, u, p, theta)
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

    problem = ODEProblem(
        model!,
        u0,
        tspan,
    )

    solution = solve(
        problem,
        Tsit5();
        saveat=0.05,
        abstol=1e-9,
        reltol=1e-9,
    )

    theta = solution.t
    M1 = getindex.(solution.u, 1)
    M2 = getindex.(solution.u, 2)

    push!(
        results,
        (
            social_factor,
            last(M1),
            last(M2),
        ),
    )

    plot!(
        plot_result,
        theta,
        M1;
        label="Firm 1, s=$(social_factor)",
        linewidth=2,
    )
end

CSV.write(
    "data/firms_parameter_scan.csv",
    results,
)

savefig(
    plot_result,
    "plots/firms_parameter_scan.png",
)

println(results)
println()
println("LAB07 PARAMETERS OK")
