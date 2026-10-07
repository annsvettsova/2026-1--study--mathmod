# # Параметрическое исследование
#
# Исследуем влияние отношения скоростей катера и лодки
# на траекторию преследования для варианта 58.

using DrWatson
@quickactivate "project"

using Plots
using DataFrames
using CSV

k = 20.2
phi = 3pi / 4

speed_ratios = [3.0, 4.0, 5.1, 6.0]

results = DataFrame(
    n = Float64[],
    start_radius = Float64[],
    interception_radius = Float64[],
    interception_x = Float64[],
    interception_y = Float64[],
)

p = plot(
    xlabel="x, км",
    ylabel="y, км",
    title="Влияние отношения скоростей на траекторию",
    aspect_ratio=:equal,
)

max_radius = 0.0

for n in speed_ratios
    x1 = k / (n + 1)

    theta0 = 0.0
    theta = range(theta0, phi; length=800)

    r = x1 .* exp.((theta .- theta0) ./ sqrt(n^2 - 1))

    xc = r .* cos.(theta)
    yc = r .* sin.(theta)

    r_int = x1 * exp((phi - theta0) / sqrt(n^2 - 1))
    px = r_int * cos(phi)
    py = r_int * sin(phi)

    push!(results, (n, x1, r_int, px, py))

    plot!(
        p,
        xc,
        yc;
        label="n = $(n)",
        linewidth=2,
    )

    global max_radius = max(max_radius, r_int)
end

boat_r = range(0.0, 1.1 * max_radius; length=400)

plot!(
    p,
    boat_r .* cos(phi),
    boat_r .* sin(phi);
    label="Лодка",
    linewidth=2,
    linestyle=:dash,
)

mkpath("data")
mkpath("plots")

CSV.write("data/pursuit_parameters.csv", results)
savefig(p, "plots/pursuit_parameters.png")

println(results)
println()
println("LAB01 PARAMETERS OK")
