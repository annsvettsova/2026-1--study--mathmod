using DrWatson
@quickactivate "project"

using Plots
using DataFrames
using CSV

k = 20.2
n = 5.1

phi = 3pi / 4

x1 = k / (n + 1)
x2 = k / (n - 1)

println("k = ", k)
println("n = ", n)
println("x1 = ", x1)
println("x2 = ", x2)

spiral_r(theta, r0, theta0, n) =
    r0 * exp((theta - theta0) / sqrt(n^2 - 1))

theta01 = 0.0
theta1 = range(theta01, phi; length=800)

r_case1 = spiral_r.(theta1, x1, theta01, n)

xc1 = r_case1 .* cos.(theta1)
yc1 = r_case1 .* sin.(theta1)

r_intersection1 = spiral_r(phi, x1, theta01, n)

p1x = r_intersection1 * cos(phi)
p1y = r_intersection1 * sin(phi)

theta02 = -pi
theta2 = range(theta02, phi; length=1200)

r_case2 = spiral_r.(theta2, x2, theta02, n)

xc2 = r_case2 .* cos.(theta2)
yc2 = r_case2 .* sin.(theta2)

r_intersection2 = spiral_r(phi, x2, theta02, n)

p2x = r_intersection2 * cos(phi)
p2y = r_intersection2 * sin(phi)

max_radius = 1.08 * max(r_intersection1, r_intersection2)

boat_r = range(0.0, max_radius; length=500)

boat_x = boat_r .* cos(phi)
boat_y = boat_r .* sin(phi)

result = DataFrame(
    case = ["Первый случай", "Второй случай"],
    start_radius = [x1, x2],
    intersection_radius = [r_intersection1, r_intersection2],
    intersection_x = [p1x, p2x],
    intersection_y = [p1y, p2y],
)

mkpath("data")
CSV.write("data/pursuit_intersections.csv", result)

println()
println(result)

p1 = plot(
    boat_x,
    boat_y;
    label="Лодка",
    linewidth=2,
    xlabel="x, км",
    ylabel="y, км",
    title="Задача о погоне — первый случай",
    aspect_ratio=:equal,
)

plot!(
    p1,
    [k, x1],
    [0.0, 0.0];
    label="Прямолинейный участок катера",
    linewidth=2,
)

plot!(
    p1,
    xc1,
    yc1;
    label="Катер",
    linewidth=2,
)

scatter!(
    p1,
    [p1x],
    [p1y];
    label="Точка пересечения",
    markersize=5,
)

mkpath("plots")
savefig(p1, "plots/pursuit_case1.png")

p2 = plot(
    boat_x,
    boat_y;
    label="Лодка",
    linewidth=2,
    xlabel="x, км",
    ylabel="y, км",
    title="Задача о погоне — второй случай",
    aspect_ratio=:equal,
)

plot!(
    p2,
    [k, -x2],
    [0.0, 0.0];
    label="Прямолинейный участок катера",
    linewidth=2,
)

plot!(
    p2,
    xc2,
    yc2;
    label="Катер",
    linewidth=2,
)

scatter!(
    p2,
    [p2x],
    [p2y];
    label="Точка пересечения",
    markersize=5,
)

savefig(p2, "plots/pursuit_case2.png")

p3 = plot(
    boat_x,
    boat_y;
    label="Лодка",
    linewidth=2,
    xlabel="x, км",
    ylabel="y, км",
    title="Сравнение двух траекторий катера",
    aspect_ratio=:equal,
)

plot!(p3, xc1, yc1; label="Катер — случай 1", linewidth=2)
plot!(p3, xc2, yc2; label="Катер — случай 2", linewidth=2)

scatter!(
    p3,
    [p1x, p2x],
    [p1y, p2y];
    label="Точки пересечения",
    markersize=5,
)

savefig(p3, "plots/pursuit_comparison.png")

println()
println("LAB01 BASE OK")
