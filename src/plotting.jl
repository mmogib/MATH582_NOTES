using CairoMakie
using Random

rng = Xoshiro(2025)

seconds = 0:0.1:2
N = length(seconds)
measurements = rand(rng, 8.2:0.1:15.5, N)

seconds = 0:0.1:2
measurements = [8.2, 8.4, 6.3, 9.5, 9.1, 10.5, 8.6, 8.2, 10.5, 8.5, 7.2,
	8.8, 9.7, 10.8, 12.5, 11.6, 12.1, 12.1, 15.1, 14.7, 13.1]

lines(seconds, measurements)
scatter(seconds, measurements)

lines(seconds, exp.(seconds) .+ 7)

scatter(seconds, measurements)
lines!(seconds, exp.(seconds) .+ 7)
current_figure()

f = Figure(size = (900, 600))

ax1 = Axis(f[1, 1],
	title = "Experimental data",
	xlabel = "Time (seconds)",
	ylabel = "Value",
)
ax2 = Axis(f[1, 2],
	title = "Exponential fit",
	xlabel = "Time (seconds)",
	ylabel = "Value",
)
ax3 = Axis(f[2, 1:2],
	title = "Experimental data and exponential fit",
	xlabel = "Time (seconds)",
	ylabel = "Value",
)

scatter!(ax1, seconds, measurements)
lines!(ax2, seconds, exp.(seconds) .+ 7)
scatter!(ax3, seconds, measurements)
lines!(ax3, seconds, exp.(seconds) .+ 7,
	color = :tomato,
	label = "f(x)")

axislegend(position = :rb)

f

save("imgs/samlpe.svg", f)
