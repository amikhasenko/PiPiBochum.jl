### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ f353bda8-4c7f-4359-be8b-7d09981fa2fd
begin
	import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	using Plots
	using LaTeXStrings
	using Plots.PlotThemes
end

# ╔═╡ ba93d210-c761-4a74-84d6-2df23fc8338f
theme(:wong2, frame=:box, grid = false, ylim = (:auto, :auto))

# ╔═╡ 3e4217f4-4228-4b7b-9831-e6c448e4c892
plot(sin, -3, 4, ylim = (-1, 1), label = L"y_1")

# ╔═╡ 697345b5-173b-4c32-a4bc-95bee15fa6c6


# ╔═╡ 2a897086-3e73-4864-9508-b7de9d83d5bb


# ╔═╡ 46dadb9f-4a35-40a0-917a-68f0d93e4e6b


# ╔═╡ f7923b0c-fa02-431b-be2b-d8e5397ec5e6


# ╔═╡ Cell order:
# ╠═f353bda8-4c7f-4359-be8b-7d09981fa2fd
# ╠═ba93d210-c761-4a74-84d6-2df23fc8338f
# ╠═3e4217f4-4228-4b7b-9831-e6c448e4c892
# ╠═697345b5-173b-4c32-a4bc-95bee15fa6c6
# ╠═2a897086-3e73-4864-9508-b7de9d83d5bb
# ╠═46dadb9f-4a35-40a0-917a-68f0d93e4e6b
# ╠═f7923b0c-fa02-431b-be2b-d8e5397ec5e6
