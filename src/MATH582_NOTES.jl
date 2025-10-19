### A Pluto.jl notebook ###
# v0.20.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 41c749c0-500a-11f0-0eb8-49496afa257e
begin
    using CommonMark
    using PlutoUI, PlutoExtras
    using Plots, PlotThemes, LaTeXStrings
    using Latexify
    using HypertextLiteral
    using Colors
    using LinearAlgebra, Random, Printf, SparseArrays, Statistics
    # using Symbolics
    # using SymPy
    using QRCoders
    using PrettyTables
	# using Primes
	using Polyhedra
    # using LinearSolve
    # using NonlinearSolve
    # using ForwardDiff
    # using Integrals
    # using OrdinaryDiffEq
	# using IntervalArithmetic
end

# ╔═╡ a0a41512-efe8-4848-ae7d-16e1606d0978
TableOfContents(title="📚 MATH582: Nonlinear Programming", indent=true, depth=4)

# ╔═╡ a6148830-0424-4d85-af87-9efa87eaa4aa
md"# Course Website"

# ╔═╡ 45572d31-300f-4e40-a755-9c099c58551a
md"# Syllabus"

# ╔═╡ ed6a2919-8df2-4907-b7cb-93e4cf4c0500
cm"""
[Please read the syllabus of the course](https://www.dropbox.com/scl/fi/jhogaom1fl083tvwi730n/T251_MATH582_Syllabus.pdf?rlkey=twh1ndt5olgqkqglkohl7lwyl&dl=0)
"""

# ╔═╡ 7bc65c80-3346-452e-b4de-e45ce3a19461
md"## Orientation (15 min)"

# ╔═╡ fcc354c2-c077-4a5a-9b84-db1d5f9f4ee7
cm"""
1. __What is Nonlinear Programming? (Motivation)__

* Optimization of nonlinear objectives with equality/inequality constraints:
  ```math
  \min f(x) \quad \text{s.t. } g_i(x) \leq 0, \; h_j(x) = 0
  ```
* Extends **linear programming** → models many real-world problems.
* Applications:
  * **Economics:** utility & production optimization.
  * **Engineering:** design of structures, energy systems, petroleum operations.
  * **AI/ML:** training models, regularization, constrained optimization.
* **Theme of the course:** Convexity → tractable optimization, strong theory, powerful algorithms.

---

2. __Course Logistics__

* **Instructor:** Dr. Mohammed Alshahrani
* **Textbook:** Bazaraa–Sherali–Shetty (*Nonlinear Programming: Theory and Algorithms*, 3rd ed.)
* **Programming:** Julia (official docs + *Think Julia*).
* **Lecture Notes:** All Pluto.jl notebooks posted at
  👉 [MATH582 Notes](https://mmogib.github.io/MATH582_NOTES/) (QR code).

---

3. __Tools & Support__
* **MATH582 TA (Custom GPT Assistant)**
  👉 [MATH582 TA GPT](https://chatgpt.com/g/g-689da632cee8819192005e8adf53b82b-math582-ta) (QR code).
  * Guides lecture content, proofs, and examples.
  * Assists with Julia programming & assignments.
  * Provides **feedback and grade access**.
  * Use responsibly: **first attempt solutions, then use TA for hints and clarification.**

---

4. __Evaluation & Projects__
* Homework & programming assignments: **10%**
* **Individual Project: 20%**
  * Week 3: proposal (1 page).
  * Week 8: progress report (1 page).
  * Week 14: 20-min presentation.
  * Week 15: final report (≤20 pages, similarity < 20%).
* Exam 1 (20%), Exam 2 (20%), Final Exam (30%).

---

5. __Expectations__
* **Mathematical rigor:** proofs, derivations, optimality conditions.
* **Computational rigor:** Julia implementation with complexity & stability checks.
* **Academic integrity:** originality required; oral defense of assignments may be used.
* **Active learning strategy:**
  * Work interactively with Pluto notebooks.
  * Use MATH582 TA as a **personal tutor**, not as a shortcut.

"""

# ╔═╡ 9d275485-9e6f-450e-8392-787ba3cda9a6
md"# Chapter 2: Convex Sets"



# ╔═╡ a889687a-ec1d-400d-901d-894e81fb5549
md"""
## 2.1 Convex Hulls

__Learning outcomes__

> 1. Define convex sets and convex combinations, with examples and counterexamples.  
> 2. Apply closure properties of convex sets (intersection, Minkowski sum, Minkowski difference).  
> 3. Describe convex hulls as both (i) all convex combinations of points, and (ii) the smallest convex set containing them.  
> 4. State Carathéodory’s theorem and explain its geometric/algorithmic significance.  
"""


# ╔═╡ f26b7a51-3583-42da-8022-76a7aa6fec5a

let 

# Axes limits (adjust if you like)
xlims = (-4.0, 4.0)
ylims = (0.0, 5.0)

# Grid for shading the feasible region x2 ≥ |x1|
nx, ny = 500, 500
xs = range(xlims[1], xlims[2], length=nx)
ys = range(ylims[1], ylims[2], length=ny)

# Mask: 1.0 where condition holds, NaN elsewhere (for contourf fill)
Z = [ (y ≥ abs(x)) ? 1.0 : NaN for y in ys, x in xs ]

p = contourf(xs, ys, Z; levels=1, fill=true, legend=false, alpha=0.5,
             color=:lightblue, xlabel="x₁", ylabel="x₂",
             title="S = {(x₁,x₂): x₂ ≥ |x₁|}", ratio=:equal, frame_style=:origin)

# Draw boundary lines x2 =  x1 and x2 = -x1
xs_n = range(xlims[1], 0, length=nx)
xs_p = range(0,xlims[2], length=nx)
plot!(xs_n,  -xs_n; lw=2, color=:blue, label="")
plot!(xs_p, xs_p; lw=2, color=:blue, label="")

p
end


# ╔═╡ 9e538d00-9783-436c-87f9-e0b9187cf5ba
let
	
	# 1) Initialize Xoshiro RNG
	rng = Xoshiro(582)              # seed for reproducibility
	
	# 2) Generate random points (40 in R²)
	pts = rand(rng, 40, 2)
	pts_= map(r->[r...],eachrow(pts)) |> collect
	p=polyhedron(convexhull(pts_...))
	# 3) Convex hull via Polyhedra (V-representation)
	# v = vrep(pts)
	# p = polyhedron(v)
	removevredundancy!(p)
	# p,polyhedron(convexhull([1, 0], [0, 1]))
	# p=polyhedron(convexhull(pts[:]))
	# removevredundancy!(p)
	# # # 4) Plot hull + points
	plot(p, ratio=:equal, alpha=0.3, label="conv(pts)")
	scatter!(pts[:,1], pts[:,2], ms=3, color=:black, label="points")
end

# ╔═╡ cbf0cb8d-0fc4-4886-87b4-a7da72210ed4
md"""
## 2.2 Closure and Interior of a Set
"""

# ╔═╡ 46f49abf-eb7c-47cf-8e27-eecaad8a2cab
cm"""
**Learning Objectives**
> 1. Define neighborhoods in ``\mathbb{R}^n``.
> 2. Define closure via neighborhoods; define closed sets by ``S = \operatorname{cl}(S)``.
> 3. Derive the **sequential definition** of closed sets as a result.
> 4. Recognize compact sets and basic consequences.
> 5. State and apply **Theorem 2.2.2** (with ``S`` convex) and its corollaries.
"""

# ╔═╡ ca2c639a-a8c5-4c31-9732-fa3bdc284347


# ╔═╡ 6bf7f29f-8c9b-46eb-ad4d-244336a02e73
md"## 2.3 Weierstrass’s Theorem"


# ╔═╡ 6a6d3f65-6cd4-4672-a007-98c92fbcb76d

cm"""
__Learning Objectives__
> By the end of this lecture, students should be able to:
> 1. **State and explain** Weierstrass’s Theorem for existence of solutions in optimization.
> 2. **Differentiate** between minimum, maximum, infimum, and supremum.
> 3. **Illustrate** cases where a minimum does not exist despite boundedness of the objective.
> 4. **Apply** Weierstrass’s Theorem to prove existence of solutions to nonlinear programs.
> 5. **Experiment** computationally with bounded and unbounded sets in Julia.
"""

# ╔═╡ 8e2102e2-5dee-4198-8326-3c83b8d07cf2
let
	f(x) = x^2
	plot(f, 0, 1, label="f(x) = x^2", linewidth=2)
	scatter!([0], [0], color=:red, label="Minimum at x=0", frame_style=:origin)
	
end

# ╔═╡ 818c7aac-3310-468c-b0c0-388bca923064
cm"""
__Computational Counterexample__
Consider ``f(x) = 1/x`` on ``S = (0,1]``.
- ``f`` is continuous on ``(0,1]``, but ``S`` is **not closed**.
- ``\inf f(x) = 1``, but there is **no minimizer**, since ``x=0 ∉ S``.
"""

# ╔═╡ 59d6f09b-2132-4437-a69f-4a51172ef6ff
let
	f2(x) = 1/x
	plot(f2, 0.001, 1, label="f(x) = 1/x", linewidth=2, ylims=(0,10), frame_style=:origin)
end

# ╔═╡ 13369f35-e62d-4f84-830d-58999af481a8
cm"""
__Formative Check__
1. Define infimum and supremum in your own words.
2. Give an example of a function on a non-compact domain where no minimum exists.
3. State Weierstrass’s Theorem formally.
4. Why is compactness essential in the theorem?
5. Suppose ``f(x) = e^x`` on ``[0,∞)``. Does a minimizer exist? Why or why not?
"""

# ╔═╡ 7958b6c5-f805-41ee-b290-c01abc51547f
cm"""
__Reading__
- **Main Textbook:** Bazaraa–Sherali–Shetty (2006), *Nonlinear Programming: Theory and Algorithms*, 3rd ed., Section 2.3, pp. 48–49.
- **Supplementary:** Rockafellar, *Convex Analysis*, Ch. 2 (for background on compactness and lower semicontinuity).

"""

# ╔═╡ 41a713fd-1e21-4034-9b77-34db9a67a8e6
md"## 2.4 Separation and Support of Sets "

# ╔═╡ dc4a1115-b3e0-45d7-ad83-9debf9e55e60
cm"""
__Learning outcomes__

> 1. Understand the concept of separating hyperplanes in convex analysis.
> 2. Derive and prove the separation theorems for convex sets.
"""

# ╔═╡ e721d5be-313a-4df9-b37e-e560c9bcaa97
md"### Hyperplanes and Separation of Two Sets "

# ╔═╡ 7559a2b5-9a9e-43db-9605-e2262f0a5c83
cm"P"

# ╔═╡ 9034672d-b285-4e44-b472-63f5d1876477
# let
# 	# θ values
# 	θ = range(0, 2π, length=500)
	
# 	# parametric circle (center (1,1), radius 1)
# 	x = 1 .+ cos.(θ)
# 	y = 1 .+ sin.(θ)
	
# 	# fill the circle
# 	plot(x, y, 
# 		 xlimits = (-4,2),
# 		 ylimits = (-2,5),
# 		 frame_style = :origin,
# 		 seriestype=:shape, c=:blue, alpha=0.3, label="")
	
# 	# outline
# 	plot!(x, y, color=:blue, linewidth=2, label="")
	
# 	# formatting
# 	plot!(aspect_ratio=:equal, xlabel="x", ylabel="y",
# 	      title="Circle centered at (1,1), radius 1", grid=true)
# 	plot!(-10:10,x->-x)
# 	quiver!([0.0], [0.0], quiver=([-4.0], [1.0]), color=:red, linewidth=2, label="y = [-4,1]")
# end

# ╔═╡ 357b4c18-6b9d-407b-b46a-e470f61e1bf4
# let
# 	p1 = -[1;1]
# 	p2 = [1;1]
# 	y = [-8, -1]
	
# 	# parameter values
# 	t = range(0, 2π, length=10)
	
# 	# circle points
# 	points = [ [1 + cos(tt), 1 + sin(tt)] for tt in t ]
	
# 	p1y, spx = p1 ⋅ y, map(x->p1 ⋅ x, points) |> maximum
# 	p2y, ipx = p2 ⋅ y, map(x->p2 ⋅ x, points) |> minimum
# 	cm"""
# 	``p^Ty=``$p1y ``\qquad \sup \{p^Tx: x
# 	\in S\}=`` $spx

# 	``p^Ty=``$p2y ``\qquad \inf \{p^Tx: x
# 	\in S\}=`` $ipx
# 	"""
# end

# ╔═╡ cc560f5c-ee2e-42a3-b0c5-b2f7e5460a74
cm"P"

# ╔═╡ 23d1eb38-8aaf-4c91-b6a6-f7eb508b88dd
md"### Support of Sets at Boundary Points"

# ╔═╡ 702156d9-521c-47ee-9965-3059684a5a8d
md"### Separation of Two Convex Sets "

# ╔═╡ b7858895-5fb6-4d58-9ff4-b51f55db32ca
md"## 2.5 Convex Cones and Polarity"

# ╔═╡ 5bdaef6d-410f-4421-b2cf-e113cb40fe5c
cm"""
__Learning outcomes__

> 1. Understand the concept of convex cones and their basic properties.
> 2. Define and explore the polarity of convex cones.
> 3. Apply polarity in the analysis of convex optimization problems.
"""

# ╔═╡ cc6206a9-e2b1-4f31-9d5c-368e5e20ba58
md"## 2.6 Polyhedral Sets, Extreme Points, and Extreme Directions"

# ╔═╡ 8e354fee-c8c6-4fa4-a6a8-f13e5e1c1ab2
cm"""
__Learning outcomes__

> 1. Understand the definition and structure of polyhedral sets.
> 2. Define and identify extreme points and extreme directions.
> 3. Characterize polyhedral sets in terms of their extreme points and directions.
"""

# ╔═╡ 3a82aafe-3d11-4743-bcfc-0eb2061b94c2
md"# Chapter 3: Convex Functions and Generalizations"

# ╔═╡ a5eea233-4510-4b7b-b7dd-51650c4a9300
md"## 3.1 Definitions and Basic Properties "

# ╔═╡ 683c53d5-7db0-49ee-b2c7-e0b87248fdd4
md"### Continuity of Convex Functions "

# ╔═╡ 71bb4d0e-8000-44ee-96c9-a356ab2afe3d
md"### Directional Derivative of Convex Functions "

# ╔═╡ cd0a647d-eeac-4798-99e0-a10549033f48
md"## 3.2 Subgradients of Convex Functions"

# ╔═╡ 2f340ca4-9fe2-4910-bec3-9fb86c7ed6a7
md"### Epigraph and Hypograph of a Function "

# ╔═╡ 788874c6-a13f-4718-86d5-bdbba1b588c8
md"## 3.3 Differentiable Convex Functions "

# ╔═╡ af8c90dd-3f32-40ed-91d1-1287919109c4
md"### Twice Differentiable Convex and Concave Functions "

# ╔═╡ b2cbaa21-30b7-4eaa-bbba-2a100b4e3f7f
let
	H = [2 1 2
			1 2 3
		 2 3 4
		]
	map(i->det(H[1:i,1:i]), 1:3)
	isposdef(H)
	M = [
    -2  -1   0;
    -1  -3  -1;
     0  -1  -2
]
	map(i->det(M[1:i,1:i]), 1:3)
	eigen(-M)
end

# ╔═╡ e8637759-f3b0-4acd-8d37-0d61098d8b16
md"## 3.4 Minima and Maxima of Convex Functions"

# ╔═╡ 6e85e225-d792-4ad7-aade-e048c78f62e2
md"### Minimizing a Convex Function "

# ╔═╡ 2adf10df-7961-4329-ad96-f031d405586f
let
	# Define the constraints:
	# -x₁ + x₂ ≤ 2  →  x₂ ≤ x₁ + 2
	# 2x₁ + 3x₂ ≤ 11  →  x₂ ≤ (11 - 2x₁)/3
	# -x₁ ≤ 0  →  x₁ ≥ 0
	# -x₂ ≤ 0  →  x₂ ≥ 0
	
	# Create a range of x₁ values
	x1_range = 0:0.1:6
	
	# Define constraint functions
	constraint1(x1) = x1 + 2          # x₂ ≤ x₁ + 2
	constraint2(x1) = (11 - 2*x1)/3   # x₂ ≤ (11 - 2x₁)/3
	
	# Create the plot
	plt = plot(xlims=(-1, 6), ylims=(-5, 6), 
	           xlabel="x₁", ylabel="x₂", 
	           title="Feasible Region", 
	           legend=:bottomright,
	           grid=true, gridwidth=1, gridcolor=:lightgray,
			   frame_style=:origin
			  )
	
	# Plot constraint lines
	plot!(plt, x1_range, constraint1.(x1_range), 
	      label="-x₁ + x₂ = 2", linewidth=2, color=:blue)
	plot!(plt, x1_range, constraint2.(x1_range), 
	      label="2x₁ + 3x₂ = 11", linewidth=2, color=:red)
	
	# Add boundary lines for non-negativity constraints
	vline!(plt, [0], label="x₁ = 0", linewidth=2, color=:green)
	hline!(plt, [0], label="x₂ = 0", linewidth=2, color=:orange)
	
	# Find intersection points to define the feasible region vertices
	vertices = []
	
	# Intersection of x₁ = 0 and x₂ = 0
	push!(vertices, (0, 0, -0.5, 0.4))
	
	# Intersection of x₁ = 0 and -x₁ + x₂ = 2
	push!(vertices, (0, 2,-0.5,0.2))
	
	# Intersection of -x₁ + x₂ = 2 and 2x₁ + 3x₂ = 11
	# Solving: -x₁ + x₂ = 2 and 2x₁ + 3x₂ = 11
	# From first: x₂ = x₁ + 2
	# Substituting: 2x₁ + 3(x₁ + 2) = 11 → 5x₁ + 6 = 11 → x₁ = 1
	# Therefore: x₂ = 3
	push!(vertices, (1, 3,0.0, 0.5))
	
	# Intersection of 2x₁ + 3x₂ = 11 and x₂ = 0
	# 2x₁ + 3(0) = 11 → x₁ = 5.5
	push!(vertices, (5.5, 0,0.0,0.5))
	
	# Extract x and y coordinates for filling
	x_coords = [v[1] for v in vertices]
	y_coords = [v[2] for v in vertices]
	
	# Fill the feasible region
	plot!(plt, x_coords, y_coords, 
	      seriestype=:shape, 
	      fillalpha=0.3, 
	      fillcolor=:lightblue,
	      linecolor=:black,
	      linewidth=2,
	      label="Feasible Region")
	
	# Mark the vertices
	scatter!(plt, x_coords, y_coords, 
	         markersize=6, 
	         markercolor=:black,
	         label="Vertices")
	
	# Add vertex labels
	for (i, (x, y,dx,dy)) in enumerate(vertices)
	    annotate!(plt, x + dx, y + dy, text("($x, $y)", 10, :black))
	end
	# Define the objective function f(x₁, x₂) = (x₁ - 3/2)² + (x₂ - 5)²
	f(x1, x2) = (x1 - 3/2)^2 + (x2 - 5)^2
	
	# Create a grid for contour plotting
	x1_grid = 0:0.1:6
	x2_grid = 0:0.1:5
	X1 = repeat(x1_grid', length(x2_grid), 1)
	X2 = repeat(x2_grid, 1, length(x1_grid))
	Z = f.(X1, X2)
	
	# Add contour lines
	contour!(plt, x1_grid, x2_grid, Z, 
	         levels=10, 
	         color=:purple, 
	         linewidth=1.5,
	         linestyle=:dash,
	         label="f(x₁,x₂) contours")
	
	# Mark the center point of the function (3/2, 5)
	scatter!(plt, [1.5], [5], 
	         markersize=8, 
	         markercolor=:purple,
	         markershape=:star,
	         label="Center (3/2, 5)")
	
	# Add annotation for the center
	annotate!(plt, 1.5 + 0.1, 5 + 0.5, text("(3/2, 5)", 10, :purple))

	# Add vector [-1, -4] starting at point (1, 3)
	start_point = (1, 3)
	vector = [-1, -4]
	vector = vector/norm(vector)
	end_point = (start_point[1] + vector[1], start_point[2] + vector[2])
	
	# Plot the vector as an arrow
	plot!(plt, [start_point[1], end_point[1]], [start_point[2], end_point[2]], 
	      arrow=true, 
	      arrowsize=0.1,
	      linewidth=3, 
	      color=:darkgreen,
	      label="Vector [-1, -4]")
	
	# Add annotation for the vector
	annotate!(plt, start_point[1] - 0.3, start_point[2] - 0.5, 
	          text("[-1, -4]", 10, :darkgreen))
	# Update the title to reflect both elements
	plot!(plt, 
		  title="Feasible Region with Contours of"*L"f(x_1,x_2) = (x_1-3/2)^2 + (x_2-5)^2",
		 titlefontsize=8)
	plt
end

# ╔═╡ 9f88e2bc-4137-42a0-bcc2-e4d377c27f00
md"### Maximizing a Convex Function "

# ╔═╡ 15af64e4-366d-464e-962a-42bf619f2e4e
md"## 3.5 Generalizations of a Convex Functions "

# ╔═╡ a874e129-7812-4515-8ffe-d875c73813e5
md"### Differentiable Quasiconvex Functions "

# ╔═╡ 3efb4f55-9325-4e02-b53f-d275dbb405b3
md"### Strictly Quasiconvex Functions "

# ╔═╡ c3f319d2-23fc-4c41-a7f2-23f831073ae0
cm"""
```math
f(x)= \begin{cases}1 & \text { if } x=0 \\ 0 & \text { if } x \neq 0\end{cases}
```

By Definition, ``f`` is strictly quasiconvex. However, ``f`` is not quasiconvex,
"""

# ╔═╡ 7deefd68-4955-4fdb-bc76-1de2ccd841c8
md"### Strongly Quasiconvex Functions"

# ╔═╡ 12e8d772-374b-438f-a5b9-df7fdab33d4a
md"### Pseudoconvex Functions "

# ╔═╡ 9cc4f6fc-b656-48d7-8c8f-da23d1e5419c
md"### Convexity at a Point "

# ╔═╡ a03f6f56-54a8-423a-a405-5f9a83e5cdb2
md"# Chapter 4: The Fritz John and Karush-Kuhn-Tucker Optimality Conditions"

# ╔═╡ 8f0be524-cabb-4b78-99a6-9ecbfc6e57a3
md"## 4.1 Unconstrained Problems "

# ╔═╡ ab11e8fc-98e8-4372-a12f-c6133dcc65e3
md"### Necessary Optimality Conditions "

# ╔═╡ 7dfbcce2-1f0f-4f43-9e01-05bd0447ed32


# ╔═╡ 3d3223ea-d105-45a5-9389-81de85a271ba
md"### Sufficient Optimality Conditions "

# ╔═╡ 807a80ef-2c20-400d-85b6-9fc71d21b5b8
md"## 4.2 Problems Having Inequality Constraints"

# ╔═╡ 72ccae91-7342-4780-92c6-b226ce13507c
cm"""
> We first develop a necessary optimality condition for the problem
> ```math
> \begin{array}{lll}
> \min & f(x) \\
> \text{subject to}\\
> & x \in S\\
> \end{array}
> ```
> for a general set ``S \subseteq \mathbb{R}^n``.
"""

# ╔═╡ ff54a546-4bd6-4d6e-9ea7-f4053635dc42
md"### Geometric Optimality Conditions"

# ╔═╡ 1391def6-0a0d-4e16-abbe-3644e56cab9b
cm"""
> We now consider the problem
> ```math
> \begin{array}{lll}
> \min & f(x) \\
> \text{subject to}\\
> & g_i(x)\leq 0 & i=1,\cdots,m\\
> & x \in X
> \end{array}
> ```
> where ``g_i: R^n \rightarrow R`` for ``i=1, \ldots, m`` and ``X`` is a nonempty open set in ``\mathbb{R}^n``.
"""

# ╔═╡ 4719b6e7-4eba-43a4-8b8d-3bdf3a46712f
begin
	x1_val_html = @bind x1_val NumberField(0.0:0.01:3.0, default=1.5)

	x2_val_html = @bind x2_val NumberField(0.0:0.01:3.0, default=1.0)
	cm"""

	``x_1=`` $(x1_val_html)
	
	``x_2=`` $(x2_val_html)
	
	"""
end

# ╔═╡ 10754f5c-e457-47a5-91d6-5f6ed203cda8
let
	# Define the objective function
	f(x1, x2) = (x1 - 3)^2 + (x2 - 2)^2
	
	# Gradient of objective function
	∇f(x1, x2) = [2*(x1 - 3), 2*(x2 - 2)]
	
	# Constraint functions (written as g(x) ≤ 0)
	g1(x1, x2) = x1^2 + x2^2 - 5  # circle constraint
	g2(x1, x2) = x1 + x2 - 3      # line constraint
	g3(x1, x2) = -x1              # x1 ≥ 0
	g4(x1, x2) = -x2              # x2 ≥ 0
	
	# Gradients of constraints
	∇g1(x1, x2) = [2*x1, 2*x2]
	∇g2(x1, x2) = [1, 1]
	∇g3(x1, x2) = [-1, 0]
	∇g4(x1, x2) = [0, -1]
	
	# Check which constraints are active (within tolerance)
	tol = 0.001
	function active_constraints(x1, x2)
		active = []
		if abs(g1(x1, x2)) < tol
			push!(active, (1, "x₁² + x₂² ≤ 5", ∇g1(x1, x2)))
		end
		if abs(g2(x1, x2)) < tol
			push!(active, (2, "x₁ + x₂ ≤ 3", ∇g2(x1, x2)))
		end
		if abs(g3(x1, x2)) < tol
			push!(active, (3, "x₁ ≥ 0", ∇g3(x1, x2)))
		end
		if abs(g4(x1, x2)) < tol
			push!(active, (4, "x₂ ≥ 0", ∇g4(x1, x2)))
		end
		return active
	end
	
	# Check if point is feasible
	function is_feasible(x1, x2)
		return g1(x1, x2) <= tol && g2(x1, x2) <= tol && 
		       g3(x1, x2) <= tol && g4(x1, x2) <= tol
	end

	# Create the plot
	x1_range = range(-0.5, 3.5, length=400)
	x2_range = range(-0.5, 3.5, length=400)
	p = plot(size=(800, 800), 
			 aspect_ratio=:equal, 
			 frame_style=:origin,
	         xlabel="x₁", ylabel="x₂", 
	         title="Gradients at ($(round(x1_val, digits=2)), $(round(x2_val, digits=2)))",
	         legend=:topright, legendfontsize=8)
	
	# Plot contour curves
	contour!(p, x1_range, x2_range, 
	         (x1, x2) -> f(x1, x2), 
	         levels=15, 
	         color=:viridis, 
	         linewidth=1.5,
	         alpha=0.5,
	         colorbar=false,
	         label="")
	
	# Plot constraint boundaries
	θ = range(0, 2π, length=200)
	circle_x1 = sqrt(5) .* cos.(θ)
	circle_x2 = sqrt(5) .* sin.(θ)
	plot!(p, circle_x1, circle_x2, 
	      linewidth=2.5, color=:red, 
	      label="x₁² + x₂² = 5")
	
	x1_line = range(0, 3, length=100)
	x2_line = 3 .- x1_line
	plot!(p, x1_line, x2_line, 
	      linewidth=2.5, color=:blue, 
	      label="x₁ + x₂ = 3")
	
	plot!(p, [0, 0], [0, 3], 
	      linewidth=2.5, color=:green, 
	      label="x₁ = 0")
	plot!(p, [0, 3], [0, 0], 
	      linewidth=2.5, color=:orange, 
	      label="x₂ = 0")
	
	# Shade feasible region
	x1_grid = range(0, 3, length=150)
	x2_grid = range(0, 3, length=150)
	feasible_x1 = Float64[]
	feasible_x2 = Float64[]
	for x1 in x1_grid
	    for x2 in x2_grid
	        if x1 >= 0 && x2 >= 0 && x1^2 + x2^2 <= 5 && x1 + x2 <= 3
	            push!(feasible_x1, x1)
	            push!(feasible_x2, x2)
	        end
	    end
	end
	scatter!(p, feasible_x1, feasible_x2, 
	         markersize=1, 
	         markerstrokewidth=0,
	         color=:lightblue, 
	         alpha=0.9,
	         label="")
	
	# Mark center of objective
	scatter!(p, [3], [2], 
	         markersize=8, 
	         color=:purple, 
	         markershape=:star,
	         label="Center (3, 2)")
	
	# Plot selected point
	point_color = is_feasible(x1_val, x2_val) ? :green : :red
	scatter!(p, [x1_val], [x2_val], 
	         markersize=10, 
	         color=point_color, 
	         markershape=:circle,
	         markerstrokewidth=3,
	         markerstrokecolor=:black,
	         label="Selected point")
	
	# Plot gradient of objective function
	grad_f = ∇f(x1_val, x2_val)
	scale = 0.3  # scale factor for gradient arrows
	quiver!(p, [x1_val], [x2_val], 
	        quiver=([scale*grad_f[1]], [scale*grad_f[2]]), 
	        color=:black, 
	        linewidth=3,
	        arrow=:closed,
	        label="∇f")
	annotate!(p, 0.1+x1_val + scale*grad_f[1], x2_val + scale*grad_f[2], 
          text(L"\nabla f = \left[%$(round(grad_f[1], digits=2)), %$(round(grad_f[2], digits=2))\right]", 
               :black, 12, :left))
	# Plot gradients of active constraints
	active = active_constraints(x1_val, x2_val)
	colors = [:red, :blue, :green, :orange]
	
	for (i, name, grad) in active
		quiver!(p, [x1_val], [x2_val], 
		        quiver=([scale*grad[1]], [scale*grad[2]]), 
		        color=colors[i], 
		        linewidth=3,
		        arrow=:closed,
		        label="∇g$i")
		# Add annotation for each active constraint gradient
		annotate!(p, 0.1+x1_val + scale*grad[1], x2_val + scale*grad[2], 
		          text(L"\nabla g_%$i = [%$(round(grad[1], digits=2)), %$(round(grad[2], digits=2))]", 
		               colors[i], 9, :left))
	end
	
	xlims!(p, -0.5, 4.5)
	ylims!(p, -0.5, 3.5)
	
	p
end

# ╔═╡ 8f24c28e-57bd-42ae-9883-1cebc0715a0c
begin
	x1_val_e2_html = @bind x1_val_e2 NumberField(0.0:0.01:3.0, default=1.5)

	x2_val_e2_html = @bind x2_val_e2 NumberField(0.0:0.01:3.0, default=1.0)
	cm"""

	``x_1=`` $(x1_val_e2_html)
	
	``x_2=`` $(x2_val_e2_html)
	
	"""
end

# ╔═╡ 93a35971-c2eb-4dd3-bd66-b3335400861d
let
	# Define the objective function
	f(x1, x2) = (x1 - 1)^2 + (x2 - 1)^2
	
	# Gradient of objective function
	∇f(x1, x2) = [2*(x1 - 1), 2*(x2 - 1)]
	
	# Constraint functions (written as g(x) ≤ 0)
	g1(x1, x2) = x1 + x2 - 1  # cubic constraint
	g2(x1, x2) = -x1              # x1 ≥ 0
	g3(x1, x2) = -x2              # x2 ≥ 0
	
	# Gradients of constraints
	∇g1(x1, x2) = [1,1]
	∇g2(x1, x2) = [-1, 0]
	∇g3(x1, x2) = [0, -1]
	
	# Check which constraints are active (within tolerance)
	tol = 0.1
	function active_constraints(x1, x2)
		active = []
		if abs(g1(x1, x2)) < tol
			push!(active, (1, "x₁ + x₂ - 1 ≤ 0", ∇g1(x1, x2)))
		end
		if abs(g2(x1, x2)) < tol
			push!(active, (2, "x₁ ≥ 0", ∇g2(x1, x2)))
		end
		if abs(g3(x1, x2)) < tol
			push!(active, (3, "x₂ ≥ 0", ∇g3(x1, x2)))
		end
		return active
	end
	
	# Check if point is feasible
	function is_feasible(x1, x2)
		return g1(x1, x2) <= tol && g2(x1, x2) <= tol && g3(x1, x2) <= tol
	end

	# Create the plot
	x1_range = range(-0.2, 2.0, length=400)
	x2_range = range(-0.2, 2.0, length=400)
	
	p = plot(size=(800, 800), 
			 aspect_ratio=:equal, 
			 frame_style=:origin,
	         xlabel="x₁", ylabel="x₂", 
	         title="Gradients at ($(round(x1_val_e2, digits=2)), $(round(x2_val_e2, digits=2)))",
	         legend=:topright, legendfontsize=8)
	
	# Plot contour curves of objective
	contour!(p, x1_range, x2_range, 
	         (x1, x2) -> f(x1, x2), 
	         levels=15, 
	         color=:viridis, 
	         linewidth=1.5,
	         alpha=0.5,
	         colorbar=false,
	         label="")
	
	# Plot constraint boundaries
	# For (x1 + x2 - 1)³ = 0, we have x1 + x2 = 1
	x1_line = range(0, 1, length=100)
	x2_line = 1 .- x1_line
	plot!(p, x1_line, x2_line, 
	      linewidth=2.5, color=:red, 
	      label="x₁ + x₂ - 1 = 0")
	
	plot!(p, [0, 0], [0, 2], 
	      linewidth=2.5, color=:green, 
	      label="x₁ = 0")
	plot!(p, [0, 2], [0, 0], 
	      linewidth=2.5, color=:orange, 
	      label="x₂ = 0")
	
	# Shade feasible region
	# Feasible: (x1 + x2 - 1)³ ≤ 0 means x1 + x2 ≤ 1, and x1, x2 ≥ 0
	x1_grid = range(0, 1.5, length=150)
	x2_grid = range(0, 1.5, length=150)
	feasible_x1 = Float64[]
	feasible_x2 = Float64[]
	for x1 in x1_grid
	    for x2 in x2_grid
	        if x1 >= 0 && x2 >= 0 && (x1 + x2 - 1)^3 <= 0
	            push!(feasible_x1, x1)
	            push!(feasible_x2, x2)
	        end
	    end
	end
	scatter!(p, feasible_x1, feasible_x2, 
	         markersize=1, 
	         markerstrokewidth=0,
	         color=:lightblue, 
	         alpha=0.8,
	         label="")
	
	# Mark center of objective
	scatter!(p, [1], [1], 
	         markersize=8, 
	         color=:purple, 
	         markershape=:star,
	         label="Center (1, 1)")
	
	# Plot selected point
	point_color = is_feasible(x1_val_e2, x2_val_e2) ? :green : :red
	scatter!(p, [x1_val_e2], [x2_val_e2], 
	         markersize=10, 
	         color=point_color, 
	         markershape=:circle,
	         markerstrokewidth=3,
	         markerstrokecolor=:black,
	         label="Selected point")
	
	# Plot gradient of objective function
	grad_f = ∇f(x1_val_e2, x2_val_e2)
	scale = 0.2  # scale factor for gradient arrows
	quiver!(p, [x1_val_e2], [x2_val_e2], 
	        quiver=([scale*grad_f[1]], [scale*grad_f[2]]), 
	        color=:black, 
	        linewidth=3,
	        arrow=:closed,
	        label="∇f")
	
	# Add annotation for ∇f
	annotate!(p, x1_val_e2 + scale*grad_f[1], x2_val_e2 + scale*grad_f[2], 
	          text(L"\nabla f = [%$(round(grad_f[1], digits=2)), %$(round(grad_f[2], digits=2))]", 
	               :black, 9, :left))
	
	# Plot gradients of active constraints
	active = active_constraints(x1_val_e2, x2_val_e2)
	colors = [:red, :green, :orange]
	
	for (i, name, grad) in active
		quiver!(p, [x1_val_e2], [x2_val_e2], 
		        quiver=([scale*grad[1]], [scale*grad[2]]), 
		        color=colors[i], 
		        linewidth=3,
		        arrow=:closed,
		        label="∇g$i")
		
		# Add annotation for each active constraint gradient
		annotate!(p, x1_val_e2 + scale*grad[1], x2_val_e2 + scale*grad[2], 
		          text(L"\nabla g_%$i = [%$(round(grad[1], digits=2)), %$(round(grad[2], digits=2))]", 
		               colors[i], 9, :left))
	end
	
	xlims!(p, -0.2, 2.0)
	ylims!(p, -0.2, 2.0)
	
	p
end

# ╔═╡ dc71a49a-1e1f-474c-bfcd-cde361115ee1
md"### Fritz John Optimality Conditions "

# ╔═╡ 6d6428eb-9442-43a3-9cad-815b7ccfc203
begin
	x1_val_e3_html = @bind x1_val_e3 NumberField(0.0:0.01:3.0, default=1.5)

	x2_val_e3_html = @bind x2_val_e3 NumberField(0.0:0.01:3.0, default=1.0)
	cm"""

	``x_1=`` $(x1_val_e3_html)
	
	``x_2=`` $(x2_val_e3_html)
	
	"""
end

# ╔═╡ e64e8a77-6b32-46ab-b126-f77c917935ad
let
	# Define the objective function
	f(x1, x2) = (x1 - 3)^2 + (x2 - 2)^2
	
	# Gradient of objective function
	∇f(x1, x2) = [2*(x1 - 3), 2*(x2 - 2)]
	
	# Constraint functions (written as g(x) ≤ 0)
	g1(x1, x2) = x1^2 + x2^2 - 5       # circle constraint
	g2(x1, x2) = x1 + 2*x2 - 4         # line constraint
	g3(x1, x2) = -x1                   # x1 ≥ 0
	g4(x1, x2) = -x2                   # x2 ≥ 0
	
	# Gradients of constraints
	∇g1(x1, x2) = [2*x1, 2*x2]
	∇g2(x1, x2) = [1, 2]
	∇g3(x1, x2) = [-1, 0]
	∇g4(x1, x2) = [0, -1]
	
	# Check which constraints are active (within tolerance)
	tol = 0.1
	function active_constraints(x1, x2)
		active = []
		if abs(g1(x1, x2)) < tol
			push!(active, (1, "x₁² + x₂² ≤ 5", ∇g1(x1, x2)))
		end
		if abs(g2(x1, x2)) < tol
			push!(active, (2, "x₁ + 2x₂ ≤ 4", ∇g2(x1, x2)))
		end
		if abs(g3(x1, x2)) < tol
			push!(active, (3, "x₁ ≥ 0", ∇g3(x1, x2)))
		end
		if abs(g4(x1, x2)) < tol
			push!(active, (4, "x₂ ≥ 0", ∇g4(x1, x2)))
		end
		return active
	end
	
	# Check if point is feasible
	function is_feasible(x1, x2)
		return g1(x1, x2) <= tol && g2(x1, x2) <= tol && 
		       g3(x1, x2) <= tol && g4(x1, x2) <= tol
	end

	# Create the plot
	x1_range = range(-0.5, 3.5, length=400)
	x2_range = range(-0.5, 3.0, length=400)
	
	p = plot(size=(800, 800), aspect_ratio=:equal, 
	         xlabel=L"x_1", ylabel=L"x_2", 
	         title="Gradients at ($(round(x1_val_e3, digits=2)), $(round(x2_val_e3, digits=2)))",
	         legend=:topright, legendfontsize=8,
	         framestyle=:origin)
	
	# Plot contour curves of objective
	contour!(p, x1_range, x2_range, 
	         (x1, x2) -> f(x1, x2), 
	         levels=15, 
	         color=:viridis, 
	         linewidth=1.5,
	         alpha=0.5,
	         colorbar=false,
	         label="")
	
	# Plot constraint boundaries
	# Circle: x₁² + x₂² = 5
	θ = range(0, 2π, length=200)
	circle_x1 = sqrt(5) .* cos.(θ)
	circle_x2 = sqrt(5) .* sin.(θ)
	plot!(p, circle_x1, circle_x2, 
	      linewidth=2.5, color=:red, 
	      label=L"x_1^2 + x_2^2 = 5", linestyle=:dash)
	
	# Line: x₁ + 2x₂ = 4
	x1_line = range(-0.5, 4, length=100)
	x2_line = (4 .- x1_line) ./ 2
	plot!(p, x1_line, x2_line, 
	      linewidth=2.5, color=:blue, 
	      label=L"x_1 + 2x_2 = 4", linestyle=:dash)
	
	# Axes constraints (already shown by framestyle=:origin, but we can highlight the positive quadrant)
	plot!(p, [0, 0], [-0.5, 3], 
	      linewidth=2.5, color=:green, 
	      label=L"x_1 = 0", linestyle=:dash, alpha=0.5)
	plot!(p, [-0.5, 3.5], [0, 0], 
	      linewidth=2.5, color=:orange, 
	      label=L"x_2 = 0", linestyle=:dash, alpha=0.5)
	
	# Shade feasible region
	x1_grid = range(0, 3, length=150)
	x2_grid = range(0, 2.5, length=150)
	feasible_x1 = Float64[]
	feasible_x2 = Float64[]
	for x1 in x1_grid
	    for x2 in x2_grid
	        if x1 >= 0 && x2 >= 0 && x1^2 + x2^2 <= 5 && x1 + 2*x2 <= 4
	            push!(feasible_x1, x1)
	            push!(feasible_x2, x2)
	        end
	    end
	end
	scatter!(p, feasible_x1, feasible_x2, 
	         markersize=1, 
	         markerstrokewidth=0,
	         color=:lightblue, 
	         alpha=0.9,
	         label="")
	
	# Mark center of objective
	scatter!(p, [3], [2], 
	         markersize=8, 
	         color=:purple, 
	         markershape=:star,
	         label="Center (3, 2)")
	
	# Plot selected point
	point_color = is_feasible(x1_val_e3, x2_val_e3) ? :green : :red
	scatter!(p, [x1_val_e3], [x2_val_e3], 
	         markersize=10, 
	         color=point_color, 
	         markershape=:circle,
	         markerstrokewidth=3,
	         markerstrokecolor=:black,
	         label="Selected point")
	
	# Plot gradient of objective function
	grad_f = ∇f(x1_val_e3, x2_val_e3)
	scale = 0.3  # scale factor for gradient arrows
	quiver!(p, [x1_val_e3], [x2_val_e3], 
	        quiver=([scale*grad_f[1]], [scale*grad_f[2]]), 
	        color=:black, 
	        linewidth=3,
	        arrow=:closed,
	        label=L"\nabla f")
	
	# Add annotation for ∇f
	annotate!(p, x1_val_e3 + scale*grad_f[1] + 0.1, x2_val_e3 + scale*grad_f[2] + 0.1, 
	          text(L"\nabla f = [%$(round(grad_f[1], digits=2)), %$(round(grad_f[2], digits=2))]", 
	               :black, 8, :left))
	
	# Plot gradients of active constraints
	active = active_constraints(x1_val_e3, x2_val_e3)
	colors = [:red, :blue, :green, :orange]
	
	for (i, name, grad) in active
		quiver!(p, [x1_val_e3], [x2_val_e3], 
		        quiver=([scale*grad[1]], [scale*grad[2]]), 
		        color=colors[i], 
		        linewidth=3,
		        arrow=:closed,
		        label=L"\nabla g_%$i")
		
		# Add annotation for each active constraint gradient
		offset_x = 0.1
		offset_y = 0.1
		annotate!(p, x1_val_e3 + scale*grad[1] + offset_x, 
		          x2_val_e3 + scale*grad[2] + offset_y, 
		          text(L"\nabla g_%$i = [%$(round(grad[1], digits=2)), %$(round(grad[2], digits=2))]", 
		               colors[i], 8, :left))
	end
	
	xlims!(p, -0.5, 4.5)
	ylims!(p, -0.5, 3.0)
	
	p
end

# ╔═╡ a60d5096-9152-44e0-b5d2-3bb789dcff5d

begin
	x1_val_e4_html = @bind x1_val_e4 NumberField(-0.5:0.01:2.0, default=0.5)
	x2_val_e4_html = @bind x2_val_e4 NumberField(0.0:0.01:2.0, default=0.5)
	cm"""
	``x_1=`` $(x1_val_e4_html)
	
	``x_2=`` $(x2_val_e4_html)
	
	"""
end


# ╔═╡ b2cdd374-e06f-4d8f-9db1-6dfb7d661a40
let
	# Define the objective function
	f_e4(x1, x2) = -x1
	
	# Gradient of objective function
	∇f_e4(x1, x2) = [-1, 0]
	
	# Constraint functions (written as g(x) ≤ 0)
	g1_e4(x1, x2) = x2 - (1 - x1)^3     # cubic constraint
	g2_e4(x1, x2) = -x2                 # x2 ≥ 0
	
	# Gradients of constraints
	∇g1_e4(x1, x2) = [3*(1 - x1)^2, 1]
	∇g2_e4(x1, x2) = [0, -1]
	
	# Check which constraints are active (within tolerance)
	tol_e4 = 0.1
	function active_constraints_e4(x1, x2)
		active = []
		if abs(g1_e4(x1, x2)) < tol_e4
			push!(active, (1, "x₂ - (1-x₁)³ ≤ 0", ∇g1_e4(x1, x2)))
		end
		if abs(g2_e4(x1, x2)) < tol_e4
			push!(active, (2, "x₂ ≥ 0", ∇g2_e4(x1, x2)))
		end
		return active
	end
	
	# Check if point is feasible
	function is_feasible_e4(x1, x2)
		return g1_e4(x1, x2) <= tol_e4 && g2_e4(x1, x2) <= tol_e4
	end

	# Create the plot
	x1_range_e4 = range(-0.5, 2.0, length=400)
	x2_range_e4 = range(-0.2, 2.0, length=400)
	
	p_e4 = plot(size=(800, 800), aspect_ratio=:equal, 
	         xlabel=L"x_1", ylabel=L"x_2", 
	         title="Gradients at ($(round(x1_val_e4, digits=2)), $(round(x2_val_e4, digits=2)))",
	         legend=:topright, legendfontsize=8,
	         framestyle=:origin)
	
	# Plot contour curves of objective (vertical lines since f = -x₁)
	for x1_contour in -0.5:0.2:2.0
		plot!(p_e4, [x1_contour, x1_contour], [-0.2, 2.0],
		      color=:gray, alpha=0.3, linewidth=1, label="")
	end
	
	# Plot constraint boundaries
	# Curve: x₂ = (1 - x₁)³
	x1_curve = range(-0.5, 2.0, length=200)
	x2_curve = (1 .- x1_curve).^3
	plot!(p_e4, x1_curve, x2_curve, 
	      linewidth=2.5, color=:red, 
	      label=L"x_2 = (1-x_1)^3", linestyle=:dash)
	
	# Axis constraint x₂ = 0
	plot!(p_e4, [-0.5, 2.0], [0, 0], 
	      linewidth=2.5, color=:blue, 
	      label=L"x_2 = 0", linestyle=:dash, alpha=0.5)
	
	# Shade feasible region
	x1_grid_e4 = range(-0.5, 2.0, length=150)
	x2_grid_e4 = range(0, 2.0, length=150)
	feasible_x1_e4 = Float64[]
	feasible_x2_e4 = Float64[]
	for x1 in x1_grid_e4
	    for x2 in x2_grid_e4
	        if x2 >= 0 && x2 <= (1 - x1)^3
	            push!(feasible_x1_e4, x1)
	            push!(feasible_x2_e4, x2)
	        end
	    end
	end
	scatter!(p_e4, feasible_x1_e4, feasible_x2_e4, 
	         markersize=1, 
	         markerstrokewidth=0,
	         color=:lightblue, 
	         alpha=0.2,
	         label="")
	
	# Plot selected point
	point_color_e4 = is_feasible_e4(x1_val_e4, x2_val_e4) ? :green : :red
	scatter!(p_e4, [x1_val_e4], [x2_val_e4], 
	         markersize=10, 
	         color=point_color_e4, 
	         markershape=:circle,
	         markerstrokewidth=3,
	         markerstrokecolor=:black,
	         label="Selected point")
	
	# Plot gradient of objective function
	grad_f_e4 = ∇f_e4(x1_val_e4, x2_val_e4)
	scale_e4 = 0.4  # scale factor for gradient arrows
	quiver!(p_e4, [x1_val_e4], [x2_val_e4], 
	        quiver=([scale_e4*grad_f_e4[1]], [scale_e4*grad_f_e4[2]]), 
	        color=:black, 
	        linewidth=3,
	        arrow=:closed,
	        label=L"\nabla f")
	
	# Add annotation for ∇f
	annotate!(p_e4, x1_val_e4 + scale_e4*grad_f_e4[1] + 0.1, 
	          x2_val_e4 + scale_e4*grad_f_e4[2] + 0.1, 
	          text(L"\nabla f = [%$(round(grad_f_e4[1], digits=2)), %$(round(grad_f_e4[2], digits=2))]", 
	               :black, 8, :left))
	
	# Plot gradients of active constraints
	active_e4 = active_constraints_e4(x1_val_e4, x2_val_e4)
	colors_e4 = [:red, :blue]
	
	for (i, name, grad) in active_e4
		quiver!(p_e4, [x1_val_e4], [x2_val_e4], 
		        quiver=([scale_e4*grad[1]], [scale_e4*grad[2]]), 
		        color=colors_e4[i], 
		        linewidth=3,
		        arrow=:closed,
		        label=L"\nabla g_%$i")
		
		# Add annotation for each active constraint gradient
		offset_x = 0.1
		offset_y = 0.1
		annotate!(p_e4, x1_val_e4 + scale_e4*grad[1] + offset_x, 
		          x2_val_e4 + scale_e4*grad[2] + offset_y, 
		          text(L"\nabla g_%$i = [%$(round(grad[1], digits=2)), %$(round(grad[2], digits=2))]", 
		               colors_e4[i], 8, :left))
	end
	
	xlims!(p_e4, -0.5, 2.0)
	ylims!(p_e4, -0.2, 2.0)
	
	p_e4
end

# ╔═╡ 9197b26f-59a7-4508-82f7-8441bea7a2e1
md"### Karush-Kuhn-Tucker Conditions "

# ╔═╡ f1b16a2c-816a-4a31-aff2-f8f10ca45cc9
cm"""
> If we require that ``u_0 > 0``, in FJ conditions, we obtain the so called __KKT__ conditions.
"""

# ╔═╡ 4699162f-ef6a-4279-9aa7-56170fa9a6ff
begin
	x1_val_e5_html = @bind x1_val_e5 NumberField(-0.5:0.01:2.0, default=0.5)
	x2_val_e5_html = @bind x2_val_e5 NumberField(0.0:0.01:2.0, default=0.5)
	cm"""
	``x_1=`` $(x1_val_e5_html)
	
	``x_2=`` $(x2_val_e5_html)
	
	"""
end



# ╔═╡ 7ed6becd-1e79-4103-950a-017d187c585c
begin
	# Define the objective function
	f_e5(x1, x2) = -x1
	
	# Gradient of objective function
	∇f_e5(x1, x2) = [-1, 0]
	
	# Constraint functions (written as g(x) ≤ 0)
	g1_e5(x1, x2) = x1 + x2 - 1     # linear constraint
	g2_e5(x1, x2) = -x2              # x2 ≥ 0
	
	# Gradients of constraints
	∇g1_e5(x1, x2) = [1, 1]
	∇g2_e5(x1, x2) = [0, -1]
	
	# Check which constraints are active (within tolerance)
	tol_e5 = 0.1
	function active_constraints_e5(x1, x2)
		active = []
		if abs(g1_e5(x1, x2)) < tol_e5
			push!(active, (1, "x₁ + x₂ - 1 ≤ 0", ∇g1_e5(x1, x2)))
		end
		if abs(g2_e5(x1, x2)) < tol_e5
			push!(active, (2, "x₂ ≥ 0", ∇g2_e5(x1, x2)))
		end
		return active
	end
	
	# Check if point is feasible
	function is_feasible_e5(x1, x2)
		return g1_e5(x1, x2) <= tol_e5 && g2_e5(x1, x2) <= tol_e5
	end

	# Create the plot
	x1_range_e5 = range(-0.5, 2.0, length=400)
	x2_range_e5 = range(-0.2, 2.0, length=400)
	
	p_e5 = plot(size=(800, 800), aspect_ratio=:equal, 
	         xlabel=L"x_1", ylabel=L"x_2", 
	         title="Gradients at ($(round(x1_val_e5, digits=2)), $(round(x2_val_e5, digits=2)))",
	         legend=:topright, legendfontsize=8,
	         framestyle=:origin)
	
	# Plot contour curves of objective (vertical lines since f = -x₁)
	for x1_contour in -0.5:0.2:2.0
		plot!(p_e5, [x1_contour, x1_contour], [-0.2, 2.0],
		      color=:gray, alpha=0.3, linewidth=1, label="")
	end
	
	# Plot constraint boundaries
	# Line: x₁ + x₂ = 1
	x1_line_e5 = range(-0.5, 2.0, length=100)
	x2_line_e5 = 1 .- x1_line_e5
	plot!(p_e5, x1_line_e5, x2_line_e5, 
	      linewidth=2.5, color=:red, 
	      label=L"x_1 + x_2 = 1", linestyle=:dash)
	
	# Axis constraint x₂ = 0
	plot!(p_e5, [-0.5, 2.0], [0, 0], 
	      linewidth=2.5, color=:blue, 
	      label=L"x_2 = 0", linestyle=:dash, alpha=0.5)
	
	# Shade feasible region
	x1_grid_e5 = range(-0.5, 2.0, length=150)
	x2_grid_e5 = range(0, 2.0, length=150)
	feasible_x1_e5 = Float64[]
	feasible_x2_e5 = Float64[]
	for x1 in x1_grid_e5
	    for x2 in x2_grid_e5
	        if x2 >= 0 && x1 + x2 <= 1
	            push!(feasible_x1_e5, x1)
	            push!(feasible_x2_e5, x2)
	        end
	    end
	end
	scatter!(p_e5, feasible_x1_e5, feasible_x2_e5, 
	         markersize=1, 
	         markerstrokewidth=0,
	         color=:lightblue, 
	         alpha=0.2,
	         label="")
	
	# Plot selected point
	point_color_e5 = is_feasible_e5(x1_val_e5, x2_val_e5) ? :green : :red
	scatter!(p_e5, [x1_val_e5], [x2_val_e5], 
	         markersize=10, 
	         color=point_color_e5, 
	         markershape=:circle,
	         markerstrokewidth=3,
	         markerstrokecolor=:black,
	         label="Selected point")
	
	# Plot gradient of objective function
	grad_f_e5 = ∇f_e5(x1_val_e5, x2_val_e5)
	scale_e5 = 0.4  # scale factor for gradient arrows
	quiver!(p_e5, [x1_val_e5], [x2_val_e5], 
	        quiver=([scale_e5*grad_f_e5[1]], [scale_e5*grad_f_e5[2]]), 
	        color=:black, 
	        linewidth=3,
	        arrow=:closed,
	        label=L"\nabla f")
	
	# Add annotation for ∇f
	annotate!(p_e5, x1_val_e5 + scale_e5*grad_f_e5[1] - 0.15, 
	          x2_val_e5 + scale_e5*grad_f_e5[2] + 0.15, 
	          text(L"\nabla f = [%$(round(grad_f_e5[1], digits=2)), %$(round(grad_f_e5[2], digits=2))]", 
	               :black, 8, :left))
	
	# Plot gradients of active constraints
	active_e5 = active_constraints_e5(x1_val_e5, x2_val_e5)
	colors_e5 = [:red, :blue]
	
	for (i, name, grad) in active_e5
		quiver!(p_e5, [x1_val_e5], [x2_val_e5], 
		        quiver=([scale_e5*grad[1]], [scale_e5*grad[2]]), 
		        color=colors_e5[i], 
		        linewidth=3,
		        arrow=:closed,
		        label=L"\nabla g_%$i")
		
		# Add annotation for each active constraint gradient
		offset_x = 0.1
		offset_y = 0.1
		annotate!(p_e5, x1_val_e5 + scale_e5*grad[1] + offset_x, 
		          x2_val_e5 + scale_e5*grad[2] + offset_y, 
		          text(L"\nabla g_%$i = [%$(round(grad[1], digits=2)), %$(round(grad[2], digits=2))]", 
		               colors_e5[i], 8, :left))
	end
	
	xlims!(p_e5, -0.5, 2.0)
	ylims!(p_e5, -1, 2.0)
	
	p_e5
end

# ╔═╡ b06c5711-3f32-4514-a44d-e3aef2e69125
cm"""
__Information about the selected point__

**Point:** ($(round(x1_val_e5, digits=3)), $(round(x2_val_e5, digits=3)))

**Feasible:** $(is_feasible_e5(x1_val_e5, x2_val_e5) ? "✓ Yes" : "✗ No")

**Objective value:** f = $(round(f_e5(x1_val_e5, x2_val_e5), digits=3))

**Gradient of f:** ∇f = [$(round(∇f_e5(x1_val_e5, x2_val_e5)[1], digits=3)), $(round(∇f_e5(x1_val_e5, x2_val_e5)[2], digits=3))]

**Constraint values:**
- g₁ (line): $(round(g1_e5(x1_val_e5, x2_val_e5), digits=3)) $(g1_e5(x1_val_e5, x2_val_e5) <= 0 ? "✓" : "✗")
- g₂ (x₂≥0): $(round(g2_e5(x1_val_e5, x2_val_e5), digits=3)) $(g2_e5(x1_val_e5, x2_val_e5) <= 0 ? "✓" : "✗")

**Active constraints:** $(length(active_constraints_e5(x1_val_e5, x2_val_e5)) > 0 ? join([name for (i, name, grad) in active_constraints_e5(x1_val_e5, x2_val_e5)], ", ") : "None")
"""


# ╔═╡ 2d7d967e-bd11-4959-978e-2325f1b78f95
md"""
__KKT Conditions Check__

For a point to be optimal, it must satisfy the KKT conditions:
1. **Stationarity:** ∇f + Σλᵢ∇gᵢ = 0 (for active constraints)
2. **Primal feasibility:** All constraints satisfied
3. **Dual feasibility:** λᵢ ≥ 0
4. **Complementary slackness:** λᵢ·gᵢ = 0

**Key Observations:**
- The objective f = -x₁ is minimized by making x₁ as large as possible
- The gradient ∇f = [-1, 0] always points left (direction of steepest decrease)
- The feasible region is a triangle bounded by x₁ + x₂ ≤ 1 and x₂ ≥ 0

**Try these interesting points:**
- **(1.0, 0.0)** - The optimal point! Both constraints are active here
- **(0.5, 0.5)** - On the line x₁ + x₂ = 1
- **(0.0, 0.0)** - At the origin
- **(0.0, 1.0)** - On the line, at the other endpoint

**At the optimal point (1, 0):**
- Both constraints are active: g₁ = 0 and g₂ = 0
- ∇f = [-1, 0]
- ∇g₁ = [1, 1]
- ∇g₂ = [0, -1]
- The KKT stationarity condition: ∇f + λ₁∇g₁ + λ₂∇g₂ = 0
- This gives: [-1, 0] + λ₁[1, 1] + λ₂[0, -1] = [0, 0]
- Solution: λ₁ = 1, λ₂ = -1
- But λ₂ < 0 violates dual feasibility! This means **the KKT conditions fail at (1,0)** even though it's the optimal point.
- However, if we only consider g₁ as active, we get: [-1, 0] + λ₁[1, 1] = [0, 0], which has no solution with λ₁ ≥ 0.
- The actual optimal point considering both constraints properly is still (1, 0), and this illustrates how constraint qualifications matter!
"""

# ╔═╡ cefc755f-2302-4213-8836-04ff8e33305d
md"### Geometric Interpretation of the KKT Conditions: Linear Programming Approximations "

# ╔═╡ 2e8796d7-06b4-467e-b295-1019d34859bf
begin
	x1_val_e12_html = @bind x1_val_e12 NumberField(0.0:0.01:2.5, default=1.0)
	x2_val_e12_html = @bind x2_val_e12 NumberField(-2.0:0.01:2.0, default=0.0)
	cm"""
	**Select a point to analyze:**
	
	``x_1=`` $(x1_val_e12_html)
	
	``x_2=`` $(x2_val_e12_html)
	"""
end

# ╠═╡ functions

# ╔═╡ 9389bae8-5492-4402-ac8a-fccadfd8351e
let
	# Objective function: f(x) = x₁
	f(x1, x2) = x1
	
	# Gradient of objective: ∇f = [1, 0]
	∇f(x1, x2) = [1.0, 0.0]
	
	# Constraints in form g(x) ≤ 0
	# g₁: (x₁-1)² + (x₂-1)² - 1 ≤ 0
	g1(x1, x2) = (x1 - 1)^2 + (x2 - 1)^2 - 1
	
	# g₂: (x₁-1)² + (x₂+1)² - 1 ≤ 0
	g2(x1, x2) = (x1 - 1)^2 + (x2 + 1)^2 - 1
	
	# Gradients of constraints
	∇g1(x1, x2) = [2*(x1 - 1), 2*(x2 - 1)]
	∇g2(x1, x2) = [2*(x1 - 1), 2*(x2 + 1)]
	
	# Tolerance for detecting active constraints
	tol = 0.1
	
	# Check which constraints are active
	function active_constraints(x1, x2)
		active = []
		if abs(g1(x1, x2)) < tol
			push!(active, (1, "(x₁-1)² + (x₂-1)² ≤ 1", ∇g1(x1, x2)))
		end
		if abs(g2(x1, x2)) < tol
			push!(active, (2, "(x₁-1)² + (x₂+1)² ≤ 1", ∇g2(x1, x2)))
		end
		return active
	end
	
	# Check feasibility
	function is_feasible(x1, x2)
		return g1(x1, x2) <= tol && g2(x1, x2) <= tol
	end

	# Setup plot ranges
	x1_range = range(-0.5, 2.5, length=400)
	x2_range = range(-2.5, 2.5, length=400)
	
	p = plot(size=(800, 800), aspect_ratio=:equal, 
	         xlabel=L"x_1", ylabel=L"x_2", 
	         title="Gradients at ($(round(x1_val_e12, digits=2)), $(round(x2_val_e12, digits=2)))",
	         legend=:topright, legendfontsize=8,
	         framestyle=:origin)
	
	# Plot contour lines of objective f = x₁
	# For linear objective, draw vertical lines (iso-cost lines)
	for x1_contour in 0.0:0.2:2.0
		plot!(p, [x1_contour, x1_contour], [-2.5, 2.5],
		      color=:gray, alpha=0.3, linewidth=1, label="")
	end
	
	# Plot first circle boundary: (x₁-1)² + (x₂-1)² = 1
	θ1 = range(0, 2π, length=200)
	circle1_x1 = 1.0 .+ cos.(θ1)
	circle1_x2 = 1.0 .+ sin.(θ1)
	plot!(p, circle1_x1, circle1_x2, 
	      linewidth=2.5, color=:red, 
	      label=L"(x_1-1)^2 + (x_2-1)^2 = 1", linestyle=:dash)
	
	# Plot second circle boundary: (x₁-1)² + (x₂+1)² = 1
	θ2 = range(0, 2π, length=200)
	circle2_x1 = 1.0 .+ cos.(θ2)
	circle2_x2 = -1.0 .+ sin.(θ2)
	plot!(p, circle2_x1, circle2_x2, 
	      linewidth=2.5, color=:blue, 
	      label=L"(x_1-1)^2 + (x_2+1)^2 = 1", linestyle=:dash)
	
	# Shade feasible region (intersection of both circles)
	x1_grid = range(-0.5, 2.5, length=150)
	x2_grid = range(-2.5, 2.5, length=150)
	feasible_x1 = Float64[]
	feasible_x2 = Float64[]
	
	for x1 in x1_grid
	    for x2 in x2_grid
	        if g1(x1, x2) <= 0 && g2(x1, x2) <= 0
	            push!(feasible_x1, x1)
	            push!(feasible_x2, x2)
	        end
	    end
	end
	
	scatter!(p, feasible_x1, feasible_x2, 
	         markersize=1, 
	         markerstrokewidth=0,
	         color=:lightblue, 
	         alpha=0.3,
	         label="Feasible region")
	
	# Mark the centers of the circles
	scatter!(p, [1], [1], 
	         markersize=6, 
	         color=:red, 
	         markershape=:x,
	         label="Center (1,1)")
	
	scatter!(p, [1], [-1], 
	         markersize=6, 
	         color=:blue, 
	         markershape=:x,
	         label="Center (1,-1)")
	
	# Mark the optimal point (1, 0)
	scatter!(p, [1], [0], 
	         markersize=8, 
	         color=:gold, 
	         markershape=:star5,
	         markerstrokewidth=2,
	         markerstrokecolor=:black,
	         label="Optimal (1,0)")
	
	# Plot selected point
	point_color = is_feasible(x1_val_e12, x2_val_e12) ? :green : :red
	scatter!(p, [x1_val_e12], [x2_val_e12], 
	         markersize=10, 
	         color=point_color, 
	         markershape=:circle,
	         markerstrokewidth=3,
	         markerstrokecolor=:black,
	         label="Selected point")
	
	# Plot gradient of objective function
	grad_f = ∇f(x1_val_e12, x2_val_e12)
	scale = 0.4
	
	quiver!(p, [x1_val_e12], [x2_val_e12], 
	        quiver=([scale*grad_f[1]], [scale*grad_f[2]]), 
	        color=:black, 
	        linewidth=3,
	        arrow=:closed,
	        label=L"\nabla f")
	
	annotate!(p, x1_val_e12 + scale*grad_f[1] + 0.15, 
	          x2_val_e12 + scale*grad_f[2] + 0.15, 
	          text(L"\nabla f = [1, 0]", 
	               :black, 8, :left))
	
	# Plot gradients of active constraints
	active = active_constraints(x1_val_e12, x2_val_e12)
	colors = [:red, :blue, :green, :orange]
	
	for (i, name, grad) in active
		# Normalize gradient for visualization if needed
		grad_norm = norm(grad)
		if grad_norm > 0.01  # Avoid division by very small numbers
			quiver!(p, [x1_val_e12], [x2_val_e12], 
			        quiver=([scale*grad[1]], [scale*grad[2]]), 
			        color=colors[i], 
			        linewidth=3,
			        arrow=:closed,
			        label=L"\nabla g_%$i")
			
			# Add annotation
			offset_x = 0.15
			offset_y = 0.15
			annotate!(p, x1_val_e12 + scale*grad[1] + offset_x, 
			          x2_val_e12 + scale*grad[2] + offset_y, 
			          text(L"\nabla g_%$i = [%$(round(grad[1], digits=2)), %$(round(grad[2], digits=2))]", 
			               colors[i], 8, :left))
		end
	end
	
	xlims!(p, -0.5, 3.5)
	ylims!(p, -2.5, 2.5)

cm"""
$p
	
**Point:** ($(round(x1_val_e12, digits=3)), $(round(x2_val_e12, digits=3)))

**Feasible:** $(is_feasible(x1_val_e12, x2_val_e12) ? "✓ Yes" : "✗ No")

**Objective value:** f = $(round(f(x1_val_e12, x2_val_e12), digits=3))

**Gradient of f:** ∇f = [1, 0] (constant, points in direction of increasing x₁)

**Constraint values:**
- g₁ (upper circle): $(round(g1(x1_val_e12, x2_val_e12), digits=3)) $(g1(x1_val_e12, x2_val_e12) <= 0 ? "✓" : "✗")
- g₂ (lower circle): $(round(g2(x1_val_e12, x2_val_e12), digits=3)) $(g2(x1_val_e12, x2_val_e12) <= 0 ? "✓" : "✗")

**Active constraints:** $(length(active_constraints(x1_val_e12, x2_val_e12)) > 0 ? join([name for (i, name, grad) in active_constraints(x1_val_e12, x2_val_e12)], ", ") : "None")

---

### Gradient Information

**∇g₁ at this point:** [$(round(∇g1(x1_val_e12, x2_val_e12)[1], digits=3)), $(round(∇g1(x1_val_e12, x2_val_e12)[2], digits=3))]

**∇g₂ at this point:** [$(round(∇g2(x1_val_e12, x2_val_e12)[1], digits=3)), $(round(∇g2(x1_val_e12, x2_val_e12)[2], digits=3))]

"""
end

# ╠═╡ info

# ╔═╡ c6554735-c50a-4c61-a9ef-95dd3956b99e
md"## 4.3 Problems Having Inequality and Equality Constraints"

# ╔═╡ e319bc66-7810-4673-bc77-c3e0d8f1eeac
md"### Fritz John Conditions"

# ╔═╡ ca7bd385-5269-494e-a2b4-7a0b779d4389
md"### Karush-Kuhn-Tucker Conditions "

# ╔═╡ ecf678cd-8702-420c-9c6c-d7010ffe42f1
begin
	x431_val_html = @bind x431_val NumberField(0.0:0.01:2.5, default=1.0)
	x432_val_html = @bind x432_val NumberField(0.0:0.01:2.5, default=1.5)
	cm"""
	__Select Point__
	
	``x_1=`` $(x431_val_html)
	
	``x_2=`` $(x432_val_html)
	"""
end

# ╔═╡ 237b661f-d295-4c50-ae22-dd3441881cc1
let
	# Objective function
	f(x1, x2) = (x1 - 3)^2 + (x2 - 2)^2
	
	# Gradient of objective
	∇f(x1, x2) = [2*(x1 - 3), 2*(x2 - 2)]
	
	# Constraints as g(x) ≤ 0
	g1(x1, x2) = x1^2 + x2^2 - 5           # x₁² + x₂² ≤ 5
	g2(x1, x2) = -x1                        # x₁ ≥ 0
	g3(x1, x2) = -x2                        # x₂ ≥ 0
	g4(x1, x2) = x1 + 2*x2 - 4             # x₁ + 2x₂ = 4 (equality)
	g5(x1, x2) = -x1 - 2*x2 + 4            # x₁ + 2x₂ = 4 (other side)
	
	# Gradients of constraints
	∇g1(x1, x2) = [2*x1, 2*x2]
	∇g2(x1, x2) = [-1, 0]
	∇g3(x1, x2) = [0, -1]
	∇g4(x1, x2) = [1, 2]
	∇g5(x1, x2) = [-1, -2]
	
	# Tolerance for active constraints
	tol = 0.1
	
	# Check active constraints
	function active_constraints(x1, x2)
		active = []
		if abs(g1(x1, x2)) < tol
			push!(active, (1, "x₁² + x₂² ≤ 5", ∇g1(x1, x2), :red))
		end
		if abs(g2(x1, x2)) < tol
			push!(active, (2, "x₁ ≥ 0", ∇g2(x1, x2), :orange))
		end
		if abs(g3(x1, x2)) < tol
			push!(active, (3, "x₂ ≥ 0", ∇g3(x1, x2), :purple))
		end
		# Equality constraint is always active
		if abs(g4(x1, x2)) < tol
			push!(active, (4, "x₁ + 2x₂ = 4", ∇g4(x1, x2), :blue))
		end
		return active
	end
	
	# Check feasibility
	function is_feasible(x1, x2)
		return g1(x1, x2) <= tol && 
		       g2(x1, x2) <= tol && 
		       g3(x1, x2) <= tol && 
		       abs(g4(x1, x2)) <= tol  # equality constraint
	end

	x1_range = range(-0.5, 2.5, length=400)
	x2_range = range(-0.5, 2.5, length=400)
	
	p = plot(size=(800, 800), aspect_ratio=:equal, 
	         xlabel=L"x_1", ylabel=L"x_2", 
	         title="Gradients at ($(round(x431_val, digits=2)), $(round(x432_val, digits=2)))",
	         legend=:topright, legendfontsize=8,
	         framestyle=:origin,
	         xlims=(-0.5, 2.5), ylims=(-0.5, 2.5))
	
	# Plot contours of objective function
	contour!(p, x1_range, x2_range, 
	         (x1, x2) -> f(x1, x2), 
	         levels=15, 
	         color=:viridis, 
	         linewidth=1.5,
	         alpha=0.5,
	         colorbar=false,
	         label="")
	
	# Plot constraint boundaries
	
	# 1. Circle constraint: x₁² + x₂² = 5
	θ = range(0, 2π, length=200)
	circle_x1 = sqrt(5) .* cos.(θ)
	circle_x2 = sqrt(5) .* sin.(θ)
	plot!(p, circle_x1, circle_x2, 
	      linewidth=2.5, color=:red, 
	      label=L"x_1^2 + x_2^2 = 5", linestyle=:dash)
	
	# 2. x₁ = 0 (vertical line)
	plot!(p, [0, 0], [-0.5, 2.5], 
	      linewidth=2.5, color=:orange, 
	      label=L"x_1 = 0", linestyle=:dash)
	
	# 3. x₂ = 0 (horizontal line)
	plot!(p, [-0.5, 2.5], [0, 0], 
	      linewidth=2.5, color=:purple, 
	      label=L"x_2 = 0", linestyle=:dash)
	
	# 4. Equality constraint: x₁ + 2x₂ = 4
	line_x1 = range(-0.5, 2.5, length=100)
	line_x2 = (4 .- line_x1) ./ 2
	plot!(p, line_x1, line_x2, 
	      linewidth=3, color=:blue, 
	      label=L"x_1 + 2x_2 = 4", linestyle=:solid)
	
	# Shade feasible region
	# Find the feasible segment along x₁ + 2x₂ = 4
	x1_line = range(0, 2.5, length=1000)
	feasible_x1 = Float64[]
	feasible_x2 = Float64[]
	
	for x1 in x1_line
		x2 = (4 - x1) / 2
		if g1(x1, x2) <= 0 && g2(x1, x2) <= 0 && g3(x1, x2) <= 0
			push!(feasible_x1, x1)
			push!(feasible_x2, x2)
		end
	end
	
	# Create a polygon ribbon around the feasible segment for shading
	if !isempty(feasible_x1)
		ribbon_width = 0.04
		upper_x1 = feasible_x1 .+ ribbon_width .* (-1/sqrt(5))
		upper_x2 = feasible_x2 .+ ribbon_width .* (2/sqrt(5))
		lower_x1 = feasible_x1 .- ribbon_width .* (-1/sqrt(5))
		lower_x2 = feasible_x2 .- ribbon_width .* (2/sqrt(5))
		
		# Create closed polygon
		poly_x1 = vcat(upper_x1, reverse(lower_x1))
		poly_x2 = vcat(upper_x2, reverse(lower_x2))
		
		plot!(p, poly_x1, poly_x2, 
		      fillrange=0, fillalpha=0.3, fillcolor=:lightgreen,
		      linewidth=0, label="Feasible Region",
		      seriestype=:shape)
	end
	
	# Plot selected point
	point_color = is_feasible(x431_val, x432_val) ? :green : :red
	scatter!(p, [x431_val], [x432_val], 
	         markersize=10, color=point_color, 
	         markerstrokewidth=2, markerstrokecolor=:black,
	         label="Selected Point")
	
	# Plot gradient of objective function (∇f)
	grad_f = ∇f(x431_val, x432_val)
	scale = 0.35
	quiver!(p, [x431_val], [x432_val], 
	        quiver=([scale*grad_f[1]], [scale*grad_f[2]]),
	        color=:black, linewidth=2.5, arrow=:closed,
	        label=L"\nabla f")
	
	# Add annotation for ∇f
	annotate!(p, x431_val + scale*grad_f[1] + 0.12, 
	          x432_val + scale*grad_f[2] + 0.12, 
	          text(L"\nabla f = [%$(round(grad_f[1], digits=2)), %$(round(grad_f[2], digits=2))]", 
	               :black, 8, :left))
	
	# Plot gradients of active constraints
	active = active_constraints(x431_val, x432_val)
	
	for (idx, name, grad, color) in active
		quiver!(p, [x431_val], [x432_val], 
		        quiver=([scale*grad[1]], [scale*grad[2]]),
		        color=color, linewidth=2.5, arrow=:closed,
		        label=L"\nabla g_%$idx")
		
		# Add annotation
		offset_x = 0.12 * (grad[1] >= 0 ? 1 : -1)
		offset_y = 0.12 * (grad[2] >= 0 ? 1 : -1)
		annotate!(p, x431_val + scale*grad[1] + offset_x, 
		          x432_val + scale*grad[2] + offset_y, 
		          text(L"\nabla g_%$idx = [%$(round(grad[1], digits=2)), %$(round(grad[2], digits=2))]", 
		               color, 8, :left))
	end
	
	p
end

# ╔═╡ 20b7851e-1cd3-4a43-9e41-80ab4ad38ccc
let
	∇f=[-2.0;-2.0]
	∇g = [4.0;2.0]
	∇h = [1.0;2.0]
	∇f+(1/3)∇g+(2/3)∇h
end

# ╔═╡ 54c9ffd1-97be-44cd-8b78-07ec3f6f6299
md"## 4.4 Second-Order Necessary and Sufficient Conditions for Constrained Problems"

# ╔═╡ ebc2e75c-21e1-452c-9e09-da88cdc488a6
md"# Chapter 5:Constraint Qualifications"

# ╔═╡ 1c78770b-4814-48ee-837a-1c0f2c99a7c8
md"## 5.1 Cone of Tangents "

# ╔═╡ bdf534ab-e2a7-456c-83dd-23fe5cb48028
md"### Abadie Constraint Qualification"

# ╔═╡ 6db4f42f-6c20-44db-a0e6-c0b061781538
cm"""
```math
T=G^{\prime}.
```
"""

# ╔═╡ af44bbd1-b9a6-4668-99a2-794fad3f6a42
md"### Linearly Constrained Problems"

# ╔═╡ 42f6c9db-97d9-4852-a4c3-f7bbcb055a0f
begin
    struct LocalImage
        filename
    end

    function Base.show(io::IO, ::MIME"image/png", w::LocalImage)
        write(io, read(w.filename))
    end
end

# ╔═╡ fc877247-39bc-4bb0-8bda-1466fcb00798
@htl("""<style>
@import url("https://mmogib.github.io/math102/custom.css");
ul {
  list-style: disc; /* or whatever default you prefer: disc, circle, square */
  padding-left: 40px; /* adjust spacing */
}

ul li:before {
  content: none; /* removes the custom emoji */
}
</style>""")

# ╔═╡ fdd3c3e3-5089-456f-adef-7ab2e311331f
begin
	function min_latex_gi()
	cm"""
	```math
	 \begin{array}{lll}
	 \min & f(x) \\
	 \text{subject to}\\
	 & g_i(x)\leq 0 & i=1,\cdots,m\\
	 & x \in X
	 \end{array}
	 ```
	"""
	end

	function eql_latex_gi()
	cm"""
	```math
	 \begin{array}{lll}
	 \min & f(x) \\
	 \text{subject to}\\
	 & g_i(x)\leq 0 & i=1,\cdots,m\\
	 & h_i(x)= 0 & i=1,\cdots,l\\
	 & x \in X
	 \end{array}
	 ```
	"""
	end
	function min_latex()
	cm"""
	```math
	\begin{array}{lll}
	\min & f(x) \\
	\text{subject to}\\
	& x \in S\\
	\end{array}
	```
	"""
	end
	function indent(txt::String, ch::Int)
	"""
	<div style="padding-left:$(ch)ch;">

		$(txt)

	</div>
	"""
	end
    function add_space(n=1)
        repeat("&nbsp;", n)
    end
    function post_img(img::String, w=500)
        res = Resource(img, :width => w)
        cm"""
      <div class="img-container">

      $(res)

      </div>"""
    end
    function poolcode()
        cm"""
      <div class="img-container">

      $(Resource("https://www.dropbox.com/s/cat9ots4ausfzyc/qrcode_itempool.com_kfupm.png?raw=1",:width=>300))

      </div>"""
    end
    function define(t="")
        beginBlock("Definition", t)
    end
    function remark(t="")
        beginBlock("Remark", t)
    end
    function remarks(t="")
        beginBlock("Remarks", t)
    end
    function bbl(t)
        beginBlock(t, "")
    end
    function bbl(t, s)
        beginBlock(t, s)
    end
    ebl() = endBlock()
    function theorem(s)
        bth(s)
    end
    function bth(s)
        beginTheorem(s)
    end
    eth() = endTheorem()
    ex(n::Int; s::String="") = ex("Example $n", s)
    ex(t::Int, s::String) = example("Example $t", s)
    ex(t, s) = example(t, s)
    function beginBlock(title, subtitle)
        """<div style="box-sizing: border-box;">
       	<div style="display: flex;flex-direction: column;border: 6px solid rgba(200,200,200,0.5);box-sizing: border-box;">
       	<div style="display: flex;">
       	<div style="background-color: #FF9733;
       	    border-left: 10px solid #df7300;
       	    padding: 5px 10px;
       	    color: #fff!important;
       	    clear: left;
       	    margin-left: 0;font-size: 112%;
       	    line-height: 1.3;
       	    font-weight: 600;">$title</div>  <div style="olor: #000!important;
       	    margin: 0 0 20px 25px;
       	    float: none;
       	    clear: none;
       	    padding: 5px 0 0 0;
       	    margin: 0 0 0 20px;
       	    background-color: transparent;
       	    border: 0;
       	    overflow: hidden;
       	    min-width: 100px;font-weight: 600;
       	    line-height: 1.5;">$subtitle</div>
       	</div>
       	<p style="padding:5px;">
       """
    end
    function beginTheorem(subtitle)
        beginBlock("Theorem", subtitle)
    end
    function endBlock()
        """</p></div></div>"""
    end
    function endTheorem()
        endBlock()
    end
    ex() = example("Example", "")
    # function example(lable, desc)
    #     """<div style="display:flex;">
    #    <div style="
    #    font-size: 112%;
    #        line-height: 1.3;
    #        font-weight: 600;
    #        color: #f9ce4e;
    #        float: left;
    #        background-color: #5c5c5c;
    #        border-left: 10px solid #474546;
    #        padding: 5px 10px;
    #        margin: 0 12px 20px 0;
    #        border-radius: 0;
    #    ">$lable:</div>
    #    <div style="flex-grow:3;
    #    line-height: 1.3;
    #        font-weight: 600;
    #        float: left;
    #        padding: 5px 10px;
    #        margin: 0 12px 20px 0;
    #        border-radius: 0;
    #    ">$desc</div>
    #    </div>"""
    # end
	function example(lable, desc)
        """<div class="example-box">
    <div class="example-header">
      $lable
    </div>
    <div class="example-title">
      $desc
    </div>
    <div class="example-content">
      
  </div>
		"""
    end
	 
    @htl("")
end


# ╔═╡ b2805e6f-a669-433f-9352-0a1f97fc2a52
cm"""
__Course website:__ (Notes, Syllabus)
$(post_img("https://www.dropbox.com/scl/fi/qxhlswxbb2cx4mgqvwuno/mshahrani_qrcode.png?rlkey=jqqkd2vo2z438dcfc776nakd3&dl=1"))

---

__ChatGPT:__ (Course AI assistant)
$(post_img("https://www.dropbox.com/scl/fi/p783o7u8qqrzgxb77qn4x/chatgpt_qrcode.png?rlkey=7wxsf0f1927loqkwij6jwxrnl&dl=1"))

"""

# ╔═╡ 98d3cb65-7c5e-49d7-89df-5b32452a7067
cm"""
$(define("Convex sets (Definition 2.1.1)"))

A set ``S\subseteq\mathbb{R}^n`` is **convex** if for all ``x^1,x^2\in S`` and all ``\lambda\in[0,1]``:  
```math
\lambda x^1 + (1-\lambda)x^2 \in S.
```

**Convex combinations:** 

Finite sums ``\sum_{i=1}^k \lambda_i x^i`` with ``\lambda_i\ge 0``, ``\sum_i \lambda_i=1``.
"""

# ╔═╡ 80e8de21-61e6-4e8f-869c-a364eb07f42d
cm"""
$(define("Linear, Affine, and Convex Combinations"))

Given points ``x^1, x^2, \dots, x^k \in \mathbb{R}^n`` and scalars ``\lambda_1, \dots, \lambda_k \in \mathbb{R}``:

- **Linear combination**  
```math
x = \sum_{i=1}^k \lambda_i x^i
````

* **Affine combination**

```math
x = \sum_{i=1}^k \lambda_i x^i, \quad \text{with } \sum_{i=1}^k \lambda_i = 1
```

* **Convex combination**

```math
x = \sum_{i=1}^k \lambda_i x^i, \quad \text{with } \sum_{i=1}^k \lambda_i = 1, \; \lambda_i \ge 0
```

"""


# ╔═╡ 660817e9-5a67-4eca-9722-d18e57bc5868
cm"""
$(bbl("Lemma","2.1.2 (closure properties)"))
Let ``S_1, S_2`` be convex subsets of ``\mathbb{R}^n``. Then:
1. ``S_1 \cap S_2`` is convex.  
2. **Minkowski sum:** ``S_1 \oplus S_2 = \{x_1+x_2 : x_1\in S_1, x_2\in S_2\}`` is convex.  
3. **Minkowski difference:** ``S_1 \ominus S_2 = \{x_1-x_2 : x_1\in S_1, x_2\in S_2\}`` is convex.
$(ebl())

**Proof ideas (1–2 lines each):**  
- *(1)* If ``y^1,y^2\in S_1\cap S_2``, then each is in both sets; convex combinations stay in each set, hence in the intersection.  
- *(2)* Write ``y^1=x_1^1+x_2^1``, ``y^2=x_1^2+x_2^2`` and use convexity of ``S_1`` and ``S_2`` componentwise:  
  ``\lambda y^1+(1-\lambda)y^2=(\lambda x_1^1+(1-\lambda)x_1^2)+(\lambda x_2^1+(1-\lambda)x_2^2)\in S_1\oplus S_2``.  
- *(3)* Same as (2) with differences.
"""

# ╔═╡ d46a3b04-cc28-4bda-90df-e95778f9bfa2
cm"""
$(define("Affine and Linear Hulls"))

**Affine hull of a set ``S\subseteq\mathbb{R}^n``**  
Smallest affine subspace containing \(S\).  Equivalently, the set of all **affine combinations** of points of \(S\):
```math
\operatorname{aff}(S)
=\left\{\;\sum_{i=1}^k \lambda_i x^i \;:\; x^i\in S,\ \sum_{i=1}^k \lambda_i = 1,\ k\in\mathbb{N}\right\}.
````

Useful identity (any ``x^0\in S``):

```math
\operatorname{aff}(S) \;=\; x^0 \;+\; \operatorname{lin}(S - x^0).
```

**Linear hull of a set ``S\subseteq\mathbb{R}^n``**
Smallest linear subspace (through the origin) containing ``S``; i.e. the **span** of ``S``.  Equivalently, all **linear combinations** of points of ``S``:

```math
\operatorname{lin}(S)
=\left\{\;\sum_{i=1}^k \lambda_i x^i \;:\; x^i\in S,\ \lambda_i\in\mathbb{R},\ k\in\mathbb{N}\right\}.
```
"""

# ╔═╡ 4b3aef5d-7009-4db1-bd29-04acf3bbaea3
cm"""
$(define("Convex hull (Definition 2.1.3)"))
For ``S\subseteq\mathbb{R}^n``, the **convex hull** ``\operatorname{conv}(S)`` is the set of all **finite** convex combinations of points of ``S``:
```math
\operatorname{conv}(S)=\Big\{\sum_{i=1}^k \lambda_i x^i:\ x^i\in S,\ \lambda_i\ge 0,\ \sum_{i=1}^k\lambda_i=1,\ k\in\mathbb{N}\Big\}.
```
"""

# ╔═╡ ed268a43-ac57-49a0-992e-8e9c16cc1d28
cm"""
$(bbl("Lemma","2.1.4 (minimality of the convex hull"))
``\operatorname{conv}(S)`` is the **smallest convex set containing** ``S``:  
```math
\operatorname{conv}(S)=\bigcap\{C\supseteq S: C\text{ convex}\}.
```
$(ebl())
**Proof idea:** Let ``T`` be all finite convex combinations from ``S``. Then ``T`` is convex and contains ``S``, so ``\operatorname{conv}(S)\subseteq T``. Conversely, any convex ``C\supseteq S`` must contain all such combinations, so ``T\subseteq C``. Hence equality.
"""

# ╔═╡ 6250662f-9589-47d4-80f6-be17341180f4
cm"""
$(define("Affinely Independent Points"))

Points ``x^1,\dots,x^{k+1}\in\mathbb{R}^n`` are **affinely independent** iff the ``k`` vectors
```math
x^2-x^1,\; x^3-x^1,\; \dots,\; x^{k+1}-x^1
```
are **linearly independent**.  Equivalently, no point is an affine combination of the others, and
``\dim\big(\operatorname{aff}\{x^1,\dots,x^{k+1}\}\big)=k``.
"""

# ╔═╡ 4c8505f9-7a53-4673-b1c2-26a82daf419d
cm"""
$(define("Polytope & simplex"))
- If ``S`` is a **finite** set, then ``\operatorname{conv}(S)`` is a **polytope**.  
- If the points of ``S`` are **affinely independent**, ``\operatorname{conv}(S)`` is a **simplex**.
"""

# ╔═╡ 91fed789-5722-4fbe-bb8f-98f44cd86a47
cm"""
$(bth("Carathéodory’s theorem"))
If ``x\in\operatorname{conv}(S)\subseteq\mathbb{R}^n``, then ``x`` is a convex combination of at most **``n+1``** points of ``S``.
$(ebl())

**Proof idea (affine dependence reduction):** Start with any convex-combination representation of ``x`` using many points. In ``\mathbb{R}^n`` they are affinely dependent, so there is a nontrivial relation among them. Adjust coefficients along this relation to set one coefficient to zero without moving ``x``. Repeat until only ``n+1`` points remain.
"""

# ╔═╡ 2010bad9-f0b6-41f9-9467-0b9e89daabaa
cm"""
$(bbl("Examples",""))
- **Two points:** ``\operatorname{conv}(\{a,b\})=\{\lambda a+(1-\lambda)b: \lambda\in[0,1]\}`` — the segment.  
- **Triangle:** ``\operatorname{conv}(\{(0,0),(1,0),(0,1)\})`` is the filled triangle.  
- **Minkowski sum:** ``[0,1]\oplus[0,2]=[0,3]``.
"""

# ╔═╡ c96c98c9-787d-4b06-a5e1-a3874b255938
cm"""
$(bbl("Julia demo — computing a 2D convex hull (Polyhedra.jl)"))
This demo builds a convex hull from random points and plots the hull polygon.
"""

# ╔═╡ fa1dec30-f740-46da-b2b3-628ea1527b5f
cm"""
$(define("Core Concepts"))

__Neighborhoods__

For ``x \in \mathbb{R}^n`` and ``\epsilon>0``, the **``\epsilon``-neighborhood** of ``x`` is
```math
N_\epsilon(x) := \{ y \in \mathbb{R}^n : \|y-x\| < \epsilon \}.
```
We use the Euclidean norm unless stated otherwise.

__Closure, Closed Sets, Interior, Boundary (Neighborhood-based)__

**Closure.** For ``S \subseteq \mathbb{R}^n``,
```math
\operatorname{cl}(S) := \{ x \in \mathbb{R}^n : N_\epsilon(x) \cap S \neq \varnothing \ \text{for all } \epsilon>0 \}.
```
**Closed set.** 
```math
S \text{  is} \textbf{ closed } \text{ iff } S = \operatorname{cl}(S).
```

**Interior.**
```math
\operatorname{int}(S) := \{ x \in S : \exists\, \epsilon>0 \text{ with } N_\epsilon(x) \subseteq S \}.
```

**Boundary.**
```math
\partial S := \operatorname{cl}(S) \setminus \operatorname{int}(S).
```
"""

# ╔═╡ 74c85703-9286-4e8b-aea2-5265987f09dd
cm"""
$(bth("Sequential characterization of closed sets"))
** A set ``S \subseteq \mathbb{R}^n`` is closed **iff** for every sequence ``{x^k} \subseteq S`` with ``x^k \to x``, one has ``x \in S``.
$(ebl())

*Sketch.* If ``S`` is closed and ``x^k \to x``, then for any ``\epsilon>0`` we have ``x^k \in N_\epsilon(x)`` for all large ``k``, and since ``x^k \in S``, the neighborhood definition of ``\operatorname{cl}(S)`` gives ``x \in S``. Conversely, if every convergent sequence in ``S`` has its limit in ``S``, then every limit point belongs to ``S``, hence ``S=\operatorname{cl}(S)``.
"""

# ╔═╡ ce80ae3e-1759-462d-bcc8-568c15bd166d
cm"""
$(define("Compact Sets"))
In ``\mathbb{R}^n``, a set ``S`` is **compact** iff it is **closed and bounded** (Heine–Borel). 
"""

# ╔═╡ 90de87f6-3707-45c6-a6a9-3cc9d4debe64
cm"""
$(ex("Examples",""))
- Interval ``[0,1]``: ``\operatorname{cl}([0,1])=[0,1]``, ``\operatorname{int}([0,1])=(0,1)``, ``\partial[0,1]=\{0,1\}``.
- Open ball ``B(0,1)``: ``\operatorname{int}(B)=B``, ``\operatorname{cl}(B)=\overline{B(0,1)}``, ``\partial B`` is the unit sphere.
- Half-space ``\{x: a^\top x \le b\}``: closed, not compact (unless bounded additionally).
- Rationals ``\mathbb{Q} \subset \mathbb{R}``: ``\operatorname{cl}(\mathbb{Q})=\mathbb{R}``, ``\operatorname{int}(\mathbb{Q})=\varnothing``.
"""

# ╔═╡ 5aa374d0-8191-40aa-a06c-d1e115d07b1d
cm"""
$(bth("Line Segment Property for Convex ``S``"))
Let ``S \subseteq \mathbb{R}^n`` is **convex**. If ``x_1 \in \operatorname{cl}(S)`` and ``x_2 \in \operatorname{int}(S)``, then the open line segment
```math
(x_1,x_2] := \{ \lambda x_1 + (1-\lambda) x_2 : 0 < \lambda \le 1 \}
```
is contained in ``\operatorname{int}(S)``.
$(ebl())

**Sketch of proof.** Since ``y \in \operatorname{int}(S)``, there exists ``\epsilon>0`` with ``N_\epsilon(y) \subseteq S``. For any ``0<\lambda\le 1``, the map ``z \mapsto \lambda x + (1-\lambda) z`` sends the ball ``N_\epsilon(y)`` to ``N_{(1-\lambda)\epsilon}(\lambda x + (1-\lambda) y)``, which lies in ``S`` by convexity. Hence each point on ``(x,y]`` has a neighborhood contained in ``S``, i.e., lies in ``\operatorname{int}(S)``.
"""

# ╔═╡ 4f44ea9c-f6de-42a2-88b2-a3c772f3e80d
cm"""
$(bbl("Corollaries",""))

1. Let ``S`` be a convex set. Then int ``S`` is convex. 
2. Let ``S`` be a convex set with a nonempty interior (``\operatorname{int}(S) \neq \varnothing``). Then ``\operatorname{cl}(S)`` is convex if .
3. Let ``S`` be a convex set with a nonempty interior (``\operatorname{int}(S) \neq \varnothing``). Then ``\operatorname{cl}(\operatorname{int}(S)) = \operatorname{cl}(S)``.
4. Let S be a convex set with a nonempty interior (``\operatorname{int}(S) \neq \varnothing``). Then ``\operatorname{int}(\operatorname{cl}(S)) = \operatorname{int}(S)``.
"""

# ╔═╡ 3cc11287-aba8-4e86-a7ef-2172d9c73fb4
cm"""
$(ex("Example","Analytical Worked Example"))

**Claim.** The hypercube ``[0,1]^n`` is compact.

*Proof.* ``[0,1]^n`` is closed as a finite product of closed intervals and bounded; by Heine–Borel it is compact. Alternatively, every sequence in ``[0,1]^n`` has a convergent subsequence (Bolzano–Weierstrass) whose limit lies in ``[0,1]^n``.
"""

# ╔═╡ 6f306138-3d90-47ef-bd46-58ea33ecfcb8
cm"""
$(bbl("Exercises",""))

1. Prove: If ``S`` is convex and ``\operatorname{int}(S) \neq \varnothing``, then ``\operatorname{int}(S)`` is dense in ``S`` iff ``\operatorname{cl}(\operatorname{int}(S))=\operatorname{cl}(S)``.
2. Give an example of a bounded, non-closed set in ``\mathbb{R}^2`` and compute its ``\operatorname{cl}``, ``\operatorname{int}``, ``\partial``.
3. True/False (justify): Every closed set in ``\mathbb{R}^n`` is compact.
"""

# ╔═╡ 2346835d-5fc8-4716-833d-0b8d54dc6cd4
cm"""
$(define(""))
- A **minimizing solution** ``x^*``
```math 
x^* ∈ S \text{ satisfies } f(x^*) ≤ f(x) \text{ for all } x ∈ S. (f(x^*)=\min f(x))
```
- The **infimum** 
```math 
\inf_{x ∈ S} f(x)\text{ is the greatest lower bound of } f \text{ over }S.
```
- The **supremum** 
```math
\sup_{x ∈ S} f(x)\text{ is the least upper bound of }f\text{ over }S.
```
- Existence of infimum/supremum does not guarantee existence of a minimizer/maximizer.
$(ebl())

$(bth("Weierstrass’s Theorem"))
 Let ``S`` be a nonempty, compact set, and let ``f : S → ℝ`` be continuous on ``S``. Then the problem
 ```math
 \min_{x ∈ S} f(x)
 ```
 attains its minimum; that is, there exists ``x^* ∈ S`` such that ``f(x^*) = \min_{x ∈ S} f(x)``.

**Interpretation:** Continuity + compactness (closed and bounded) guarantees that a minimizer exists.
"""

# ╔═╡ 3cefbf49-7946-4978-802c-9c46ff83835e
cm"""
$(ex("","Analytical Example"))
Consider the problem:
```math
\min_{x ∈ [0,1]} f(x) = x^2
```
- ``S = [0,1]`` is compact.
- ``f(x) = x^2`` is continuous.
- By Weierstrass, a minimizer exists.
- Direct evaluation: ``f(0) = 0`` is the minimum.
"""

# ╔═╡ fb465e6a-2a6d-4ba3-9624-19acc0cb236c
cm"""
$(bbl("Parallelogram Law",""))

---

For any ``a, b \in \mathbb{R}^n``, the parallelogram law states
```math
\|a+b\|_2^2 + \|a-b\|_2^2 = 2\|a\|_2^2 + 2\|b\|_2^2.
```

This identity characterizes the Euclidean norm among all norms.
"""

# ╔═╡ 5b4a35bb-89a9-46bb-910a-790af826a6d3
cm"""
$(bth("Theorem 2.4.1 (Minimum Distance to a Closed Convex Set)"))

Let ``S \subset \mathbb{R}^n`` be a nonempty **closed convex** set and let ``y \in \mathbb{R}^n``. Then there exists a **unique** point ``\bar{x} \in S`` such that
```math
\|y-\bar{x}\| = \inf_{x \in S} \|y-x\|.
```
Moreover, ``\bar{x}`` satisfies the **obtuse-angle condition**
```math
(y-\bar{x})^{\top}(x-\bar{x}) \le 0, \quad \forall\, x \in S,
```
and hence ``S`` lies in the half-space
```math
(y-\bar{x})^{\top}(x-\bar{x}) \le 0
```
relative to the hyperplane
```math
(y-\bar{x})^{\top}(x-\bar{x}) = 0
```
passing through ``\bar{x}`` with normal ``a = y-\bar{x}``.

$(ebl())
"""

# ╔═╡ 13796a55-9f61-4aea-9349-f4c3498ec209
cm"""
$(ex("Example", "Equation of a Hyperplane in R^4"))

Consider the hyperplane
```math
H = \{ (x_1, x_2, x_3, x_4) : x_1 + x_2 - x_3 + 2x_4 = 4 \}.
```
The normal vector is
```math
p = (1, 1, -1, 2)^T.
```

Alternatively, the hyperplane can be written in reference to any point in ``H``. For example, take
```math
\bar{x} = (0, 6, 0, -1)^T \in H.
```
Then we can write
```math
H = \{ (x_1, x_2, x_3, x_4) : x_1 + (x_2 - 6) - x_3 + 2(x_4 + 1) = 0 \}.
```
"""


# ╔═╡ 7b542f58-3b2b-4bf5-8183-2ac3588fe464
cm"""
$(define("Separating Hyperplane"))

A **hyperplane** ``H`` in ``\mathbb{R}^n`` is a collection of points of the form
```math
H = \{ x \in \mathbb{R}^n : p^{\top}x = \alpha \},
```
where ``p`` is a nonzero vector in ``\mathbb{R}^n`` and ``\alpha`` is a scalar. The vector ``p`` is called the **normal vector** of the hyperplane.

A hyperplane ``H`` defines two **closed half-spaces**
```math
H^+ = \{ x : p^{\top}x \geq \alpha \}, \quad H^- = \{ x : p^{\top}x \leq \alpha \},
```
and two **open half-spaces**
```math
\{ x : p^{\top}x > \alpha \}, \quad \{ x : p^{\top}x < \alpha \}.
```

Let ``S_1`` and ``S_2`` be nonempty sets in ``\mathbb{R}^n``. Then:

- ``H`` is said to **separate** ``S_1`` and ``S_2`` if
```math
p^{\top}x \geq \alpha \quad \forall x \in S_1, \\
p^{\top}x \leq \alpha \quad \forall x \in S_2.
```
- If, in addition, ``S_1 \cup S_2 \not\subset H``, then ``H`` is said to **properly separate** ``S_1`` and ``S_2``.
- ``H`` is said to **strictly separate** ``S_1`` and ``S_2`` if
```math
p^{\top}x > \alpha \quad \forall x \in S_1, \\
p^{\top}x < \alpha \quad \forall x \in S_2.
```
- ``H`` is said to **strongly separate** ``S_1`` and ``S_2`` if there exists ``\varepsilon > 0`` such that
```math
p^{\top}x \geq \alpha + \varepsilon \quad \forall x \in S_1, \\
p^{\top}x \leq \alpha \quad \forall x \in S_2.
```
"""


# ╔═╡ 43d02e39-037d-445f-9f8e-96866c34e858
cm"""
__Various types of separation.__
$(post_img("https://www.dropbox.com/scl/fi/5hij066cg197299a7b6ib/fig2.9.png?rlkey=b0uy09hylqlsasrtmby0ry9u9&dl=1"))
"""

# ╔═╡ 68aa1d60-0777-4066-9037-3513a705d3f3
cm"""
$(bth("2.4.4 (Separation of a Convex Set and a Point)"))

Let ``S`` be a nonempty closed convex set in ``\mathbb{R}^n`` and ``y \notin S``. Then there exists a nonzero vector ``p`` and a scalar ``\alpha`` such that
```math
p^{\top} y > \alpha, \quad p^{\top} x \leq \alpha \quad \forall x \in S.
```

$(ebl())
"""


# ╔═╡ cfc46ffc-5b02-4b61-bfec-9c313a8d60c4
cm"""
$(bbl("Corollary 1",""))

Let ``S`` be a closed convex set in ``\mathbb{R}^n``. Then ``S`` is the intersection of all half-spaces containing ``S``.
$(ebl())

$(bbl("Corollary 2",""))

Let ``S`` be a nonempty set, and let ``y \notin \operatorname{cl} \operatorname{conv}(S)``, the closure of the convex hull of ``S``. Then there exists a strongly separating hyperplane for ``S`` and ``y``.
"""


# ╔═╡ ea540681-5897-4745-bd43-f1671bbe2d94
cm"""
$(bbl("Remark",""))

The conclusion of Theorem 2.4.4 is equivalent to the following statements:

1. There exists a hyperplane that **strictly** separates ``S`` and ``y``.
2. There exists a hyperplane that **strongly** separates ``S`` and ``y``.
3. There exists a vector ``p`` such that ``p^{\top}y > \sup\{p^{\top}x : x \in S\}``.
4. There exists a vector ``p`` such that ``p^{\top}y < \inf\{p^{\top}x : x \in S\}``.
"""


# ╔═╡ a8cf2766-7784-421c-8a2c-d2fcb26a7344
cm"""
$(bth("2.4.5 (Farkas's Theorem)"))

Let ``A`` be an ``m \times n`` matrix and ``c`` be an ``n``-vector. Then exactly one of the following two systems has a solution:

**System 1:**
```math
Ax \leq 0, \quad c^{\top}x > 0 \quad \text{for some } x \in \mathbb{R}^n.
```

**System 2:**
```math
A^{\top}y = c, \quad y \geq 0 \quad \text{for some } y \in \mathbb{R}^m.
```

$(ebl())
"""


# ╔═╡ 97af3754-7d2c-4db1-8124-c78cacaa4d3f
cm"""
__Farkas's Theorem__
$(post_img("https://www.dropbox.com/scl/fi/s6srf62hlewmc9q6xenz0/fig2.10.png?rlkey=wns7uqrm9x7essqxbbst4bo5l&dl=1"))
"""

# ╔═╡ 7ba2f271-9064-4f72-b51c-3834a3deec56
cm"""
$(bbl("Corollary 1 (Gordan's Theorem)",""))

Let ``A`` be an ``m \times n`` matrix. Then exactly one of the following two systems has a solution:

**System 1:**
```math
Ax < 0 \quad \text{for some } x \in \mathbb{R}^n.
```

**System 2:**
```math
A^{\top}y = 0, \quad y \geq 0 \quad \text{for some nonzero } y \in \mathbb{R}^m.
```

$(bbl("Corollary 2",""))

Let ``A`` be an ``m \times n`` matrix and ``c`` be an ``n``-vector. Then exactly one of the following two systems has a solution:

**System 1:**
```math
Ax \leq 0, \; x \geq 0, \; c^{\top}x > 0 \quad \text{for some } x \in \mathbb{R}^n.
```

**System 2:**
```math
A^{\top}y \geq c, \; y \geq 0 \quad \text{for some } y \in \mathbb{R}^m.
```

$(bbl("Corollary 3",""))

Let ``A`` be an ``m \times n`` matrix, ``B`` be an ``\ell \times n`` matrix, and ``c`` be an ``n``-vector. Then exactly one of the following two systems has a solution:

**System 1:**
```math
Ax \leq 0, \; Bx = 0, \; c^{\top}x > 0 \quad \text{for some } x \in \mathbb{R}^n.
```

**System 2:**
```math
A^{\top}y + B^{\top}z = c, \; y \geq 0 \quad \text{for some } y \in \mathbb{R}^m, z \in \mathbb{R}^{\ell}.
```
"""


# ╔═╡ adb20e86-d93b-4b4f-89a9-1005af9df039
cm"""
$(define("Supporting Hyperplane"))

Let ``S`` be a nonempty set in ``\mathbb{R}^n``, and let ``\bar{x} \in \partial S``. A hyperplane
```math
H = \{ x : p^{\top}(x - \bar{x}) = 0 \}
```
is called a **supporting hyperplane** of ``S`` at ``\bar{x}`` if either 
```math
S \subseteq H^+, \text{ that is, } \qquad p^{\top}(x - \bar{x}) \geq 0 \quad \forall x \in S,
```
or else 
```math
S \subseteq H^-, \text{ that is, } \qquad p^{\top}(x - \bar{x}) \leq 0 \quad \forall x \in S.
```

If, in addition, ``S \not\subset H``, then ``H`` is called a **proper supporting hyperplane** of ``S`` at ``\bar{x}``.
"""


# ╔═╡ 0fffa7d0-74af-4d71-bc18-dd64baf1eb13
cm"""
__Supporting hyperplanes.__
$(post_img("https://www.dropbox.com/scl/fi/1ixb2sca0368s79yvijes/fig2.11.png?rlkey=21yb2gv1apvavyhotnmq5eogq&dl=1"))

$(post_img("https://www.dropbox.com/scl/fi/2cd8mt4trcy0iv42hk3co/fig2.11_2.png?rlkey=d7zbg3lvn2k97b2rerpfj8nz5&dl=1", 300))
"""

# ╔═╡ 887631f8-c063-44eb-a795-be849ed6b284
cm"""
$(bth("2.4.7"))

Let ``S`` be a nonempty convex set in ``\mathbb{R}^n``, and let ``\bar{x} \in \partial S``. Then there exists a hyperplane that supports ``S`` at ``\bar{x}``; that is, there exists a nonzero vector ``p`` such that
```math
p^{\top}(x - \bar{x}) \leq 0 \quad \forall x \in \operatorname{cl} S.
```

$(ebl())
"""


# ╔═╡ a4dcf88f-3df2-4971-a130-c575623d7bc0
cm"""
$(bbl("Corollary",""))

- **(1)** Let ``S`` be a nonempty convex set in ``\mathbb{R}^n`` and let ``\bar{x} \notin \operatorname{int} S``. Then there exists a nonzero vector ``p`` such that
```math
p^{\top}(x - \bar{x}) \leq 0 \quad \forall x \in \operatorname{cl} S.
```

- **(2)** Let ``S`` be a nonempty set in ``\mathbb{R}^n``, and let ``y \notin \operatorname{int} \operatorname{conv}(S)``. Then there exists a hyperplane that separates ``S`` and ``y``.

- **(3)** Let ``S`` be a nonempty set in ``\mathbb{R}^n``, and let ``\bar{x} \in \partial S \cap \partial \operatorname{conv}(S)``. Then there exists a hyperplane that supports ``S`` at ``\bar{x}``.
"""


# ╔═╡ 98c2e6e3-5669-49a8-ab90-ed787f738700
cm"""
$(bth("2.4.8"))

Let ``S_1`` and ``S_2`` be nonempty convex sets in ``\mathbb{R}^n`` and suppose that ``S_1 \cap S_2 = \emptyset``. Then there exists a hyperplane that separates ``S_1`` and ``S_2``; that is, there exists a nonzero vector ``p`` in ``\mathbb{R}^n`` such that
```math
\inf\{p^{\top}x : x \in S_1\} \geq \sup\{p^{\top}x : x \in S_2\}.
```

$(ebl())
"""


# ╔═╡ 17dfacb3-49e6-471b-843a-b3fad927cb26
cm"""$(post_img("https://www.dropbox.com/scl/fi/rph7asup0ip9583fxdnw7/fig_th_2.4.8.png?rlkey=6hjt4s3z97wpvbgy451dmflox&dl=1", 300))"""

# ╔═╡ 1a7107ab-143c-4dc5-bb1b-ec4493915682
cm"""
$(bbl("Corollary 1",""))

Let ``S_1`` and ``S_2`` be nonempty convex sets in ``\mathbb{R}^n``. Suppose that ``\operatorname{int} S_2`` is not empty and that ``S_1 \cap \operatorname{int} S_2 = \emptyset``. Then there exists a hyperplane that separates ``S_1`` and ``S_2``; that is, there exists a nonzero ``p`` such that
```math
\inf\{p^{\top}x : x \in S_1\} \geq \sup\{p^{\top}x : x \in S_2\}.
```
"""


# ╔═╡ df3497ce-8187-4f19-8351-c3cab6fc3e04
cm"""
$(bbl("Corollary 2",""))

Let ``S_1`` and ``S_2`` be nonempty sets in ``\mathbb{R}^n`` such that ``\operatorname{int}(\operatorname{conv}(S_i)) \neq \emptyset`` for ``i=1,2``, but
```math
\operatorname{int}(\operatorname{conv}(S_1)) \cap \operatorname{int}(\operatorname{conv}(S_2)) = \emptyset.
```
Then there exists a hyperplane that separates ``S_1`` and ``S_2``.
"""


# ╔═╡ 55ee2a12-4187-40c0-b832-89de362cca82
cm"""
$(bth("2.4.10 (Strong Separation)"))

Let ``S_1`` and ``S_2`` be closed convex sets, and suppose that ``S_1`` is bounded. If ``S_1 \cap S_2 = \emptyset``, then there exists a hyperplane that strongly separates ``S_1`` and ``S_2``; that is, there exists a nonzero vector ``p`` and ``\varepsilon > 0`` such that
```math
\inf\{p^{\top}x : x \in S_1\} \geq \varepsilon + \sup\{p^{\top}x : x \in S_2\}.
```

$(ebl())
"""


# ╔═╡ be931167-e339-4392-9245-145a8aa6df53
cm"""
__Nonexistence of a strongly separating hyperplane.__
$(post_img("https://www.dropbox.com/scl/fi/93jwy5c99ic9y0db39l5b/fig2.13.png?rlkey=aivrz6xrrpwaitw555jw3h84n&dl=1"))
"""

# ╔═╡ ff77a554-3b51-4592-98dd-a7a7b12efcaa
cm"""
$(bbl("Corollary 1",""))

Let ``S_1`` and ``S_2`` be nonempty sets in ``\mathbb{R}^n``, and suppose that ``S_1`` is bounded. If
```math
\operatorname{cl}(\operatorname{conv}(S_1)) \cap \operatorname{cl}(\operatorname{conv}(S_2)) = \emptyset,
```
then there exists a hyperplane that strongly separates ``S_1`` and ``S_2``.
"""


# ╔═╡ 2dbd94c1-6303-4453-930e-4caee827ea02
cm"""
$(define("Convex Cone"))

A set ``\emptyset\ne C \subseteq \mathbb{R}^n`` is called a **cone** if for every ``x \in C`` and every scalar ``\alpha \geq 0``, we have ``\alpha x \in C``.  

A cone ``C`` is called a **convex cone** if it is also convex, i.e., if for all ``x, y \in C`` and ``\lambda_1, \lambda_2 \geq 0``:  
```math
\lambda_1 x + \lambda_2 y \in C.
```
"""

# ╔═╡ af6e8b0b-efe3-4168-b2ea-035c08969ba8
cm"""
$(define("Polar Cone"))


Let ``S`` be a nonempty set in ``\mathbb{R}^n``. The **polar cone** of ``S``, denoted by ``S^*``, is defined as
```math
S^* = \{ p : p^T x \leq 0, \; \forall x \in S \}.
```
If ``S`` is empty, then ``S^*`` is interpreted as ``\mathbb{R}^n``.
"""

# ╔═╡ 5f1d27f3-3808-4c00-850f-7f2da641e03e
cm"""
$(bbl("Lemma", "2.5.3"))

Let ``S, S_1, S_2`` be nonempty sets in ``\mathbb{R}^n``. Then the following statements hold true:

1. ``S^*`` is a closed convex cone.  
2. ``S \subseteq S^{**}``, where ``S^{**}`` is the polar cone of ``S^*``.
3. ``S_1 \subseteq S_2``, implies that ``S_2^* \subseteq S_1^{*}``.
"""


# ╔═╡ 9b459391-83c7-47a0-a414-92f934da112c
cm"""
$(bth("Theorem 2.5.4"))

Let ``C`` be a nonempty closed convex cone. Then  
```math
C = C^{**}.
```
"""

# ╔═╡ 0d19c555-5d0b-46c4-a765-748094298d75
cm"""
$(bbl("Remark","Farkas's Theorem as a Consequence of Theorem 2.5.4"))

Let ``A`` be an ``m \times n`` matrix, and let  
```math
C = \{ A^t y : y \geq 0 \}.
```
Note that ``C`` is a closed convex cone. It can be verified that  
```math
C^* = \{ x : A x \leq 0 \}.
```
By Theorem 2.5.4, ``c \in C^{**}`` if and only if ``c \in C``. But ``c \in C^{**}`` means that whenever ``x \in C^*``, we have  
```math
c^t x \leq 0.
```
Equivalently, ``A x \leq 0`` implies ``c^t x \leq 0``. By the definition of ``C``, ``c \in C`` means that  
```math
c = A^t y, \quad y \geq 0.
```
Thus, the result ``C = C^{**}`` can be stated as follows:

_System 1:_  
```math
A x \leq 0 \; \Rightarrow \; c^t x \leq 0.
```

_System 2:_  
```math
A^t y = c, \quad y \geq 0.
```

System 1 is consistent if and only if System 2 has a solution ``y``.

This statement can be put in the more usual and equivalent form of Farkas’s theorem. Exactly one of the following two systems has a solution:


_System 1:_
```math
A x \leq 0, \; c^t x > 0 \quad (i.e., \; c \notin C^{**} = C).
```


_System 2:_
```math
A^t y = c, \; y \geq 0 \quad (i.e., \; c \in C).
```
"""

# ╔═╡ 6c72b325-3abc-4ede-aff5-2f5bf192649d
cm"""
$(define("Polyhedral Set"))

A set ``S`` in ``\mathbb{R}^n`` is called a **polyhedral set** if it is the intersection of a finite number of closed half-spaces:
```math
S = \{ x : p_i^T x \leq \alpha_i, \; i=1,\ldots,m \},
```
where ``p_i`` is a nonzero vector and ``\alpha_i`` is a scalar.

Examples:
```math
S = \{ x : A x \leq b \},
```
```math
S = \{ x : A x = b, x \geq 0 \},
```
```math
S = \{ x : A x \geq b, x \geq 0 \}.
```
"""

# ╔═╡ e8cbb8c4-f137-453d-8e8e-87199c79a562
cm"""
$(ex("Example","Polyhedral Set in R2"))

Consider the polyhedral set
```math
S = \{ (x_1, x_2) : -x_1 + x_2 \leq 2, \; x_2 \leq 4, \; x_1 \geq 0, \; x_2 \geq 0 \}.
```

The feasible region is a polygon in the ``(x_1,x_2)``-plane. Its extreme points are located at the intersections of the boundary lines and are illustrated below.
"""

# ╔═╡ 6a62cb74-f9bc-465f-9cdc-9381232b5b3c
cm"""
$(define("Extreme Point"))

Let ``S`` be a nonempty convex set in ``\mathbb{R}^n``. A vector ``x \in S`` is called an **extreme point** of ``S`` if whenever ``x = \lambda x_1 + (1-\lambda) x_2`` with ``x_1, x_2 \in S`` and ``\lambda \in (0,1)``, then ``x = x_1 = x_2``.

Examples:
1. ``S = \{ (x_1,x_2): x_1^2 + x_2^2 \leq 1 \}`` → extreme points are on the boundary circle.
2. ``S = \{ (x_1,x_2): x_1 + x_2 \leq 2, -x_1 + 2x_2 \leq 2, x_1, x_2 \geq 0 \}`` → extreme points are the vertices.
3. Polytopes have finitely many extreme points.
"""

# ╔═╡ c102db9c-b4f3-4398-a055-c0b024a46446
cm"""
$(define("Extreme Direction"))

Let ``S`` be a nonempty, closed convex set in ``\mathbb{R}^n``. A nonzero vector ``d`` is called a **direction**,  or a __recession direction__,  of ``S`` if for each ``x \in S``, ``x + \lambda d \in S`` for all ``\lambda \geq 0``.

It is called an **extreme direction** if it cannot be written as a positive combination of two distinct directions.
"""

# ╔═╡ 95183221-8558-470b-aa94-1a40359c2562
cm"""
$(post_img("https://www.dropbox.com/scl/fi/lcgbda4vbf9eug8z8hwgk/fig2.18.png?rlkey=12u0h5bj1u1hw3pdgbiv3dv41&dl=1"))
"""

# ╔═╡ a9601a07-00a8-4a5c-a1eb-060b47a98912
cm"""
$(bth("Characterization of Extreme Points"))

Let ``S = \{ x : A x = b, x \geq 0 \}``, where ``A`` has rank ``m``. Then ``x`` is an extreme point of ``S`` if and only if ``A`` can be decomposed into ``[B,N]`` such that
```math
x = \begin{bmatrix} x_B \\ x_N \end{bmatrix} = \begin{bmatrix} B^{-1} b \\ 0 \end{bmatrix}, \quad B^{-1} b \geq 0.
```
Such a solution is called a **basic feasible solution (BFS)**.
"""

# ╔═╡ 3225e74a-4f0e-4cdd-a687-b3cad3349823
cm"""
$(bbl("Corollary","Finiteness of Extreme Points"))

The number of extreme points of ``S`` is finite and less than or equal to
```math
\binom{n}{m} = \frac{n!}{m!(n-m)!}.
```
"""

# ╔═╡ 9bd6deb3-425b-451b-9ec1-43ca55b51445
cm"""
$(bth("Existence of Extreme Points"))

If ``S = \{ x : A x = b, x \geq 0 \}`` is nonempty, then ``S`` has at least one extreme point.
"""

# ╔═╡ 965ab318-1528-4d9e-a586-e7ea4382ac57
cm"""
$(bth("Characterization of Extreme Directions"))

Let ``S = \{ x : A x = b, x \geq 0 \} \neq \emptyset``. A vector ``\bar{d}`` is an extreme direction of ``S`` if and only if ``A`` can be decomposed into ``[B, N]`` such that ``B^{-1} a_j \leq 0`` for some column ``a_j`` of ``N`` and ``\bar{d}`` is a positive multiple of
```math
d = \begin{bmatrix} -B^{-1} a_j \\ e_j \end{bmatrix},
```
where ``e_j`` is an  ``n - m`` unit vector.
"""

# ╔═╡ 64cedf28-216b-4b03-924b-c07c11b33a51
cm"""
$(bbl("Corollary",""))
The number of extreme directions of ``S`` is finite.
"""

# ╔═╡ c38ad4af-2245-4482-900e-7691577bcff2
cm"""
$(bth("Representation Theorem"))

Let ``S`` be a nonempty polyhedral set ``S = \{ x : A x = b, x \geq 0 \}``. Let ``x_1, \ldots, x_k`` be the extreme points and ``d_1, \ldots, d_\ell`` be the extreme directions. Then ``x \in S`` if and only if
```math
x = \sum_{j=1}^k \lambda_j x_j + \sum_{j=1}^\ell \mu_j d_j,
```
with
```math
\sum_{j=1}^k \lambda_j = 1, \quad \lambda_j \geq 0, \; \mu_j \geq 0.
```
This gives the inner representation of ``S``.
"""

# ╔═╡ 78422f45-8847-45a9-a269-3f0ee5918076
cm"""
$(bbl("Corollary","Existence of Extreme Directions"))

A nonempty polyhedral set ``S = \{ x : A x = b, x \geq 0 \}`` has at least one extreme direction if and only if it is unbounded.
"""

# ╔═╡ 69c4460b-93f7-440e-98c0-9dd9c66eb8fb
cm"""
$(define("Convex and Concave Functions"))

Let ``f : S \to \mathbb{R}``, where ``S`` is a nonempty convex set in ``\mathbb{R}^n``.  
- The function ``f`` is said to be **convex** on ``S`` if

```math
f(\lambda x_1 + (1 - \lambda)x_2) \leq \lambda f(x_1) + (1 - \lambda) f(x_2)
```

$(add_space(15))for each ``x_1, x_2 \in S`` and for each ``\lambda \in (0,1)``.

- The function ``f`` is called **strictly convex** on ``S`` if the above inequality is true as a **strict inequality** for each distinct ``x_1`` and ``x_2`` in ``S`` and for each ``\lambda \in (0,1)``.

- The function ``f : S \to \mathbb{R}`` is called **concave** (resp. **strictly concave**) on ``S`` if ``-f`` is convex (resp. strictly convex) on ``S``.
"""


# ╔═╡ 9120ebc6-f004-40e5-a2bb-69f5c7b6f74b
cm"""
$(bbl("Remarks", ""))
1. Let ``f_1, f_2, \ldots, f_k : \mathbb{R}^n \to \mathbb{R}`` be convex functions. Then:

   (a) ``f(x) = \sum_{j=1}^k \alpha_j f_j(x)``, where ``\alpha_j > 0`` for ``j = 1,2,\ldots,k``, is a convex function (see Exercise 3.8).

   (b) ``f(x) = \max \{ f_1(x), f_2(x), \ldots, f_k(x) \}`` is a convex function (see Exercise 3.9).

2. Suppose that ``g : \mathbb{R}^n \to \mathbb{R}`` is a concave function. Let ``S = \{ x : g(x) > 0 \}``, and define ``f : S \to \mathbb{R}`` as ``f(x) = 1 / g(x)``. Then ``f`` is convex over ``S`` (see Exercise 3.11).

3. Let ``g : \mathbb{R} \to \mathbb{R}`` be a nondecreasing, univariate, convex function, and let ``h : \mathbb{R}^n \to \mathbb{R}`` be a convex function. Then the composite function ``f : \mathbb{R}^n \to \mathbb{R}`` defined as ``f(x) = g[h(x)]`` is a convex function (see Exercise 3.10).

4. Let ``g : \mathbb{R}^m \to \mathbb{R}`` be a convex function, and let ``h : \mathbb{R}^n \to \mathbb{R}^m`` be an affine function of the form ``h(x) = A x + b``, where ``A`` is an ``m \times n`` matrix and ``b`` is an ``m \times 1`` vector. Then the composite function ``f : \mathbb{R}^n \to \mathbb{R}`` defined as ``f(x) = g[h(x)]`` is a convex function (see Exercise 3.16).

"""


# ╔═╡ c6bd9fa2-e9e8-4bb3-a866-1b450d3bf8d8
cm"""
$(define("Level Set"))

Let ``f``  be a convex function. The set  

```math
S_\alpha = \{ x \in S : f(x) \leq \alpha \}, \quad \alpha \in \mathbb{R},
```  
is called  a **level set**. Sometimes this set is called a **lower-level set**, to differentiate it from the **upper-level set**  

```math
\{ x \in S : f(x) \geq \alpha \}.
```
"""


# ╔═╡ e949636a-b55d-4ae3-b605-26ad244a5be2
cm"""
$(bbl("Lemma", "Convexity of Level Sets"))

Let ``S`` be a nonempty convex set in ``\mathbb{R}^n``, and let ``f : S \to \mathbb{R}`` be a convex function.  
Then the level set  

```math
S_\alpha = \{ x \in S : f(x) \leq \alpha \}, \quad \alpha \in \mathbb{R},
```
is convex.

"""

# ╔═╡ 248cef95-c012-446c-8e4d-258fb8f06410
cm"""
$(bth("3.1.3"))

Let ``S`` be a nonempty convex set in ``\mathbb{R}^n``, and let ``f: S \to \mathbb{R}`` be convex. Then ``f`` is continuous on the interior of ``S``.
"""


# ╔═╡ 5617b3db-ed53-4c81-b7a5-073c8f34fd9f
cm"""
$(bbl("Remark",""))

- Discontinuity only are allowed at the boundary of ``S``, as illustrated by the following convex function defined on ``S = \{ x : -1 \leq x \leq 1 \}``:

```math
f(x) =
\begin{cases}
x^2 & \text{for } |x| < 1, \\
2   & \text{for } |x| = 1.
\end{cases}
```
"""

# ╔═╡ 0d650b2c-6ec8-4822-b875-43a6a5b52879
cm"""
$(define("Directional Derivative "))

Let ``S`` be a nonempty set in ``\mathbb{R}^n``, and let ``f: S \to \mathbb{R}``. Let ``\bar{x} \in S`` and ``d`` be a nonzero vector such that ``\bar{x} + \lambda d \in S`` for ``\lambda > 0`` and sufficiently small.  

The **directional derivative** of ``f`` at ``\bar{x}`` along the vector ``d``, denoted by ``f'(\bar{x}; d)``, is given by the following limit if it exists:

```math
f'(\bar{x}; d) = \lim_{\lambda \to 0^+} \frac{f(\bar{x} + \lambda d) - f(\bar{x})}{\lambda}.
```
"""


# ╔═╡ 6e464ce7-f3d2-43dd-b529-a6171e9dd898
cm"""
$(bbl("Lemma",""))

Let ``f: \mathbb{R}^n \to \mathbb{R}`` be a convex function.  Consider any point ``\bar{x} \in \mathbb{R}^n`` and a nonzero direction ``d \in \mathbb{R}^n``.  
Then the directional derivative ``f'(\bar{x}; d)`` of ``f`` at ``\bar{x}`` in the direction ``d`` __exists__.
"""

# ╔═╡ ae88db5b-7923-4216-98af-570d1ecf39dc
cm"""
$(define("Epigraph"))
Let ``S`` be a nonempty set in ``R^n``, and let ``f: S \rightarrow R``. The epigraph of ``f``, denoted by epi ``f``, is a subset of ``R^{n+1}`` defined by
```math
\{(\mathbf{x}, y): \mathbf{x} \in S, y \in R, y \geq f(\mathbf{x})\}
```

The hypograph of ``f``, denoted by hyp ``f``, is a subset of ``R^{n+1}`` defined by
```math
\{(\mathbf{x}, y): \mathbf{x} \in S, y \in R, y \leq f(\mathbf{x})\}
```
"""

# ╔═╡ 525e28ee-6af9-4e5b-9ce1-1f87881ff681
cm"""
$(bth("3.2.2"))

Let ``S`` be a nonempty convex set in ``R^n``, and let ``f: S \rightarrow R``. Then ``f`` is convex if and only if epi ``f`` is a convex set.
"""

# ╔═╡ 69aecb16-725a-49dd-aba4-9775d797aaae
cm"""
$(define("Subgradients"))
Let ``S`` be a nonempty convex set in ``R^n``, and let ``f: S \rightarrow R`` be convex. Then ``\xi`` is called a subgradient of ``f`` at ``\overline{\mathbf{x}} \in S`` if
```math
f(\mathbf{x}) \geq f(\overline{\mathbf{x}})+\xi^t(\mathbf{x}-\overline{\mathbf{x}}) \quad \text { for all } \mathbf{x} \in S
```

Similarly, let ``f: S \rightarrow R`` be concave. Then ``\xi`` is called a subgradient of ``f`` at ``\overline{\mathbf{x}} \in S`` if
```math
f(\mathbf{x}) \leq f(\overline{\mathbf{x}})+\xi^t(\mathbf{x}-\overline{\mathbf{x}}) \quad \text { for all } \mathbf{x} \in S
```
$(ebl())
$(bbl("Remarks",""))
The collection of subgradients of ``f`` at ``\overline{\mathbf{x}}`` (known as the subdifferential of ``f`` at ``\overline{\mathbf{x}}`` ) is a convex set. 
We will denote this set by 
```math
\partial f(\overline{\mathbf{x}})
```
"""

# ╔═╡ 017ac677-b8d3-40a9-90c3-ec7d4c463f0f
cm"""
$(bth("3.2.5"))

Let ``S`` be a nonempty convex set in ``R^n``, and let ``f: S \rightarrow R`` be convex. Then for ``\overline{\mathbf{x}} \in`` int ``S``, there exists a vector ``\xi`` such that the hyperplane
```math
H=\left\{(\mathbf{x}, y): y=f(\overline{\mathbf{x}})+\xi^t(\mathbf{x}-\overline{\mathbf{x}})\right\}
```
supports epi ``f`` at ``[\overline{\mathbf{x}}, f(\overline{\mathbf{x}})]``. In particular,
```math
f(\mathbf{x}) \geq f(\overline{\mathbf{x}})+\xi^t(\mathbf{x}-\overline{\mathbf{x}}) \quad \text { for each } \mathbf{x} \in S ;
```
that is, ``\boldsymbol{\xi}`` is a subgradient of ``f`` at ``\overline{\mathbf{x}}``.
"""

# ╔═╡ 44c72d9b-0adc-49e7-a878-184951cefe0d
cm"""
$(bbl("Corollary",""))
Let ``S`` be a nonempty convex set in ``R^n``, and let ``f: S \rightarrow R`` be strictly convex. Then for ``\overline{\mathbf{x}} \in \operatorname{int} S`` there exists a vector ``\xi`` such that
```math
f(\mathbf{x})>f(\overline{\mathbf{x}})+\xi^t(\mathbf{x}-\overline{\mathbf{x}}) \quad \text { for all } \mathbf{x} \in S, \mathbf{x} \neq \overline{\mathbf{x}} .
```
"""

# ╔═╡ 3c14bb03-4fff-47d0-8fc8-643561950b2a
cm"""
$(bth("3.2.6"))
Let ``S`` be a nonempty convex set in ``R^n``, and let ``f: S \rightarrow R``. Suppose that for each point ``\overline{\mathbf{x}} \in`` int ``S`` there exists a subgradient vector ``\xi`` such that
```math
f(\mathbf{x}) \geq f(\overline{\mathbf{x}})+\xi^t(\mathbf{x}-\overline{\mathbf{x}}) \quad\text{ for each } \mathbf{x} \in S.
```
Then, ``f`` is convex on ``\operatorname{int} S``.
"""

# ╔═╡ 657034bf-4034-4b37-995b-cc3e22a6ff19
cm"""
$(define("Differentiable Function"))
Let ``S`` be a nonempty set in ``R^n``, and let ``f: S \rightarrow R``. Then ``f`` is said to be differentiable at ``\overline{\mathbf{x}} \in \operatorname{int} S`` if there exist a vector ``\nabla f(\overline{\mathbf{x}})``, called the gradient vector, and a function ``\alpha: R^n \rightarrow R`` such that
```math
f(\mathbf{x})=f(\overline{\mathbf{x}})+\nabla f(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}})+\|\mathbf{x}-\overline{\mathbf{x}}\| \alpha(\overline{\mathbf{x}} ; \mathbf{x}-\overline{\mathbf{x}}) \quad \text { for each } \mathbf{x} \in S
```
where ``\lim _{\mathbf{x} \rightarrow \overline{\mathbf{x}}} \alpha(\overline{\mathbf{x}} ; \mathbf{x}-\overline{\mathbf{x}})=0``. 

- The function ``f`` is said to be differentiable on the open set ``S^{\prime} \subseteq S`` if it is differentiable at each point in ``S^{\prime}``. 
- The representation of ``f`` above is called a first-order (Taylor series) expansion of ``f`` at (or about) the point ``\overline{\mathbf{x}}``; and without the implicitly defined remainder term involving the function ``\alpha``, the resulting representation is called a first-order (Taylor series) approximation of ``f`` at (or about) the point ``\overline{\mathbf{x}}``.
$(ebl())

"""

# ╔═╡ f22887ef-aabb-4a52-9f38-afb2f82ed16b
cm"""
$(bbl("Lemma","3.3.2"))
Let ``S`` be a nonempty convex set in ``R^n``, and let ``f: S \rightarrow R`` be convex. Suppose that ``f`` is differentiable at ``\overline{\mathbf{x}} \in`` int ``S``. Then the collection of subgradients of ``f`` at ``\overline{\mathbf{x}}`` is the singleton set ``\{\nabla f(\overline{\mathbf{x}})\}``.
"""

# ╔═╡ 4a15eb36-dd71-4ad4-9512-6f189a448bb8
cm"""
$(bth("3.3.3"))

Let ``S`` be a nonempty open convex set in ``R^n``, and let ``f: S \rightarrow R`` be differentiable on ``S``. Then ``f`` is convex if and only if for any ``\overline{\mathbf{x}} \in S``, we have
```math
f(\mathbf{x}) \geq f(\overline{\mathbf{x}})+\nabla f(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}}) \quad \text { for each } \mathbf{x} \in S
```

Similarly, ``f`` is strictly convex if and only if for each ``\overline{\mathbf{x}} \in S``, we have
```math
f(\mathbf{x})>f(\overline{\mathbf{x}})+\nabla f(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}}) \quad \text { for each } \mathbf{x} \neq \overline{\mathbf{x}} \text { in } S
```
"""

# ╔═╡ de381e37-2c8e-4e42-a4b8-292cb349ed3d
cm"""
$(bth("3.3.4"))

Let ``S`` be a nonempty open convex set in ``R^n`` and let ``f: S \rightarrow R`` be differentiable on ``S``. Then ``f`` is convex if and only if for each ``\mathbf{x}_1, \mathbf{x}_2 \in S`` we have
```math
\left[\nabla f\left(\mathbf{x}_2\right)-\nabla f\left(\mathbf{x}_1\right)\right]^t\left(\mathbf{x}_2-\mathbf{x}_1\right) \geq 0 .
```

Similarly, ``f`` is strictly convex if and only if, for each distinct ``\mathbf{x}_1, \mathbf{x}_2 \in S``, we have
```math
\left[\nabla f\left(\mathbf{x}_2\right)-\nabla f\left(\mathbf{x}_1\right)\right]^t\left(\mathbf{x}_2-\mathbf{x}_1\right)>0 .
```
"""

# ╔═╡ 3e9400d7-48f7-4b1c-8bdf-441b874e99a7
cm"""
$(define("Twice Differentiability"))
Let ``S`` be a nonempty set in ``R^n``, and let ``f: S \rightarrow R``. Then ``f`` is said to be twice differentiable at ``\overline{\mathbf{x}} \in \operatorname{int} S`` if there exist a vector ``\nabla f(\overline{\mathbf{x}})``, and an ``n \times n`` symmetric matrix ``\mathbf{H}(\overline{\mathbf{x}})``, called the Hessian matrix, and a function ``\alpha: R^n \rightarrow R`` such that
```math
f(\mathbf{x})=f(\overline{\mathbf{x}})+\nabla f(\mathbf{x})^t(\mathbf{x}-\overline{\mathbf{x}})+\frac{1}{2}(\mathbf{x}-\overline{\mathbf{x}})^t \mathbf{H}(\overline{\mathbf{x}})(\mathbf{x}-\overline{\mathbf{x}})+\|\mathbf{x}-\overline{\mathbf{x}}\|^2 \alpha(\overline{\mathbf{x}} ; \mathbf{x}-\overline{\mathbf{x}})
```
for each ``\mathbf{x} \in S``, where ``\lim _{\mathbf{x} \rightarrow \overline{\mathbf{x}}} \alpha(\overline{\mathbf{x}} ; \mathbf{x}-\overline{\mathbf{x}})=0``. 
- The function ``f`` is said to be twice differentiable on the open set ``S^{\prime} \subseteq S`` if it is twice differentiable at each point in ``S^{\prime}``.
"""

# ╔═╡ c78f468a-75c8-41dc-896c-1021765adf83
cm"""
$(ex(1))
Let ``f\left(x_1, x_2\right)=2 x_1+6 x_2-2 x_1^2-3 x_2^2+4 x_1 x_2``.
"""

# ╔═╡ 57d1b602-1b3f-4de6-85d3-ac158a01bcc8
cm"""
$(bth("3.3.7"))

Let ``S`` be a nonempty open convex set in ``R^n``, and let ``f: S \rightarrow R`` be twice differentiable on ``S``. Then ``f`` is convex if and only if the Hessian matrix is __positive semidefinite__ at each point in ``S``.
"""

# ╔═╡ 20fdb0be-e603-4414-bb8e-d9df5b3d4666
cm"""
$(bth("3.3.8"))

Let ``S`` be a nonempty open convex set in ``R^n``, and let ``f: S \rightarrow R`` be twice differentiable on ``S``. If the Hessian matrix is __positive definite__ at each point in ``S, f`` is strictly convex. Conversely, if ``f`` is strictly convex, the Hessian matrix is __positive semidefinite__ at each point in ``S``. However, if ``f`` is strictly convex and quadratic, its Hessian is positive definite.
"""

# ╔═╡ 1e87dde2-48a2-4158-b20a-af94fc2f3308
cm"""
$(bth("3.3.9 (R version)"))

Let ``S`` be a nonempty open convex set in ``R``, and let ``f: S \rightarrow R`` be infinitely differentiable. Then ``f`` is strictly convex on ``S`` if and only if for each ``\bar{x} \in S``, there exists an even ``n`` such that ``f^{(n)}(\bar{x}) > 0``, while ``f^{(j)}(\bar{x})=0`` for any ``1 < j < n``, where ``f^{(j)}`` denotes the ``j`` th-order derivative of ``f``.
"""

# ╔═╡ 5aaa0e88-cdd3-4688-adfa-11487ab512ac
cm"""
$(bth("3.3.10"))

Consider a function ``f: R^n \rightarrow R``, and for any point ``\overline{\mathbf{x}} \in R^n`` and a nonzero direction ``\mathbf{d} \in R^n``, define ``F_{(\overline{\mathbf{x}} ; \mathbf{d})}(\lambda)=f(\overline{\mathbf{x}}+\lambda \mathbf{d})`` as a function of ``\lambda \in R``. Then ``f`` is (strictly) convex if and only if ``F_{(\overline{\mathbf{x}} ; \mathbf{d})}`` is (strictly) convex for all ``\overline{\mathbf{x}}`` and ``\mathbf{d} \neq \mathbf{0}`` in ``R^n``.
"""

# ╔═╡ 53023197-6e74-4968-8fa4-ce0fee8bf6d9
cm"""
$(bth("3.3.12 (Checking for PSD/PD)"))

Let ``\mathbf{H}`` be a symmetric ``n \times n`` matrix with elements ``h_{i j}``.
- (a) If ``h_{i i} \leq 0`` for any ``i \in\{1, \ldots, n\}, \mathbf{H}`` is not positive definite; and if ``h_{i i}<`` 0 for any ``i \in\{1, \ldots, n\}, \mathbf{H}`` is not positive semidefinite.

- (b) If ``h_{i i}=0`` for any ``i \in\{1, \ldots, n\}``, we must have ``h_{i j}=h_{j i}=0`` for all ``j= 1, \ldots, n`` as well, or else ``\mathbf{H}`` is not positive semidefinite.

- (c) If ``n=1, \mathbf{H}`` is positive semidefinite (positive definite) if and only if ``h_{11} \geq 0(>0)``. Otherwise, if ``n \geq 2``, let
```math
\mathbf{H}=\left[\begin{array}{cc}
h_{11} & \mathbf{q}^t \\
\mathbf{q} & \mathbf{G}
\end{array}\right]
```
$(add_space(10))in partitioned form, where ``\mathbf{q}=\mathbf{0}`` if ``h_{11}=0``, and otherwise, ``h_{11}>0``. 
Perform elementary Gauss-Jordan operations using the first row of ``\mathbf{H}`` to reduce it to the following matrix in either case:
```math
\mathbf{H}=\left[\begin{array}{cc}
h_{11} & \mathbf{q}^t \\
\mathbf{0} & \mathbf{G}_{\text {new }}
\end{array}\right]
```
$(add_space(10))Then ``\mathbf{G}_{\text {new }}`` is a symmetric ``(n-1) \times(n-1)`` matrix, and ``\mathbf{H}`` is positive semidefinite if and only if ``\mathbf{G}_{\text {new }}`` is positive semidefinite. Moreover, if ``h_{11}>0, \mathbf{H}`` is positive definite if and only if ``\mathbf{G}_{\text {new }}`` is positive definite.
"""

# ╔═╡ 08a79da6-d7be-4749-8c7a-960dd85d9404
cm"""
$(bbl("Corollary"))
Let ``\mathbf{H}`` be an ``n \times n`` symmetric matrix. Then ``\mathbf{H}`` is positive definite if and only if it is positive semidefinite and nonsingular.
"""

# ╔═╡ 9d66f924-30a8-428e-9f48-06b0be9b9687
cm"""
$(ex(3)) 
Consider the matrix
```math
\mathbf{H}=\left[\begin{array}{lll}
2 & 1 & 2 \\
1 & 2 & 3 \\
2 & 3 & 4
\end{array}\right]
```
"""

# ╔═╡ 997295b0-bc83-4a2e-a81f-52212c041152
cm"""
$(define("Solutions"))
Let ``f: R^n \rightarrow R`` and consider the problem 
```math
\min f(\mathbf{x})
\quad \text{subject to}\quad  \mathbf{x} \in S.
```
- A point ``\mathbf{x} \in S`` is called a __feasible solution__ to the problem. 
- A feasible solution ``\overline{\mathbf{x}}`` is called an __optimal solution__, a __global optimal solution__, or simply a __solution__ to the problem if  
```math
f(\overline{\mathbf{x}}) \leq f(\mathbf{x}) \quad \text{for each} \quad \mathbf{x} \in S.
``` 

- The collection of optimal solutions are called __alternative optimal solutions__.

- A feasible solution ``\overline{\mathbf{x}}`` is called a __local optimal solution__ if there exists an ``\varepsilon``-neighborhood ``N_{\varepsilon}(\overline{\mathbf{x}})`` around ``\overline{\mathbf{x}}`` such that 

```math
f(\overline{\mathbf{x}})\leq f(\mathbf{x}) \quad \text{for each }\quad \mathbf{x} \in S \cap N_{\varepsilon}(\overline{\mathbf{x}}).
```

- A feasible solution ``\overline{\mathbf{x}}``  is called a __strict local optimal solution__ if there exists ``\varepsilon>0`` such that 
```math
f(\overline{\mathbf{x}})< f(\mathbf{x}) \quad \text{for all}\quad \mathbf{x} \in S \cap N_{\varepsilon}(\overline{\mathbf{x}}), \mathbf{x} \neq \overline{\mathbf{x}}.
```

-  A feasible solution ``\overline{\mathbf{x}}`` is called a __strong__ or __isolated__ local optimal solution if it is the only local minimum in ``S \cap N_{\varepsilon}(\overline{\mathbf{x}})`` for some ``\varepsilon`` neighborhood ``N_{\varepsilon}(\overline{\mathbf{x}})`` around ``\overline{\mathbf{x}}``.

- All these types (``\varepsilon``) of local optima or minima are sometimes also referred to as __relative minima__. 

Figure below illustrates instances of local and global minima for the problem of minimizing ``f(\mathbf{x})`` subject to ``\mathbf{x} \in S``, where ``f`` and ``S`` are shown in the figure.

$(post_img("https://www.dropbox.com/scl/fi/op414xeqtdnz5o5qsiogt/fig3.6.png?rlkey=xopzhhjcubqvnuowdpthu2ju3&dl=1"))
"""

# ╔═╡ 9f633dba-e6b5-4d9f-b7b0-2f505c4642ab
cm"""
$(bbl("Remarks",""))
- Any isolated local minimum is always strict.
- A **strict local minimum** does **not** always have to be **isolated**. That means you can have a strict local minimum point, but there could still be other local minima arbitrarily close to it.
"""

# ╔═╡ da04715d-9713-4230-9d25-df069c19c9d4
cm"""
$(ex("Example",""))

Consider
```math
  f(x) = \begin{cases} 
  1 & x = 1, \\ 
  2 & \text{otherwise}.
  \end{cases}
```

* Domain: ``S = \mathbb{R}``.
* ``x=1`` is a **strict local minimum**.
* But it’s **not isolated**: any neighborhood around ``x=1`` contains other points (all with ``f(x)=2``) that are also local minima.

---

$(ex("Example",""))

* Function: ``f(x) = x^2``.
* Domain: ``S = \{ 1/2^k : k=0,1,2,\ldots \} \cup \{0\}``. (A nonconvex set made of discrete points approaching 0.)
* For any ``k \ge 0``, ``x = 1/2^k`` is a **strict local minimum** because it’s an isolated point of the domain.
* At ``x=0``, we also have a **strict local minimum** (in fact, the **global minimum**) with value ``f(0)=0``.
* But ``x=0`` is **not isolated** because points like ``1/2, 1/4, 1/8, \ldots`` approach it.



"""

# ╔═╡ 1c9a2d2f-39d5-453e-8a82-bda18570e762
cm"""
$(bth("3.4.2"))
Let ``S`` be a nonempty convex set in ``R^n``, and let ``f: S \rightarrow R`` be convex on ``S``. Consider the problem 
```math 
\min f(\mathbf{x}) \quad \text{subject to}\quad \mathbf{x} \in S.
```
Suppose that ``\overline{\mathbf{x}} \in S`` is a local optimal solution to the problem.
1. Then ``\overline{\mathbf{x}}`` is a global optimal solution.
2. If either ``\overline{\mathbf{x}}`` is a strict local minimum or ``f`` is strictly convex, ``\overline{\mathbf{x}}`` is the unique global optimal solution and is also a strong local minimum.
"""

# ╔═╡ 61e3e889-4444-4ecb-91a1-3d8f91d0054a
cm"""
$(bth("3.4.3"))

Let ``f: R^n \rightarrow R`` be a convex function, and let ``S`` be a nonempty convex set in ``R^n``. Consider the problem 
```math 
\min f(\mathbf{x}) \quad \text{subject to}\quad \mathbf{x} \in S.
```
The point ``\overline{\mathbf{x}} \in S`` is an optimal solution to this problem if and only if ``f`` has a subgradient ``\boldsymbol{\xi}`` at ``\overline{\mathbf{x}}`` such that ``\xi^t(\mathbf{x}-\overline{\mathbf{x}}) \geq 0`` for all ``\mathbf{x} \in S``.
"""

# ╔═╡ 44ecea1f-be97-4fb5-a62e-0177b98d5404
cm"""
$(bbl("Corollary","1"))
Under the assumptions of Theorem 3.4.3, if ``S`` is __open__, ``\overline{\mathbf{x}}`` is an optimal solution to the problem if and only if there exists a zero subgradient of ``f`` at ``\overline{\mathbf{x}}``. In other words,
```math
0\in \partial f(\overline{\bf x}).
```
In particular, if ``S=R^n, \overline{\mathbf{x}}`` is a global minimum if and only if there exists a zero subgradient of ``f`` at ``\overline{\mathbf{x}}``.
"""

# ╔═╡ c3cfdf7e-0621-4add-9d85-945508a44eec
cm"""
$(bbl("Corollary","2"))
In addition to the assumptions of the theorem, suppose that ``f`` is differentiable. Then ``\overline{\mathbf{x}}`` is an optimal solution if and only if ``\nabla f(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}}) \geq 0`` for all ``\mathbf{x} \in S``. Furthermore, if ``S`` is open, ``\overline{\mathbf{x}}`` is an optimal solution if and only if ``\nabla f(\overline{\mathbf{x}})=0``.
"""

# ╔═╡ 8c011096-3474-4d80-b40d-d72fd50e621b
cm"""
$(post_img("https://www.dropbox.com/scl/fi/nl863u4ywfqq9u4sa11ds/fig3.8.png?rlkey=1zximze8kzaoy6rpqepavm5r7&dl=1"))
"""

# ╔═╡ b283d6eb-0d60-4eb8-af1b-2a30f4f5e596
cm"""
$(bth("3.4.4"))
Consider the problem 
```math
\min\quad f(\mathbf{x}) \text{ subject to } \mathbf{x} \in S,
```
where ``f`` is a convex and twice differentiable function and ``S`` is a convex set, and suppose that there exists an optimal solution ``\overline{\mathbf{x}}``. Then the set of alternative optimal solutions is characterized by the set
```math
S^*=\left\{\mathbf{x} \in S: \nabla f(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}}) \leq 0 \text { and } \nabla f(\mathbf{x})=\nabla f(\overline{\mathbf{x}})\right\}
```
"""

# ╔═╡ f6c2cbee-b331-4132-946e-bbe9fd3b5881
cm"""

$(bbl("Corollary","1"))

The set ``S^*`` of alternative optimal solutions can equivalently be defined as
```math
S^*=\left\{\mathbf{x} \in S: \nabla f(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}})=0 \text { and } \nabla f(\mathbf{x})=\nabla f(\overline{\mathbf{x}})\right\}
```
"""

# ╔═╡ 8d3beee4-64ef-4942-adc2-739122182dc9
cm"""
$(bbl("Corollary", "2"))
Suppose that ``f`` is a quadratic function given by ``f(\mathbf{x})=\mathbf{c}^t \mathbf{x}+(1 / 2) \mathbf{x}^t \mathbf{H} \mathbf{x}`` and that ``S`` is polyhedral. Then ``S^*`` is a polyhedral set given by
```math
\begin{gathered}
S^*=\left\{\mathbf{x} \in S: \mathbf{c}^t(\mathbf{x}-\overline{\mathbf{x}}) \leq 0, \mathbf{H}(\mathbf{x}-\overline{\mathbf{x}})=\mathbf{0}\right\}\\ =\left\{\mathbf{x} \in S: \mathbf{c}^t(\mathbf{x}-\overline{\mathbf{x}})=0\right. ,
\mathbf{H}(\mathbf{x}-\overline{\mathbf{x}})=\mathbf{0}\}
\end{gathered}
```
"""

# ╔═╡ f6180dbc-4c0d-4b06-9670-9a4122e7681c
cm"""
$(ex("Example",""))
```math
\begin{aligned}
\operatorname{Minimize}\left(x_1-\frac{3}{2}\right)^2 & +\left(x_2-5\right)^2 \\
\text { subject to }-x_1+x_2 & \leq 2 \\
2 x_1+3 x_2 & \leq 11 \\
-x_1 & \leq 0 \\
-x_2 & \leq 0
\end{aligned}
```
"""

# ╔═╡ d3bcf477-8bf6-4193-a6b1-d0f924c32bf5
cm"""
$(bth("3.4.6"))

Let ``f: R^n \rightarrow R`` be a convex function, and let ``S`` be a nonempty convex set in ``R^n``. Consider the problem to maximize ``f(\mathbf{x})`` subject to ``\mathbf{x} \in S``. If ``\overline{\mathbf{x}} \in S`` is a local optimal solution, ``\boldsymbol{\xi}^t(\mathbf{x}-\overline{\mathbf{x}}) \leq 0`` for each ``\mathbf{x} \in S``, where ``\boldsymbol{\xi}`` is any subgradient of ``f`` at ``\overline{\mathbf{x}}``.
"""

# ╔═╡ 8b0eaffa-18c3-4174-a4db-8cb514c728bf
cm"""
$(bbl("Corollary",""))
In addition to the assumptions of the theorem, suppose that ``f`` is differentiable. If ``\overline{\mathbf{x}} \in S`` is a local optimal solution, ``\nabla f(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}}) \leq 0`` for all ``\mathbf{x} \in S``.
"""

# ╔═╡ a5ccf185-1a61-47ea-9fe2-b33b7bba8e6c
cm"""
$(bth("3.4.7"))

Let ``f: R^n \rightarrow R`` be a convex function, and let ``S`` be a nonempty compact polyhedral set in ``R^n``. Consider the problem to maximize ``f(\mathbf{x})`` subject to ``\mathbf{x} \in S``. An optimal solution ``\overline{\mathbf{x}}`` to the problem then exists, where ``\overline{\mathbf{x}}`` is an extreme point of ``S``.
"""

# ╔═╡ a0e01d7c-9ea7-46b1-a980-9d4e43315c13
cm"""
$(define("Quasiconvex Functions"))
Let ``f: S \rightarrow R``, where ``S`` is a nonempty convex set in ``R^n``. The function ``f`` is said to be quasiconvex if for each ``\mathbf{x}_1`` and ``\mathbf{x}_2 \in S``, the following inequality is true:
```math
f\left[\lambda \mathbf{x}_1+(1-\lambda) \mathbf{x}_2\right] \leq \max \left\{f\left(\mathbf{x}_1\right), f\left(\mathbf{x}_2\right)\right\} \text { for each } \lambda \in(0,1) .
```

The function ``f`` is said to be quasiconcave if ``-f`` is quasiconvex.
"""

# ╔═╡ 46187f5f-d96b-45bd-94a3-f9f1d7446960
cm"""
$(bth("3.5.2"))

Let ``f: S \rightarrow R`` where ``S`` is a nonempty convex set in ``R^n``. The function ``f`` is quasiconvex if and only if 
```math
S_\alpha=\{\mathbf{x} \in S: f(\mathbf{x}) \leq \alpha\}
``` 
is convex for each real number ``\alpha``.
"""

# ╔═╡ ff143a15-6e2e-4eb3-8faa-d683c417bf64
cm"""
$(bth("3.5.3 "))

Let ``S`` be a nonempty compact polyhedral set in ``R^n``, and let ``f: R^n \rightarrow R`` be quasiconvex and continuous on ``S``. Consider the problem 
```math
\max \quad f(\mathbf{x})\quad \text{subject to } \mathbf{x} \in S.
```
Then an optimal solution ``\overline{\mathbf{x}}`` to the problem exists, where ``\overline{\mathbf{x}}`` is an extreme point of ``S``.
"""

# ╔═╡ f26fed75-941b-4ab4-9720-405281d00170
cm"""
$(bth("3.5.4"))

Let ``S`` be a nonempty open convex set in ``\mathbb{R}^n``, and let ``f: S \rightarrow R`` be differentiable on ``S``. Then ``f`` is __quasiconvex__ if and only if either one of the following equivalent statements holds true:
1. If ``\mathbf{x}_1, \mathbf{x}_2 \in S`` and ``f\left(\mathbf{x}_1\right) \leq f\left(\mathbf{x}_2\right), \nabla f\left(\mathbf{x}_2\right)^t\left(\mathbf{x}_1-\mathbf{x}_2\right) \leq 0``.
2. If ``\mathbf{x}_1, \mathbf{x}_2 \in S`` and ``\nabla f\left(\mathbf{x}_2\right)^t\left(\mathbf{x}_1-\mathbf{x}_2\right)>0, f\left(\mathbf{x}_1\right)>f\left(\mathbf{x}_2\right)``.
"""

# ╔═╡ 2268c221-9e2e-403e-9729-53611b5bc63c
cm"""
$(define("Strictly Quasiconvex Functions"))
Let ``f: S \rightarrow R``, where ``S`` is a nonempty convex set in ``R^n``. The function ``f`` is said to be __strictly quasiconvex__ if for each ``\mathbf{x}_1, \mathbf{x}_2 \in S`` with ``f\left(\mathbf{x}_1\right) \neq f\left(\mathbf{x}_2\right)``, we have
```math
f\left[\lambda \mathbf{x}_1+(1-\lambda) \mathbf{x}_2\right]<\max \left\{f\left(\mathbf{x}_1\right), f\left(\mathbf{x}_2\right)\right\} \quad \text { for each } \lambda \in(0,1) .
```

The function ``f`` is called strictly quasiconcave if ``-f`` is strictly quasiconvex. Strictly quasiconvex functions are also sometimes referred to as semi-strictly quasiconvex, functionally convex, or explicitly quasiconvex.
"""

# ╔═╡ 6629edc2-72c3-4aa9-9e19-56b9a30adeb4
cm"""
$(post_img("https://www.dropbox.com/scl/fi/vpq7nl1ljuwo1bw2iokkh/fig3.11.png?rlkey=72axcy4a5hvd7pyyrj9sdk2u8&dl=1"))
"""

# ╔═╡ 375e3f9e-45fb-42da-b177-54c77c1081db
cm"""
$(bth("3.5.6"))

Let ``f: R^n \rightarrow R`` be strictly quasiconvex. Consider the problem to 
```math
\min f(\mathbf{x}) \text{ subject to }\mathbf{x} \in S,
```
where ``S`` is a nonempty convex set in ``R^n``. If ``\overline{\mathbf{x}}`` is a local optimal solution, ``\overline{\mathbf{x}}`` is also a global optimal solution.
"""

# ╔═╡ 7a0aa831-6dd4-4b77-891d-aa56892f7759
cm"""
$(bbl("Lemma","3.5.7"))

Let ``S`` be a nonempty convex set in ``R^n`` and let ``f: S \rightarrow R`` be strictly quasiconvex and lower semicontinuous. Then ``f`` is quasiconvex.
"""

# ╔═╡ 819fbc46-b696-4522-8fa4-ee68c89058b4
cm"""
$(define("Strongly Quasiconvex"))
Let ``S`` be a nonempty convex set in ``R^n``, and let ``f: S \rightarrow R``. The function ``f`` is said to be __strongly quasiconvex__ if for each ``\mathbf{x}_1, \mathbf{x}_2 \in S``, with ``\mathbf{x}_1 \neq \mathbf{x}_2``, we have
```math
f\left[\lambda \mathbf{x}_1+(1-\lambda) \mathbf{x}_2\right]<\max \left\{f\left(\mathbf{x}_1\right), f\left(\mathbf{x}_2\right)\right\}
```
for each ``\lambda \in(0,1)``. The function ``f`` is said to be __strongly quasiconcave__ if ``-f`` is strongly quasiconvex. 


The following statements hold true:
1. __strictly convex__ ``\Rightarrow`` __strongly quasiconvex__.
2. __strongly quasiconvex__ ``\Rightarrow`` __strictly quasiconvex__.
3. __strongly quasiconvex__ ``\Rightarrow`` __quasiconvex__ (even in the absence of any semicontinuity assumption.)
"""

# ╔═╡ 7a73b652-eae0-4ebe-b472-2ac1984607cf
cm"""
$(bth("3.5.9"))

Let ``f: R^n \rightarrow R`` be __strongly quasiconvex__. 
Consider the problem to 
```math
\min \; f(\mathbf{x}) \quad \text{subject to } \mathbf{x} \in S,
```
where ``S`` is a nonempty convex set in ``\mathbb{R}^n``. If ``\overline{\mathbf{x}}`` is a local optimal solution, ``\overline{\mathbf{x}}`` is the __unique__ global optimal solution.
"""

# ╔═╡ 434b0a5a-6937-4ed4-ac9a-8961d9578145
cm"""
$(define("Pseudoconvex Functions"))

Let ``S`` be a nonempty open set in ``R^n``, and let ``f: S \rightarrow R`` be differentiable on ``S``. The function ``f`` is said to be __pseudoconvex__ if for each ``\mathbf{x}_1, \mathbf{x}_2 \in S``, then
```math
\nabla f\left(\mathbf{x}_1\right)^t \left(\mathbf{x}_2-\mathbf{x}_1\right) \geq 0\quad \Rightarrow \quad f\left(\mathbf{x}_2\right) \geq f\left(\mathbf{x}_1\right);
```
or equivalently, 
```math
f\left(\mathbf{x}_2\right) <  f\left(\mathbf{x}_1\right)\quad\Rightarrow\quad  \nabla f\left(\mathbf{x}_1\right)^t\left(\mathbf{x}_2-\mathbf{x}_1\right)<0.
```

The function ``f`` is said to be __pseudoconcave__ if ``-f`` is pseudoconvex.



The function ``f`` is said to be __strictly pseudoconvex__ if for each distinct ``\mathbf{x}_1``, ``\mathbf{x}_2 \in S`` satisfying 
```math
\nabla f\left(\mathbf{x}_1\right)^t\left(\mathbf{x}_2-\mathbf{x}_1\right) \geq 0\quad\Rightarrow \quad f\left(\mathbf{x}_2\right) \geq f\left(\mathbf{x}_1\right);
```
or equivalently, if for each distinct ``\mathbf{x}_1, \mathbf{x}_2 \in S``,
```math
f\left(\mathbf{x}_2\right) \leq f\left(\mathbf{x}_1\right)\quad\Rightarrow\quad \nabla f\left(\mathbf{x}_1\right)^t\left(\mathbf{x}_2-\mathbf{x}_1\right)<0.
```
The function ``f`` is said to be __strictly pseudoconcave__ if ``-f`` is strictly pseudoconvex.
"""

# ╔═╡ f902fff5-6175-48d2-bebc-27c7e0f72d10
cm"""
$(post_img("https://www.dropbox.com/scl/fi/5q4mnjipmwm6s22i2aop5/fig3.12.png?rlkey=wrzzu5hdqlsgo5vy14avg5j9o&dl=1"))
"""

# ╔═╡ 960eff24-3070-4051-a4e0-f25a00b935b7
cm"""
$(bth("3.5.11"))

Let ``S`` be a nonempty open convex set in ``R^n``, and let ``f: S \rightarrow R`` be a differentiable pseudoconvex function on ``S``. Then ``f`` is both strictly quasiconvex and quasiconvex.
"""

# ╔═╡ e87f8c80-c756-4ff8-9e7d-ef35d8941afd
cm"""
$(bth("3.5.12"))

Let ``S`` be a nonempty open convex set in ``R^n``, and let ``f: S \rightarrow R`` be a differentiable strictly pseudoconvex function. Then ``f`` is strongly quasiconvex.

"""

# ╔═╡ 22b3e9fa-cf5d-43c3-b964-758b26e33468
cm"""
$(post_img("https://www.dropbox.com/scl/fi/x11nrxhomjrjsdvhvbakp/fig3.13.png?rlkey=8aow6soltthdbjipkx2w1d3jk&dl=1"))
"""

# ╔═╡ 4cb4d21c-464e-4675-b0b8-3872195ecc76
cm"""
$(define("Convexity at a Point"))
Let ``S`` be a nonempty convex set in ``R^n``, and let ``f: S \rightarrow R``. The following are relaxations of various forms of convexity presented in this chapter:

- Convexity at ``\overline{\mathbf{x}}``. The function ``f`` is said to be convex at ``\overline{\mathbf{x}} \in S`` if
```math
f[\lambda \overline{\mathbf{x}}+(1-\lambda) \mathbf{x}] \leq \lambda f(\overline{\mathbf{x}})+(1-\lambda) f(\mathbf{x})
```
for each ``\lambda \in(0,1)`` and each ``\mathbf{x} \in S``.

- Strict convexity at ``\overline{\mathbf{x}}``. The function ``f`` is said to be strictly convex at ``\overline{\mathbf{x}} \in S`` if
```math
f[\lambda \overline{\mathbf{x}}+(1-\lambda) \mathbf{x}]<\lambda f(\overline{\mathbf{x}})+(1-\lambda) f(\mathbf{x})
```
for each ``\lambda \in(0,1)`` and for each ``\mathbf{x} \in S, \mathbf{x} \neq \overline{\mathbf{x}}``.

- Quasiconvexity at ``\overline{\mathbf{x}}``. The function ``f`` is said to be quasiconvex at ``\overline{\mathbf{x}} \in S`` if
```math
f[\lambda \overline{\mathbf{x}}+(1-\lambda) \mathbf{x}] \leq \max \{f(\mathbf{x}), f(\overline{\mathbf{x}})\}
```
for each ``\lambda \in(0,1)`` and each ``\mathbf{x} \in S``.

- Strict quasiconvexity at ``\overline{\mathbf{x}}``. The function is said to be strictly quasiconvex at ``\overline{\mathbf{x}} \in S`` if
```math
f[\lambda \overline{\mathbf{x}}+(1-\lambda) \mathbf{x}]<\max \{f(\mathbf{x}), f(\overline{\mathbf{x}})\}
```
for each ``\lambda \in(0,1)`` and each ``\mathbf{x} \in S`` such that ``f(\mathbf{x}) \neq f(\overline{\mathbf{x}})``.

- Strong quasiconvexity at ``\overline{\mathbf{x}}``. The function ``f`` is said to be strongly quasiconvex at ``\overline{\mathbf{x}} \in S`` if
```math
f[\lambda \overline{\mathbf{x}}+(1-\lambda) \mathbf{x}]<\max \{f(\mathbf{x}), f(\overline{\mathbf{x}})\}
```
for each ``\lambda \in(0,1)`` and each ``\mathbf{x} \in S, \mathbf{x} \neq \overline{\mathbf{x}}``.

- Pseudoconvexity at ``\overline{\mathbf{x}}``. The function ``f`` is said to be pseudoconvex at ``\overline{\mathbf{x}} \in S`` if ``\nabla f(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}}) \geq 0`` for ``\mathbf{x} \in S`` implies that ``f(\mathbf{x}) \geq f(\overline{\mathbf{x}})``.

- Strict pseudoconvexity at ``\overline{\mathbf{x}}``. The function ``f`` is said to be strictly pseudoconvex at ``\overline{\mathbf{x}} \in S`` if ``\nabla f(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}}) \geq 0`` for ``\mathbf{x} \in S, \mathbf{x} \neq \overline{\mathbf{x}}``, implies that ``f(\mathbf{x})>f(\overline{\mathbf{x}})``.
"""

# ╔═╡ fc4c567e-3cc2-4b85-9294-eae1e2da69fd
cm"""
$(define(""))
Consider the problem of minimizing 
```math
f(\mathbf{x})\quad \text{over}\quad  \mathbb{R}^n,
```
and let ``\overline{\mathbf{x}} \in R^n``. 

⚫ If ``f(\overline{\mathbf{x}}) \leq f(\mathbf{x})`` for all ``\mathbf{x} \in R^n, \overline{\mathbf{x}}`` is called a __global minimum__. 

⚫ If there exists an ``\varepsilon`` neighborhood ``N_{\varepsilon}(\overline{\mathbf{x}})`` around ``\overline{\mathbf{x}}`` such that 
```math
f(\overline{\mathbf{x}}) \leq f(\mathbf{x})\quad \text{for each}\quad \mathbf{x} \in N_{\varepsilon}(\overline{\mathbf{x}}), \overline{\mathbf{x}}
``` 
is called a __local minimum__, while if ``f(\overline{\mathbf{x}}) < f(\mathbf{x})`` for all ``\mathbf{x} \in N_{\varepsilon}(\overline{\mathbf{x}}), \mathbf{x} \neq \overline{\mathbf{x}}``, for some ``\varepsilon > 0, \overline{\mathbf{x}}`` is called a __strict local minimum__. Clearly, a global minimum is also a local minimum.
"""

# ╔═╡ 17154407-9b5d-46be-92b9-2d02005e5c9c
cm"""

$(bth("4.1.2"))

Suppose that ``f: R^n \rightarrow R`` is differentiable at ``\overline{\mathbf{x}}``. If there is a vector ``\mathbf{d}`` such that ``\nabla f(\overline{\mathbf{x}})^t \mathbf{d} < 0``, there exists a ``\delta > 0`` such that ``f(\overline{\mathbf{x}}+\lambda \mathbf{d}) < f(\overline{\mathbf{x}})`` for each ``\lambda \in(0, \delta)``, so that ``\mathbf{d}`` is a descent direction of ``f`` at ``\overline{\mathbf{x}}``.
"""

# ╔═╡ 8eba2357-e825-4232-875a-f18be55bd38a
cm"""
$(bbl("Corollary",""))
Suppose that ``f: R^n \rightarrow R`` is differentiable at ``\overline{\mathbf{x}}``. If ``\overline{\mathbf{x}}`` is a local minimum, ``\nabla f(\overline{\mathbf{x}})=\mathbf{0}``.
"""

# ╔═╡ 4eebca6a-7e54-4d1c-84b2-cd893d4d2a3f
cm"""
$(bth("4.1.3 "))

Suppose that ``f: R^n \rightarrow R`` is twice differentiable at ``\overline{\mathbf{x}}``. If ``\overline{\mathbf{x}}`` is a local minimum, ``\nabla f(\overline{\mathbf{x}})=\mathbf{0}`` and ``\mathbf{H}(\overline{\mathbf{x}})`` is positive semidefinite.
"""

# ╔═╡ b6809792-4885-4bdd-9dea-d4e93ebe68c6
cm"""
$(bth("4.1.4"))

Suppose that ``f: R^n \rightarrow R`` is twice differentiable at ``\overline{\mathbf{x}}``. If ``\nabla f(\overline{\mathbf{x}})=0`` and ``\mathbf{H}(\overline{\mathbf{x}})`` is positive definite, ``\overline{\mathbf{x}}`` is a strict local minimum.
"""

# ╔═╡ 563d0e84-e2f3-49fe-b92e-6c0716ed8523
cm"""
$(bth("4.1.5"))

Let ``f: R^n \rightarrow R`` be pseudoconvex at ``\overline{\mathbf{x}}``. Then ``\overline{\mathbf{x}}`` is a global minimum if and only if ``\nabla f(\overline{\mathbf{x}})=\mathbf{0}``.
"""

# ╔═╡ 010e228c-7590-4178-99b6-bb50819dfe1c
cm"""
$(bth("4.1.6"))

Let ``f: R \rightarrow R`` be an infinitely differentiable univariate function. Then ``\bar{x} \in R`` is a local minimum if and only if either ``f^{(j)}(\bar{x})=0`` for all ``j=1,2, \ldots``, or else there exists an even ``n \geq 2`` such that ``f^{(n)}(\bar{x})>0`` while ``f^{(j)}(\bar{x})=0`` for all ``1 \leq j < n``, where ``f^{(j)}`` denotes the ``j`` th-order derivative of ``f``.
"""

# ╔═╡ d478679e-e6a0-48cc-8119-7e822f26dd02
min_P = min_latex();"";

# ╔═╡ ee88be92-5954-4f41-9b6d-c1db938b7368
cm"""
$(define("Feasible Solutions"))
Let ``S`` be a nonempty set in ``R^n``, and let ``\overline{\mathbf{x}} \in \mathrm{cl} S``. 

The __cone of feasible directions__ of ``S`` at ``\overline{\mathbf{x}}``, denoted by ``D``, is given by
```math
D=\{\mathbf{d}: \mathbf{d} \neq \mathbf{0}, \text { and } \overline{\mathbf{x}}+\lambda \mathbf{d} \in S \text { for all } \lambda \in(0, \delta) \text { for some } \delta >0\} .
```

Each nonzero vector ``\mathbf{d} \in D`` is called a __feasible direction__. Moreover, given a function ``f: R^n \rightarrow R``, the __cone of improving directions__ at ``\overline{\mathbf{x}}``, denoted by ``F``, is given by
```math
F=\{\mathbf{d}: f(\overline{\mathbf{x}}+\lambda \mathbf{d}) < f(\overline{\mathbf{x}}) \text { for all } \lambda \in(0, \delta) \text { for some } \delta >0\} .
```

Each direction ``\mathbf{d} \in F`` is called an __improving direction__, or a __descent direction__, of ``f``
"""

# ╔═╡ 6a7ff171-88af-431f-9ee8-4724f8eff61d
cm"""
$(bth("4.2.2"))

Consider the problem 

$(min_P)

where ``f: R^n \rightarrow R`` and ``S`` is a nonempty set in ``R^n``. Suppose that ``f`` is __differentiable__ at a point ``\overline{\mathbf{x}} \in S``. 

If ``\overline{\mathbf{x}}`` is a local optimal solution, ``F_0 \cap D=\varnothing``, where 
```math
F_0=\left\{\mathbf{d}: \nabla f(\overline{\mathbf{x}})^t \mathbf{d}<0\right\}
``` 
and ``D`` is the cone of feasible directions of ``S`` at ``\overline{\mathbf{x}}``. 

__Conversely__, suppose that ``F_0 \cap D=\varnothing``, ``f`` is pseudoconvex at ``\overline{\mathbf{x}}``, and that there exists an ``\varepsilon``-neighborhood ``N_{\varepsilon}(\overline{\mathbf{x}}), \varepsilon>0``, such that ``\mathbf{d}=(\mathbf{x}-\overline{\mathbf{x}}) \in D`` for any ``\mathbf{x} \in S \cap N_{\varepsilon}(\overline{\mathbf{x}})``. Then ``\overline{\mathbf{x}}`` is a local minimum of ``f``.
"""

# ╔═╡ b14b538e-e980-4e49-b9e0-74b43fe4620b
cm"""
$(post_img("https://www.dropbox.com/scl/fi/rxxwly6y6btoagymgyvda/fig4.3.png?rlkey=022v1odr7s2m4nblnmko57bqw&dl=1"))
"""

# ╔═╡ 2eea535f-8ab5-434b-b685-71539f73f062
cm"""
$(bbl("Remarks",""))
- ``F_0\subseteq F`` (Th 4.1.2).
- If ``d\in F``, then we have to have 
```math
\nabla f(\mathbf{x})^td \leq 0.
```
- ``F_0 \subseteq F \subseteq F_0^{\prime}=\left\{\mathbf{d} \neq \mathbf{0}: \nabla f(\overline{\mathbf{x}})^t \mathbf{d} \leq 0\right\}``. 
- The above inclusion can be __strict__. Meaning we may have ``\nabla f(\mathbf{x})^td=0``, and there might exist 
directions of motion that give descent or ascent, or even hold the value of ``f`` 
constant as we move away from ``\mathbf{x}``. 
"""

# ╔═╡ 8e1b7954-4159-4fd4-8c5a-2f1ac529e463
cm"""
$(bbl("Lemma","4.2.3"))

Consider a differentiable function ``f: R^n \rightarrow R``, and let ``F, F_0, F_0^{\prime}`` be as defined above. 

Then we have ``F_0 \subseteq F \subseteq F_0^{\prime}``. Moreover, if ``f`` is pseudoconvex at ``\overline{\mathbf{x}}, F=F_0``, and if ``f`` is strictly pseudoconcave at ``\overline{\mathbf{x}}, F=F_0^{\prime}``.
"""

# ╔═╡ 2f98c41a-3854-4a4a-b496-596d6f953d33
cm"""
$(bbl("Lemma","4.2.4"))

Consider the feasible region 
```math 
S=\left\{\mathbf{x} \in X: g_i(\mathbf{x}) \leq 0\text{ for } i=1, \ldots, m\right\}, \text{ where } X \text{ is a nonempty open set in }\mathbb{R}^n
``` 
and where ``g_i: R^n \rightarrow R`` for ``i=1, \ldots, m``. 

Given a feasible point ``\overline{\mathbf{x}} \in S``, let ``I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}`` be the index set for the __binding, or active, or tight constraints__, and assume that 
- ``g_i`` for ``i \in I`` are differentiable at ``\overline{\mathbf{x}}`` and that
- ``g_i`` for ``i \notin I`` are continuous at ``\overline{\mathbf{x}}``. 

Define the sets
```math
\begin{gathered}
G_0=\left\{\mathbf{d}: \nabla g_i(\overline{\mathbf{x}})^t \mathbf{d}<0 \text { for each } i \in I\right\} \\
G_0^{\prime}=\left\{\mathbf{d} \neq \mathbf{0}: \nabla g_i(\overline{\mathbf{x}})^t \mathbf{d} \leq 0 \text { for each } i \in I\right\}
\end{gathered}
```
Then we have
```math
G_0 \subseteq D \subseteq G_0^{\prime}
```

Moreover, 
- if ``g_i, i \in I``, are strictly pseudoconvex at ``\overline{\mathbf{x}}``, then ``D=G_0``; and 
- if ``g_i, i \in I``, are pseudoconcave at ``\overline{\mathbf{x}}``, then ``D=G_0^{\prime}``.
"""

# ╔═╡ e4c10925-07a1-4aa8-a921-bba1f49976ba
cm"""
$(bth("4.2.5"))

Consider Problem 

$(min_latex_gi())

where ``X`` is a nonempty open set in ``R^n, f: R^n \rightarrow R``, and ``g_i: R^n \rightarrow R``, for ``i =1, \ldots, m``. 

Let ``\overline{\mathbf{x}}`` be a feasible point, and denote 
```math
I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}.
```
Furthermore, suppose that 
- ``f`` and ``g_i`` for ``i \in I`` are differentiable at ``\overline{\mathbf{x}}`` and that 
- ``g_i`` for ``i \notin I`` are continuous at ``\overline{\mathbf{x}}``. 

If ``\overline{\mathbf{x}}`` is a __local optimal solution__, then
```math 
F_0 \cap G_0=\varnothing,
```
where
```math
\begin{array}{lcl}
F_0=\left\{\mathbf{d}: \nabla f(\overline{\mathbf{x}})^t \mathbf{d}<0\right\}
&\text{and}&G_0=\left\{\mathbf{d}: \nabla g_i(\overline{\mathbf{x}})^t \mathbf{d}<0\right.\quad  \text{for each } \left.i \in I\right\}.
\end{array}
```

__Conversely__, if 
- ``F_0 \cap G_0=\varnothing``, and if 
- ``f`` is pseudoconvex at ``\overline{\mathbf{x}}`` and 
- ``g_i, i \in I``, are strictly pseudoconvex over some ``\varepsilon``-neighborhood of ``\overline{\mathbf{x}}``.

Then ``\overline{\mathbf{x}}`` is a local minimum.
"""

# ╔═╡ 524c222e-f9e1-4421-88b0-072cb6b74d95
cm"""
$(bbl("Remarks",""))
- Observe that under the __converse hypothesis__ of Theorem 4.2.5, and assuming that ``g_i, i \notin I``, are continuous at ``\overline{\mathbf{x}}``, we have, noting (4.9a),
```math
\overline{\mathbf{x}} \text { is a local minimum } \Leftrightarrow F_0 \cap D=\varnothing \Leftrightarrow F_0 \cap G_0=\varnothing .
```
- If ``\overline{\mathbf{x}}`` is a local minimum, then clearly we must have 
```math 
F \cap D=\varnothing.
```
- However, the __converse__ is not necessarily true. That is, if ``F \cap D=\varnothing``, this does not necessarily imply that ``\overline{\mathbf{x}}`` is a local minimum. 
"""

# ╔═╡ f68e53bc-5c85-4b4b-95b7-d94af1a7f124
cm"""
$(ex("Example","4.2.6"))
```math
\begin{aligned}
\operatorname{Minimize}\left(x_1-3\right)^2 & +\left(x_2-2\right)^2 \\
\text { subject to } x_1^2+x_2^2 & \leq 5 \\
x_1+x_2 & \leq 3 \\
x_1 & \geq 0 \\
x_2 & \geq 0 .
\end{aligned}
```
"""

# ╔═╡ ee994ad9-5fc9-4722-9201-536996e19b4c
cm"""
$(bth("4.2.8(Fritz John Necessary Conditions)"))

Let ``X`` be a nonempty open set in ``\mathbb{R}^n`` and let ``f: \mathbb{R}^n \rightarrow R``, and ``g_i: \mathbb{R}^n \rightarrow R`` for ``i= 1, \ldots, m``.

Consider Problem P 

$(min_latex_gi())

Let ``\overline{\mathbf{x}}`` be a __feasible solution__, and denote 
```math
I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}.
```
Furthermore, suppose that 
- ``f`` and ``g_i`` for ``i \in I`` are differentiable at ``\overline{\mathbf{x}}`` and that 
- ``g_i`` for ``i \notin I`` are continuous at ``\overline{\mathbf{x}}``. 

If ``\overline{\mathbf{x}}`` solves Problem P locally, there exist scalars ``u_0`` and ``u_i`` for ``i \in I`` such that
```math
\begin{aligned}
u_0 \nabla f(\overline{\mathbf{x}})+\sum_{i \in I} u_i \nabla g_i(\overline{\mathbf{x}}) & =0 \\
u_0, u_i & \geq 0, \quad \text { for } i \in I \\
\left(u_0, \mathbf{u}_I\right) & \neq(0, \mathbf{0}), 
\end{aligned}
```
where ``\mathbf{u}_I`` is the vector whose components are ``u_i`` for ``i \in I``. Furthermore, if ``g_i`` for ``i \notin I`` are also differentiable at ``\overline{\mathbf{x}}``, the foregoing conditions can be written in the following equivalent form:
```math
\begin{aligned}
u_0 \nabla f(\overline{\mathbf{x}})+\sum_{i=1}^m u_i \nabla g_i(\overline{\mathbf{x}}) & =\mathbf{0} & & \\
u_i g_i(\overline{\mathbf{x}}) & =0 & & \text { for } i=1, \ldots, m \\
u_0, u_i & \geq 0 & & \text { for } i=1, \ldots, m \\
\left(u_0, \mathbf{u}\right) & \neq(0, \mathbf{0}), & &
\end{aligned}
```
where ``\mathbf{u}`` is the vector whose components are ``u_i`` for ``i=1, \ldots, m``.
"""

# ╔═╡ 5196f462-29fe-4f16-8e8d-c9b233a7dca8
cm"""
$(ex("Example","4.2.9"))
```math
\begin{aligned}
\operatorname{Minimize}\left(x_1-3\right)^2 & +\left(x_2-2\right)^2 \\
\text { subject to } x_1^2+x_2^2 & \leq 5 \\
x_1+2 x_2 & \leq 4 \\
-x_1 & \leq 0 \\
-x_2 & \leq 0
\end{aligned}
```
"""

# ╔═╡ b7f2b3d3-8aa5-4a1e-97e5-a10f17aebbaf
cm"""
$(ex("Example","4.2.10 "))
```math
\begin{aligned} & \operatorname{Minimize}-x_1 \\ & \text { subject to } x_2-\left(1-x_1\right)^3 \leq 0 \\ & -x_2 \leq 0\end{aligned}
```
"""

# ╔═╡ 10f03ebd-0da6-4ecf-afdd-e524a9ddf530
cm"""
$(bth("4.2.12 (Fritz John Sufficient Conditions)"))

Let ``X`` be a nonempty open set in ``R^n``, and let ``f: R^n \rightarrow R`` and ``g_i: R^n \rightarrow R`` for ``i= 1, \ldots, m``. 
Consider Problem P

$(min_latex_gi())


Let ``\overline{\mathbf{x}}`` be a FJ solution and denote 
```math 
I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}.
```

Define ``S`` as the relaxed feasible region for Problem ``P`` in which the nonbinding constraints are dropped.

- __(a)__ If there exists an ``\varepsilon``-neighborhood ``N_{\varepsilon}(\overline{\mathbf{x}}), \varepsilon>0``, such that 
   - ``f`` is pseudoconvex over ``N_{\mathcal{E}}(\overline{\mathbf{x}}) \cap S``, and 
   - ``g_i, i \in I``, are strictly pseudoconvex over ``N_{\varepsilon}(\overline{\mathbf{x}}) \cap S``. 

$(add_space(10)) Then ``\overline{\mathbf{x}}`` is a __local minimum for Problem P__.

- __(b)__ If 
   - ``f`` is pseudoconvex at ``\overline{\mathbf{x}}``, and if 
   - ``g_i, i \in I``, are both strictly pseudoconvex and quasiconvex at ``\overline{\mathbf{x}}``. 

$(add_space(10))Then ``\overline{\mathbf{x}}`` is a __global optimal solution for Problem P__. 

<div style="padding-left:4ch;">

In particular, if these generalized convexity assumptions hold true only by restricting the domain of ``f`` to ``N_{\varepsilon}(\overline{\mathbf{x}})`` for some ``\varepsilon> 0, \overline{\mathbf{x}}`` is a local minimum for Problem P .
</div>
"""

# ╔═╡ 75d8b9f6-8d89-49b9-b8a4-e7cbbd066a31
cm"""
$(bbl("Remarks","Issues"))

- ``\overline{x}`` is an FJ point ``\Longleftrightarrow`` ``F_0\cap G_0 = \varnothing``.

$(post_img("https://www.dropbox.com/scl/fi/bk367fwcohea2rm1pk1bl/fig4.9.png?rlkey=8d1eny6oybrjfr6zsfakt1tb4&dl=1"))
"""

# ╔═╡ fd688b17-447f-4aa6-8584-20b04ef1f822
cm"""
$(ex("Example", "LP"))

Minimize -x₁

Subject to:
- x₁ + x₂ - 1 ≤ 0
- x₂ ≥ 0

Enter values to select a point and see the gradients!

## Select a point (x₁, x₂)

"""

# ╔═╡ 1791db8c-eb23-450c-b79d-82124a027cee
cm"""
$(bbl("Constraint Qualification (CQ)"))

A __(CQ)__ is an assumption made about constraint functions (equality and inequality) that, when satisfied at a local minimizer ``\overline{x}``, ensures that ``\overline{x}`` is a KKT point.

For example ``G_0\neq \varnothing`` is (CQ).
"""

# ╔═╡ 2b18aae5-24c1-4f3b-8926-9993c3f47db2
cm"""
$(bth("4.2.13 (Karush-Kuhn-Tucker Necessary Conditions)"))

Let ``X`` be a nonempty open set in ``R^n``, and let ``f: R^n \rightarrow R`` and ``g_i: R^n \rightarrow R`` for ``i= 1, \ldots, m``. Consider the Problem P 
$(min_latex_gi())

Let ``\overline{\mathbf{x}}`` be a __feasible solution__, and denote 
```math
I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}.
```
Suppose that 
- ``f`` and ``g_i`` for ``i \in I`` are differentiable at ``\overline{\mathbf{x}}`` and that 
- ``g_i`` for ``i \notin I`` are continuous at ``\overline{\mathbf{x}}``.

Furthermore, suppose that 
```math
\nabla g_i(\overline{\mathbf{x}}) \text{ for } i \in I \text{ are linearly independent}.\qquad (CQ')
```

If ``\overline{\mathbf{x}}`` solves Problem P locally, there exist scalars ``u_i`` for ``i \in I`` such that
```math
\begin{aligned}
\nabla f(\overline{\mathbf{x}})+\sum_{i \in I} u_i \nabla g_i(\overline{\mathbf{x}}) & =0 \\
u_i & \geq 0 \quad \text { for } i \in I .
\end{aligned}
```

In addition to the above assumptions, if ``g_i`` for each ``i \notin I`` is also differentiable at ``\overline{\mathbf{x}}``, the foregoing conditions can be written in the following equivalent form:
```math
\begin{aligned}
\nabla f(\overline{\mathbf{x}})+\sum_{i=1}^m u_i \nabla g_i(\overline{\mathbf{x}}) & =0 & & \\
u_i g_i(\overline{\mathbf{x}}) & =0 & & \text { for } i=1, \ldots, m \\
u_i & \geq 0 & & \text { for } i=1, \ldots, m
\end{aligned}
```
"""

# ╔═╡ 55c127e9-9647-4f9f-b9fd-4d9752071a90
cm"""
$(bbl("Remarks",""))
- Together, these PF, DF, and CS conditions are called the __KKT conditions__. 
- Any point ``\overline{x}`` for which there exist Lagrangian (or Lagrange) multipliers ``\overline{u}`` such that ``(\overline{x},\overline{u})`` 
satisfies the KKT conditions is called a __KKT point__.
"""

# ╔═╡ 38143cd0-3552-4083-954a-bb10efb8c090
cm"""
$(post_img("https://www.dropbox.com/scl/fi/way294o6eux2ak3qrduuy/fig4.10.png?rlkey=jvrcim4ig2tsr25w001119rh0&dl=1"))
"""

# ╔═╡ 1552ddcb-1e1d-43de-a9e3-f4a7c186970b
cm"""
$(bbl("4.2.15 (KKT Conditions and First-Order LP Approximations)"))

Let ``X`` be a nonempty open set in ``R^n``, and let ``f: R^n \rightarrow R`` and ``g_i: R^n \rightarrow R, i=1, \ldots``, ``m`` be differentiable functions. 
Consider Problem P , 
$(min_latex_gi())

Let ``\overline{\mathbf{x}}`` be a feasible solution, and denote ``I =\left\{i: g_i(\overline{\mathbf{x}})=0\right\}``. 

Define 
```math
F_0=\left\{\mathbf{d}: \nabla f(\overline{\mathbf{x}})^t \mathbf{d}<0\right\}
``` 
and 
```math
G_0^{\prime}=\left\{\mathbf{d} \neq \mathbf{0}: \nabla g_i(\overline{\mathbf{x}})^t \mathbf{d} \leq 0\right., \text{ for each  }i \in I\}
``` 
as before, and let 
```math
G^{\prime}=\left\{\mathbf{d}: \nabla g_i(\overline{\mathbf{x}})^t \mathbf{d} \leq 0\right. \text{ for each } \left.i \in I\right\}= G_0^{\prime} \cup\{\mathbf{0}\}.
```
Then ``\overline{\mathbf{x}}`` is a KKT solution if and only if ``F_0 \cap G^{\prime}=\varnothing``, which is equivalent to ``F_0 \cap G_0^{\prime}=\varnothing``. Furthermore, consider the first-order linear programming approximation to Problem P:
```math
\begin{gathered}
\operatorname{LP}(\overline{\mathbf{x}}): \operatorname{Minimize}\left\{f(\overline{\mathbf{x}})+\nabla f(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}}): g_i(\overline{\mathbf{x}})+\nabla g_i(\overline{\mathbf{x}})^t(\mathbf{x}-\overline{\mathbf{x}}) \leq 0\right. \\
\text { for } i=1, \ldots, m\} .
\end{gathered}
```

Then, ``\overline{\mathbf{x}}`` is a KKT solution if and only if ``\overline{\mathbf{x}}`` solves ``\mathrm{LP}(\overline{\mathbf{x}})``.
"""

# ╔═╡ bccbbf3a-1cf3-41e2-bc0d-b32e122421ec
cm"""
$(bbl("Proof","First Part"))

```math
\text{A feasible solution } x \text{ is a KKT point iff there exist multipliers } 
\{u_i,\, i \in J\} 
\text{ satisfying:}
```
```math
\sum_{i \in J} u_i \nabla g_i(x) = -\nabla f(x), \quad u_i \ge 0 \ \forall i \in J.
```

By Farkas’s Lemma (see, e.g., Corollary 2 to Theorem 2.7.3), this holds true if and only if there does **not** exist a vector \( d \) such that:
```math
\nabla g_i(x)^{T} d < 0 \quad \forall i \in J, 
\quad \text{and} \quad 
\nabla f(x)^{T} d < 0.
```

Hence, \( x \) is a KKT point if and only if 
```math
F_y \cap G' = \emptyset.
```
Clearly, this is equivalent to 
```math
F_y \cap G_0 = \emptyset.
```


"""

# ╔═╡ f4969ec2-4f5b-4893-9385-d74b038f4738
cm"""
$(bth("4.2.16 (Karush–Kuhn–Tucker Sufficient Conditions)"))

Let `` X `` be a nonempty open set in `` \mathbb{R}^n ``, and let 
`` f: \mathbb{R}^n \to \mathbb{R} `` and `` g_i: \mathbb{R}^n \to \mathbb{R} `` 
for `` i = 1, \ldots, m ``. 

Consider Problem P

$(min_latex_gi())

Let `` \bar{x} `` be a KKT solution, and denote 
`` I = \{ i : g_i(\bar{x}) = 0 \} ``.  
Define `` S `` as the relaxed feasible region for Problem `` P `` 
in which the constraints that are not binding at `` \bar{x} `` are dropped.  

Then:
<ul>

<li> 

**(a)**
If there exists an `` \varepsilon ``-neighborhood `` N_{\varepsilon}(\bar{x}) `` about 
`` \bar{x} ``, `` \varepsilon > 0 ``, such that `` f `` is pseudoconvex over 
`` N_{\varepsilon}(\bar{x}) \cap S `` and `` g_i,\, i \in I, `` are differentiable at `` \bar{x} `` 
and are quasiconvex over `` N_{\varepsilon}(\bar{x}) \cap S ``,  
then `` \bar{x} `` is a local minimum for Problem `` P ``.

</li>

<li> 

**(b)**  If `` f `` is pseudoconvex at `` \bar{x} ``, and if `` g_i,\, i \in I, `` are differentiable 
and quasiconvex at `` \bar{x} ``, then `` \bar{x} `` is a global optimal solution to Problem `` P ``.  
In particular, if this assumption holds true with the domain of the feasible restriction 
to `` N_{\varepsilon}(\bar{x}) ``, for some `` \varepsilon > 0 ``, then `` \bar{x} `` is a local minimum 
for `` P ``.

</li>

</ul>

"""

# ╔═╡ 6e386883-6678-4ed6-9c09-bd9a41201783
cm"""
$(bbl("Remark","💣"))
__KKT conditions are not necessary for optimality for convex programming problems.__
"""

# ╔═╡ 995242c4-a629-4491-9db5-645730d6b9bb
cm"""
$(ex("Example",""))

**Problem:**

Minimize ``x_1``

subject to:
- ``(x_1 - 1)^2 + (x_2 - 1)^2 \leq 1`` (Circle centered at (1,1))
- ``(x_1 - 1)^2 + (x_2 + 1)^2 \leq 1`` (Circle centered at (1,-1))

The feasible region is the **lens-shaped intersection** of two unit circles.
"""

# ╠═╡ input

# ╔═╡ 09a738b0-d156-4123-a5c1-a3fc349d680a
cm"""
We study the Problem

$(eql_latex_gi())

"""

# ╔═╡ 402f4f37-fc08-4cb7-9b64-ffd028711de4
cm"""
$(bth("4.3.1"))

Let ``X`` be a nonempty open set in ``\mathbb{R}^n`` and let ``f: \mathbb{R}^n \rightarrow R``, ``g_i: \mathbb{R}^n \rightarrow R`` for ``i= 1, \ldots, m`` and ``h_i: \mathbb{R}^n \rightarrow R`` for ``i= 1, \ldots, l``.

Consider Problem P 

$(eql_latex_gi())

Let ``\overline{\mathbf{x}}`` be a __local optimal solution__, and denote 
```math
I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}.
```
Furthermore, suppose that 
- ``f`` and ``g_i`` for ``i \in I`` are differentiable at ``\overline{\mathbf{x}}`` and that 
- ``g_i`` for ``i \notin I`` are continuous at ``\overline{\mathbf{x}}``. 
- ``h_i`` for all ``i=1,\cdots l`` are __continuously differentiable__ at ``\overline{\mathbf{x}}``. 
- ``\{\nabla h_i(\overline{\mathbf{x}})\}`` for ``i=1,\cdots l`` is __linearly independent__. 

``\underline{\large{\text{THEN}}}``
```math
F_0 \cap G_0 \cap H_0 = \varnothing.
```
where
```math
\begin{array}{lcl}
F_0 &=& \left\{d: \nabla f(\overline{x})^td < 0\right\}\\
G_0 &=& \left\{d: \nabla g_i(\overline{x})^td < 0 \text{ for } i\in I\right\}\\
H_0 &=& \left\{d: \nabla h_i(\overline{x})^td = 0 \text{ for } i =1,\cdots, l\right\}\\
\end{array}
```
Conversely, suppose that 
```math
F_0 \cap G_0 \cap H_0=\varnothing.
```
and if 
- ``f`` is __pseudoconvex__ at ``\overline{\mathbf{x}},`` 
- ``g_i`` for ``i \in I`` are __strictly pseudoconvex__ over some ``\varepsilon``-neighborhood of ``\overline{\mathbf{x}}``; and if 
- ``h_i`` for ``i=1, \ldots, \ell`` are __affine__. 


``\underline{\large{\text{THEN}}}``

``\overline{\mathbf{x}}`` is a local optimal solution.

"""

# ╔═╡ a8133296-0ec7-4504-96bb-39e443193c7f
cm"""
$(bth("4.3.1  (Fritz John Necessary Conditions)"))

Let ``X`` be a nonempty open set in ``\mathbb{R}^n`` and let ``f: \mathbb{R}^n \rightarrow R``, ``g_i: \mathbb{R}^n \rightarrow R`` for ``i= 1, \ldots, m`` and ``h_i: \mathbb{R}^n \rightarrow R`` for ``i= 1, \ldots, l``.

Consider Problem P 

$(eql_latex_gi())

Let ``\overline{\mathbf{x}}`` be a __local optimal solution__, and denote 
```math
I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}.
```
Furthermore, suppose that 
- ``f`` and ``g_i`` for ``i \in I`` are differentiable at ``\overline{\mathbf{x}}`` and that 
- ``g_i`` for ``i \notin I`` are continuous at ``\overline{\mathbf{x}}``. 
- ``h_i`` for all ``i=1,\cdots l`` are __continuously differentiable__ at ``\overline{\mathbf{x}}``. 
- ``\{\nabla h_i(\overline{\mathbf{x}})\}`` for ``i=1,\cdots l`` is __linearly independent__. 

``\underline{\large{\text{THEN}}}``
<div style="text-align: center;font-weight:700;margin-bottom:40px;">

There exist scalars ``u_0``, ``u_i, \text{ for } i \in I``,and  ``v_i \text{ for } i=1,\cdots, l`` such that
</div>

```math
\begin{aligned}
u_0 \nabla f(\overline{\mathbf{x}})+\sum_{i \in I} u_i \nabla g_i(\overline{\mathbf{x}}) +\sum_{i=1}^{l} v_i \nabla h_i(\overline{\mathbf{x}})& =0 \\
u_0, u_i & \geq 0, \quad \text { for } i \in I \\
\left(u_0, \mathbf{u}_I, \mathbf{v}\right) & \neq(0, \mathbf{0}, \mathbf{0}), 
\end{aligned}
```
where ``\mathbf{u}_I`` is the vector whose components are ``u_i`` for ``i \in I`` and ``\mathbf{v}=(v_1,\cdots,v_l)``. Furthermore, if ``g_i`` for ``i \notin I`` are also differentiable at ``\overline{\mathbf{x}}``, the foregoing conditions can be written in the following equivalent form:
```math
\begin{aligned}
u_0 \nabla f(\overline{\mathbf{x}})+\sum_{i=1}^m u_i \nabla g_i(\overline{\mathbf{x}}) +\sum_{i=1}^{l} v_i \nabla h_i(\overline{\mathbf{x}})& =\mathbf{0} & & \\
u_i g_i(\overline{\mathbf{x}}) & =0 & & \text { for } i=1, \ldots, m \\
u_0, u_i & \geq 0 & & \text { for } i=1, \ldots, m \\
\left(u_0, \mathbf{u}, \mathbf{v}\right) & \neq(0, \mathbf{0}, \mathbf{0}), & &
\end{aligned}
```
where ``\mathbf{u}^t=(u_1,\cdots,u_m)`` and ``\mathbf{v}^t=(v_1,\cdots,v_l)``.

"""

# ╔═╡ 3758e93e-ed91-478e-b431-caccf55bdda0
cm"""
$(bth("4.3.6 (Fritz John Sufficient Conditions)"))

Let ``X`` be a nonempty open set in ``R^n``, and let ``f: R^n \rightarrow R``, ``g_i: R^n \rightarrow R`` for ``i= 1, \ldots, m``, and ``h_i: R^n \rightarrow R`` for ``i= 1, \ldots, l``
Consider Problem P

$(eql_latex_gi())


Let ``\overline{\mathbf{x}}`` be a FJ solution and denote 
```math 
I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}, \quad S=\{x\in X: g_i(x)\leq 0\text{ for } i\in I, h_i(x)=0 \text{ for } i =1,\cdots, l\}
```
If
- ``h_i`` for ``i=1,\cdots,l`` are affine and ``\nabla h_i(\overline{x})``, ``i=1,\cdots l`` are linearly independent,
- ``f`` is pseudoconvex on ``S\cap N_{\varepsilon}(\overline{\mathbf{x}})`` for some nieghbothood ``N_{\varepsilon}(\overline{x})`` of ``\overline{x}`` and ``\varepsilon>0``, and if 
- ``g_i, i \in I``, are both strictly pseudoconvex on ``S\cap N_{\varepsilon}(\overline{\mathbf{x}})``. 

Then ``\overline{\mathbf{x}}`` is a __local optimal solution for Problem P__. 

"""

# ╔═╡ 0f4d581b-9e6b-47d7-b014-d0d161cf84a9
cm"""
$(bth("4.3.7 (Karush-Kuhn-Tucker Necessary Conditions)"))

Let ``X`` be a nonempty open set in ``R^n``, and let ``f: R^n \rightarrow R``, ``g_i: R^n \rightarrow R`` for ``i= 1, \ldots, m``, and ``h_i: R^n \rightarrow R`` for ``i= 1, \ldots, l``. Consider the Problem P 
$(eql_latex_gi())

Let ``\overline{\mathbf{x}}`` be a __feasible solution__, and denote 
```math
I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}.
```
Suppose that 
- ``f`` and ``g_i`` for ``i \in I`` are differentiable at ``\overline{\mathbf{x}}`` and that 
- ``g_i`` for ``i \notin I`` are continuous at ``\overline{\mathbf{x}}``.
- ``h_i`` for ``i=1,\cdots l`` are continuously differentiable at ``\overline{\mathbf{x}}``.

Furthermore, suppose that ``\overline{x}`` is a __regular point``, i.e.
```math
\nabla g_i(\overline{\mathbf{x}}) \text{ for } i \in I \text{ and } \nabla h_i(\overline{\mathbf{x}}) \text{ for } i=1, \cdots, l \text{ are linearly independent}.\qquad (CQ')
```

If ``\overline{\mathbf{x}}`` solves Problem P locally, 

``\underline{\large{\text{THEN}}}``

there exist __unique__ scalars ``u_i`` for ``i \in I`` and ``v_i`` for ``i=1,\cdots,l`` such that
```math
\begin{aligned}
\nabla f(\overline{\mathbf{x}})+\sum_{i \in I} u_i \nabla g_i(\overline{\mathbf{x}}) +\sum_{i=1}^{l} v_i \nabla h_i(\overline{\mathbf{x}}) & =0 \\
u_i & \geq 0 \quad \text { for } i \in I .
\end{aligned}
```

In addition to the above assumptions, if ``g_i`` for each ``i \notin I`` is also differentiable at ``\overline{\mathbf{x}}``, the foregoing conditions can be written in the following equivalent form:
```math
\begin{aligned}
\nabla f(\overline{\mathbf{x}})+\sum_{i=1}^m u_i \nabla g_i(\overline{\mathbf{x}})
+\sum_{i=1}^m v_i \nabla h_i(\overline{\mathbf{x}})& =0 & & \\
u_i g_i(\overline{\mathbf{x}}) & =0 & & \text { for } i=1, \ldots, m \\
u_i & \geq 0 & & \text { for } i=1, \ldots, m
\end{aligned}
```
"""

# ╔═╡ fa6c7bd2-e304-4555-833d-989d6414964a
cm"""
$(bth("4.3.8 (Karush-Kuhn-Tucker Necessary Conditions)"))

Let ``X`` be a nonempty open set in ``R^n``, and let ``f: R^n \rightarrow R``, ``g_i: R^n \rightarrow R`` for ``i= 1, \ldots, m`` and ``h_i: R^n \rightarrow R`` for ``i= 1, \ldots, l``. Consider the Problem P 
$(eql_latex_gi())

Let ``\overline{\mathbf{x}}`` be a __feasible solution__, and denote 
```math
I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}.
```
Suppose that the __KKT__ conditions hold true at ``\overline{\mathbf{x}}``; that is, there exist scalars ``\bar{u}_i \geq 0`` for ``i \in I`` and ``\bar{v}_i`` for ``i= 1, \ldots, \ell`` such that
```math
\nabla f(\overline{\mathbf{x}})+\sum_{i \in I} \bar{u}_i \nabla g_i(\overline{\mathbf{x}})+\sum_{i=1}^{\ell} \bar{v}_i \nabla h_i(\overline{\mathbf{x}})=\mathbf{0} .
```

Let ``J=\left\{i: \bar{v}_i>0\right\}`` and ``K=\left\{i: \widetilde{v}_i<0\right\}``. Further, suppose that 
- ``f`` is pseudoconvex at ``\overline{\mathbf{x}}``, 
- ``g_i`` is quasiconvex at ``\overline{\mathbf{x}}`` for ``i \in I,`` 
- ``h_i`` is quasiconvex at ``\overline{\mathbf{x}}`` for ``i \in J``, and ``h_i`` is quasiconcave at ``\overline{\mathbf{x}}`` for ``i \in K``. 

``\underline{\large{\text{THEN}}}``

``\overline{\mathbf{x}}`` is a __global optimal solution__ to Problem P. 

In particular, if the generalized convexity assumptions on the objective and constraint functions are restricted to the domain ``N_{\varepsilon}(\overline{\mathbf{x}})`` for some ``\varepsilon>0, \overline{\mathbf{x}}`` is a __local minimum__ for P .
"""

# ╔═╡ 07cb2ece-065f-4790-b8e2-ef164907d69d
cm"""
$(ex("Example","4.3.4"))
```math
\begin{aligned}
\operatorname{Minimize}\left(x_1-3\right)^2 & +\left(x_2-2\right)^2 \\
\text { subject to } x_1^2+x_2^2 & \leq 5 \\
-x_1 & \leq 0 \\
-x_2 & \leq 0 \\
x_1+2 x_2 & =4
\end{aligned}
```
"""

# ╔═╡ e15eb170-02e9-42e9-b8d5-35d9284fd30e
cm"""
Consider the problem (P)

$(eql_latex_gi())
Let 
```math
S=\left\{\mathbf{x}: g_i(\mathbf{x}) \leq 0\right. \text{ for }i=1, \ldots, m, h_i(\mathbf{x})=0\text{ for }i=1, \ldots, \ell\text{ and }\left.\mathbf{x} \in X\right\}.
```


and let's define the Lagrangiun function for this problem.
```math
\phi(\mathbf{x}, \mathbf{u}, \mathbf{v})=f(\mathbf{x})+\sum_{i=1}^m u_i g_i(\mathbf{x})+\sum_{i=1}^{\ell} v_i h_i(\mathbf{x}).
```

- Now, let ``\overline{\mathbf{x}}`` be a KKT point for Problem P , with associated Lagrangian multipliers ``\overline{\mathbf{u}}`` and ``\overline{\mathbf{v}}`` corresponding to the inequality and equality constraints, respectively. 
- Conditioned on ``\overline{\mathbf{u}}`` and ``\overline{\mathbf{v}}``, define the restricted Lagrangian function
```math
L(\mathbf{x}) \equiv \phi(\mathbf{x}, \overline{\mathbf{u}}, \overline{\mathbf{v}})=f(\mathbf{x})+\sum_{i \in I} \bar{u}_i g_i(\mathbf{x})+\sum_{i=1}^{\ell} \bar{v}_i h_i(\mathbf{x}),
```
where ``I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}`` is the index set of the binding inequality constraints at ``\overline{\mathbf{x}}``.
"""

# ╔═╡ 4c0fffcb-327a-4a45-a8b6-cea81005ab2b
cm"""
$(bbl("Lemma","4.4.1"))

Consider Problem P 
$(eql_latex_gi())

where the objective and constraint defining functions are all __twice differentiable__, and where ``X`` is a nonempty, open set in ``\mathbb{R}^n``.

Suppose that ``\overline{\mathbf{x}}`` is a __KKT point__ for Problem P with Lagrangian multipliers ``\overline{\mathbf{u}}`` and ``\overline{\mathbf{v}}`` associated with the inequality and equality constraints, respectively. 

Define the restricted Lagrangian function ``L`` as

```math
L(\mathbf{x}) \equiv \phi(\mathbf{x}, \overline{\mathbf{u}}, \overline{\mathbf{v}})=f(\mathbf{x})+\sum_{i \in I} \bar{u}_i g_i(\mathbf{x})+\sum_{i=1}^{\ell} \bar{v}_i h_i(\mathbf{x}),
```

and denote its Hessian by ``\nabla^2 L``.

- (a.) If ``\nabla^2 L`` is __positive semidefinite__ for all ``\mathbf{x} \in S``. __THEN__ ``  \overline{\mathbf{x}}`` is a global minimum for Problem P . On the other hand, if ``\nabla^2 L`` is positive semidefinite for all ``\mathbf{x} \in S \cap N_{\varepsilon}(\overline{\mathbf{x}})`` for some ``\varepsilon``-neighborhood ``N_{\varepsilon}(\overline{\mathbf{x}})`` about ``\overline{\mathbf{x}}, \varepsilon>0, \overline{\mathbf{x}}`` is a local minimum for Problem P.

- (b.) If ``\nabla^2 L(\overline{\mathbf{x}})`` is positive definite, ``\overline{\mathbf{x}}`` is a strict local minimum for Problem P.
"""

# ╔═╡ d19630de-5c1d-47b5-b91e-359ae2898453
cm"""
$(bth("4.4.2 Theorem (KKT Second-Order Sufficient Conditions)"))

Consider Problem P 
$(eql_latex_gi())

where the objective and constraint defining functions are all __twice differentiable, and where ``X`` is a nonempty, open set in ``\mathbb{R}^n``__. 

Let ``\overline{\mathbf{x}}`` be a __KKT point__ for Problem P , with Lagrangian multipliers ``\overline{\mathbf{u}}`` and ``\overline{\mathbf{v}}`` associated with the inequality and equality constraints, respectively. 

Let 
```math
I =\left\{i: g_i(\overline{\mathbf{x}})=0\right\},\quad \text{and denote } I^{+}=\left\{i \in I: \bar{u}_i>0\right\}\text { and } I^0=\left\{i \in I: \bar{u}_i=0\right\}.
```
( ``I^{+}`` and ``I^0`` are sometimes referred to as the set of strongly active and weakly active constraints, respectively.) 

Define the restricted Lagrangian function ``L(\mathbf{x})`` as 

```math
L(\mathbf{x}) \equiv \phi(\mathbf{x}, \overline{\mathbf{u}}, \overline{\mathbf{v}})=f(\mathbf{x})+\sum_{i \in I} \bar{u}_i g_i(\mathbf{x})+\sum_{i=1}^{\ell} \bar{v}_i h_i(\mathbf{x}),
```

and denote its Hessian at ``\overline{\mathbf{x}}`` by
```math
\nabla^2 L(\overline{\mathbf{x}}) \equiv \nabla^2 f(\overline{\mathbf{x}})+\sum_{i \in I} \bar{u}_i \nabla^2 g_i(\overline{\mathbf{x}})+\sum_{i=1}^{\ell} \bar{v}_i \nabla^2 h_i(\overline{\mathbf{x}}),
```

where ``\nabla^2 f(\overline{\mathbf{x}}), \nabla^2 g_i(\overline{\mathbf{x}})`` for ``i \in I``, and ``\nabla^2 h_i(\overline{\mathbf{x}})`` for ``i=1, \ldots, \ell``, are the Hessians of ``f, g_i`` for ``i \in I``, and ``h_i`` for ``i=1, \ldots, \ell``, respectively, all evaluated at ``\overline{\mathbf{x}}``. Define the cone
```math
\begin{aligned}
C=\left\{\mathbf{d} \neq \mathbf{0}: \nabla g_i(\overline{\mathbf{x}})^t \mathbf{d}=0\right. & \text { for } i \in I^{+}, \nabla g_i(\overline{\mathbf{x}})^t \mathbf{d} \leq 0 \text { for } i \in I^0, \\
\nabla h_i(\overline{\mathbf{x}})^t \mathbf{d}=0 & \text { for } i=1, \ldots, \ell\} .
\end{aligned}
```

Then if ``\mathbf{d}^t \nabla^2 L(\overline{\mathbf{x}}) \mathbf{d}>0`` for all ``\mathbf{d} \in C``, we have that ``\overline{\mathbf{x}}`` is a strict local minimum for P.
"""

# ╔═╡ fdcba702-afd5-4b24-860a-03f7b04d6f1e
cm"""
$(bbl("Corollary"))
Consider Problem P as defined in the theorem, and let ``\overline{\mathbf{x}}`` be a KKT point with associated Lagrangian multipliers ``\overline{\mathbf{u}}`` and ``\overline{\mathbf{v}}`` corresponding to the inequality and equality constraints, respectively. 

Furthermore, suppose that the collection of vectors ``\nabla g_i(\overline{\mathbf{x}})`` for ``i \in I^{+}=\left\{i \in I: \bar{u}_i>0\right\}`` and ``\nabla h_i(\overline{\mathbf{x}})`` for ``i=1, \ldots, \ell`` contains a set of ``n`` __linearly independent vectors__. 

Then ``\overline{\mathbf{x}}`` is a strict local minimum for P .
"""

# ╔═╡ 1511e16d-6573-46dd-8f71-930117000963
cm"""
$(bth("4.4.3 Theorem (KKT Second-Order Necessary Conditions)"))

Consider Problem P

$(eql_latex_gi())

where the objective and constraint defining functions are all __twice differentiable, and where ``X`` is a nonempty, open set in ``\mathbb{R}^n``__. 

Let ``\overline{\mathbf{x}}`` be a local minimum for Problem P , and denote 
```math 
I=\left\{i: g_i(\overline{\mathbf{x}})=\right. 0\}.
```

Define the restricted Lagrangian function ``L(\mathbf{x})`` as,
```math
L(\mathbf{x}) \equiv \phi(\mathbf{x}, \overline{\mathbf{u}}, \overline{\mathbf{v}})=f(\mathbf{x})+\sum_{i \in I} \bar{u}_i g_i(\mathbf{x})+\sum_{i=1}^{\ell} \bar{v}_i h_i(\mathbf{x}),
```

and denote its Hessian at ``\overline{\mathbf{x}}`` by
```math
\nabla^2 L(\overline{\mathbf{x}}) \equiv \nabla^2 f(\overline{\mathbf{x}})+\sum_{i \in I} \bar{u}_i \nabla^2 g_i(\overline{\mathbf{x}})+\sum_{i=1}^{\ell} \bar{v}_i \nabla^2 h_i(\overline{\mathbf{x}}),
```
where ``\nabla^2 f(\overline{\mathbf{x}}), \nabla^2 g_i(\overline{\mathbf{x}})`` for ``i \in I``, and ``\nabla^2 h_i(\overline{\mathbf{x}})`` for ``i=1, \ldots, \ell`` are the Hessians of ``f, g_i`` for ``i \in I``, and ``h_i`` for ``i=1, \ldots, \ell``, respectively, all evaluated at ``\overline{\mathbf{x}}``. 

Assume that ``\nabla g_i(\overline{\mathbf{x}})`` for ``i \in I``, and ``\nabla h_i(\overline{\mathbf{x}})`` for ``i=1 \ldots, \ell``, are __linearly independent__.

Then ``\overline{\mathbf{x}}`` is a __KKT point__ having Lagrange multipliers ``\overline{\mathbf{u}} \geq \mathbf{0}`` and ``\overline{\mathbf{v}}`` associated with the inequality and the equality constraints, respectively. 

Moreover, ``\mathbf{d}^t \nabla^2 L(\overline{\mathbf{x}}) \mathbf{d} \geq 0`` for all ``\mathbf{d} \in C`` where
```math 
\begin{array}{lcl}
C&=\left\{\right.&\mathbf{d} \neq \mathbf{0}: \nabla g_i(\overline{\mathbf{x}})^t \mathbf{d}=0\text{ for }i \in I^{+}, \\
&& \nabla g_i(\overline{\mathbf{x}})^t \mathbf{d} \leq 0\text{ for }i \in I^0,\\
&& \nabla h_i(\overline{\mathbf{x}})^t \mathbf{d}=0\text{ for all }i=1, \ldots, \ell\\
&&\},
\end{array}
```
where ``I^{+}=\left\{i \in I: \bar{u}_i>0\right\}`` and ``I^0=\left\{i \in I: \bar{u}_i=0\right\}``.
"""

# ╔═╡ 727d41be-4f75-4dca-b117-434ea4517aae
cm"""
$(ex("Example",""))
```math
\mathrm{P}: \min\left\{\left(x_1-1\right)^2+x_2^2: g_1(\mathbf{x})=2 k x_1-x_2^2 \leq 0\right\},\qquad k>0.
```

"""

# ╔═╡ 1214ea8b-3295-4ba6-b713-0f7db0ece34c
cm"""
$(define("Cone of tangent"))

Let ``S`` be a nonempty set in ``R^n``, and let ``\overline{\mathbf{x}} \in \mathrm{cl} S``. 

The __cone of tangents of ``S`` at ``\overline{\mathbf{x}}``__, denoted by ``T``, is 

```math
T =  \left\{\mathbf{d} \;|\; \mathbf{d}=\lim _{k \rightarrow \infty} \lambda_k\left(\mathbf{x}_k-\overline{\mathbf{x}}\right),
 \lambda_k > 0, \mathbf{x}_k \in S \text{ for each }k, \text{ and } \mathbf{x}_k \rightarrow \overrightarrow{\mathbf{x}}.
\right\}
```
"""

# ╔═╡ 2bb7e096-9311-4f92-8c25-74b0da8993d2
cm"""
$(bbl("Remarks",""))

From the above definition, it is clear that ``\mathbf{d}`` belongs to the cone of tangents if there is a feasible sequence ``\left\{\mathbf{x}_k\right\}`` converging to ``\overline{\mathbf{x}}`` such that the directions ``\mathbf{x}_k-\overline{\mathbf{x}}`` converge to ``\mathbf{d}``. (See Exercise 5.1),

$(post_img("https://www.dropbox.com/scl/fi/i7d73wt0riqwlrhbug05q/fig5.1.png?rlkey=zvgeopv5pcfon73h1774unw45&dl=1"))

"""

# ╔═╡ 64fe0bf7-8de4-4d43-ba35-c8d473970bff
cm"""
$(bth("5.1.2"))

Let ``S`` be a nonempty set in ``R^n``, and let ``\overline{\mathbf{x}} \in S``. Furthermore, suppose that ``f: R^n \rightarrow R`` is differentiable at ``\overline{\mathbf{x}}``. 

If ``\overline{\mathbf{x}}`` locally solves the problem to minimize ``f(\mathbf{x})``
subject to ``\mathbf{x} \in S``.

THEN

```math
F_0 \cap T=\varnothing,
```

where ``F_0=\left\{\mathbf{d}: \nabla f(\overline{\mathbf{x}})^t \mathbf{d}<0\right\}`` and ``T`` is the cone of tangents of ``S`` at ``\overline{\mathbf{x}}``.
"""

# ╔═╡ ea9e38bc-63f1-44e6-939d-7f9fc17e44c3
cm"""
$(bth("5.1.3 Theorem (Karush-Kuhn-Tucker Necessary Conditions)"))

Let ``X`` be a nonempty set in ``R^n``, and let ``f: R^n \rightarrow R`` and ``g_i: R^n \rightarrow R`` for ``i=1, \ldots``, ``m``. 

Consider the problem: 
$(min_latex_gi())

Let ``\overline{\mathbf{x}}`` be a feasible solution, and let ``I=\left\{i: g_i(\overline{\mathbf{x}})=0\right\}``. Suppose that ``f`` and ``g_i`` for ``i \in I`` are differentiable at ``\overline{\mathbf{x}}``. 

Furthermore, suppose that the constraint qualification ``T=G^{\prime}`` holds true, where ``T`` is the cone of tangents of the feasible region at ``\overline{\mathbf{x}}`` and ``G^{\prime}=\left\{\mathbf{d}: \nabla g_i(\overline{\mathbf{x}})^t \mathbf{d} \leq 0\right.`` for ``\left.i \in I\right\}``. If ``\overline{\mathbf{x}}`` is a local optimal solution, there exist nonnegative scalars ``u_i`` for ``i \in I`` such that
```math
\nabla f(\overline{\mathbf{x}})+\sum_{i \in I} u_i \nabla g_i(\overline{\mathbf{x}})=0
```
"""

# ╔═╡ 81eb374a-68ee-4bc7-ad02-b352536ad8a8
cm"""
$(bbl(" Lemma","5.1.4"))

Let ``\mathbf{A}`` be an ``m \times n`` matrix, let ``\mathbf{b}`` be an ``m``-vector, and let ``S=\{\mathbf{x}: \mathbf{A x} \leq \mathbf{b}\}``. 

Suppose that ``\overline{\mathbf{x}} \in S`` is such that ``\mathbf{A}_1 \overline{\mathbf{x}}=\mathbf{b}_1`` and ``\mathbf{A}_2 \overline{\mathbf{x}}<\mathbf{b}_2``, where ``\mathbf{A}^t=\left(\mathbf{A}_1^t, \mathbf{A}_2^t\right)`` and ``\mathbf{b}^t=\left(\mathbf{b}_1^t, \mathbf{b}_2^t\right)``. Then ``T=G^{\prime}``, where ``T`` is the cone of tangents of ``S`` at ``\overline{\mathbf{x}}`` and ``G^{\prime} =\left\{\mathbf{d}: \mathbf{A}_1 \mathbf{d} \leq \mathbf{0}\right\}``.
"""

# ╔═╡ 944649de-8db0-44a7-98b9-7de1855917e7
cm"""
$(post_img("https://www.dropbox.com/scl/fi/15x0oz0xjvpuyb0k0hacp/summary_51.png?rlkey=9fyktiu6f4cveams2kanbtzqc&dl=1"))
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlotThemes = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoExtras = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polyhedra = "67491407-f73d-577b-9b50-8179a7c68029"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
QRCoders = "f42e9828-16f3-11ed-2883-9126170b272d"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Colors = "~0.12.11"
CommonMark = "~0.9.1"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.4.0"
Latexify = "~0.16.8"
PlotThemes = "~3.3.0"
Plots = "~1.40.14"
PlutoExtras = "~0.7.15"
PlutoUI = "~0.7.65"
Polyhedra = "~0.8.1"
PrettyTables = "~2.4.0"
QRCoders = "~1.4.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "d9dac5ca2be7c023aa422d3b7b2ed3f4451f10dc"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BenchmarkTools]]
deps = ["Compat", "JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "e38fbc49a620f5d0b660d7f543db1009fe0f8336"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.6.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "TranscodingStreams"]
git-tree-sha1 = "84990fa864b7f2b4901901ca12736e45ee79068c"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.8.5"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonMark]]
deps = ["PrecompileTools"]
git-tree-sha1 = "351d6f4eaf273b753001b2de4dffb8279b100769"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.9.1"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "910febccb28d493032495b7009dce7d7f7aee554"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.0.1"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "b5c7fe9cea653443736d264b85466bad8c574f4a"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.9"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "4424dca1462cc3f19a0e6f07b809ad948ac1d62b"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.16"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d7ecfaca1ad1886de4f9053b5b8aef34f36ede7f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.16+0"

[[deps.GenericLinearAlgebra]]
deps = ["LinearAlgebra", "Printf", "Random", "libblastrampoline_jll"]
git-tree-sha1 = "54ee4866eb8c982ee23cf79230ca0aaf916c382b"
uuid = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
version = "0.3.15"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "8e233d5167e63d708d41f87597433f59a0f213fe"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.4"

[[deps.GeoInterface]]
deps = ["DataAPI", "Extents", "GeoFormatTypes"]
git-tree-sha1 = "294e99f19869d0b0cb71aef92f19d03649d028d5"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.4.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "b62f2b2d76cee0d61a2ef2b3118cd2a3215d3134"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.11"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43ba3d3c82c18d88471cfd2924931658838c9d8f"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+4"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "fee60557e4f19d0fe5cd169211fdda80e494f4e8"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.84.0+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed5e9c58612c4e081aecdb6e1a479e18462e041e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "437abb322a41d527c197fa800455f79d414f0a3c"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.8"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d65554bad8b16d9562050c67e7223abf91eaba2f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.13+0"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "196b41e5a854b387d99e5ede2de3fcb4d0422aae"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.2"

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

    [deps.JSON3.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "4f34eaabe49ecb3fb0d58d6015e32fd31a733199"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.8"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "fa7fd067dca76cadd880f1ca937b4f387975a9f5"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.16.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.MarchingCubes]]
deps = ["PrecompileTools", "StaticArrays"]
git-tree-sha1 = "0e893025924b6becbae4109f8020ac0e12674b01"
uuid = "299715c1-40a9-479a-aaf9-4a633d36f717"
version = "0.1.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON3", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays", "SpecialFunctions", "Test"]
git-tree-sha1 = "c2514ede436071529470010540bc3a1441e8140c"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.39.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "491bdcdc943fcbc4c005900d7463c9f216aabf4c"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

    [deps.OffsetArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "libpng_jll"]
git-tree-sha1 = "7dc7028a10d1408e9103c0a77da19fdedce4de6c"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.5.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9216a80ff3682833ac4b733caa8c00390620ba5d"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.0+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "275a9a6d85dc86c24d03d1837a0010226a96f540"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.3+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "28ea788b78009c695eb0d637587c81d26bdf0e36"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.14"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoExtras]]
deps = ["AbstractPlutoDingetjes", "DocStringExtensions", "HypertextLiteral", "InteractiveUtils", "Markdown", "PlutoUI", "REPL", "Random"]
git-tree-sha1 = "91d3820f5910572fd9c6077f177ba375e06f7a0e"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.7.15"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3151a0c8061cc3f887019beebf359e6c4b3daa08"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.65"

[[deps.Polyhedra]]
deps = ["GenericLinearAlgebra", "LinearAlgebra", "MathOptInterface", "MutableArithmetics", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "754bc39995daff07ed01d7ebdc8c9cf6681d241e"
uuid = "67491407-f73d-577b-9b50-8179a7c68029"
version = "0.8.1"

    [deps.Polyhedra.extensions]
    PolyhedraGeometryBasicsExt = "GeometryBasics"
    PolyhedraJuMPExt = "JuMP"
    PolyhedraRecipesBaseExt = "RecipesBase"

    [deps.Polyhedra.weakdeps]
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "13c5103482a8ed1536a54c08d0e742ae3dca2d42"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.4"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.QRCoders]]
deps = ["FileIO", "ImageCore", "ImageIO", "ImageMagick", "StatsBase", "UnicodePlots"]
git-tree-sha1 = "b3e5fcc7a7ade2d43f0ffd178c299b7a264c268a"
uuid = "f42e9828-16f3-11ed-2883-9126170b272d"
version = "1.4.5"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "eb38d376097f47316fe089fc62cb7c6d85383a52"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "2766344a35a1a5ec1147305c4b343055d7c22c90"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "9537ef82c42cdd8c5d443cbc359110cbb36bae10"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.21"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "38f139cc4abf345dd4f22286ec000728d5e8e097"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.2"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "24c1c558881564e2217dcf7840a8b2e10caeb0f9"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.UnicodePlots]]
deps = ["ColorTypes", "Contour", "Crayons", "Dates", "FileIO", "FreeTypeAbstraction", "LazyModules", "LinearAlgebra", "MarchingCubes", "NaNMath", "Printf", "SparseArrays", "StaticArrays", "StatsBase", "Unitful"]
git-tree-sha1 = "ae67ab0505b9453655f7d5ea65183a1cd1b3cfa0"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "2.12.4"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d2282232f8a4d71f79e85dc4dd45e5b12a6297fb"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.23.1"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "af305cc62419f9bd61b6644d19170a4d258c7967"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.7.0"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "XML2_jll"]
git-tree-sha1 = "49be0be57db8f863a902d59c0083d73281ecae8e"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.23.1+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "54b8a029ac145ebe8299463447fd1590b2b1d92f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.44.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "9caba99d38404b285db8801d5c45ef4f4f425a6d"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.1+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522c1df09d05a71785765d19c9524661234738e9"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.11.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "cd155272a3738da6db765745b89e466fa64d0830"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.49+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "c950ae0a3577aec97bfccf3381f66666bc416729"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.8.1+0"
"""

# ╔═╡ Cell order:
# ╟─a0a41512-efe8-4848-ae7d-16e1606d0978
# ╟─a6148830-0424-4d85-af87-9efa87eaa4aa
# ╟─b2805e6f-a669-433f-9352-0a1f97fc2a52
# ╟─45572d31-300f-4e40-a755-9c099c58551a
# ╟─ed6a2919-8df2-4907-b7cb-93e4cf4c0500
# ╟─7bc65c80-3346-452e-b4de-e45ce3a19461
# ╟─fcc354c2-c077-4a5a-9b84-db1d5f9f4ee7
# ╟─9d275485-9e6f-450e-8392-787ba3cda9a6
# ╟─a889687a-ec1d-400d-901d-894e81fb5549
# ╟─98d3cb65-7c5e-49d7-89df-5b32452a7067
# ╟─80e8de21-61e6-4e8f-869c-a364eb07f42d
# ╟─f26b7a51-3583-42da-8022-76a7aa6fec5a
# ╟─660817e9-5a67-4eca-9722-d18e57bc5868
# ╟─d46a3b04-cc28-4bda-90df-e95778f9bfa2
# ╟─4b3aef5d-7009-4db1-bd29-04acf3bbaea3
# ╟─ed268a43-ac57-49a0-992e-8e9c16cc1d28
# ╟─6250662f-9589-47d4-80f6-be17341180f4
# ╟─4c8505f9-7a53-4673-b1c2-26a82daf419d
# ╟─91fed789-5722-4fbe-bb8f-98f44cd86a47
# ╟─2010bad9-f0b6-41f9-9467-0b9e89daabaa
# ╟─c96c98c9-787d-4b06-a5e1-a3874b255938
# ╠═9e538d00-9783-436c-87f9-e0b9187cf5ba
# ╟─cbf0cb8d-0fc4-4886-87b4-a7da72210ed4
# ╟─46f49abf-eb7c-47cf-8e27-eecaad8a2cab
# ╟─fa1dec30-f740-46da-b2b3-628ea1527b5f
# ╟─74c85703-9286-4e8b-aea2-5265987f09dd
# ╟─ce80ae3e-1759-462d-bcc8-568c15bd166d
# ╟─90de87f6-3707-45c6-a6a9-3cc9d4debe64
# ╟─5aa374d0-8191-40aa-a06c-d1e115d07b1d
# ╟─4f44ea9c-f6de-42a2-88b2-a3c772f3e80d
# ╟─3cc11287-aba8-4e86-a7ef-2172d9c73fb4
# ╠═ca2c639a-a8c5-4c31-9732-fa3bdc284347
# ╟─6f306138-3d90-47ef-bd46-58ea33ecfcb8
# ╟─6bf7f29f-8c9b-46eb-ad4d-244336a02e73
# ╟─6a6d3f65-6cd4-4672-a007-98c92fbcb76d
# ╟─2346835d-5fc8-4716-833d-0b8d54dc6cd4
# ╟─3cefbf49-7946-4978-802c-9c46ff83835e
# ╟─8e2102e2-5dee-4198-8326-3c83b8d07cf2
# ╟─818c7aac-3310-468c-b0c0-388bca923064
# ╟─59d6f09b-2132-4437-a69f-4a51172ef6ff
# ╟─13369f35-e62d-4f84-830d-58999af481a8
# ╟─7958b6c5-f805-41ee-b290-c01abc51547f
# ╟─41a713fd-1e21-4034-9b77-34db9a67a8e6
# ╟─dc4a1115-b3e0-45d7-ad83-9debf9e55e60
# ╟─fb465e6a-2a6d-4ba3-9624-19acc0cb236c
# ╟─5b4a35bb-89a9-46bb-910a-790af826a6d3
# ╟─13796a55-9f61-4aea-9349-f4c3498ec209
# ╟─e721d5be-313a-4df9-b37e-e560c9bcaa97
# ╟─7b542f58-3b2b-4bf5-8183-2ac3588fe464
# ╟─43d02e39-037d-445f-9f8e-96866c34e858
# ╟─68aa1d60-0777-4066-9037-3513a705d3f3
# ╟─7559a2b5-9a9e-43db-9605-e2262f0a5c83
# ╟─cfc46ffc-5b02-4b61-bfec-9c313a8d60c4
# ╟─ea540681-5897-4745-bd43-f1671bbe2d94
# ╟─9034672d-b285-4e44-b472-63f5d1876477
# ╟─357b4c18-6b9d-407b-b46a-e470f61e1bf4
# ╟─a8cf2766-7784-421c-8a2c-d2fcb26a7344
# ╟─cc560f5c-ee2e-42a3-b0c5-b2f7e5460a74
# ╠═97af3754-7d2c-4db1-8124-c78cacaa4d3f
# ╟─7ba2f271-9064-4f72-b51c-3834a3deec56
# ╟─23d1eb38-8aaf-4c91-b6a6-f7eb508b88dd
# ╟─adb20e86-d93b-4b4f-89a9-1005af9df039
# ╟─0fffa7d0-74af-4d71-bc18-dd64baf1eb13
# ╟─887631f8-c063-44eb-a795-be849ed6b284
# ╟─a4dcf88f-3df2-4971-a130-c575623d7bc0
# ╟─702156d9-521c-47ee-9965-3059684a5a8d
# ╟─98c2e6e3-5669-49a8-ab90-ed787f738700
# ╟─17dfacb3-49e6-471b-843a-b3fad927cb26
# ╟─1a7107ab-143c-4dc5-bb1b-ec4493915682
# ╟─df3497ce-8187-4f19-8351-c3cab6fc3e04
# ╟─55ee2a12-4187-40c0-b832-89de362cca82
# ╟─be931167-e339-4392-9245-145a8aa6df53
# ╟─ff77a554-3b51-4592-98dd-a7a7b12efcaa
# ╟─b7858895-5fb6-4d58-9ff4-b51f55db32ca
# ╟─5bdaef6d-410f-4421-b2cf-e113cb40fe5c
# ╟─2dbd94c1-6303-4453-930e-4caee827ea02
# ╟─af6e8b0b-efe3-4168-b2ea-035c08969ba8
# ╟─5f1d27f3-3808-4c00-850f-7f2da641e03e
# ╟─9b459391-83c7-47a0-a414-92f934da112c
# ╟─0d19c555-5d0b-46c4-a765-748094298d75
# ╟─cc6206a9-e2b1-4f31-9d5c-368e5e20ba58
# ╟─8e354fee-c8c6-4fa4-a6a8-f13e5e1c1ab2
# ╟─6c72b325-3abc-4ede-aff5-2f5bf192649d
# ╟─e8cbb8c4-f137-453d-8e8e-87199c79a562
# ╟─6a62cb74-f9bc-465f-9cdc-9381232b5b3c
# ╟─c102db9c-b4f3-4398-a055-c0b024a46446
# ╟─95183221-8558-470b-aa94-1a40359c2562
# ╟─a9601a07-00a8-4a5c-a1eb-060b47a98912
# ╟─3225e74a-4f0e-4cdd-a687-b3cad3349823
# ╟─9bd6deb3-425b-451b-9ec1-43ca55b51445
# ╟─965ab318-1528-4d9e-a586-e7ea4382ac57
# ╟─64cedf28-216b-4b03-924b-c07c11b33a51
# ╟─c38ad4af-2245-4482-900e-7691577bcff2
# ╟─78422f45-8847-45a9-a269-3f0ee5918076
# ╟─3a82aafe-3d11-4743-bcfc-0eb2061b94c2
# ╟─a5eea233-4510-4b7b-b7dd-51650c4a9300
# ╟─69c4460b-93f7-440e-98c0-9dd9c66eb8fb
# ╟─9120ebc6-f004-40e5-a2bb-69f5c7b6f74b
# ╟─c6bd9fa2-e9e8-4bb3-a866-1b450d3bf8d8
# ╟─e949636a-b55d-4ae3-b605-26ad244a5be2
# ╟─683c53d5-7db0-49ee-b2c7-e0b87248fdd4
# ╟─248cef95-c012-446c-8e4d-258fb8f06410
# ╟─5617b3db-ed53-4c81-b7a5-073c8f34fd9f
# ╟─71bb4d0e-8000-44ee-96c9-a356ab2afe3d
# ╟─0d650b2c-6ec8-4822-b875-43a6a5b52879
# ╟─6e464ce7-f3d2-43dd-b529-a6171e9dd898
# ╟─cd0a647d-eeac-4798-99e0-a10549033f48
# ╟─2f340ca4-9fe2-4910-bec3-9fb86c7ed6a7
# ╟─ae88db5b-7923-4216-98af-570d1ecf39dc
# ╟─525e28ee-6af9-4e5b-9ce1-1f87881ff681
# ╟─69aecb16-725a-49dd-aba4-9775d797aaae
# ╟─017ac677-b8d3-40a9-90c3-ec7d4c463f0f
# ╟─44c72d9b-0adc-49e7-a878-184951cefe0d
# ╟─3c14bb03-4fff-47d0-8fc8-643561950b2a
# ╟─788874c6-a13f-4718-86d5-bdbba1b588c8
# ╟─657034bf-4034-4b37-995b-cc3e22a6ff19
# ╟─f22887ef-aabb-4a52-9f38-afb2f82ed16b
# ╟─4a15eb36-dd71-4ad4-9512-6f189a448bb8
# ╟─de381e37-2c8e-4e42-a4b8-292cb349ed3d
# ╟─af8c90dd-3f32-40ed-91d1-1287919109c4
# ╟─3e9400d7-48f7-4b1c-8bdf-441b874e99a7
# ╟─c78f468a-75c8-41dc-896c-1021765adf83
# ╟─57d1b602-1b3f-4de6-85d3-ac158a01bcc8
# ╟─20fdb0be-e603-4414-bb8e-d9df5b3d4666
# ╟─1e87dde2-48a2-4158-b20a-af94fc2f3308
# ╟─5aaa0e88-cdd3-4688-adfa-11487ab512ac
# ╟─53023197-6e74-4968-8fa4-ce0fee8bf6d9
# ╟─08a79da6-d7be-4749-8c7a-960dd85d9404
# ╟─9d66f924-30a8-428e-9f48-06b0be9b9687
# ╠═b2cbaa21-30b7-4eaa-bbba-2a100b4e3f7f
# ╟─e8637759-f3b0-4acd-8d37-0d61098d8b16
# ╟─6e85e225-d792-4ad7-aade-e048c78f62e2
# ╟─997295b0-bc83-4a2e-a81f-52212c041152
# ╟─9f633dba-e6b5-4d9f-b7b0-2f505c4642ab
# ╟─da04715d-9713-4230-9d25-df069c19c9d4
# ╟─1c9a2d2f-39d5-453e-8a82-bda18570e762
# ╟─61e3e889-4444-4ecb-91a1-3d8f91d0054a
# ╟─44ecea1f-be97-4fb5-a62e-0177b98d5404
# ╟─c3cfdf7e-0621-4add-9d85-945508a44eec
# ╟─8c011096-3474-4d80-b40d-d72fd50e621b
# ╟─b283d6eb-0d60-4eb8-af1b-2a30f4f5e596
# ╟─f6c2cbee-b331-4132-946e-bbe9fd3b5881
# ╟─8d3beee4-64ef-4942-adc2-739122182dc9
# ╟─f6180dbc-4c0d-4b06-9670-9a4122e7681c
# ╟─2adf10df-7961-4329-ad96-f031d405586f
# ╟─9f88e2bc-4137-42a0-bcc2-e4d377c27f00
# ╟─d3bcf477-8bf6-4193-a6b1-d0f924c32bf5
# ╟─8b0eaffa-18c3-4174-a4db-8cb514c728bf
# ╟─a5ccf185-1a61-47ea-9fe2-b33b7bba8e6c
# ╟─15af64e4-366d-464e-962a-42bf619f2e4e
# ╟─a0e01d7c-9ea7-46b1-a980-9d4e43315c13
# ╟─46187f5f-d96b-45bd-94a3-f9f1d7446960
# ╟─ff143a15-6e2e-4eb3-8faa-d683c417bf64
# ╟─a874e129-7812-4515-8ffe-d875c73813e5
# ╟─f26fed75-941b-4ab4-9720-405281d00170
# ╟─3efb4f55-9325-4e02-b53f-d275dbb405b3
# ╟─2268c221-9e2e-403e-9729-53611b5bc63c
# ╟─6629edc2-72c3-4aa9-9e19-56b9a30adeb4
# ╟─375e3f9e-45fb-42da-b177-54c77c1081db
# ╟─c3f319d2-23fc-4c41-a7f2-23f831073ae0
# ╟─7a0aa831-6dd4-4b77-891d-aa56892f7759
# ╟─7deefd68-4955-4fdb-bc76-1de2ccd841c8
# ╟─819fbc46-b696-4522-8fa4-ee68c89058b4
# ╟─7a73b652-eae0-4ebe-b472-2ac1984607cf
# ╟─12e8d772-374b-438f-a5b9-df7fdab33d4a
# ╟─434b0a5a-6937-4ed4-ac9a-8961d9578145
# ╟─f902fff5-6175-48d2-bebc-27c7e0f72d10
# ╟─960eff24-3070-4051-a4e0-f25a00b935b7
# ╟─e87f8c80-c756-4ff8-9e7d-ef35d8941afd
# ╟─22b3e9fa-cf5d-43c3-b964-758b26e33468
# ╟─9cc4f6fc-b656-48d7-8c8f-da23d1e5419c
# ╟─4cb4d21c-464e-4675-b0b8-3872195ecc76
# ╟─a03f6f56-54a8-423a-a405-5f9a83e5cdb2
# ╟─8f0be524-cabb-4b78-99a6-9ecbfc6e57a3
# ╟─fc4c567e-3cc2-4b85-9294-eae1e2da69fd
# ╟─ab11e8fc-98e8-4372-a12f-c6133dcc65e3
# ╟─17154407-9b5d-46be-92b9-2d02005e5c9c
# ╠═7dfbcce2-1f0f-4f43-9e01-05bd0447ed32
# ╟─8eba2357-e825-4232-875a-f18be55bd38a
# ╟─4eebca6a-7e54-4d1c-84b2-cd893d4d2a3f
# ╟─3d3223ea-d105-45a5-9389-81de85a271ba
# ╟─b6809792-4885-4bdd-9dea-d4e93ebe68c6
# ╟─563d0e84-e2f3-49fe-b92e-6c0716ed8523
# ╟─010e228c-7590-4178-99b6-bb50819dfe1c
# ╟─807a80ef-2c20-400d-85b6-9fc71d21b5b8
# ╟─d478679e-e6a0-48cc-8119-7e822f26dd02
# ╟─72ccae91-7342-4780-92c6-b226ce13507c
# ╟─ff54a546-4bd6-4d6e-9ea7-f4053635dc42
# ╟─ee88be92-5954-4f41-9b6d-c1db938b7368
# ╟─6a7ff171-88af-431f-9ee8-4724f8eff61d
# ╟─b14b538e-e980-4e49-b9e0-74b43fe4620b
# ╟─2eea535f-8ab5-434b-b685-71539f73f062
# ╟─8e1b7954-4159-4fd4-8c5a-2f1ac529e463
# ╟─1391def6-0a0d-4e16-abbe-3644e56cab9b
# ╟─2f98c41a-3854-4a4a-b496-596d6f953d33
# ╟─e4c10925-07a1-4aa8-a921-bba1f49976ba
# ╟─524c222e-f9e1-4421-88b0-072cb6b74d95
# ╟─f68e53bc-5c85-4b4b-95b7-d94af1a7f124
# ╟─4719b6e7-4eba-43a4-8b8d-3bdf3a46712f
# ╟─10754f5c-e457-47a5-91d6-5f6ed203cda8
# ╟─8f24c28e-57bd-42ae-9883-1cebc0715a0c
# ╟─93a35971-c2eb-4dd3-bd66-b3335400861d
# ╟─dc71a49a-1e1f-474c-bfcd-cde361115ee1
# ╟─ee994ad9-5fc9-4722-9201-536996e19b4c
# ╟─5196f462-29fe-4f16-8e8d-c9b233a7dca8
# ╟─6d6428eb-9442-43a3-9cad-815b7ccfc203
# ╟─e64e8a77-6b32-46ab-b126-f77c917935ad
# ╟─b7f2b3d3-8aa5-4a1e-97e5-a10f17aebbaf
# ╟─a60d5096-9152-44e0-b5d2-3bb789dcff5d
# ╟─b2cdd374-e06f-4d8f-9db1-6dfb7d661a40
# ╟─10f03ebd-0da6-4ecf-afdd-e524a9ddf530
# ╟─75d8b9f6-8d89-49b9-b8a4-e7cbbd066a31
# ╟─9197b26f-59a7-4508-82f7-8441bea7a2e1
# ╟─f1b16a2c-816a-4a31-aff2-f8f10ca45cc9
# ╟─fd688b17-447f-4aa6-8584-20b04ef1f822
# ╟─4699162f-ef6a-4279-9aa7-56170fa9a6ff
# ╟─7ed6becd-1e79-4103-950a-017d187c585c
# ╟─b06c5711-3f32-4514-a44d-e3aef2e69125
# ╟─2d7d967e-bd11-4959-978e-2325f1b78f95
# ╟─1791db8c-eb23-450c-b79d-82124a027cee
# ╟─2b18aae5-24c1-4f3b-8926-9993c3f47db2
# ╟─55c127e9-9647-4f9f-b9fd-4d9752071a90
# ╟─cefc755f-2302-4213-8836-04ff8e33305d
# ╟─38143cd0-3552-4083-954a-bb10efb8c090
# ╟─1552ddcb-1e1d-43de-a9e3-f4a7c186970b
# ╟─bccbbf3a-1cf3-41e2-bc0d-b32e122421ec
# ╟─f4969ec2-4f5b-4893-9385-d74b038f4738
# ╟─6e386883-6678-4ed6-9c09-bd9a41201783
# ╟─995242c4-a629-4491-9db5-645730d6b9bb
# ╟─2e8796d7-06b4-467e-b295-1019d34859bf
# ╟─9389bae8-5492-4402-ac8a-fccadfd8351e
# ╟─c6554735-c50a-4c61-a9ef-95dd3956b99e
# ╟─09a738b0-d156-4123-a5c1-a3fc349d680a
# ╟─402f4f37-fc08-4cb7-9b64-ffd028711de4
# ╟─e319bc66-7810-4673-bc77-c3e0d8f1eeac
# ╟─a8133296-0ec7-4504-96bb-39e443193c7f
# ╟─3758e93e-ed91-478e-b431-caccf55bdda0
# ╟─ca7bd385-5269-494e-a2b4-7a0b779d4389
# ╟─0f4d581b-9e6b-47d7-b014-d0d161cf84a9
# ╟─fa6c7bd2-e304-4555-833d-989d6414964a
# ╟─07cb2ece-065f-4790-b8e2-ef164907d69d
# ╟─ecf678cd-8702-420c-9c6c-d7010ffe42f1
# ╟─237b661f-d295-4c50-ae22-dd3441881cc1
# ╠═20b7851e-1cd3-4a43-9e41-80ab4ad38ccc
# ╟─54c9ffd1-97be-44cd-8b78-07ec3f6f6299
# ╟─e15eb170-02e9-42e9-b8d5-35d9284fd30e
# ╟─4c0fffcb-327a-4a45-a8b6-cea81005ab2b
# ╟─d19630de-5c1d-47b5-b91e-359ae2898453
# ╟─fdcba702-afd5-4b24-860a-03f7b04d6f1e
# ╟─1511e16d-6573-46dd-8f71-930117000963
# ╟─727d41be-4f75-4dca-b117-434ea4517aae
# ╟─ebc2e75c-21e1-452c-9e09-da88cdc488a6
# ╟─1c78770b-4814-48ee-837a-1c0f2c99a7c8
# ╟─1214ea8b-3295-4ba6-b713-0f7db0ece34c
# ╟─2bb7e096-9311-4f92-8c25-74b0da8993d2
# ╟─64fe0bf7-8de4-4d43-ba35-c8d473970bff
# ╟─bdf534ab-e2a7-456c-83dd-23fe5cb48028
# ╟─6db4f42f-6c20-44db-a0e6-c0b061781538
# ╟─ea9e38bc-63f1-44e6-939d-7f9fc17e44c3
# ╟─af44bbd1-b9a6-4668-99a2-794fad3f6a42
# ╟─81eb374a-68ee-4bc7-ad02-b352536ad8a8
# ╟─944649de-8db0-44a7-98b9-7de1855917e7
# ╠═41c749c0-500a-11f0-0eb8-49496afa257e
# ╟─42f6c9db-97d9-4852-a4c3-f7bbcb055a0f
# ╟─fc877247-39bc-4bb0-8bda-1466fcb00798
# ╟─fdd3c3e3-5089-456f-adef-7ab2e311331f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
