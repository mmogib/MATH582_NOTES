# Interactive Optimization Visualization with Pluto.jl - Complete Guide

## Quick Start Prompt

Copy and paste this to a new chat to recreate interactive optimization visualizations:

---

I need help creating interactive Pluto.jl notebooks for visualizing constrained optimization problems. Here's what I need:

### Requirements

**Structure:** Create a Pluto notebook with:
1. Interactive input using `NumberField` from PlutoUI (not sliders)
2. Plot showing:
   - Contour lines of the objective function
   - Constraint boundaries
   - Shaded feasible region
   - Selected point (green if feasible, red if not)
   - Gradient arrows of objective function (∇f) in black
   - Gradient arrows of active constraints (∇g₁, ∇g₂, etc.) in different colors
   - LaTeX annotations showing gradient values next to arrows
3. Information section displaying point details, feasibility, constraint values
4. KKT conditions explanation

**Packages needed:**
```julia
using Plots
using PlutoUI
using LinearAlgebra
using LaTeXStrings
using CommonMark
```

**Key features:**
- Use `framestyle=:origin` to show axes through origin
- Constraints written as g(x) ≤ 0
- Tolerance of 0.1 for detecting active constraints
- Use `cm"""` with CommonMark for LaTeX in input labels
- Scale factor around 0.3-0.4 for gradient arrows
- Format: `@bind x1_val NumberField(0.0:0.01:3.0, default=1.5)`

**Input format:**
```julia
begin
	x1_val_html = @bind x1_val NumberField(min:step:max, default=val)
	x2_val_html = @bind x2_val NumberField(min:step:max, default=val)
	cm"""
	``x_1=`` $(x1_val_html)
	``x_2=`` $(x2_val_html)
	"""
end
```

**Gradient annotations format:**
```julia
annotate!(p, x1_val + scale*grad[1] + 0.1, 
          x2_val + scale*grad[2] + 0.1, 
          text(L"\nabla f = [%$(round(grad[1], digits=2)), %$(round(grad[2], digits=2))]", 
               :black, 8, :left))
```

**Active constraint checking:**
```julia
function active_constraints(x1, x2)
	active = []
	if abs(g1(x1, x2)) < tol
		push!(active, (1, "constraint name", ∇g1(x1, x2)))
	end
	return active
end
```

### Problem to Solve

[Insert your specific optimization problem here]

Example format:
```
Minimize (x₁ - 3)² + (x₂ - 2)²
subject to x₁² + x₂² ≤ 5
           x₁ + 2x₂ ≤ 4
           x₁ ≥ 0
           x₂ ≥ 0
```

Please create the complete Pluto notebook code with all sections.

---

## Detailed Technical Reference

### 1. Package Setup
```julia
using Plots
using PlutoUI
using LinearAlgebra
using LaTeXStrings
using CommonMark
```

### 2. Interactive Input with NumberField
```julia
begin
	x1_val_html = @bind x1_val NumberField(0.0:0.01:3.0, default=1.5)
	x2_val_html = @bind x2_val NumberField(0.0:0.01:3.0, default=1.0)
	cm"""
	``x_1=`` $(x1_val_html)
	``x_2=`` $(x2_val_html)
	"""
end
```

### 3. Function Definitions Structure
```julia
begin
	# Objective function
	f(x1, x2) = (x1 - 3)^2 + (x2 - 2)^2
	
	# Gradient of objective
	∇f(x1, x2) = [2*(x1 - 3), 2*(x2 - 2)]
	
	# Constraints as g(x) ≤ 0
	g1(x1, x2) = x1^2 + x2^2 - 5
	g2(x1, x2) = x1 + x2 - 3
	g3(x1, x2) = -x1  # x1 ≥ 0
	g4(x1, x2) = -x2  # x2 ≥ 0
	
	# Gradients of constraints
	∇g1(x1, x2) = [2*x1, 2*x2]
	∇g2(x1, x2) = [1, 1]
	∇g3(x1, x2) = [-1, 0]
	∇g4(x1, x2) = [0, -1]
	
	# Tolerance for active constraints
	tol = 0.1
	
	# Check active constraints
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
	
	# Check feasibility
	function is_feasible(x1, x2)
		return g1(x1, x2) <= tol && g2(x1, x2) <= tol && 
		       g3(x1, x2) <= tol && g4(x1, x2) <= tol
	end
end
```

### 4. Plotting Setup
```julia
begin
	x1_range = range(-0.5, 3.5, length=400)
	x2_range = range(-0.5, 3.0, length=400)
	
	p = plot(size=(800, 800), aspect_ratio=:equal, 
	         xlabel=L"x_1", ylabel=L"x_2", 
	         title="Gradients at ($(round(x1_val, digits=2)), $(round(x2_val, digits=2)))",
	         legend=:topright, legendfontsize=8,
	         framestyle=:origin)  # Shows axes through origin
```

### 5. Plotting Contours

**For nonlinear objectives:**
```julia
contour!(p, x1_range, x2_range, 
         (x1, x2) -> f(x1, x2), 
         levels=15, 
         color=:viridis, 
         linewidth=1.5,
         alpha=0.5,
         colorbar=false,
         label="")
```

**For linear objectives (e.g., f = -x₁):**
```julia
for x1_contour in -0.5:0.2:2.0
	plot!(p, [x1_contour, x1_contour], [ymin, ymax],
	      color=:gray, alpha=0.3, linewidth=1, label="")
end
```

### 6. Plotting Constraint Boundaries

**Circle constraint (x₁² + x₂² = r²):**
```julia
θ = range(0, 2π, length=200)
circle_x1 = sqrt(5) .* cos.(θ)
circle_x2 = sqrt(5) .* sin.(θ)
plot!(p, circle_x1, circle_x2, 
      linewidth=2.5, color=:red, 
      label=L"x_1^2 + x_2^2 = 5", linestyle=:dash)
```

**Line constraint (ax₁ + bx₂ = c):**
```julia
x1_line = range(xmin, xmax, length=100)
x2_line = (c .- a .* x1_line) ./ b  # Solve for x2
plot!(p, x1_line, x2_line, 
      linewidth=2.5, color=:blue, 
      label=L"x_1 + x_2 = 3", linestyle=:dash)
```

**Cubic/polynomial constraint:**
```julia
x1_curve = range(xmin, xmax, length=200)
x2_curve = (1 .- x1_curve).^3
plot!(p, x1_curve, x2_curve, 
      linewidth=2.5, color=:red, 
      label=L"x_2 = (1-x_1)^3", linestyle=:dash)
```

**Axis constraints:**
```julia
plot!(p, [0, 0], [ymin, ymax], 
      linewidth=2.5, color=:green, 
      label=L"x_1 = 0", linestyle=:dash, alpha=0.5)
plot!(p, [xmin, xmax], [0, 0], 
      linewidth=2.5, color=:orange, 
      label=L"x_2 = 0", linestyle=:dash, alpha=0.5)
```

### 7. Shading Feasible Region
```julia
x1_grid = range(xmin, xmax, length=150)
x2_grid = range(ymin, ymax, length=150)
feasible_x1 = Float64[]
feasible_x2 = Float64[]

for x1 in x1_grid
    for x2 in x2_grid
        if x1 >= 0 && x2 >= 0 && g1(x1, x2) <= 0 && g2(x1, x2) <= 0
            push!(feasible_x1, x1)
            push!(feasible_x2, x2)
        end
    end
end

scatter!(p, feasible_x1, feasible_x2, 
         markersize=1, 
         markerstrokewidth=0,
         color=:lightblue, 
         alpha=0.2,
         label="")
```

### 8. Plotting Selected Point
```julia
# Mark center of objective (if applicable)
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
```

### 9. Plotting Gradients with Annotations

**Objective function gradient:**
```julia
grad_f = ∇f(x1_val, x2_val)
scale = 0.3  # Scale factor for arrow length

quiver!(p, [x1_val], [x2_val], 
        quiver=([scale*grad_f[1]], [scale*grad_f[2]]), 
        color=:black, 
        linewidth=3,
        arrow=:closed,
        label=L"\nabla f")

# Annotation with LaTeX
annotate!(p, x1_val + scale*grad_f[1] + 0.1, 
          x2_val + scale*grad_f[2] + 0.1, 
          text(L"\nabla f = [%$(round(grad_f[1], digits=2)), %$(round(grad_f[2], digits=2))]", 
               :black, 8, :left))
```

**Active constraint gradients:**
```julia
active = active_constraints(x1_val, x2_val)
colors = [:red, :blue, :green, :orange]

for (i, name, grad) in active
	quiver!(p, [x1_val], [x2_val], 
	        quiver=([scale*grad[1]], [scale*grad[2]]), 
	        color=colors[i], 
	        linewidth=3,
	        arrow=:closed,
	        label=L"\nabla g_%$i")
	
	# Annotation
	offset_x = 0.1
	offset_y = 0.1
	annotate!(p, x1_val + scale*grad[1] + offset_x, 
	          x2_val + scale*grad[2] + offset_y, 
	          text(L"\nabla g_%$i = [%$(round(grad[1], digits=2)), %$(round(grad[2], digits=2))]", 
	               colors[i], 8, :left))
end

xlims!(p, xmin, xmax)
ylims!(p, ymin, ymax)

p  # Display plot
end
```

### 10. Information Display Section
```julia
md"""
## Information about the selected point

**Point:** ($(round(x1_val, digits=3)), $(round(x2_val, digits=3)))

**Feasible:** $(is_feasible(x1_val, x2_val) ? "✓ Yes" : "✗ No")

**Objective value:** f = $(round(f(x1_val, x2_val), digits=3))

**Gradient of f:** ∇f = [$(round(∇f(x1_val, x2_val)[1], digits=3)), $(round(∇f(x1_val, x2_val)[2], digits=3))]

**Constraint values:**
- g₁ (circle): $(round(g1(x1_val, x2_val), digits=3)) $(g1(x1_val, x2_val) <= 0 ? "✓" : "✗")
- g₂ (line): $(round(g2(x1_val, x2_val), digits=3)) $(g2(x1_val, x2_val) <= 0 ? "✓" : "✗")
- g₃ (x₁≥0): $(round(g3(x1_val, x2_val), digits=3)) $(g3(x1_val, x2_val) <= 0 ? "✓" : "✗")
- g₄ (x₂≥0): $(round(g4(x1_val, x2_val), digits=3)) $(g4(x1_val, x2_val) <= 0 ? "✓" : "✗")

**Active constraints:** $(length(active_constraints(x1_val, x2_val)) > 0 ? join([name for (i, name, grad) in active_constraints(x1_val, x2_val)], ", ") : "None")
"""
```

### 11. KKT Conditions Section
```julia
md"""
## KKT Conditions Check

For a point to be optimal, it must satisfy the KKT conditions:
1. **Stationarity:** ∇f + Σλᵢ∇gᵢ = 0 (for active constraints)
2. **Primal feasibility:** All constraints satisfied
3. **Dual feasibility:** λᵢ ≥ 0
4. **Complementary slackness:** λᵢ·gᵢ = 0

**Try these interesting points:**
- **(2.0, 1.0)** - Description
- **(1.5, 1.0)** - Description

[Add problem-specific insights here]
"""
```

## Common Constraint Types and Gradients

### Constraint Types

| Constraint | Form g(x) ≤ 0 | Gradient ∇g |
|------------|---------------|-------------|
| Circle | x₁² + x₂² - r² | [2x₁, 2x₂] |
| Line | ax₁ + bx₂ - c | [a, b] |
| Cubic | (x₁ + x₂ - c)³ | [3(x₁+x₂-c)², 3(x₁+x₂-c)²] |
| Cubic (alt) | x₂ - (1-x₁)³ | [3(1-x₁)², 1] |
| Non-negative x₁ | -x₁ | [-1, 0] |
| Non-negative x₂ | -x₂ | [0, -1] |

### Objective Function Gradients

| Objective | Gradient ∇f |
|-----------|-------------|
| (x₁ - a)² + (x₂ - b)² | [2(x₁-a), 2(x₂-b)] |
| -x₁ | [-1, 0] |
| x₁² + x₂² | [2x₁, 2x₂] |
| ax₁ + bx₂ | [a, b] |

## Design Guidelines

### Visual Design
1. **Color scheme:** 
   - Black for ∇f
   - Red, blue, green, orange for ∇g₁, ∇g₂, ∇g₃, ∇g₄
   - Green for feasible points, red for infeasible
   - Light blue for feasible region shading

2. **Arrow scaling:** 
   - Use scale = 0.3-0.4 for gradient arrows
   - Adjust based on gradient magnitude if needed

3. **Annotation placement:**
   - Add small offsets (0.1) to avoid overlap with arrows
   - Use `:left` alignment for most annotations

### Numerical Settings
- **Tolerance:** 0.1 for detecting active constraints
- **Grid resolution:** 150x150 for feasible region
- **Curve resolution:** 100-200 points for smooth boundaries
- **Number field step:** 0.01 for precise input

### Plot Settings
- **Size:** (800, 800) for square plots
- **Aspect ratio:** `:equal` for undistorted visualization
- **Frame style:** `:origin` to show coordinate axes
- **Legend:** `:topright` with fontsize 8

## Example Problems

### Problem 1: Quadratic with Circle and Line
```
Minimize (x₁ - 3)² + (x₂ - 2)²
subject to x₁² + x₂² ≤ 5
           x₁ + 2x₂ ≤ 4
           x₁ ≥ 0
           x₂ ≥ 0
```

### Problem 2: Linear with Cubic Constraint
```
Minimize -x₁
subject to x₂ - (1 - x₁)³ ≤ 0
           x₂ ≥ 0
```

### Problem 3: Linear with Linear Constraints
```
Minimize -x₁
subject to x₁ + x₂ - 1 ≤ 0
           x₂ ≥ 0
```

### Problem 4: Quadratic with Cubic (Degenerate)
```
Minimize (x₁ - 1)² + (x₂ - 1)²
subject to (x₁ + x₂ - 1)³ ≤ 0
           x₁ ≥ 0
           x₂ ≥ 0
```

## Special Cases and Notes

### Constraint Qualifications
- **LICQ (Linear Independence):** When gradients of active constraints are linearly independent
- **Degenerate cases:** When constraint gradient becomes zero (e.g., cubic at boundary)
- Always mention when constraint qualifications fail

### Educational Notes to Include
1. Explain KKT conditions clearly
2. Suggest specific test points to try
3. Note special behaviors (degenerate cases, multiple active constraints)
4. Explain why certain points are or aren't optimal
5. Discuss the geometric interpretation of gradients

### Common Issues
- **Zero gradients:** Handle cubic constraints carefully
- **Multiple active constraints:** Show how gradients combine
- **Constraint qualification failures:** Explain when/why they occur
- **Feasibility:** Color-code points appropriately

## Template Structure

Every notebook should have:
1. **Title and problem description**
2. **Interactive input section**
3. **Function definitions block**
4. **Main plotting block**
5. **Information display**
6. **KKT conditions explanation**
7. **Suggested test points**

## Additional Tips

- Use unique variable names (e.g., `x1_val_e3`) if creating multiple examples in one notebook
- Test with corner cases (origin, boundary intersections, etc.)
- Verify constraint gradients are correct before implementing
- Check that feasible region appears correctly
- Ensure annotations don't overlap with plot elements

## Resources

- [Pluto.jl Documentation](https://plutojl.org/)
- [PlutoUI.jl](https://github.com/JuliaPluto/PlutoUI.jl)
- [Plots.jl Documentation](https://docs.juliaplots.org/)
- [LaTeXStrings.jl](https://github.com/JuliaStrings/LaTeXStrings.jl)

---

## Quick Checklist

Before finalizing a notebook, verify:
- [ ] All packages imported
- [ ] NumberField inputs with appropriate ranges
- [ ] Constraints written as g(x) ≤ 0
- [ ] All gradients calculated correctly
- [ ] Active constraint detection works
- [ ] Feasible region shaded correctly
- [ ] All constraint boundaries plotted
- [ ] Gradient arrows display with annotations
- [ ] Information section shows all details
- [ ] KKT conditions explained
- [ ] Test points suggested
- [ ] Plot limits set appropriately
- [ ] Legend not cluttered
- [ ] Colors consistent and meaningful

---

**Version:** 1.0  
**Date:** 2025  
**Compatible with:** Pluto.jl, Julia 1.6+