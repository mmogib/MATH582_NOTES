param(
    [string]$Message
)

# Run Julia
julia --project=. src/export.jl

# Git commands
git add .
git commit -m "$Message"
git push
