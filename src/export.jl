using PlutoSliderServer

notebook_name = "MATH582_NOTES"
notebook_dir = "src/$(notebook_name).jl"
docs = "docs"
source = "$docs/$(notebook_name).html"
dest = "$docs/index.html"

PlutoSliderServer.export_notebook(notebook_dir; Export_output_dir = docs)

mv(source, dest, force = true)
