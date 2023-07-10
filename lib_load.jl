include("constants.jl")
include("Harmonic.jl")

# How to find the path of a package in Julia?
# https://stackoverflow.com/questions/61440940/how-to-find-the-path-of-a-package-in-julia

# Function module file looks up the location of a module from the method table
module_file(modu) = String(first(methods(getfield(modu, :eval))).file)

println("\nFile lib_load.jl")
println("Include file: ", module_file(constants))
println("Include file: ", module_file(Harmonic))
