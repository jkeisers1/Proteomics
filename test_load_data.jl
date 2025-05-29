# test_load_data.jl

# Load utility and data functions
include("src/utils.jl")
include("src/load_data.jl")

println("🔍 Loading all data from /data folder...")
data = load_all_data()

println("✅ Loaded datasets:")
for key in keys(data)
    println("  • ", key)
end

println("\n🔢 Preview of sheet 2 in 'samples2':")
if haskey(data["mass_fractions2"], 2)
    display(first(data["mass_fractions2"][2], 5))
else
    println("Sheet 2 not found.")
end
