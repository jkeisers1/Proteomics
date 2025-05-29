# test_sector_by_length.jl

# Ensure project directory is correct
cd(@__DIR__)

# Load modules
include("src/preprocessing.jl")
include("src/sector_analysis.jl")
include("src/sector_by_length.jl")

println("📦 Loading processed data...")
data = load_processed_data()

gene_ranges = [
    ("short", 0:1000),
    ("medium", 1001:1500),
    ("long", 1501:3000)
]

limitations = [:C_lim , :A_lim]
df_g_length = df_by_gene_length(gene_ranges, limitations, :mass, data)

save_gene_length_results(df_g_length, gene_ranges)
