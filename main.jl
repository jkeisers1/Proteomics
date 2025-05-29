cd(@__DIR__)
include("src/preprocessing.jl")
include("src/sector_analysis.jl")      # provides `load_processed_data`
include("src/sector_by_length.jl")     # provides `df_by_gene_length`
include("src/plotting.jl")


using .Plotting

# ── load (or regenerate) the processed data bundle ───────────────────────────
data = load_data()   # wrapper around load_processed_data()

# ── sector analysis: total *mass* fraction of rmf under C limitation ─────────
Plotting.plot_sector_analysis(
    ["hpf"],
    limitation = :C_lim,
    mode       = :concentration,
    data       = data)

# ── gene-length bins example – summed *concentration* under C & A limitation ─
gene_length = [
    ("0:500",       1:499),
    ("500:999",   500:999),
    ("1000:2000",   1000:2000),
    ("2000:5000", 2000:5000)
    ]

Plotting.plot_gene_length_analysis(
    gene_length,
    limitations = [:R_lim],
    mode        = :mass,
    ylims = nothing,
    xlims = nothing,
    data        = data)




