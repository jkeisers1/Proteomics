
using JLD2, DataFrames, Plots
using GLM  
using Pkg
Pkg.add("GLM")

# -----------------------------------------------------------------------------
# Local dependencies (assumes the existing project structure)
# -----------------------------------------------------------------------------

include("sector_analysis.jl")      # provides `load_processed_data`
include("sector_by_length.jl")     # provides `df_by_gene_length`

println("🔬 Analyzing gene list under carbon limitation...")
result_df = analyze_gene_list(["rmf"], data, :C_lim; mode = :mass)

g_len_result = df_by_gene_length(gene_ranges, [:C_lim, :A_lim], data; mode = :concentration)


plot_growth_vs_abundance(g_len_result["short"] ; mode = :concentration)

"""
    plot_growth_vs_abundance(path_to_jld2::AbstractString;
                             abundance_col::Symbol = :concentration,
                             savepath::Union{Nothing,String} = nothing) -> Plots.Plot

Load a 3-column DataFrame from `path_to_jld2` and create a scatter plot of
`abundance_col` (y-axis) against `:growth_rate` (x-axis).  
A simple linear regression line is overlaid.  
Pass `savepath = \"figure.pdf\"` (or .png, .svg) to write the figure to disk.
"""
function plot_growth_vs_abundance(df; mode::Symbol = :concentration, savepath::Union{Nothing,String}=nothing)
      # the variable name saved in the file, e.g. `df`

    # ------------------------------------------------------------------
    # Ensure the growth-rate column is called :growth_rate
    # ------------------------------------------------------------------
    long_name = Symbol("Growthrate(1/h)")
    if hasproperty(df, long_name)
        rename!(df, long_name => :growth_rate)
    end

    # ------------------------------------------------------------------
    # Basic scatter
    # ------------------------------------------------------------------
    p = scatter(df.growth_rate, df[!, 2];
                xlabel = "Growth rate (1/h)",
                ylabel = string(mode),
                markersize = 6,
                legend = false)

    # ------------------------------------------------------------------
    # Linear fit overlay
    # ------------------------------------------------------------------
    if savepath !== nothing
        savefig(p, savepath)
        @info "📊 figure written to $(abspath(savepath))"
    end
    return p
end