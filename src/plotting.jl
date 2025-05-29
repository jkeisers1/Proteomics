module Plotting

using DataFrames, JLD2, Plots

# ── Local project logic ------------------------------------------------------
include("sector_analysis.jl")      # load_processed_data, analyze_gene_list
include("sector_by_length.jl")     # df_by_gene_length

# ── Public entry: load bundled data ----------------------------------------
load_data() = load_processed_data()

# ── Utility: canonical abundance column name --------------------------------
_default_abundance_sym(mode::Symbol) = mode == :concentration ? :concentration : Symbol(string(mode, "_fraction"))

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Core helper: scatter of abundance vs. growth rate                       ║
# ╚══════════════════════════════════════════════════════════════════════════╝

"""
    plot_growth_vs_abundance(df::DataFrame;
                              mode::Symbol = :mass,
                              abundance_col::Union{Symbol,Nothing} = nothing,
                              xlims = nothing,
                              ylims = nothing,
                              logscale::Bool = false,
                              label::Union{Nothing,String} = nothing,
                              savepath::Union{Nothing,String} = nothing,
                              kwargs...) -> Plots.Plot

Create a scatter‑plot of *abundance* (y‑axis) against *growth_rate* (x‑axis).
The abundance column is inferred from `mode` unless you pass `abundance_col`.
No regression fit is drawn (by request).
"""
function plot_growth_vs_abundance(df::DataFrame;
                                  mode::Symbol = :mass,
                                  abundance_col::Union{Symbol,Nothing} = nothing,
                                  xlims = nothing,
                                  ylims = nothing,
                                  logscale::Bool = false,
                                  label = nothing,
                                  savepath = nothing,
                                  kwargs...)

    abundance_col = abundance_col === nothing ? _default_abundance_sym(mode) : abundance_col


    p = scatter(df.growth_rate, df[!, abundance_col];
                xlabel = "Growth rate (1/h)",
                ylabel = string(abundance_col),
                label  = label === nothing ? "" : label,
                kwargs...)

    xlims !== nothing && xlims!(p, xlims)
    ylims !== nothing && ylims!(p, ylims)
    logscale && yscale!(p, :log10)
    savepath !== nothing && savefig(p, savepath)
    return p
end

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Gene‑length bins plot                                                   ║
# ╚══════════════════════════════════════════════════════════════════════════╝

function plot_gene_length_analysis(
    gene_ranges;
    limitations   = [:C_lim],
    mode::Symbol  = :mass,
    data          = load_data(),
    xlims = (0,1),
    ylims = (0,1),
    logscale::Bool = false,
    savepath::Union{Nothing,String} = nothing,
    kwargs...)

    results = df_by_gene_length(gene_ranges, limitations, data,mode = mode)

    abundance_col = _default_abundance_sym(mode)

    for (range, dfs) in results
        p = scatter(
            dfs.growth_rate, dfs[!, abundance_col];
            title = "$range",
            xlabel = "Growth rate (1/h)",
            ylabel = string(abundance_col),
            legend = :topright)
        xlims !== nothing && xlims!(p, xlims)
        ylims !== nothing && ylims!(p, ylims)
        logscale && yscale!(p, :log10)
        savepath !== nothing && savefig(p, savepath)
        display(p)
    end
end

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Sector plot                                                             ║
# ╚══════════════════════════════════════════════════════════════════════════╝

function plot_sector_analysis(gene_list::Vector{String};
                              limitation::Symbol = :C_lim,
                              mode::Symbol       = :mass,
                              data               = load_data(),
                              xlims = nothing,
                              ylims = nothing,
                              logscale::Bool = false,
                              savepath::Union{Nothing,String} = nothing,
                              kwargs...)

    df = analyze_gene_list(gene_list, data, limitation; mode = mode)

    abundance_sym = _default_abundance_sym(mode)


    long = Symbol("Growthrate(1/h)")
    hasproperty(df, long) && rename!(df, long => :growth_rate)

    label_str = join(gene_list, ", ")

    return plot_growth_vs_abundance(df;
            mode = mode,
            abundance_col = abundance_sym,
            xlims = xlims,
            ylims = ylims,
            logscale = logscale,
            label = label_str,
            savepath = savepath,
            kwargs...)
end

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Quick helper to plot straight from a saved JLD2 file                    ║
# ╚══════════════════════════════════════════════════════════════════════════╝

function plot_from_file(path::AbstractString; kwargs...)
    @load path df
    return plot_growth_vs_abundance(df; kwargs...)
end

end # module Plotting
