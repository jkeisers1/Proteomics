# src/sector_by_length.jl

using DataFrames
using JLD2
using Dates
# Load dependencies
include("preprocessing.jl")    # provides mass/number/concentration data frames
include("sector_analysis.jl")  # provides analyze_gene_list, prepare_limitation_data, etc.

"""
    combine_limitations(
        limitations::Vector{Symbol},
        data::Dict{Symbol,Any}
    ) => DataFrame

Concatenate multiple growth-rate DataFrames (e.g., :C_lim, :A_lim) into one table.
Returned DataFrame has columns :SampleID and :Growthrate(1/h).
"""
function combine_limitations(limitations::Vector{Symbol}, data::Dict{Symbol,Any})
    growth_tables = [data[:growth][lim] for lim in limitations]
    return vcat(growth_tables...)
end

"""
    prepare_limitation_data(
        lim::Symbol,
        data::Dict{Symbol,Any};
        mode::Symbol = :mass
    ) => DataFrame

For a single limitation (e.g. :C_lim), pivot the abundance DataFrame into long form
and attach the growth rate. Columns:
- :Genename, :gene_length
- :SampleID, :value (abundance)
- :Growthrate(1/h)
"""
function prepare_limitation_data(lim::Symbol, data::Dict{Symbol,Any}; mode::Symbol = :mass)
    # 1) Extract growth rates for this limitation
    growth_df = data[:growth][lim]

    # 2) Select corresponding columns in the abundance matrix
    ab = data[mode]  # wide: rows=genes, cols=samples
    sample_cols = Symbol.(growth_df.SampleID)
    wide = select(ab, [:Genename, :gene_length, sample_cols...])

    # 3) Melt wide → long format: one row per gene per sample
    long = stack(wide, sample_cols; variable_name = :SampleID, value_name = :value)

    # 4) Join growth rates by SampleID
    return leftjoin(long, growth_df, on = :SampleID)
end

"""
    summarize_gene_list(
        gene_list::Vector{<:AbstractString},
        df_abund::DataFrame,
        df_growth::DataFrame
    ) => DataFrame

Take a list of genes, a wide abundance table, and a growth-rate table, then:
1) Subset the abundance table to those genes & samples in the growth table.
2) Sum abundance per sample.
3) Join to growth rates.
Result columns: :SampleID, :fraction, :Growthrate(1/h).
"""
function summarize_gene_list(gene_list, df_abund::DataFrame, df_growth::DataFrame; mode= mode)
    # Convert growth SampleIDs to symbols for column matching
    sample_syms = Symbol.(df_growth.SampleID)

    # Determine intersecting sample columns
    cols = intersect(Symbol.(names(df_abund)), sample_syms)

    # Subset to genes and those sample columns
    sub = select(df_abund, [:Genename, :gene_length, cols...])
    sub = filter(row -> row.Genename in gene_list, sub)

    # Sum each sample column yields total fraction per sample
    fractions = [sum(sub[!, c]) for c in cols]

    if mode == :concentration
        mode_col_name = String(mode)
    else
        mode_col_name = String(mode) * "_fraction"
    end

    # Build fraction table
    df_frac = DataFrame(;SampleID = String.(cols), Symbol(mode_col_name) => fractions)
    
    # Join with growth rates and return
    return leftjoin(df_frac, df_growth, on = :SampleID)
end

"""
    group_genes_by_length(
        df::DataFrame,
        length_ranges::Vector{Tuple{String, UnitRange}}
    ) => Dict{String, Vector{String}}

Bins all genes in `df` by their lengths according to `length_ranges`,
returning a Dict mapping each range label to the vector of gene names.
"""
function group_genes_by_length(df::DataFrame, length_ranges)
    groups = Dict{String, Vector{String}}()
    for (label, _) in length_ranges
        groups[label] = String[]
    end
    for row in eachrow(df)
        g, L = row.Genename, row.gene_length
        for (label, rng) in length_ranges
            if L in rng
                push!(groups[label], g)
                break
            end
        end
    end
    return groups
end

"""
    df_by_gene_length(
        gene_ranges,
        limitations::Vector{Symbol},
        measure::Symbol,
        data::Dict{Symbol,Any}
    ) => Dict{String, DataFrame}

For each gene-length bin in `gene_ranges`, under specified `limitations` and `measure` (:mass/:number/:concentration),
returns a Dict where keys are length-range labels and values are DataFrames with columns
:SampleID, :fraction, and :Growthrate(1/h).
"""
function df_by_gene_length(gene_ranges, limitations::Vector{Symbol}, data::Dict{Symbol,Any}; mode::Symbol = :mass)
    
    # Combine growth tables for all limitations
    growth_df = combine_limitations(limitations, data)
    # Abundance wide-form for the chosen measure
    df_abund = data[mode]

    # Bin genes by length
    length_groups = group_genes_by_length(df_abund, gene_ranges)

    results = Dict{Any, Any}()
    for (label, gene_list) in length_groups
        # For each bin, produce summarized DataFrame
        df = summarize_gene_list(gene_list, df_abund, growth_df, mode = mode)

        results[label] = df
    end
    return results 
end


"""
    save_gene_length_results(
        results::Dict{String, DataFrame};
        output_dir::String = "output",
        filename_prefix::String = "df_",
        filename_suffix::String = "_gene_length.jld2"
    )

Saves each DataFrame in `results` to `output_dir` as JLD2 files named
`<filename_prefix><label><filename_suffix>`.
"""
function save_gene_length_results(
    results::Dict{String, DataFrame},
    gene_ranges::Vector{Tuple{String, UnitRange{Int64}}};
    output_dir::String = "output",
    filename_prefix::String = "df_",
    filename_suffix::String = "_gene_length.jld2"
)
    final_dict = Dict("ranges" => gene_ranges,"data" => results)
    # Ensure output directory exists
    isdir(output_dir) || mkpath(output_dir)
    # Insert current date (YYYY-MM-DD) into the filename
    date_tag = Dates.format(Dates.now(), "yyyy-mm-dd")
    date_tag = replace(date_tag, "-" => "_")

    filename = joinpath(
        output_dir,
        filename_prefix * date_tag * filename_suffix
    )

    @save filename final_dict

end


"""
    print_available_protein_sectors(df::DataFrame)

Prints the unique sector labels present in a Proteinsector column of `df`.
"""
function print_available_protein_sectors(df::DataFrame)
    println("Available protein sectors: ", unique(df.Proteinsector))
end

"""
    filter_df_by_genelist(
        df::DataFrame,
        gene_list::Vector{<:AbstractString}
    ) => DataFrame

Return only the rows of `df` whose :Genename is in `gene_list`.
"""
function filter_df_by_genelist(df::DataFrame, gene_list::Vector{<:AbstractString})
    return filter(row -> row.Genename in gene_list, df)
end
