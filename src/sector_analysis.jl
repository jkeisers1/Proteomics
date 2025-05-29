# sector_analysis.jl

using DataFrames
using JLD2

include("utils.jl")
include("load_data.jl")
include("preprocessing.jl")


"""
    analyze_gene_list(gene_list::Vector{String}, data::Dict{Symbol, Any},
                      condition::Symbol, mode::Symbol = :mass) => DataFrame

Returns a DataFrame with total fraction (or concentration) of a gene list under a specific condition,
joined with growth rate.
"""
function analyze_gene_list(gene_list::AbstractVector{<:AbstractString},
                           data::Dict{Symbol, Any},
                           condition::Symbol; 
                           mode::Symbol = :mass)

    df = data[mode]   # :mass, :number, or :concentration
    growth_df = data[:growth][condition]  # :C_lim, :A_lim, etc.
    
    filtered_df = df[in.(df.Genename, Ref(gene_list)), :]
    sample_ids = growth_df.SampleID
    matching_cols = Symbol.(intersect(names(filtered_df), sample_ids))

    total = [sum(filtered_df[!, col]) for col in matching_cols]

    if mode == :concentration
        mode_col_name = String(mode)
    else
        mode_col_name = String(mode) * "_fraction"
    end

    result_df = DataFrame(;SampleID = String.(matching_cols), Symbol(mode_col_name) => total)


    return leftjoin(result_df, growth_df, on = :SampleID)
end

"""
    load_gene_list_from_txt(path::String) => Vector{String}

Loads a gene list from a plain text file, one gene per line.
Strips whitespace and skips empty lines.
"""
function load_gene_list_from_txt(path::String)
    file_text = read(path, String)
    return [strip(line) for line in split(file_text, '\n') if !isempty(strip(line))]
end