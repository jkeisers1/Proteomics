# src/utils.jl

using DataFrames

"""
    standardise_growth_rate!(df) -> DataFrame

If a column called `\"Growthrate(1/h)\"` is present, rename it to `:growth_rate`.
Returns the same DataFrame (mutated in place).
"""
function standardise_growth_rate!(df::DataFrame)
    long = Symbol("Growthrate(1/h)")
    hasproperty(df, long) && rename!(df, long => :growth_rate)
    return df
end

"""
    remove_spaces(s::AbstractString) -> String
Remove all whitespace from a string.
"""
remove_spaces(s::AbstractString) = replace(s, r"\s+" => "")

"""
    change_minus_symbol(s::AbstractString) -> String
Replaces hyphens with underscores.
"""
change_minus_symbol(s::AbstractString) = replace(s, "-" => "_")

"""
    clean_up_df(df::DataFrame) -> DataFrame
Clean up DataFrame:
- Drops columns beyond the first contiguous non-missing header
- Sets cleaned header names
- Replaces missing with 0
- Removes whitespace and hyphens in string columns
"""
function clean_up_df(df::DataFrame)
    first_missing_col = length(findall(x -> ismissing(x) == 0, first(df)))
    df = df[!, 1:first_missing_col]

    # Clean column headers
    new_header = first(df)
    string_header = [String(h) for h in new_header]
    clean_header = [Symbol(change_minus_symbol(remove_spaces(s))) for s in string_header]
    rename!(df, clean_header, makeunique=true)

    df = df[2:end, :]
    df = coalesce.(df, 0)

    # Clean string columns
    string_columns = [col for col in names(df) if eltype(df[!, col]) == String]
    for col in string_columns
        df[!, col] = replace.(df[!, col], r"\s+" => "")
        df[!, col] = replace.(df[!, col], "-" => "_")
    end

    return df
end

"""
    strip_non_empty(lines::Vector{String}) -> Vector{String}
Strips whitespace from each line and drops empty ones.
"""
strip_non_empty(lines::Vector{String}) = [strip(l) for l in lines if !isempty(strip(l))]