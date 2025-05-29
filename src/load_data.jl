# src/load_data.jl

using CSV
using DataFrames
using XLSX
using DelimitedFiles

include(joinpath(@__DIR__, "utils.jl"))

"""
    read_sheets(xlsx_file_path::String) => Dict{Int, DataFrame}

Reads all sheets from an Excel file, cleans them, and returns a dictionary of DataFrames.
"""
function read_sheets(xlsx_file_path::String)
    xlsx_file = XLSX.readxlsx(xlsx_file_path)
    sheets = Dict{Int, DataFrame}()
    i = 1

    for sheet_name in XLSX.sheetnames(xlsx_file)
        sheet_data = xlsx_file[sheet_name]
        sheets[i] = DataFrame(sheet_data[:], :auto)
        i += 1
    end

    cleaned_sheets = Dict(i => clean_up_df(sheet) for (i, sheet) in sheets)
    return cleaned_sheets
end

"""
    load_all_data(data_path::String="data/") => Dict{String, Any}

Loads and prepares all relevant datasets from the data folder.
Returns a dictionary with meaningful keys.
"""
function load_all_data(data_path::String = "data/")
    return Dict(
        "mass_fractions1" => read_sheets(data_path * "absolute_mass_fractions1.xlsx"),
        "mass_fractions2" => read_sheets(data_path * "absolute_mass_fractions2.xlsx"),
        "samples1" => read_sheets(data_path * "samples1.xlsx"),
        "samples2" => read_sheets(data_path * "samples2.xlsx"),
        "protein_sectors" => read_sheets(data_path * "protein_sectors.xlsx"),
        "gene_length" => read_sheets(data_path * "gene_length.xlsx"),
        "constitutive_genes" => CSV.read(data_path * "constitutively_expressed_genes.csv", DataFrame),
        "mRNA_fractions" => read_sheets(data_path * "mRNA_concentration.xlsx"),
        "ribosomal_genes" => strip_non_empty(readlines(data_path * "ribosomal_genes_e_coli.txt"))
    )
end
