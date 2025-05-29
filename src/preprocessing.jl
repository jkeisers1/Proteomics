# src/preprocessing.jl

using DataFrames
using JLD2

"""
    create_combined_df(gene_length_df, mass_fractions_df) => DataFrame

Joins gene length and mass fraction data on :Genename.
"""
function create_combined_df(gene_length_df::Dict{Int, DataFrame}, mass_fractions_df::Dict{Int, DataFrame})
    return innerjoin(gene_length_df[1][!, [:Genename, :gene_length]], mass_fractions_df[2], on = :Genename)
end

"""
    convert_to_number_fraction(combined_df::DataFrame) => DataFrame

Converts protein mass fractions to number fractions.
"""
function convert_to_number_fraction(combined_df::DataFrame)
    number_df = deepcopy(combined_df)
    num_cols = names(number_df)[5:end]
    aa_length = number_df.gene_length ./ 3

    for col in num_cols
        L_tot = sum(combined_df[!, col] .* aa_length)
        number_df[!, col] .= combined_df[!, col] .* (L_tot ./ aa_length)
    end

    for col in num_cols
        total = sum(number_df[!, col])
        number_df[!, col] .= number_df[!, col] ./ total
    end

    return number_df
end

"""
    conversion_factor_to_conc(gene_name::String, df::DataFrame) => Float64

Calculates the conversion factor from mass fraction to concentration for a gene.
"""
function conversion_factor_to_conc(gene_name::String, df::DataFrame)
    row = findfirst(x -> x == gene_name, df[!, :Genename])
    gene_len = df[row, :gene_length]

    nt_to_aa = 3               # 3 nucleotides to one amino acid
    Da_to_fg = 1.66e-9         # Conversion between Dalton and femtogram
    ρ_cell = 314               # Cell density in fg/μm³
    m_aa = 110                 # Average mass of an amino acid in Da

    return ρ_cell / (m_aa * Da_to_fg * gene_len / nt_to_aa)
end

"""
    convert_to_concentrations(df::DataFrame) => DataFrame

Converts mass fractions to concentrations using gene-specific conversion factors.
"""
function convert_to_concentrations(df::DataFrame)
    conc_df = deepcopy(df)
    sample_cols = names(df)[5:end]

    for row in eachrow(conc_df)
        conv_factor = conversion_factor_to_conc(row.Genename, df)
        for col in sample_cols
            row[col] *= conv_factor
        end
    end

    return conc_df
end

"""
    generate_mass_number_conc_dfs(gene_length_df, mass_fractions_df)
        => (mass_df, number_df, concentration_df)

Creates and returns the mass, number, and concentration dataframes.
Also saves them as JLD2 files to the output folder.
"""
function generate_mass_number_conc_dfs(
    gene_length_df::Dict{Int, DataFrame},
    mass_fractions_df::Dict{Int, DataFrame};
    savepath::String = "output/"
)
    mass_df = create_combined_df(gene_length_df, mass_fractions_df)
    number_df = convert_to_number_fraction(mass_df)
    concentration_df = convert_to_concentrations(mass_df)

    mkpath(savepath)
    @save joinpath(savepath, "mass_df.jld2") mass_df
    @save joinpath(savepath, "number_df.jld2") number_df
    @save joinpath(savepath, "concentration_df.jld2") concentration_df

    return mass_df, number_df, concentration_df
end

"""
    get_data_frames(samples_df::DataFrame; savepath::String="output/")
        => (df_Cali, df_C_lim, df_A_lim, df_R_lim)

Splits the sample info by limitation group and saves the growth-rate mappings.
"""
function get_data_frames(samples_df::DataFrame; savepath::String = "output/")
    sampleID_growth_rate_df = samples_df[:, [1,2,3]]
    groups = sampleID_growth_rate_df[!, :Group]

    R_lim_ID, C_lim_ID, A_lim_ID, Cali_ID = Int[], Int[], Int[], Int[]
    limitations = unique(groups)

    for i in eachindex(groups)
        if groups[i] == limitations[1]
            push!(Cali_ID, i)
        elseif groups[i] == limitations[2]
            push!(C_lim_ID, i)
        elseif groups[i] == limitations[3]
            push!(A_lim_ID, i)
        elseif groups[i] == limitations[4]
            push!(R_lim_ID, i)
        end
    end

    df_Cali = sampleID_growth_rate_df[collect(Cali_ID[1]:1:Cali_ID[end]), [1,3]]
    df_R_lim = sampleID_growth_rate_df[collect(R_lim_ID[1]:1:R_lim_ID[end]), [1,3]]
    df_A_lim = sampleID_growth_rate_df[collect(A_lim_ID[1]:1:A_lim_ID[end]), [1,3]]
    df_C_lim = sampleID_growth_rate_df[collect(C_lim_ID[1]:1:C_lim_ID[end]), [1,3]]

    sort!(df_Cali, :SampleID)
    sort!(df_R_lim, :SampleID)
    sort!(df_C_lim, :SampleID)
    sort!(df_A_lim, :SampleID)

    mkpath(savepath)
    @save joinpath(savepath, "df_Cali.jld2") df_Cali
    @save joinpath(savepath, "df_C_lim.jld2") df_C_lim
    @save joinpath(savepath, "df_A_lim.jld2") df_A_lim
    @save joinpath(savepath, "df_R_lim.jld2") df_R_lim

    return df_Cali, df_C_lim, df_A_lim, df_R_lim
end

"""
    load_processed_data(output_dir::String="output/")

Loads preprocessed mass, number, concentration, and growth data from JLD2 files.
If not all files exist, runs preprocessing.
Returns a dictionary with keys :mass, :number, :concentration, and :growth (subkeys).
"""
function load_processed_data(output_dir::String = "output/")
    mass_path = joinpath(output_dir, "mass_df.jld2")
    number_path = joinpath(output_dir, "number_df.jld2")
    conc_path = joinpath(output_dir, "concentration_df.jld2")
    samples_path = joinpath(output_dir, "df_C_lim.jld2")

    if !(isfile(mass_path) && isfile(number_path) && isfile(conc_path) && isfile(samples_path))
        println("⚙️ Preprocessed files missing. Running preprocessing...")
        data_raw = load_all_data()
        generate_mass_number_conc_dfs(data_raw["gene_length"], data_raw["mass_fractions2"])
        get_data_frames(data_raw["samples2"][2])
    end


    data = Dict{Symbol, Any}()

    @load mass_path mass_df
    @load number_path number_df
    @load conc_path concentration_df

    @load joinpath(output_dir, "df_Cali.jld2") df_Cali
    @load joinpath(output_dir, "df_C_lim.jld2") df_C_lim
    @load joinpath(output_dir, "df_A_lim.jld2") df_A_lim
    @load joinpath(output_dir, "df_R_lim.jld2") df_R_lim
    
            # ── standardise the column right away ───────────────────────────────
    for df in (df_Cali, df_C_lim, df_A_lim, df_R_lim, mass_df, number_df, concentration_df)
        standardise_growth_rate!(df)
    end
    
    data[:mass] = mass_df
    data[:number] = number_df
    data[:concentration] = concentration_df
    data[:growth] = Dict(
        :Cali => df_Cali,
        :C_lim => df_C_lim,
        :A_lim => df_A_lim,
        :R_lim => df_R_lim,
    )

    return data
end
