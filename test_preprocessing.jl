# test_preprocessing.jl

include("src/utils.jl")
include("src/load_data.jl")
include("src/preprocessing.jl")

println("📦 Loading data...")
data = load_all_data()

println("🔬 Generating mass, number, and concentration DataFrames...")
mass_df, number_df, concentration_df = generate_mass_number_conc_dfs(
    data["gene_length"],
    data["mass_fractions2"]
)


df_Cali, df_C_lim, df_A_lim, df_R_lim = get_data_frames(data["samples2"][2])


println("✅ DataFrames created:")
println("  • mass_df:           ", size(mass_df))
println("  • number_df:         ", size(number_df))
println("  • concentration_df:  ", size(concentration_df))
println("  • df_Cali:           ", size(df_Cali))
println("  • df_C_lim:         ", size(df_C_lim))
println("  • df_A_lim:  ", size(df_A_lim))
println("  • df_R_lim:  ", size(df_R_lim))


println("\n🔢 Preview (first 5 rows of mass_df):")
display(first(mass_df, 5))


display(keys(data["protein_sectors"]))

data