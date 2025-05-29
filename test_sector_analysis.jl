# test_sector_analysis.jl

include("src/sector_analysis.jl")

println("📦 Loading processed data or triggering preprocessing...")
data = load_processed_data()

println("📄 Loading gene list from file...")
gene_list = load_gene_list_from_txt("data/ribosomal_genes_e_coli.txt")

println("🔬 Analyzing gene list under carbon limitation...")
result_df = analyze_gene_list(["rmf"], data, :C_lim; mode = :number)

println("✅ Analysis complete. Preview:")
display(result_df)
