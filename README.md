# 🧬 Gene Sector Analysis — Julia Pipeline

This repository provides a Julia workflow to analyze **gene abundance, protein sectors, and growth-rate relationships** under various nutrient limitations.  
It loads raw omics data, preprocesses it into usable forms (**mass**, **number**, and **concentration**), and generates interpretable plots across **gene sectors** and **gene-length bins**.

---

## 📁 Project Structure

```text
.
├── main.jl
├── src/
│   ├── load_data.jl
│   ├── preprocessing.jl
│   ├── sector_analysis.jl
│   ├── sector_by_length.jl
│   ├── plotting.jl
│   └── utils.jl
├── plotting_test.jl
├── data/
│   ├── absolute_mass_fractions1.xlsx
│   ├── absolute_mass_fractions2.xlsx
│   ├── samples1.xlsx
│   ├── samples2.xlsx
│   ├── protein_sectors.xlsx
│   ├── gene_length.xlsx
│   ├── constitutively_expressed_genes.csv
│   ├── mRNA_concentration.xlsx
│   └── ribosomal_genes_e_coli.txt
└── output/
    ├── mass_df.jld2
    ├── number_df.jld2
    ├── concentration_df.jld2
    └── df_*_gene_length.jld2
```

---

## ⚙️ Installation

1. Install **Julia ≥ 1.9**.
2. In the Julia REPL:
```julia
using Pkg
Pkg.add(["DataFrames", "CSV", "XLSX", "DelimitedFiles", "JLD2", "Plots", "GLM"])
```

---

## ▶️ Usage

### Run the pipeline
```bash
julia main.jl
```

This will:
- Load raw input files from `data/`  
- Generate preprocessed datasets for **mass**, **number**, **concentration**  
- Save outputs in `output/` (as `.jld2`)  
- Produce example plots

### Example — analyze a gene/sector
```julia
using .Plotting
data = load_processed_data()
Plotting.plot_sector_analysis(["hpf"], limitation = :C_lim, mode = :concentration, data = data)
```

### Example — group by gene length
```julia
gene_length = [
    ("0:500", 1:499),
    ("500:999", 500:999),
    ("1000:2000", 1000:2000),
    ("2000:5000", 2000:5000)
]
Plotting.plot_gene_length_analysis(gene_length, limitations = [:R_lim], mode = :mass, data = data)
```

---

## 🧩 Modules Overview

| File | Purpose |
|------|---------|
| `src/utils.jl` | Data cleaning and growth‑rate standardization |
| `src/load_data.jl` | Reads and cleans Excel/CSV inputs |
| `src/preprocessing.jl` | Builds mass/number/concentration tables and saves `.jld2` |
| `src/sector_analysis.jl` | Sums abundance for gene lists under conditions |
| `src/sector_by_length.jl` | Bins genes by length and aggregates abundance |
| `src/plotting.jl` | Plotting utilities (growth vs abundance, length bins) |
| `main.jl` | Example workflow entry point |

---

## 📊 Outputs

- `output/mass_df.jld2`, `output/number_df.jld2`, `output/concentration_df.jld2`  
- `output/df_<date>_gene_length.jld2` (summaries by gene-length bin)  
- Figures saved if `savepath` is provided in plotting functions

---

## 🧾 License

MIT © 2025
