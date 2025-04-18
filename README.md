# RNA-seq Pipeline Assistant

A comprehensive tool for analyzing RNA-seq data, identifying differentially expressed genes (DEGs), and performing pathway enrichment analysis.

## Overview

This pipeline takes normalized RNA-seq data and performs a complete analysis workflow, generating various visualizations and identifying biologically significant patterns. The tool features both a command-line interface and a user-friendly Streamlit web application, making it accessible to both bioinformaticians and researchers with limited programming experience.

## Features

- **Data Processing**: Creates standardized metadata from normalized counts
- **Quality Control**: Generates boxplots, correlation heatmaps, and PCA plots
- **Differential Expression Analysis**: Processes DESeq2 files to identify DEGs
- **Visualization**: Creates MA plots, volcano plots, and DEG heatmaps
- **Ortholog Mapping**: Maps genes to human orthologs for cross-species analysis
- **Enrichment Analysis**: Performs GO and KEGG pathway enrichment analysis
- **Interactive UI**: Streamlit-based web interface for easy parameter configuration
- **AI Assistant**: Integrated AI functionality to help with parameter selection

## Installation

### Prerequisites
- Python 3.8+
- Required Python packages (install via pip):
  ```
  pandas
  numpy
  matplotlib
  seaborn
  scipy
  statsmodels
  scikit-learn
  streamlit
  click
  requests
  gseapy
  mygene
  pydeseq2
  openai
  ```

### Setup
1. Clone the repository
   ```bash
   git clone https://github.com/yourusername/rnaseq-pipeline-assistant.git
   ```
2. Navigate to the project directory
   ```bash
   cd rnaseq-pipeline-assistant
   ```
3. Install dependencies
   ```bash
   pip install -r requirements.txt
   ```
4. Set up OpenAI API key (for AI assistant functionality)
   ```bash
   export OPENAI_API_KEY="your-api-key"
   ```

## Usage

### Command Line Interface

The pipeline can be run using the CLI with two main commands:

1. **Setup**: Initialize the directory structure
   ```bash
   python pipeline_CLI.py --norm-file "path/to/normalized_data.csv" --pw-data "path/to/pairwise_data" --base-dir "output_directory" --num-replicates 3 --conditions "Control" --conditions "Treatment1" --conditions "Treatment2" setup
   ```

2. **Run Analysis**: Execute the full analysis pipeline
   ```bash
   python pipeline_CLI.py --norm-file "path/to/normalized_data.csv" --pw-data "path/to/pairwise_data" --base-dir "output_directory" --num-replicates 3 --conditions "Control" --conditions "Treatment1" --conditions "Treatment2" run-all --model-organism "drerio" --pw-interest "Treatment1_vs_Control" --pw-interest "Treatment2_vs_Control" --log2fc-threshold 1.0 --padj-threshold 0.05 --enrich-sig-cutoff 0.05
   ```

### Web Interface

For a more user-friendly experience, run the Streamlit app:

```bash
streamlit run app_v2.py
```

This will open a web interface where you can:
- Input all parameters through form fields
- Generate and review pipeline commands
- Execute commands directly from the interface
- View real-time command output

## Parameters

### Setup Parameters
- **norm-file**: Path to the normalized data file (CSV format)
- **pw-data**: Directory containing pairwise comparison files
- **base-dir**: Base output directory for results
- **conditions**: Experimental conditions (specify multiple with repeated flags)
- **num-replicates**: Number of replicates per condition

### Analysis Parameters
- **model-organism**: Model organism code (e.g., "drerio" for zebrafish)
- **pw-interest**: Pairwise comparisons of interest (e.g., "Treatment_vs_Control")
- **log2fc-threshold**: Log2 fold change threshold for DEG identification (default: 1.0)
- **padj-threshold**: Adjusted p-value threshold (default: 0.05)
- **enrich-sig-cutoff**: Significance cutoff for enrichment analysis (default: 0.05)

### Optional Flags
Skip specific analysis steps with these flags:
- **--skip-boxplots**: Skip boxplot generation
- **--skip-correlation-heatmap**: Skip correlation heatmap
- **--skip-pca**: Skip PCA analysis
- **--skip-ma-plots**: Skip MA plot generation
- **--skip-volcano-plots**: Skip volcano plot generation
- **--skip-heatmap**: Skip DEG heatmap
- **--skip-go**: Skip GO enrichment analysis
- **--skip-kegg**: Skip KEGG enrichment analysis

## Output Structure

The pipeline creates a standardized directory structure:

```
base_dir/
├── data/
│   ├── DESeq2/
│   ├── ortholog_mapping/
│   └── GO_and_KEGG_enrichments/
├── plots/
│   ├── boxplots/
│   ├── correlation_heatmap/
│   ├── PCA/
│   ├── MA_plots/
│   ├── volcano_plots/
│   ├── DEG_heatmap/
│   ├── GO_enrichment/
│   └── KEGG_enrichment/
└── results/
    ├── DEGs/
    ├── GO_enrichment/
    └── KEGG_enrichment/
```

## File Descriptions

- **app_v2.py**: Streamlit web application interface that  
- **ai_assistant.py**: AI functionality that collects parameters from the user in natural language and transforms them into valid commands 
- **pipeline_CLI.py**: Command-line interface for the pipeline
- **DEGpipeline_v6.py**: Core pipeline functionality and analysis modules

## Example Workflow

1. Prepare normalized RNA-seq data in CSV format
2. Prepare pairwise comparison files 
3. Run the setup command to create directory structure
4. Run the analysis command with appropriate parameters
5. Examine the generated plots and results files
6. Interpret biological significance of DEGs and enriched pathways

## License

This project hsa no special license. Thus, the default copyright laws apply, meaning that the owner retain all rights to the source code and no one may reproduce, distribute, or create derivative works from this work unless explicit permission is granted. Please contact me if you would like to use this pipeline. 

## Contact

For questions or support, please contact [i.banwait7@gmail.com].

## Acknowledgments

- This pipeline uses several open-source libraries including pandas, matplotlib, seaborn, and gseapy
- Ortholog mapping is performed using the g:Profiler API
- Pathway enrichment analysis uses the Enrichr API through gseapy
