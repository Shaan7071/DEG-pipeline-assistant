# %% [markdown]
# # DEG Pipeline

# %% [markdown]
# This pipeline takes normalized RNA-seq data and outputs diffrentially expressed genes and pathways. 
# 
# 
# Version 7: updated figure generations to improve readability of outputs. 
# 
# Author: Ishaan Banwait
# 
# Last Updated: 2025-09-14

# %%
# Import libraries
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as scipy_stats
import statsmodels.stats.multitest as multi
import requests
import gseapy as gp
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import json
import gseapy as gp
import mygene
from matplotlib.colors import LinearSegmentedColormap
import itertools
import statsmodels.stats.multitest as smm
import textwrap

# %% [markdown]
# ## Metadata Creation

# %%
'''
(create_metadata) takes normalized RNA-seq data and extracts the relevant columns to create a metadata DataFrame with standardized columns: Sample_ID, Condition, and Replicate.

'''

def create_metadata(norm_file_name, conditions, num_replicates, base_dir):
    # Load the normalized data including the first row as header
    norm_data = pd.read_csv(norm_file_name, header=None)

    # Extract sample names from the first row
    sample_names = norm_data.iloc[0,1:].tolist()

    # Generate metadata based on provided conditions and replicates
    metadata = pd.DataFrame({
        'Sample_ID': sample_names,
        'Condition': [condition for condition in conditions for i in range(num_replicates)],
        'Replicate': [rep for i in conditions for rep in range(1, num_replicates + 1)]
    })
    
    # export the metadata to a CSV file in the data folder
    metadata.to_csv(f"{base_dir}/data/metadata.csv", index=False)

# %% [markdown]
# ## Quality Control: Boxplots of Normalized Counts

# %%
'''
(create_boxplots) takes the normalized RNA-seq data and generates boxplots of normalized counts for each set of replicates.
'''

def create_boxplots(norm_file_name, conditions, num_replicates, base_dir):
    # Load the normalized data, excluding the first row as header
    norm_data = pd.read_csv(norm_file_name)
    # drop the first column (gene names) 
    norm_data = norm_data.drop(columns=norm_data.columns[0])

    set = 0
    multiplier = 1

    for i in conditions:
        plt.figure(figsize=(12, 10))
        sns.boxplot(data=norm_data.iloc[:,multiplier*num_replicates-num_replicates:multiplier*num_replicates])
        plt.title('Boxplots of Normalized Counts Across Replicates for ' + conditions[set] + ' Condition', fontsize=18, fontweight='bold')
        plt.xlabel('Samples', fontsize=16, fontweight='bold')
        plt.ylabel('Normalized Counts', fontsize=16, fontweight='bold')
        plt.xticks(rotation=45, fontsize=14)
        plt.yticks(fontsize=14)

        plt.tight_layout()

        # save to appropriate sub
        plt.savefig(f"{base_dir}/plots/boxplots/boxplot_normalized_counts_{conditions[set].replace(" ", "_")}.png")
        plt.close()

        set +=1
        multiplier += 1
        


# %% [markdown]
# ## Sample Correlation Heatmap

# %%
'''
(create_correlation_heatmap) takes the normalized RNA-seq data and generates a correlation heatmap of the samples.'''

def create_correlation_heatmap(norm_file_name, base_dir):
    # Load the normalized data assuming header is in the first row
    norm_data = pd.read_csv(norm_file_name)
    # drop the first column (gene names) 
    norm_data = norm_data.drop(columns=norm_data.columns[0])

    # Calculate the correlation matrix between samples (columns)
    corr_matrix = norm_data.corr(method='pearson')

    # Create the heatmap with hierarchical clustering
    plt.figure(figsize=(12, 10))
    clustered_heatmap = sns.clustermap(
        corr_matrix,
        cmap="RdBu",  # Blue for similar (high correlation), red for dissimilar (low correlation)
        figsize=(12, 10),
        xticklabels=True,
        yticklabels=True,
        annot=False,  # Set to True if you want to see correlation values in cells
        cbar_kws={'label': 'Pearson Correlation'},
        method='average',  # Clustering method
        metric='euclidean'  # Distance metric for clustering
    )

    # Adjust the plot
    plt.title('Sample Correlation Heatmap', pad=50)

    # save to appropriate sub
    plt.savefig(f'{base_dir}/plots/correlation_heatmap/correlation_heatmap.png')
    plt.close()



# %% [markdown]
# ## Principal Component Analysis (Dimensionality Reduction)

# %%
"""
(create_pca) takes the normalized RNA-seq data and generates a PCA plot to visualize the relationship between samples.
"""

def create_pca(norm_file_name, base_dir):
    # Load the normalized data assuming header is in the first row
    norm_data = pd.read_csv(norm_file_name)
    # drop the first column (gene names) 
    norm_data = norm_data.drop(columns=norm_data.columns[0])

    # Transpose if needed (PCA expects samples as rows, genes as columns)
    if norm_data.shape[0] > norm_data.shape[1]:  # More rows than columns suggests genes as rows
        counts_data = norm_data.T

    # Standardize the data to ensure that each feature (gene) contributes equally to the analysis
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(counts_data)

    # Perform PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(scaled_data)

    # Create a DataFrame with the principal components
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

    # Add metadata information
    metadata = pd.read_csv(f"{base_dir}/data/metadata.csv")
    pca_df['condition'] = metadata['Condition']
    pca_df['replicate'] = metadata['Replicate']

    # Create a plot
    plt.figure(figsize=(10, 8))

    # Dynamically create markers dictionary from unique conditions in metadata
    unique_conditions = metadata['Condition'].unique()
    marker_options = ['o', '^', '*', 'v', 'X', 'D', 's', 'p', 'h',]
    markers = {condition: marker_options[i % len(marker_options)] 
            for i, condition in enumerate(unique_conditions)}

    # Dynamically create colors dictionary from unique replicates
    unique_replicates = metadata['Replicate'].unique()
    color_options = ['red', 'gold', 'green', 'cyan', 'blue', 'magenta', 'orange', 'purple', 'brown', 'gray']
    colors = {replicate: color_options[i % len(color_options)] 
            for i, replicate in enumerate(unique_replicates)}

    # Plot each point
    for condition in markers.keys():
        for replicate in colors.keys():
            subset = pca_df[(pca_df['condition'] == condition) & (pca_df['replicate'] == replicate)]
            if not subset.empty:
                plt.scatter(subset['PC1'], subset['PC2'], 
                        marker=markers[condition], 
                        c=colors[replicate],
                        s=80)

    # Add legends for condition and replicate with headers
    from matplotlib.lines import Line2D

    # Create a list for all legend elements
    legend_elements = []

    # Add condition header
    legend_elements.append(Line2D([0], [0], color='w', marker='', linestyle='none', 
                                label='condition', markerfacecolor='w'))

    # Add condition elements
    for condition, marker in markers.items():
        legend_elements.append(Line2D([0], [0], marker=marker, color='w', 
                                    markerfacecolor='black', markersize=10, label=condition))
        
    # Add empty line for spacing
    legend_elements.append(Line2D([0], [0], color='w', marker='', linestyle='none', 
                                label='', markerfacecolor='w'))

    # Add replicate header
    legend_elements.append(Line2D([0], [0], color='w', marker='', linestyle='none', 
                                label='replicate', markerfacecolor='w'))

    # Add replicate elements
    for replicate, color in colors.items():
        legend_elements.append(Line2D([0], [0], marker='o', color=color, 
                                    markersize=10, label=replicate))

    # Create the legend with all elements and format headers
    legend = plt.legend(handles=legend_elements, 
                       loc='center left', 
                       bbox_to_anchor=(1, 0.5))
    
    # Make headers bold and size 14
    legend_texts = legend.get_texts()
    for text in legend_texts:
        if text.get_text() in ['condition', 'replicate']:
            text.set_fontsize(14)
            text.set_fontweight('bold')
        else:
            text.set_fontsize(12)  # Keep other items readable

    # Show the percentage of variance explained by each component
    explained_variance = pca.explained_variance_ratio_ * 100
    plt.xlabel(f'PC1 ({explained_variance[0]:.2f}%)', fontsize=16, fontweight='bold')
    plt.ylabel(f'PC2 ({explained_variance[1]:.2f}%)', fontsize=16, fontweight='bold')
    plt.title('PCA of Gene Expression Data', fontsize=18, fontweight='bold')
    plt.grid(True)

    plt.tight_layout()
    # save to appropriate sub
    plt.savefig(f'{base_dir}/plots/PCA/pca_plot.png')
    plt.close()


# %% [markdown]
# ## Pairwise Comparisons

# %% [markdown]
# Create the function to read the pairwise comparison files. Statistical testing is already completed. DEG filtering is done as files are read in. A summary of results is generated. 

# %%

'''
(process_deseq2_files) takes a list of DESeq2 result files and processes them to categorize genes as upregulated, downregulated, or not significant based on specified thresholds. It also generates a summary of the results for each comparison and saves the processed files to a new directory.
'''

def process_deseq2_files(input_dir_pw_data, base_dir, log2fc_threshold, padj_threshold):
    
    results = {}
    summary = {}
    
    pw_file_names = os.listdir(input_dir_pw_data)
    for file_name in pw_file_names:
        comparison = file_name.replace(".csv", "")
        file_path = os.path.join(input_dir_pw_data, file_name)
        
        try:
            # Check if file exists
            if Path(file_path).is_file():
                # Read the file into a DataFrame
                df = pd.read_csv(file_path)
                
                # Verify that the file has the expected columns
                expected_columns = ['ensembl_gene_id', 'baseMean', 'log2FoldChange', 'pvalue', 'padj']
                missing_columns = [col for col in expected_columns if col not in df.columns]
                
                if missing_columns:
                    print(f"Warning: Missing expected columns in {file_name}: {missing_columns}")
                
                # Add 'Expression' column to categorize genes as upregulated, downregulated, or not significant
                df['Expression'] = 'Not Significant'
                significant_mask = (abs(df['log2FoldChange']) >= log2fc_threshold) & (df['padj'] < padj_threshold)
                df.loc[significant_mask & (df['log2FoldChange'] > 0), 'Expression'] = 'Upregulated'
                df.loc[significant_mask & (df['log2FoldChange'] < 0), 'Expression'] = 'Downregulated'
                
                # Count significant genes
                sig_up = df[(df['Expression'] == 'Upregulated')].shape[0]
                sig_down = df[(df['Expression'] == 'Downregulated')].shape[0]
                
                # Save processed file to appropriate sub
                processed_file_name = file_name.replace(".csv", "_processed.csv")
                df.to_csv(f'{base_dir}/data/DESeq2/{processed_file_name}', index=False)
                
                # Store results
                results[comparison] = df
                summary[comparison] = {
                    'total_genes': len(df),
                    'significant_genes': sig_up + sig_down,
                    'upregulated': sig_up,
                    'downregulated': sig_down
                }
                
                print(f"Successfully processed: {file_name}")
                print(f"  - Total genes: {len(df)}")
                print(f"  - Significant genes: {sig_up + sig_down} (Up: {sig_up}, Down: {sig_down})")
                
            else:
                print(f"File not found: {file_path}")
                results[comparison] = None
                summary[comparison] = None
                
        except Exception as e:
            print(f"Error processing {file_name}: {str(e)}")
            results[comparison] = None
            summary[comparison] = None
    
    # Save summary as CSV
    summary_df = pd.DataFrame.from_dict(summary, orient='index')
    summary_df.index.name = 'Comparison'
    summary_df.to_csv(f'{base_dir}/results/DEGs/pairwise_comparisons_summary.csv')
    print(f"Created summary table with {len(summary_df)} rows")

    # Print overall summary
    print("\nOverall Summary:")
    for comparison, stats in summary.items():
        if stats:
            print(f"{comparison}: {stats['significant_genes']} significant genes " +
                    f"({stats['upregulated']} up, {stats['downregulated']} down)")
        else:
            print(f"{comparison}: Failed to process")


# %% [markdown]
# ## MA Plots

# %% [markdown]
# Create the function that will generate the MA plots.

# %%
# Create MA plots from the processed DESeq2 files

'''
(create_ma_plots) reads processed DESeq2 files and generates an MA plot for it.
'''

def create_ma_plots(input_dir_pw_data, base_dir):

    pw_file_names = os.listdir(input_dir_pw_data)

    # Process each file
    for file_name in pw_file_names:
        file_name = file_name.replace('.csv', '_processed.csv')
        
        # Read the processed data file
        df = pd.read_csv(f'{base_dir}/data/DESeq2/{file_name}')
        
        # Extract necessary columns
        # Assuming the processed files contain 'baseMean' and 'log2FoldChange'
        A = np.log2(df['baseMean'])  # Average expression (log2 scale)
        M = df['log2FoldChange']     # Log2 fold change
        
        # Create the MA plot
        plt.figure(figsize=(10, 8))
        
        # Plot all points in gray
        plt.scatter(A, M, color='gray', alpha=0.5, s=10)
        
        # Highlight significant genes (assuming padj < 0.05 and |log2FC| > 1)
        significant = (df['padj'] < 0.05) & (abs(df['log2FoldChange']) > 1)
        plt.scatter(A[significant], M[significant], color='red', alpha=0.7, s=15)
        
        # Add horizontal line at M=0
        plt.axhline(y=0, color='blue', linestyle='--', alpha=0.7)
        
        # Add labels and title
        plt.xlabel('A (Average log2 expression)', fontsize=16, fontweight='bold')
        plt.ylabel('M (Log2 fold change)', fontsize=16, fontweight='bold')
        plt.title(f'MA Plot: {os.path.basename(file_name).replace("_processed.csv", "")}', fontsize=18, fontweight='bold')

        # Set tick font size
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        
        # Add a legend
        plt.legend(['Not significant', 'Significant (padj<0.05, |log2FC|>1)', 'M=0'], fontsize=14)
        
        # Save the plot
        plt.tight_layout()
        plt.savefig(f'{base_dir}/plots/MA_plots/{(file_name).replace("_processed.csv", "")}_MA_plot.png')
        plt.close()


# %% [markdown]
# ## Volcano Plots

# %% [markdown]
# Create the function that will generate the volcano plots.

# %%
# Create volcano plots from the processed DESeq2 files

'''
(create_volcano_plot) reads processed DESeq2 files and generates a volcano plot for it.
'''


def create_volcano_plots(input_dir_pw_data, base_dir, log2fc_threshold, padj_threshold):
    # Process each file

    pw_file_names = os.listdir(input_dir_pw_data)
    
    for file_name in pw_file_names:
        file_name = file_name.replace('.csv', '_processed.csv')
        
        # Read the processed data file
        df = pd.read_csv(f'{base_dir}/data/DESeq2/{file_name}')

        # Before plotting, filter out rows with NaN values
        df = df.dropna(subset=['padj', 'log2FoldChange'])

        # Create the volcano plot
        plt.figure(figsize=(10, 8))
        
        # Plot points with different colors based on expression status
        colors = {'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'gray'}
        
        # Plot each category separately for better control
        for category, color in colors.items():
            subset = df[df['Expression'] == category]
            plt.scatter(subset['log2FoldChange'], 
                    -np.log10(subset['padj']), 
                    c=color, 
                    alpha=0.6, 
                    s=25,
                    label=category)
            
        # Add threshold lines
        plt.axhline(y=-np.log10(padj_threshold), color='gray', linestyle='--', linewidth=1)
        plt.axvline(x=log2fc_threshold, color='gray', linestyle='--', linewidth=1)
        plt.axvline(x=-log2fc_threshold, color='gray', linestyle='--', linewidth=1)

        # Customize the plot
        plt.title(f'Volcano Plot: {os.path.basename(file_name).replace("_processed.csv", "")}', fontsize=18, fontweight='bold')
        plt.xlabel('Log2 Fold Change', fontsize=16, fontweight='bold')
        plt.ylabel('-Log10 Adjusted P-value', fontsize=16, fontweight='bold')
        
        # Add legend
        plt.legend(loc='upper right', fontsize=14)
        
        # Add grid for better readability
        plt.grid(True, linestyle='--', alpha=0.3)
        plt.tick_params(axis='both', which='major', labelsize=14)

        # matplotlib is automatically determining the axis limits based on the data range of the values being plotted.
        
        # Save the plot
        plt.tight_layout()
        plt.savefig(f'{base_dir}/plots/volcano_plots/{(file_name).replace("_processed.csv", "")}_volcano_plot.png')
        plt.close()
    

# %% [markdown]
# ## Results Table of Top DEGs

# %% [markdown]
# Create the function to query gProfiler for otholog mapping

# %%
# Map zebrafish genes to human orthologs using gProfiler g:Orth tool

'''
(map_to_human_orthologs) maps a list of zebrafish genes to human orthologs using the g:Orth tool from g:Profiler.
'''

def map_to_human_orthologs(pw_interest, model_organism, base_dir):

    for file in pw_interest:
        df = pd.read_csv(os.path.join(base_dir, f'data/DESeq2/{file}_processed.csv'))

        # Sort the data by log2FoldChange
        df = df.sort_values('log2FoldChange', ascending=False)
        
        # Extract Ensembl gene IDs 
        gene_list = df['ensembl_gene_id'].tolist()

        # gProfiler API endpoint for ortholog mapping (g:Orth)
        url = "https://biit.cs.ut.ee/gprofiler/api/orth/orth/"
        
        # Prepare the request payload
        payload = {
            'organism': model_organism,
            'target': 'hsapiens',
            'query': gene_list
        }
        
        # Make the API request
        response = requests.post(url, json=payload)
        
        if response.status_code == 200:
            result = response.json()
            # Convert to DataFrame
            if result and 'result' in result:
                results_df = pd.DataFrame(result['result'])
                # Rename columns to match desired output format
                results_df = results_df.rename(columns={
                    'incoming': 'input',
                    'name': 'ortholog_name',
                    'ensg': 'ortholog_ensg',
                    'description': 'description'
                })
                
            else:
                return pd.DataFrame()
        else:
            print(f"Error in API request: {response.status_code}")
        
        # Merge with original data, handling cases of multiple orthologs for a single model organism gene
        result_df = pd.merge(
            results_df,
            df,
            left_on='input',  # 'input' column from orthologs df
            right_on='ensembl_gene_id',  # 'ensemblgeneid' column from original df
            how='left'
            )

        # Drop rows with missing ortholog Ensembl ID
        result_df = result_df[result_df['ortholog_ensg'] != "N/A"]

        # Save only significant genes to a seperate df
        result_df_sig = result_df[result_df['Expression'] != 'Not Significant']

        # Select and reorder columns
        final_df = result_df[[
            'input',  # zebrafish Ensembl ID
            'ortholog_name', 
            'ortholog_ensg', 
            'description', 
            'log2FoldChange',
            'Expression', 
            'padj'
            ]]  
        
        final_df_sig = result_df_sig[[
            'input',  # zebrafish Ensembl ID
            'ortholog_name', 
            'ortholog_ensg', 
            'description',
            'log2FoldChange', 
            'Expression', 
            'padj'
            ]]  
        
        # Save to a CSV file
        final_df.to_csv(f'{base_dir}/data/ortholog_mapping/{file}_orthologs.csv', index=False)
        final_df_sig.to_csv(f'{base_dir}/data/ortholog_mapping/{file}_orthologs_sig.csv', index=False)

    # Create one summary table with all the orthologs

    # Initialize a dictionary to store dataframes
    comparison_dfs = {}

    # Read each comparison file
    for file in pw_interest:
        df = pd.read_csv(os.path.join(base_dir, f'data/ortholog_mapping/{file}_orthologs_sig.csv'))
        comp_name = file
        comparison_dfs[comp_name] = df
        print(f"Loaded {file} with {len(df)} rows")

    # Create the summary table
    summary_data = []

    # Process each comparison file and add its data to the summary table
    for comp_name, df in comparison_dfs.items():
        for _, row in df.iterrows():
            # Create a dictionary for each ortholog entry
            entry = {
                'Zebrafish Gene ID': row['input'],
                'Ortholog Name': row['ortholog_name'],
                f"{comp_name} Expression": row['Expression'],
                f"{comp_name} L2FC": row['log2FoldChange']
            }
            
            # Check if this entry already exists in summary_data
            existing_entry = None
            for i, existing in enumerate(summary_data):
                if (existing['Zebrafish Gene ID'] == entry['Zebrafish Gene ID'] and 
                    existing['Ortholog Name'] == entry['Ortholog Name']):
                    existing_entry = i
                    break
            
            # If entry exists, update it with this comparison's data
            if existing_entry is not None:
                summary_data[existing_entry].update({
                    f"{comp_name} Expression": entry[f"{comp_name} Expression"],
                    f"{comp_name} L2FC": entry[f"{comp_name} L2FC"]
                })
            # Otherwise add as new entry
            else:
                summary_data.append(entry)

    # Create the summary dataframe
    summary_df = pd.DataFrame(summary_data)

    # Reorder columns
    column_order = ['Zebrafish Gene ID', 'Ortholog Name']
    for comp in comparison_dfs.keys():
        column_order.extend([f"{comp} L2FC", f"{comp} Expression"])

    # Select only columns that exist in the dataframe
    valid_columns = [col for col in column_order if col in summary_df.columns]
    summary_df = summary_df[valid_columns]

    # Sort by average log2 fold change across comparisons
    summary_df['Average L2FC'] = summary_df[[col for col in summary_df.columns if 'L2FC' in col]].mean(axis=1)
    summary_df = summary_df.sort_values('Average L2FC', ascending=False)

    # Save the summary table
    summary_df.to_csv(f'{base_dir}/results/DEGs/ortholog_summary_table.csv', index=False)
    print(f"Created summary table with {len(summary_df)} rows")

# %% [markdown]
# ## Heatmap of Top DEGs

# %%

'''
(create_heatmap_topDEGs) reads the DESeq2 files, extracts the DEGs to a list, extracts the normalized expression data for these genes, and plots a heatmap across the samples
'''

def create_heatmap_topDEGs (pw_interest, norm_file_name, base_dir):

    all_degs = pd.DataFrame()

    # Process each comparison file 
    for file in pw_interest:
        df = pd.read_csv(os.path.join(base_dir, f'data/DESeq2/{file}_processed.csv'))
        # Filter to keep only DEGs (Upregulated or Downregulated)
        df_degs = df[df['Expression'].isin(['Upregulated', 'Downregulated'])]
        # Add to collection
        all_degs = pd.concat([all_degs, df_degs])
        
    # Sort by absolute log2FoldChange and get top genes across all comparisons
    top_degs = all_degs.sort_values('log2FoldChange', key=abs, ascending=False)
    top_gene_ids = top_degs['ensembl_gene_id'].unique().tolist()[:30]  # Get top 30 unique genes

    # Load normalized counts and extract data for top genes
    norm_data = pd.read_csv(norm_file_name, index_col=0)
    expression_data = norm_data.loc[norm_data.index.isin(top_gene_ids)]

    # Create the heatmap
    plt.figure(figsize=(12, 10))
    sns.clustermap(
        expression_data,
        cmap="viridis",
        z_score=0,
        col_cluster=True,
        row_cluster=True,
        dendrogram_ratio=(.1, .2),
        cbar_pos=(0.02, .32, .03, .2)
    )
    
    plt.savefig(f'{base_dir}/plots/DEG_heatmap/TopDEG_heatmap.png')
    plt.close()




# %% [markdown]
# ## Enrichment Analysis

# Function to wrap GO term labels
def wrap_terms(terms, max_chars=35):
    wrapped_terms = []
    for term in terms:
        if len(term) > max_chars:
            # Use textwrap for clean line breaks
            wrapped_term = '\n'.join(textwrap.wrap(term, width=max_chars))
        else:
            wrapped_term = term
        wrapped_terms.append(wrapped_term)
    return wrapped_terms

# %% [markdown]
# Perform GO enrichment analysis on the significant genes for each comparison of interest

# %%
'''
(perform_GO_enrichment) reads the DESeq2 files with significant orthologs, extracts them to a list, performs GO enrichment analysis, and saves the results as bar plots.
'''

def perform_GO_enrichment (pw_interest, base_dir, enrich_sig_cutoff):

    # Initialize a dictionary to store dataframes
    comparison_dfs = {}

    # Process each comparison file
    for file in pw_interest:
        df = pd.read_csv(os.path.join(base_dir, f'data/ortholog_mapping/{file}_orthologs_sig.csv'))
        comp_name = file
        comparison_dfs[comp_name] = df
        print(f"Loaded {file} with {len(df)} rows")

        # Convert Ensembl IDs to gene symbols
        mg = mygene.MyGeneInfo()
        ensembl_ids = df['ortholog_ensg'].tolist()

        # I need to supress a terminal output here for user experience
        original_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')
        gene_info = mg.getgenes(ensembl_ids, fields='symbol', species='human')
        sys.stderr = original_stderr

        # Extract symbols and convert to uppercase
        gene_list = [gene.get('symbol', '').upper() for gene in gene_info if gene.get('symbol')]
        
        # Save gene list to file
        gene_list_file = f"{comp_name}_significant_genes.txt"
        with open(f"{base_dir}/data/GO_and_KEGG_enrichments/{gene_list_file}", "w") as f:
            f.write("\n".join(gene_list))
        
        print(f"Found {len(gene_list)} significant genes")
        
        # Perform GO enrichment analysis
        if len(gene_list) > 0:
            
            # Run enrichment analysis
            enrichr_results = gp.enrichr(
            gene_list=gene_list,
            gene_sets='GO_Biological_Process_2025',
            organism='human',
            cutoff=enrich_sig_cutoff
            )
            # Get results as DataFrames
            go_bp_results = enrichr_results.results

            # Sort the results by adjusted p-value (ascending)
            go_bp_results = go_bp_results.sort_values('Adjusted P-value')

            # Extract gene count from the Overlap column (first number before the slash)
            go_bp_results['Gene Count'] = go_bp_results['Overlap'].str.split('/').str[0].astype(int)
            
            # Save results to CSV
            go_bp_results.to_csv(f"{base_dir}/results/GO_enrichment/{comp_name}_GO_BP_enrichment.csv", index=False)

            # Visualize the results
            plt.figure(figsize=(14, 10))

            # Select top 10 most significant terms
            top10_terms = go_bp_results[:10].copy()

            # Wrap the GO term labels
            top10_terms['Wrapped_Term'] = wrap_terms(top10_terms['Term'].tolist())

            # Create a colormap based on p-values of only the top 10 terms
            norm = plt.Normalize(vmin=top10_terms['Adjusted P-value'].min(), vmax=top10_terms['Adjusted P-value'].max())

            # Create a monochromatic green colormap
            from matplotlib.colors import LinearSegmentedColormap
            green_cmap = LinearSegmentedColormap.from_list('green_cmap', ['#c8e6c9', '#1b5e20'])

            # Plot bars with colors based on p-values
            plt.barh(top10_terms['Wrapped_Term'], top10_terms['Gene Count'],
                color=green_cmap(norm(top10_terms['Adjusted P-value'])))
            plt.yticks(fontsize=14)
            plt.xticks(fontsize=14)
            plt.xlabel('Gene Count', fontsize=16, fontweight='bold')
            plt.ylabel('GO Term Description', fontsize=16, fontweight='bold')
            plt.title(f'BP GO Terms: {comp_name}', fontsize=18, fontweight='bold')
            plt.gca().invert_yaxis()
            
            # Adjust layout to accommodate wrapped labels
            plt.subplots_adjust(left=0.45)  # Make more room for longer labels

            # Add color bar for p-values
            sm = plt.cm.ScalarMappable(cmap=green_cmap, norm=norm)
            sm.set_array([])
            cax = plt.colorbar(sm, ax=plt.gca(), label='P-adjusted')
            cax.ax.tick_params(labelsize=14)  # Colorbar tick labels
            cax.set_label('P-adjusted', fontsize=16)  # Colorbar title
            
            plt.tight_layout()
            plt.savefig(f'{base_dir}/plots/GO_enrichment/{comp_name}_GO_BP_enrichment_plot.png')
            plt.close()

            print(f"GO enrichment analysis completed. Results saved for {comp_name}")
        else:
            print(f"No significant GO terms found for {comp_name}")
    else:
        print(f"No significant genes found for {comp_name}")

    print("All comparisons processed.")


# %% [markdown]
# Perform KEGG enrichment analysis on the significant genes for each comparison of interest

# %%
'''
(perform_KEGG_enrichment) reads the DESeq2 files with significant orthologs, extracts them to a list, performs KEGG enrichment analysis, and saves the results as bar plots.
'''

def perform_KEGG_enrichment (pw_interest, base_dir, enrich_sig_cutoff):

    # Initialize a dictionary to store dataframes
    comparison_dfs = {}

    # Process each comparison file
    for file in pw_interest:
        df = pd.read_csv(f'{base_dir}/data/ortholog_mapping/{file}_orthologs_sig.csv')
        comp_name = file
        comparison_dfs[comp_name] = df
        print(f"Loaded {file} with {len(df)} rows")

        # Convert Ensembl IDs to gene symbols
        mg = mygene.MyGeneInfo()
        ensembl_ids = df['ortholog_ensg'].tolist()

        # I need to supress a terminal output here for user experience
        original_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')
        gene_info = mg.getgenes(ensembl_ids, fields='symbol', species='human')
        sys.stderr = original_stderr
        
        # Extract symbols and convert to uppercase
        gene_list = [gene.get('symbol', '').upper() for gene in gene_info if gene.get('symbol')]
        
        # Save gene list to file
        gene_list_file = f"{comp_name}_significant_genes.txt"
        with open(f"{base_dir}/data/GO_and_KEGG_enrichments/{gene_list_file}", "w") as f:
            f.write("\n".join(gene_list))
        
        print(f"Found {len(gene_list)} significant genes")
        
        # Perform KEGG enrichment analysis
        if len(gene_list) > 0:
        
            # Run enrichment analysis
            enrichr_results = gp.enrichr(
            gene_list=gene_list,
            gene_sets='KEGG_2021_Human',
            organism='human',
            cutoff=enrich_sig_cutoff
            )
            # Get results as DataFrames
            go_bp_results = enrichr_results.results

            # Sort the results by adjusted p-value (ascending)
            go_bp_results = go_bp_results.sort_values('Adjusted P-value')

            # Extract gene count from the Overlap column (first number before the slash)
            go_bp_results['Gene Count'] = go_bp_results['Overlap'].str.split('/').str[0].astype(int)
            
            # Save results to CSV
            go_bp_results.to_csv(f"{base_dir}/results/KEGG_enrichment/{comp_name}_KEGG_enrichment.csv", index=False)

            # Visualize the results
            plt.figure(figsize=(14, 10))
            
            # Select top 10 most significant terms
            top10_terms = go_bp_results[:10]
            
            # Wrap the KEGG term labels
            top10_terms['Wrapped_Term'] = wrap_terms(top10_terms['Term'].tolist())
            
            # Create a colormap based on p-values of only the top 10 terms
            norm = plt.Normalize(vmin=top10_terms['Adjusted P-value'].min(), vmax=top10_terms['Adjusted P-value'].max())

            # Create a monochromatic green colormap
            from matplotlib.colors import LinearSegmentedColormap
            green_cmap = LinearSegmentedColormap.from_list('green_cmap', ['#c8e6c9', '#1b5e20'])

            # Plot bars with colors based on p-values
            plt.barh(top10_terms['Wrapped_Term'], top10_terms['Gene Count'],
                color=green_cmap(norm(top10_terms['Adjusted P-value'])))
            plt.yticks(fontsize=14)
            plt.xticks(fontsize=14)
            plt.xlabel('Gene Count', fontsize=16, fontweight='bold')
            plt.ylabel('KEGG Term Description', fontsize=16, fontweight='bold')
            plt.title(f'KEGG Terms: {comp_name}', fontsize=18, fontweight='bold')
            plt.gca().invert_yaxis()

            # Adjust layout to accommodate wrapped labels
            plt.subplots_adjust(left=0.45)  # Make more room for longer labels

            # Add color bar for p-values with larger font sizes
            sm = plt.cm.ScalarMappable(cmap=green_cmap, norm=norm)
            sm.set_array([])
            cax = plt.colorbar(sm, ax=plt.gca())
            cax.ax.tick_params(labelsize=14)  # Colorbar tick labels
            cax.set_label('P-adjusted', fontsize=16, fontweight='bold')  # Colorbar title
            
            plt.tight_layout()
            plt.savefig(f'{base_dir}/plots/KEGG_enrichment/{comp_name}_KEGG_enrichment_plot.png')
            plt.close()

            print(f"KEGG enrichment analysis completed. Results saved for {comp_name}")
        else:
            print(f"No significant KEGG terms found for {comp_name}")


