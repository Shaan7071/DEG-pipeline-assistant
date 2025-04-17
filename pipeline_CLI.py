import os
import click
import sys
from pathlib import Path


# Import functions from your pipeline module
from DEGpipeline_v6 import (
    create_metadata,
    create_boxplots,
    create_correlation_heatmap,
    create_pca,
    process_deseq2_files,
    create_ma_plots,
    create_volcano_plots,
    map_to_human_orthologs,
    create_heatmap_topDEGs,
    perform_GO_enrichment,
    perform_KEGG_enrichment
)

@click.group()
@click.option('--norm-file', type=click.Path(exists=True), required=True, help='Normalized data file')
@click.option('--pw-data', type=click.Path(exists=True), help='Directory containing pairwise comparison files')
@click.option('--base-dir', type=click.Path(), required=True, help='Base directory for outputs')
@click.option('--conditions', '-c', multiple=True, required=True, help='Experimental conditions')
@click.option('--num-replicates', '-n', type=int, required=True, help='Number of replicates per condition')
@click.pass_context
def cli(ctx, norm_file, pw_data, base_dir, conditions, num_replicates):
    """DEG Pipeline CLI - Process RNA-seq data to find differentially expressed genes and pathways."""
    ctx.ensure_object(dict)
    # Create objects from the command line arguments
    ctx.obj['NORM_FILE'] = norm_file
    ctx.obj['PW_DATA'] = pw_data
    ctx.obj['BASE_DIR'] = base_dir
    ctx.obj['NUM_REPLICATES'] = num_replicates
    ctx.obj['CONDITIONS'] = list(conditions)

@cli.command()
@click.pass_context
def setup(ctx):
    """Initialize the pipeline and create directory structure."""
    base_dir=ctx.obj['BASE_DIR']

    # Create directories
    directories = [f"{base_dir}/data", f"{base_dir}/results", f"{base_dir}/plots"]
    
    # Create subdirectories within 'data'
    subdirectories = [f"{base_dir}/data/DESeq2", f"{base_dir}/data/ortholog_mapping", 
                      f"{base_dir}/data/GO_and_KEGG_enrichments"]
    directories.extend(subdirectories)
    
    # Create subdirectories within 'plots'
    subdirectories = [f"{base_dir}/plots/boxplots", f"{base_dir}/plots/correlation_heatmap", 
                      f"{base_dir}/plots/PCA", f"{base_dir}/plots/MA_plots", 
                      f"{base_dir}/plots/volcano_plots", f"{base_dir}/plots/DEG_heatmap", 
                      f"{base_dir}/plots/GO_enrichment", f"{base_dir}/plots/KEGG_enrichment"]
    directories.extend(subdirectories)
    
    # Create subdirectories within 'results'
    subdirectories = [f"{base_dir}/results/DEGs", f"{base_dir}/results/GO_enrichment", 
                      f"{base_dir}/results/KEGG_enrichment"]
    directories.extend(subdirectories)
    
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)

    click.echo("Directories created successfully")

    
@cli.command()
@click.pass_context
def metadata(ctx):
    """Create metadata file from normalized RNA-seq data."""
    # Pass context values to your function
    create_metadata(
        norm_file_name=ctx.obj['NORM_FILE'],
        conditions=ctx.obj.get('CONDITIONS', []),
        num_replicates=ctx.obj['NUM_REPLICATES'],
        base_dir=ctx.obj['BASE_DIR']
    )
    click.echo("Metadata created successfully")

@cli.command()
@click.pass_context
def boxplots(ctx):
    """Create boxplots of samples for quality control."""
    create_boxplots(
        norm_file_name=ctx.obj['NORM_FILE'],
        conditions=ctx.obj.get('CONDITIONS', []),
        num_replicates=ctx.obj['NUM_REPLICATES'],
        base_dir=ctx.obj['BASE_DIR']
    )
    
    click.echo("Boxplots created")

@cli.command()
@click.pass_context
def correlation_heatmap(ctx):
    """Create correlation heatmap of samples for quality control."""
    create_correlation_heatmap(
        norm_file_name=ctx.obj['NORM_FILE'],
        base_dir=ctx.obj['BASE_DIR']
    )
    
    click.echo("Sample correlation heatmap created")

@cli.command()
@click.pass_context
def pca(ctx):
    """Perform principal component analysis."""
    create_pca(
        norm_file_name=ctx.obj['NORM_FILE'],
        base_dir=ctx.obj['BASE_DIR']
    )

    click.echo("PCA plot created")

@cli.command()
@click.pass_context
def process_deseq2(ctx):
    """Process DESeq2 result files and categorize genes."""
    
    process_deseq2_files(
        input_dir_pw_data=ctx.obj['PW_DATA'],
        base_dir=ctx.obj['BASE_DIR'],
        log2fc_threshold=ctx.obj['LOG2FC_THRESHOLD'],
        padj_threshold=ctx.obj['PADJ_THRESHOLD']
    )
    click.echo("Pairwise comparison files processed successfully")

@cli.command()
@click.pass_context
def ma_plots(ctx):
    """Generate MA plots from processed DESeq2 files."""
    create_ma_plots(
        input_dir_pw_data=ctx.obj['PW_DATA'],
        base_dir=ctx.obj['BASE_DIR']
    )
    click.echo("MA plots created successfully")

@cli.command()
@click.pass_context
def volcano_plots(ctx):
    """Generate volcano plots from processed DESeq2 files."""
    create_volcano_plots(
        input_dir_pw_data=ctx.obj['PW_DATA'],
        base_dir=ctx.obj['BASE_DIR'],
        log2fc_threshold=ctx.obj['LOG2FC_THRESHOLD'],
        padj_threshold=ctx.obj['PADJ_THRESHOLD']
    )
    
    click.echo("Volcano plots created successfully")

@cli.command()
@click.pass_context
def ortholog_mapping(ctx):
    """Map model organism genes to human orthologs."""
    
    map_to_human_orthologs(
        pw_interest=ctx.obj.get('PW_INTEREST', []),
        model_organism=ctx.obj['MODEL_ORGANISM'],
        base_dir=ctx.obj['BASE_DIR']
    )
    click.echo("Ortholog mapping completed successfully")

@cli.command()
@click.pass_context
def deg_heatmap(ctx):
    """Create heatmap of top differentially expressed genes."""
    create_heatmap_topDEGs(
        pw_interest=ctx.obj.get('PW_INTEREST', []),
        norm_file_name=ctx.obj['NORM_FILE'],
        base_dir=ctx.obj['BASE_DIR']
    )
    click.echo("DEG heatmap created successfully")

@cli.command()
@click.pass_context
def go_enrichment(ctx):
    """Perform GO enrichment analysis on significant genes."""
    perform_GO_enrichment(
        pw_interest=ctx.obj.get('PW_INTEREST', []),
        base_dir=ctx.obj['BASE_DIR'],
        enrich_sig_cutoff=ctx.obj['ENRICH_SIG_CUTOFF']
    )
    click.echo("GO enrichment analysis completed successfully")

@cli.command()
@click.pass_context
def kegg_enrichment(ctx):
    """Perform KEGG enrichment analysis on significant genes."""
    perform_KEGG_enrichment(
        pw_interest=ctx.obj.get('PW_INTEREST', []),
        base_dir=ctx.obj['BASE_DIR'],
        enrich_sig_cutoff=ctx.obj['ENRICH_SIG_CUTOFF']
    )
    click.echo("KEGG enrichment analysis completed successfully")

@cli.command()
# @click.option('--run-metadata/--skip-metadata', default=True, help='Run metadata creation step')
@click.option('--run-boxplots/--skip-boxplots', default=True, help='Create boxplots of samples for quality control')
@click.option('--run-correlation-heatmap/--skip-correlation-heatmap', default=True, help='Create correlation heatmap for the samples')
@click.option('--run-pca/--skip-pca', default=True, help='Perform principal component analysis on the data')
@click.option('--run-ma-plots/--skip-ma-plots', default=True, help='Generate MA plots')
@click.option('--run-volcano-plots/--skip-volcano-plots', default=True, help='Generate volcano plots')
@click.option('--run-heatmap/--skip-heatmap', default=True, help='Create DEG heatmap')
@click.option('--run-go/--skip-go', default=True, help='Perform GO enrichment analysis')
@click.option('--run-kegg/--skip-kegg', default=True, help='Perform KEGG enrichment analysis')


@click.option('--log2fc-threshold', type=float, default=1, help='Log2 fold change threshold for DESeq2')
@click.option('--padj-threshold', type=float, default=0.05, help='Adjusted p-value threshold for DESeq2')
@click.option('--enrich-sig-cutoff', type=float, default=0.05, help='Significance cutoff for enrichment analysis')
@click.option('--model-organism', type=str, required=True, help='Model organism for ortholog mapping')
@click.option('--pw-interest', multiple=True, required=True, help='Comparisons of interest )')

@click.pass_context
def run_all(ctx, log2fc_threshold, padj_threshold, enrich_sig_cutoff, model_organism, pw_interest,
            run_boxplots, run_correlation_heatmap, run_pca, run_ma_plots, run_volcano_plots, run_heatmap, run_go, run_kegg):
    """Run selected components of the pipeline in sequence."""

    ctx.obj['MODEL_ORGANISM'] = model_organism
    ctx.obj['PW_INTEREST'] = list(pw_interest)


    ctx.obj['LOG2FC_THRESHOLD'] = log2fc_threshold
    ctx.obj['PADJ_THRESHOLD'] = padj_threshold
    ctx.obj['ENRICH_SIG_CUTOFF'] = enrich_sig_cutoff


    click.echo("Running DEG pipeline with selected components...")

    # Run mandatory metadata step
    ctx.invoke(metadata)
    
    # Run optional steps between metadata and DESeq2 processing based on flags
    if run_boxplots:
        ctx.invoke(boxplots)

    if run_correlation_heatmap:
        ctx.invoke(correlation_heatmap)

    if run_pca:
        ctx.invoke(pca)
    
    # Run mandatory DESeq2 processing 
    ctx.invoke(process_deseq2)
    
    # Run optional plots from DESeq2 data 
    if run_ma_plots:
        ctx.invoke(ma_plots)
    
    if run_volcano_plots:
        ctx.invoke(volcano_plots)
    
    # Run mandatory ortholog mapping 
    ctx.invoke(ortholog_mapping)
    
    # Run optional plots from otholog data
    if run_heatmap:
        ctx.invoke(deg_heatmap)
    
    if run_go:
        ctx.invoke(go_enrichment)
    
    if run_kegg:
        ctx.invoke(kegg_enrichment)
    
    click.echo("Pipeline completed successfully!")

# entry point
if __name__ == '__main__':
    cli()
