# A Streamlit-based interface that simplifies pipeline usage by allowing users to upload files and set parameters in a browser.
# Author: Ishaan Banwait
# Last Updated: 2025-04-20

import streamlit as st
import subprocess
import sys
import os
import time
from pathlib import Path
import glob
import shlex

# Import your AI assistant functionality
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from ai_assistant import get_pipeline_parameters

# Function to run a command and display its output in real-time
def run_and_display_stdout(cmd_string):
    print(f"Running command: {cmd_string}")
    
    # Use shlex.split to properly handle quoted arguments
    parts = shlex.split(cmd_string)
    
    # Replace 'python' with sys.executable
    if parts[0] == 'python':
        parts[0] = sys.executable
        
    st.session_state.command_output = []
    st.session_state.command_output.append(f"Running command: {' '.join(parts)}")
    
    # Create a placeholder for live output
    output_placeholder = st.empty()
    
    process = subprocess.Popen(
        parts,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env={**os.environ, 'PYTHONUNBUFFERED': '1'}
    )
    
    # Display output in real-time
    for line in iter(lambda: process.stdout.readline(), ''):
        if line:
            st.session_state.command_output.append(line.strip())
            output_placeholder.code('\n'.join(st.session_state.command_output))
            time.sleep(0.05)  # Small delay to allow UI to update
    
    # Check process return code
    return_code = process.wait()
    if return_code == 0:
        st.session_state.command_output.append("Command completed successfully!")
    else:
        st.session_state.command_output.append(f"Command failed with return code {return_code}.")
    
    output_placeholder.code('\n'.join(st.session_state.command_output))

st.set_page_config(page_title="RNA-seq Pipeline Assistant", layout="wide")

# Initialize session state for conversation context
if 'context' not in st.session_state:
    st.session_state.context = {'messages': []}

if 'setup_command' not in st.session_state:
    st.session_state.setup_command = ""

if 'runall_command' not in st.session_state:
    st.session_state.runall_command = ""

if 'command_output' not in st.session_state:
    st.session_state.command_output = []

st.title("DEG Pipeline Assistant")

st.markdown("""
## About the Pipeline

This pipeline takes normalized RNA-seq data and outputs differentially expressed genes and pathways. It performs several steps:

- Creates metadata from normalized counts
- Generates quality control plots (boxplots, correlation heatmaps, PCA)
- Processes DESeq2 files to identify differentially expressed genes
- Creates visualization plots (MA plots, volcano plots)
- Maps genes to human orthologs
- Performs GO and KEGG enrichment analysis
""")

# Create two columns for the form
col1, col2 = st.columns(2)

with col1:
    st.header("Setup Parameters")
    
    uploaded_norm_file = st.file_uploader("Upload Normalized Data File", type=["csv"])
    norm_file_path = None
    if uploaded_norm_file is not None:
        norm_file_path = f"/tmp/{uploaded_norm_file.name}"
        with open(norm_file_path, "wb") as f:
            f.write(uploaded_norm_file.getbuffer())
    
    # Create a dedicated directory for pairwise files
    pw_data_dir = "/tmp/pairwise_files/"
    os.makedirs(pw_data_dir, exist_ok=True)
    uploaded_pw_files = st.file_uploader("Upload Pairwise Data Files", type=["csv"], accept_multiple_files=True)
    pw_data_paths = []
    if uploaded_pw_files:
        for uploaded_file in uploaded_pw_files:
            file_path = f"{pw_data_dir}{uploaded_file.name}"
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            pw_data_paths.append(file_path)
            
    conditions = st.text_input("Conditions", "e.g., Control, Treatment1, Treatment2")
    num_replicates = st.number_input("Number of Replicates", min_value=1, value=3)

with col2:
    st.header("Analysis Parameters")
    
    model_organism = st.text_input("Model Organism", "e.g., drerio")
    pw_interest = st.text_input("Pairwise Comparisons of Interest (comma-separated)", "e.g., ConditionA_vs_ConditionB")
    log2fc_threshold = st.number_input("Log2 Fold Change Threshold", value=1.0, step=0.1)
    padj_threshold = st.number_input("Adjusted p-value Threshold", value=0.05, step=0.01)
    enrich_sig_cutoff = st.number_input("Enrichment Significance Cutoff", value=0.05, step=0.01)
    skip_steps = st.multiselect(
        "Steps to Skip",
        options=["boxplots", "sample correlation heatmap", "PCA", "ma plots", "volcano plots", "DEG heatmap", "GO enrichment", "KEGG enrichment"],
        default=[],
        help="Select the steps you want to skip in the pipeline."
    )

# Single button to generate pipeline commands
if st.button("Generate Pipeline Commands"):
    # Collect natural language parameters
    param_inputs = {
        "conditions": conditions,
        "model_organism": model_organism,
        "pw_interest": pw_interest,
    }

    # Initialize param_final with direct values that don't need processing
    param_final = {
        'norm_file': norm_file_path,
        'pw_data': "/tmp/pairwise_files/",
        'base_dir': "/tmp/output/",
        'num_replicates': num_replicates,
        'log2fc_threshold': log2fc_threshold,
        'padj_threshold': padj_threshold,
        'enrich_sig_cutoff': enrich_sig_cutoff,
        'skip_steps': skip_steps
    }

    if any(not value for key, value in param_inputs.items()) or any(not value for key, value in param_final.items() if key != "skip_steps"):
        st.error("Please fill in all required parameters.")
    else:
        # Make a single call to the AI assistant to generate commands
        with st.spinner("Generating commands..."):
            # This would call your existing get_pipeline_commands function
            command_suggestion = get_pipeline_parameters(param_inputs, st.session_state.context)

            # Update param_final with AI-generated values
            for param, value in command_suggestion.items():
                param_final[param] = value
            
            # Generate base command
            base_command = (
                f"python pipeline_CLI.py "
                f"--norm-file {param_final['norm_file']} "
                f"--pw-data {param_final['pw_data']} "
                f"--base-dir {param_final['base_dir']} "
                f"--num-replicates {param_final['num_replicates']} "
            )
            for condition in param_final['conditions']:
                base_command += f"--conditions '{condition}' "

            # Generate setup command
            setup_command = (
                f"{base_command}setup "
            )

            # Generate run-all command
            runall_command = (
                f"{base_command}run-all "
                f"--model-organism '{param_final['model_organism']}' "
            )
            # Add each pairwise interest as a separate --pw-interest argument
            for pw in param_final['pw_interest']:
                runall_command += f"--pw-interest '{pw}' "

            if param_final['log2fc_threshold'] != "1":
                runall_command += f" --log2fc-threshold {param_final['log2fc_threshold']}"
            
            if param_final['padj_threshold'] != "0.05":
                runall_command += f" --padj-threshold {param_final['padj_threshold']}"
            
            if param_final['enrich_sig_cutoff'] != "0.05":
                runall_command += f" --enrich-sig-cutoff {param_final['enrich_sig_cutoff']}"
            
            for step in param_final['skip_steps']:
                if step:
                    runall_command += f" --skip-{step}"
            
            
            st.session_state.setup_command = setup_command
            st.session_state.runall_command = runall_command
            st.success("Commands generated successfully!")

# Display generated commands
if st.session_state.setup_command:
    st.header("Generated Commands")
    
    st.subheader("Setup Command")
    st.code(st.session_state.setup_command, language="bash")

    if st.button("Run Setup Command"):
        run_and_display_stdout(st.session_state.setup_command)

    st.subheader("Analysis Command")
    st.code(st.session_state.runall_command, language="bash")
    
    if st.button("Run Analysis Command"):
        run_and_display_stdout(st.session_state.runall_command)
    
    output_dir = "/tmp/output"
    files = glob.glob(os.path.join(output_dir, '**'), recursive=True)
    files = [f for f in files if os.path.isfile(f)]
    if files:  # Only show if there are files to download
        st.header("Pipeline Output Files", divider="gray")
        st.markdown("Download the generated analysis files below")
        
        # Create categories for different file types
        file_categories = {
            "Plots": [f for f in files if f.endswith(('.png', '.jpg', '.pdf', '.svg'))],
            "Data Tables": [f for f in files if f.endswith(('.csv', '.tsv', '.xlsx', '.txt'))],
            "Reports": [f for f in files if f.endswith(('.html', '.md', '.json'))],
            "Other Files": [f for f in files if not any(f.endswith(ext) for ext in 
                            ('.png', '.jpg', '.pdf', '.svg', '.csv', '.tsv', '.xlsx', '.txt', '.html', '.md', '.json'))]
        }
        
        # Display files by category in a tabbed interface
        tabs = st.tabs([f"{category} ({len(files)})" for category, files in file_categories.items() if files])
        
        tab_index = 0
        for category, category_files in file_categories.items():
            if not category_files:
                continue
                
            with tabs[tab_index]:
                # Create a grid layout for download buttons
                cols_per_row = 3
                for i in range(0, len(category_files), cols_per_row):
                    cols = st.columns(cols_per_row)
                    for j in range(cols_per_row):
                        if i + j < len(category_files):
                            filepath = category_files[i + j]
                            filename = os.path.relpath(filepath, output_dir)
                            with cols[j]:
                                with open(filepath, "rb") as f:
                                    st.download_button(
                                        label=os.path.basename(filename),
                                        data=f,
                                        file_name=os.path.basename(filename),
                                        key=f"download_{i}_{j}_{category}",
                                        use_container_width=True
                                    )
                                # Add a small caption with file info
                                file_size = os.path.getsize(filepath) / 1024  # KB
                                st.caption(f"{file_size:.1f} KB | {os.path.splitext(filename)[1][1:].upper()}")
            
            tab_index += 1
        
        # Add a button to download all files as a zip
        st.markdown("### Download All Files")
        
        if st.button("Prepare All Files as ZIP Archive", use_container_width=True):
            # Create a temporary zip file
            import zipfile
            import io
            
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for file in files:
                    filename = os.path.relpath(file, output_dir)
                    zip_file.write(file, arcname=filename)
            
            zip_buffer.seek(0)
            st.download_button(
                label="Download All Files (ZIP)",
                data=zip_buffer,
                file_name="pipeline_output.zip",
                mime="application/zip",
                use_container_width=True
            )

