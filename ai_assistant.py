# Author: Ishaan Banwait 
# Last Updated: 2025-04-16

import os
import openai
from openai import OpenAI
import subprocess
import json

# Initialize the OpenAI client
client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))

def get_pipeline_parameters(user_input_dict, pipeline_context=None):
    """
    Extract pipeline parameters from a dictionary of user inputs in natural language using OpenAI API
    
    Args:
        user_input_dict (dict): The user's natural language request for the command where each key is a param and each value is the user input for that param
        pipeline_context (dict, optional): Previous conversation context
        
    Returns:
        dict: Extracted parameters for the pipeline
    """
    # Create system message for parameter extraction
    system_message = """
    You are an assistant that helps researchers extract parameters for a DEG (Differentially Expressed Genes) pipeline for RNA-seq analysis.
    
    Your job is to extract these parameters from natural language descriptions:
    
    Setup parameters:
    - norm_file: Path to the normalized data file (string)
    - pw_data: Directory containing pairwise data files (string)
    - base_dir: Directory to store results (string)
    - conditions: List of experimental conditions being tested (list of strings)
    - num_replicates: Number of replicates per condition (int)
    
    Analysis parameters:
    - model_organism: Model organism code. You may need to convert into the g:Profiler organism id if it isn't already. The dict is here if needed: https://biit.cs.ut.ee/gprofiler/page/organism-list (output should be a string e.g., "drerio" for zebrafish.)
    - pw_interest: Pairwise comparisons of interest (list of strings)
    - log2fc_threshold: Log2 fold change threshold (int)
    - padj_threshold: Adjusted p-value threshold (int)
    - enrich_sig_cutoff: Enrichment significance cutoff (int)
    - skip_steps: List of steps to skip (list of boxplots, correlation-heatmap, pca, ma-plots, volcano-plots, heatmap, go, kegg)
    
    FORMAT YOUR RESPONSE AS A VALID JSON OBJECT with parameter names as keys.
    Only include parameters that are explicitly mentioned or clearly implied in the text.
    For file paths and directory paths, preserve the exact formatting including quotes and backslashes.
    Ensure that skip_steps match the given options above (eg. skip--go, not skip-GO
    """
    # Format the user input dictionary into a structured text
    formatted_input = ""
    for param_name, param_value in user_input_dict.items():
        formatted_input += f"Parameter: {param_name}\nUser Input: {param_value}\n\n"
    
    # Initialize messages list
    messages = [
        {"role": "system", "content": system_message},
    ]
    
    # Add previous conversation context if available
    if pipeline_context and 'messages' in pipeline_context:
        messages.extend(pipeline_context['messages'])
    
    # Add the current user input
    messages.append({"role": "user", "content": formatted_input})
    
    # Call the OpenAI API
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=messages,
            response_format={"type": "json_object"},  # Request JSON response
            max_tokens=500,
            temperature=0.2
        )
        
        # Extract the parameters as JSON
        parameters_json = response.choices[0].message.content.strip()
        parameters = json.loads(parameters_json)
        
        # Update conversation context
        if pipeline_context is not None:
            pipeline_context['messages'] = messages + [{"role": "assistant", "content": parameters_json}]
        
        return parameters
    except json.JSONDecodeError:
        return {"error": "Failed to parse parameters from response"}
    except Exception as e:
        return {"error": f"Error extracting parameters: {str(e)}"}


def interactive_command_builder():
    """
    Interactive session to build the pipeline commands through conversation.
    """
    print("\n===== DEG Pipeline AI Assistant =====")
    print("I'll help you set up your RNA-seq analysis. Describe your experiment, and I'll create the appropriate command.")
    print("Type 'exit' to quit at any time.\n")

    # Initialize conversation context
    context = {'messages': []}
    
    # Track collected information for setup
    setup_info = {
        'norm_file': None,
        'pw_data': None,
        'base_dir': None,
        'conditions': None,
        'num_replicates': None
    }

    setup_inputs = {
        'norm_file': None,
        'pw_data': None,
        'base_dir': None,
        'conditions': None,
        'num_replicates': None
    }
    
    # Track collected information for run-all
    runall_info = {
        'model_organism': None,
        'pw_interest': None,
        'log2fc_threshold': None,
        'padj_threshold': None,
        'enrich_sig_cutoff': None,
        'skip_steps': None
    }

    runall_inputs = {
        'model_organism': None,
        'pw_interest': None,
        'log2fc_threshold': None,
        'padj_threshold': None,
        'enrich_sig_cutoff': None,
        'skip_steps': None
    }
    
    # Stage tracking
    current_stage = "setup"
    setup_command = None
    
    # Function to check if all required parameters are collected
    def all_parameters_collected(info_dict, required_keys):
        return all(info_dict[key] is not None for key in required_keys)
    
    # Start interactive session
    print("Assistant: First, we'll construct an output directory for all the results of the analysis. Then we'll proceed with the actual rna-seq analysis. Let's start with the setup parameters.")
    
    while True:
        
        # Check for exit command
        
        user_input = ""
        if user_input.lower() == 'exit':
            print("Assistant: Goodbye!")
            return None
        

        # Setup command stage
        if current_stage == "setup":
            # print the values of setup_inputs
            print(f"\nAssistant: Current setup parameters: {setup_inputs}")

            # Collect setup parameters
            if not setup_inputs['norm_file']:
                print("Assistant: What is the path to your normalized data file?")
                user_input = input("You: ")
                setup_inputs['norm_file'] = user_input
                continue
            
            if not setup_inputs['pw_data']:
                print("Assistant: What is the directory containing your pairwise data files?")
                user_input = input("You: ")
                setup_inputs['pw_data'] = user_input
                continue
            
            if not setup_inputs['base_dir']:
                print("Assistant: What is the base directory you would like the outputs to be saved in?")
                user_input = input("You: ")
                setup_inputs['base_dir'] = user_input
                continue
            
            if not setup_inputs['conditions']:
                print("Assistant: What are your experimental conditions?")
                user_input = input("You: ")
                setup_inputs['conditions'] = user_input
                continue
            
            if not setup_inputs['num_replicates']:
                print("Assistant: How many replicates do you have per condition?")
                user_input = input("You: ")
                setup_inputs['num_replicates'] = user_input
                continue

            # Process natural language input
            command_suggestion = get_pipeline_parameters(setup_inputs, context)
            # Update setup_info with extracted parameters
            for param, value in command_suggestion.items():
                if param in setup_info:
                    setup_info[param] = value
            
            # Check if all setup parameters are collected
            if all_parameters_collected(setup_info, ['norm_file', 'pw_data', 'base_dir', 'conditions', 'num_replicates']):
                # Generate base command
                base_command = (
                f"python pipeline_CLI.py "
                f"--norm-file \"{setup_info['norm_file']}\" "
                f"--pw-data \"{setup_info['pw_data']}\" "
                f"--base-dir \"{setup_info['base_dir']}\" "
                f"--num-replicates {setup_info['num_replicates']} "
            )

                # Add each condition as a separate --conditions argument
                for condition in setup_info['conditions']:
                    base_command += f"--conditions \"{condition}\" "

                
                # Generate setup command
                setup_command = (
                    f"{base_command}setup "
                )
                print("\nAssistant: Great, let's move on:")
                
                # Move to run-all stage
                current_stage = "run-all"
                continue

            
        # Run-all command stage
        elif current_stage == "run-all":
            # print the values of runall_inputs
            print(f"\nAssistant: Current analysis parameters: {runall_inputs}")

            # Prompt to collect run-all parameters
            if not runall_inputs['model_organism']:
                print("Assistant: What was your model organism for this experiment?")
                user_input = input("You: ")
                runall_inputs['model_organism'] = user_input
                continue

            if not runall_inputs['pw_interest']:
                print("Assistant: Which pairwise comparisons are you interested in? (e.g., 'ConditionA vs ConditionB, ConditionA vs ConditionC')")
                user_input = input("You: ")
                runall_inputs['pw_interest'] = user_input
                continue

            if runall_inputs['log2fc_threshold'] is None:
                print("Assistant: Would you like to specify a log2 fold change threshold? (default is 1, press Enter to use default)")
                user_input = input("You: ")
                runall_inputs['log2fc_threshold'] = user_input if user_input else "1"
                continue
            
            if runall_inputs['padj_threshold'] is None:
                print("Assistant: Would you like to specify an adjusted p-value threshold? (default is 0.05, press Enter to use default)")
                user_input = input("You: ")
                runall_inputs['padj_threshold'] = user_input if user_input else "0.05"
                continue
            
            if runall_inputs['enrich_sig_cutoff'] is None:
                print("Assistant: Would you like to specify an enrichment significance cutoff? (default is 0.05, press Enter to use default)")
                user_input = input("You: ")
                runall_inputs['enrich_sig_cutoff'] = user_input if user_input else "0.05"
                continue
                
            if runall_inputs['skip_steps'] is None:
                print("Assistant: Would you like to skip any steps? Enter the steps to skip or press Enter for none.")
                print("Assistant: Available steps to skip: boxplots, correlation-heatmap, pca, ma-plots, volcano-plots, heatmap, go, kegg")
                user_input = input("You: ")
                runall_inputs['skip_steps'] = [user_input] if user_input else []
                continue

            # Process natural language input
            command_suggestion = get_pipeline_parameters(runall_inputs, context)
            # Update runall_info with extracted parameters
            for param, value in command_suggestion.items():
                if param in runall_info:
                    runall_info[param] = value
            
            # Check if all run-all parameters are collected
            if all_parameters_collected(runall_info, ['model_organism', 'pw_interest', 'log2fc_threshold', 'padj_threshold', 'enrich_sig_cutoff']):
                # Generate run-all command
                runall_command = (
                    f"{base_command}run-all "
                    f"--model-organism \"{runall_info['model_organism']}\" "
                )
                
                # Add each pairwise interest as a separate --pw-interest argument
                for pw in runall_info['pw_interest']:
                    runall_command += f"--pw-interest \"{pw}\" "
                
                # Add optional parameters if not default
                if runall_info['log2fc_threshold'] != "1":
                    runall_command += f" --log2fc-threshold {runall_info['log2fc_threshold']}"
                
                if runall_info['padj_threshold'] != "0.05":
                    runall_command += f" --padj-threshold {runall_info['padj_threshold']}"
                
                if runall_info['enrich_sig_cutoff'] != "0.05":
                    runall_command += f" --enrich-sig-cutoff {runall_info['enrich_sig_cutoff']}"
                
                # Add skip flags
                for step in runall_info['skip_steps']:
                    if step:
                        runall_command += f" --skip-{step}"
                
                print("\nAssistant: Based on your inputs, here is the setup command:")
                print(f"\n{setup_command}\n")
                print("\nAssistant: And here is the run-all command:")
                print(f"\n{runall_command}\n")
    
                current_stage = "confirmation"
                continue
        
        # CONFIRMATION STAGE
        elif current_stage == "confirmation":
            # Ask if user wants to run commands
            print("Assistant: Would you like to run these commands? (yes/no/modify)")
            user_input = input("You: ")

            if user_input.lower() == "yes":
                print("Assistant: Running the commands...")
                try:
                    # Execute the setup command
                    print(f"\nExecuting: {setup_command}")
                    subprocess.run(setup_command, shell=True, check=True)
                    print(f"\nSetup command executed successfully! You should see the standardized directory structure in {setup_info['base_dir']}.")

                    # Execute the run-all command
                    print(f"\nExecuting: {runall_command}")
                    subprocess.run(runall_command, shell=True, check=True)
                    print('\nAnalysis completed successfully! You can find the results in the specified base directory.')
                    
                    return [setup_command, runall_command]
                except subprocess.CalledProcessError as e:
                    print(f"Assistant: Error executing command: {e}")
            elif user_input.lower() == "no":
                print("Assistant: Commands generated but not executed. You can copy and run them manually.")
                return [setup_command, runall_command]
            else:
                current_stage = "modify"
                continue
        
        # MODIFICATION STAGE
        elif current_stage == "modify":
            # Ask user if they want to modify setup or run-all parameters
            print("Assistant: What would you like to modify? Type 'setup' to modify setup arguments or 'run-all' to modify run-all arguments.")
            user_input = input("You: ")

            if user_input.lower() == "setup":
                # Show current setup parameters
                print("\nAssistant: Current setup parameters:")
                for param, value in setup_info.items():
                    print(f"  - {param}: {value}")
                
                current_stage = "select_setup_params"
                continue
                
            elif user_input.lower() == "run-all":
                # Show current run-all parameters
                print("\nAssistant: Current run-all parameters:")
                for param, value in runall_info.items():
                    print(f"  - {param}: {value}")
                
                current_stage = "select_runall_params"
                continue
                
            else:
                print("Assistant: Invalid selection. Type 'setup' to modify setup parameters or 'run-all' to modify run-all parameters.")
                current_stage = "modify"
                continue

        # PARAMETER SELECTION STAGE FOR SETUP
        elif current_stage == "select_setup_params":
            # Ask which parameters to modify
            print("\nAssistant: Which parameters would you like to modify? (Enter parameter names separated by commas, or 'all' for all parameters)")
            user_input = input("You: ")

            if user_input.lower() == "all":
                # Reset all setup parameters
                setup_info = {key: None for key in setup_info}
                print("Assistant: All setup parameters will be reset. Let's start again.")
            else:
                # Parse the parameter names to reset
                params_to_reset = [param.strip() for param in user_input.split(",")]
                valid_params = []
                
                # Reset only the specified parameters
                for param in params_to_reset:
                    if param in setup_info:
                        setup_info[param] = None
                        valid_params.append(param)
                    else:
                        print(f"Assistant: Warning: '{param}' is not a valid setup parameter and will be ignored.")
                
                if not valid_params:
                    print("Assistant: No valid parameters were selected. Please try again.")
                    current_stage = "select_setup_params"
                    continue
                    
                print(f"Assistant: The following parameters will be reset: {', '.join(valid_params)}.")
            
            # Move to setup stage to collect new values
            current_stage = "setup"
            continue

        # PARAMETER SELECTION STAGE FOR RUN-ALL
        elif current_stage == "select_runall_params":
            # Ask which parameters to modify
            print("\nAssistant: Which parameters would you like to modify? (Enter parameter names separated by commas, or 'all' for all parameters)")
            user_input = input("You: ")

            if user_input.lower() == "all":
                # Reset all run-all parameters
                runall_info = {key: None for key in runall_info}
                print("Assistant: All run-all parameters will be reset. Let's start again.")
            else:
                # Parse the parameter names to reset
                params_to_reset = [param.strip() for param in user_input.split(",")]
                valid_params = []
                
                # Reset only the specified parameters
                for param in params_to_reset:
                    if param in runall_info:
                        runall_info[param] = None
                        valid_params.append(param)
                    else:
                        print(f"Assistant: Warning: '{param}' is not a valid run-all parameter and will be ignored.")
                
                if not valid_params:
                    print("Assistant: No valid parameters were selected. Please try again.")
                    current_stage = "select_runall_params"
                    
                print(f"Assistant: The following parameters will be reset: {', '.join(valid_params)}.")
            
            # Move to run-all stage to collect new values
            current_stage = "run-all"
            continue
            


if __name__ == "__main__":
    interactive_command_builder()

