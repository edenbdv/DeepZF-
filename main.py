from helper_functions_script import run_deepzf_for_protein

print("enter sequence")
seq = input().strip()
num_protein = 1  # Set the num_protein value as needed

predictions_file_path = run_deepzf_for_protein(seq, num_protein)
print(f"Predictions file created at: {predictions_file_path}")