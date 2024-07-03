import re
import csv
import numpy as np
import os
import pandas as pd
from shutil import copyfile
from  PWMpredictor.code  import main_PWMpredictor
from BindZF_predictor.code import main_bindzfpredictor_predict


def find_zf_binding_domains(protein_seq):
    # Define the correct regular expression pattern with character classes for amino acids
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    protein_seq = protein_seq.replace(" ", "")
    zf_pattern = re.compile(
        r'(?P<zf_binding>[{0}]{{2}}C[{0}]{{2,4}}C(?P<zf_center>[{0}]{{12}})[H[{0}]{{3,5}}H)'.format(
            amino_acids))

    # Find all matches in the protein sequence
    matches = zf_pattern.finditer(protein_seq)
    zf_seq_list = []

    # Extract and print the matched ZF binding domains
    for match in matches:
        zf_domain = match.group('zf_binding')
        start, end = match.span()

        zf_center_bd = match.group('zf_center')
        start_center , end_center = match.span('zf_center')  # Span of the {{12}} part

        # Extract 40 amino acids on each side of the ZF binding domain
        left_context = protein_seq[max(0, start_center - 40):start_center]
        right_context = protein_seq[end_center:min(end_center + 40, len(protein_seq))]

        # Combine ZF domain and its context
        zf_with_context = left_context + zf_center_bd + right_context
        # list of tuples, where each tuple is the 12
        # center aa and the 92 length aa (center with neighbors)

        zf_seq_list.append((zf_center_bd, zf_with_context))

    return  zf_seq_list


def create_zf_dataset(zf_lst,num_protein):
    print("Creating ZF dataset...")
    print(zf_lst)
    fieldnames = ['zf', 'label', 'seq', 'group']
    file_path = f'/content/DeepZF/Results/zf_40_dataset{num_protein}.csv'

    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for zf_sequence in zf_lst:
            # Fill 'label' and 'group' with default values of 0
            writer.writerow({'zf': zf_sequence[0], 'label': 0.0,
                             'seq': zf_sequence[1], 'group': 0})

    print(f"ZF dataset created successfully. Saved to {file_path}")


def create_pwm_dataset(num_protein):
    # Read TSV file into a pandas DataFrame without header
    tsv_file_path = f'/content/DeepZF/Results/results{num_protein}.tsv'
    tsv_df = pd.read_csv(tsv_file_path, delimiter='\t', header=None)

    # Add a column name to t+he TSV file
    tsv_df.columns = ['probabilities'] + list(tsv_df.columns[1:])

    # Read CSV file into a pandas DataFrame
    csv_file_path = f'/content/DeepZF/Results/zf_40_dataset{num_protein}.csv'
    csv_df = pd.read_csv(csv_file_path)

    # Concatenate the first column of TSV with the first column of CSV
    result_df = pd.concat([csv_df.iloc[:, 0], tsv_df.iloc[:, 0]], axis=1)

    # Write the result DataFrame to a new CSV file
    result_csv_file_path = 'Data/PWMpredictor/PWM_perdictor_input.csv'
    result_df.to_csv(result_csv_file_path, index=False)

    print(f"Concatenation completed. Result saved to {result_csv_file_path}.")

     # Load the CSV file into a pandas DataFrame
    df = pd.read_csv(result_csv_file_path)

    # Define the condition to filter rows based on the 'prob' column
    condition = df['probabilities'] > 0.5

    # Use the condition to filter rows and create a new DataFrame without those rows
    filtered_df = df[condition]

    # Save the filtered DataFrame back to a new CSV file or overwrite the original file
    filtered_df.to_csv(result_csv_file_path, index=False)

    print("the filter was done")
     # Load the updated CSV file into a pandas DataFrame
    df = pd.read_csv(result_csv_file_path)
     # Split the values in the "zf" column and create new columns 'AA1' through 'AA12'
    df[['AA1', 'AA2', 'AA3', 'AA4', 'AA5', 'AA6', 'AA7', 'AA8', 'AA9', 'AA10',
        'AA11', 'AA12']] = df['zf'].apply(lambda x: pd.Series(list(x)))
    column_order = ['AA1', 'AA2', 'AA3',  'AA4', 'AA5',
                    'AA6', 'AA7', 'AA8', 'AA9', 'AA10',
                     'AA11', 'AA12']
    df =df[column_order]
    # Save the modified DataFrame back to the original CSV file or create a new CSV file
    df.to_csv(result_csv_file_path ,sep=' ' ,index=False)
    print("split to aa was done")

def clear_aa_columns():
  # Load the CSV file into a DataFrame
  df = pd.read_csv('c_rc_df_copy.csv', sep=' ')  # Assuming the values are separated by space


  # Define the columns to clear
  columns_to_clear = ['AA1', 'AA2', 'AA3', 'AA4', 'AA5', 'AA6', 'AA7', 'AA8', 'AA9', 'AA10', 'AA11', 'AA12','res_12']

  # Set all values in the specified columns to NaN
  df[columns_to_clear] = np.nan

  # Save the modified DataFrame back to the original CSV file (override)
  df.to_csv('c_rc_df_copy.csv', index=False, sep=' ')

  def override_aa_columns():
      # Load the original DataFrame
      original_df = pd.read_csv('c_rc_df_copy.csv', sep=' ')

      # Load the DataFrame from another CSV file
      second_df = pd.read_csv('PWM_perdictor_input.csv', sep=' ')

      # Define the columns to replace
      columns_to_replace = ['AA1', 'AA2', 'AA3', 'AA4', 'AA5', 'AA6', 'AA7',
                            'AA8', 'AA9', 'AA10', 'AA11', 'AA12']

      # Replace NaN values in the specified columns with values from the second DataFrame
      original_df[columns_to_replace] = original_df[columns_to_replace].fillna(
          second_df[columns_to_replace])

      # Concatenate 'AA1' to 'AA12' and store the result in 'res_12', replacing NaN values with an empty string
      original_df['res_12'] = original_df[
          ['AA1', 'AA2', 'AA3', 'AA4', 'AA5', 'AA6', 'AA7', 'AA8', 'AA9',
           'AA10', 'AA11', 'AA12']].apply(
          lambda row: ''.join(row.dropna().astype(str)), axis=1).replace('',
                                                                         'NA')

      # Save the modified DataFrame back to the original CSV file (override)
      original_df.to_csv('c_rc_df_copy.csv', index=False, sep=' ')

def override_aa_columns():
  # Load the original DataFrame
  original_df = pd.read_csv('c_rc_df_copy.csv', sep=' ')

  # Load the DataFrame from another CSV file
  second_df = pd.read_csv('PWM_perdictor_input.csv', sep=' ')

  # Define the columns to replace
  columns_to_replace = ['AA1', 'AA2', 'AA3', 'AA4', 'AA5', 'AA6', 'AA7', 'AA8', 'AA9', 'AA10', 'AA11', 'AA12']

  # Replace NaN values in the specified columns with values from the second DataFrame
  original_df[columns_to_replace] = original_df[columns_to_replace].fillna(second_df[columns_to_replace])


  # Concatenate 'AA1' to 'AA12' and store the result in 'res_12', replacing NaN values with an empty string
  original_df['res_12'] = original_df[['AA1', 'AA2', 'AA3', 'AA4', 'AA5', 'AA6', 'AA7', 'AA8', 'AA9', 'AA10', 'AA11', 'AA12']].apply(lambda row: ''.join(row.dropna().astype(str)), axis=1).replace('', 'NA')

  # Save the modified DataFrame back to the original CSV file (override)
  original_df.to_csv('c_rc_df_copy.csv', index=False, sep=' ')

def filter_csv():
    pwm_input_file = "../../Data/PWMpredictor/PWM_perdictor_input.csv"
    c_rc_df_file = "../../Data/PWMpredictor/c_rc_df_copy.csv"
    output_file = "../../Data/PWMpredictor/c_rc_df_filtered.csv"

    # Read the number of lines in PWM_perdictor_input.csv
    with open(pwm_input_file, 'r') as f:
        num_lines = sum(1 for line in f)
    print("num of lines:", num_lines)

    # Read the content of c_rc_df_copy.csv up to the specified number of lines
    with open(c_rc_df_file, 'r') as f:
        content = ''.join([next(f) for _ in range(num_lines)])

    # Write the content to c_rc_df_filtered.csv
    with open(output_file, 'w') as f:
        f.write(content)


def pwm_format():

  with open("predictions.txt", "r") as file:
      probabilities = [float(line.strip()) for line in file]

  # Create a list of lists, where each inner list contains four probabilities
  formatted_probs = [probabilities[i:i+4] for i in range(0, len(probabilities), 4)]

  with open("predictions.txt", "w") as output_file:
      for row in formatted_probs:
          output_file.write("\t".join(map(str, row)) + "\n")



def run_deepzf_for_protein(protein_seq, num_protein):
    # First model:
    os.chdir('/DeepZF')

    zf_lst = find_zf_binding_domains(protein_seq)
    create_zf_dataset(zf_lst, num_protein)

    input_file_path = f'/DeepZF/Results/zf_40_dataset{num_protein}.csv'  # Dynamic input file name based on num_protein
    print(input_file_path)
    print(os.getcwd())

    with open(input_file_path, 'r') as file:
        print(file.read())

    os.chdir('/DeepZF')
    with open('BindZF_predictor/code/model.p', 'wb') as model_file:
        for part in ['BindZF_predictor/code/x01',
                     'BindZF_predictor/code/x02']:  # Assuming x?? means x01, x02, etc.
            with open(part, 'rb') as part_file:
                model_file.write(part_file.read())

    output_file_path = f"/DeepZF/Results/results{num_protein}.tsv"
    print(output_file_path)
    print("trying to get into model1")

    main_bindzfpredictor_predict(input_file=input_file_path,
                          output_file=output_file_path,
                          model='BindZF_predictor/code/model.p',
                          encoder='BindZF_predictor/code/encoder.p', r=1)

    with open(output_file_path, 'r') as file:
        print(file.read())

    # Second model:
    create_pwm_dataset(num_protein)
    os.chdir('/DeepZF/Data/PWMpredictor')

    copyfile('c_rc_df.csv', 'c_rc_df_copy.csv')
    clear_aa_columns()
    override_aa_columns()

    filter_csv()
    os.chdir('/DeepZF/PWMpredictor/code')

    main_PWMpredictor(
        input_file='../../Data/PWMpredictor/c_rc_df_filtered.csv',
        output_file='predictions.txt', model='/DeepZF/transfer_model100.h5')

    pwm_format()

    with open('predictions.txt', 'r') as file:
        print(file.read(1000))  # Print the first 1000 characters
        print(sum(1 for line in file))  # Count the number of lines

    predictions_file_path = "/DeepZF/PWMpredictor/code/predictions.txt"
    return predictions_file_path
