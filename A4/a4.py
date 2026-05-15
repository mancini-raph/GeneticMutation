import math

def count_differences(sequence_a, sequence_b):
    """Count the differences between two sequences."""
    differences_count = 0
    for index in range(min(len(sequence_a), len(sequence_b))):
        if sequence_a[index] != sequence_b[index]:
            differences_count += 1
    return differences_count

def jukes_cantor_distance(normalized_differences):
    """Calculate evolutionary distance using the Jukes-Cantor model."""
    distance = -0.75 * math.log(1 - (4.0 / 3.0) * normalized_differences)
    return distance

def calculate_distance_matrix(sequences):
    """Calculate the distance matrix of sequences."""
    num_sequences = len(sequences)
    distance_matrix = [[0] * num_sequences for _ in range(num_sequences)]  # Create an n x n matrix filled with zeros

    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            differences = count_differences(sequences[i][1], sequences[j][1])
            normalized_differences = differences / len(sequences[i][1])
            evolutionary_distance = jukes_cantor_distance(normalized_differences)
            final_distance = evolutionary_distance * len(sequences[i][1])
            # Update the matrix symmetrically
            distance_matrix[i][j] = distance_matrix[j][i] = final_distance

    return distance_matrix

def parse_sequences_from_file(file_path):
    """Parse sequences from a file in FASTA-like format, along with their headers."""
    parsed_sequences = []
    current_sequence = ''
    header = ''

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if current_sequence:
                    # Append the previous sequence before starting a new one
                    parsed_sequences.append((header, current_sequence))
                    current_sequence = ''
                # Extract and abbreviate the header
                header = line.strip()[1:].split('_')[0]
            else:
                # Add line to the current sequence
                current_sequence += line.strip()

        # Add the last sequence if the file does not end with '>'
        if current_sequence:
            parsed_sequences.append((header, current_sequence))

    return parsed_sequences

# File path to the sequences file
file_path = 'sequence.txt'

# Parse sequences from the file
parsed_sequences = parse_sequences_from_file(file_path)

# Calculate the distance matrix
distance_matrix = calculate_distance_matrix(parsed_sequences)

# Print the distance matrix with headers and aligned values
# Extract headers from sequences for labeling the matrix
headers = [sequence[0] for sequence in parsed_sequences]

# Create the header line for the matrix
header_line = '\t\t' + '\t'.join(f'{h:<10}' for h in headers)
print(header_line)

# Print each row of the matrix
for i, row in enumerate(distance_matrix):
    # Format each row for aligned printing
    formatted_row = '\t'.join(f'{cell:<10.3f}' for cell in row)
    # Print the header and the row
    print(f'{headers[i]:<7}\t{formatted_row}')
