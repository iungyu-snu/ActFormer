import sys

# Get the sequence and position from command-line arguments
sequence = sys.argv[1]
position = int(sys.argv[2])

# Get the residue at the specified position
residue = sequence[position - 1]  # Convert 1-based index to 0-based index

# Output the residue
print("The residue at position", position, "is:", residue)

