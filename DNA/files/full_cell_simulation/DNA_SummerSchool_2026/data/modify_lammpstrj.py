def modify_lammpstrj(input_file, output_file, center_index, width_per_frame, max_index=54338):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        frame = 0
        line = fin.readline()
        while line:
            if line.startswith("ITEM: TIMESTEP"):
                fout.write(line)
                fout.write(fin.readline())
                fout.write(fin.readline())  # ITEM: NUMBER OF ATOMS
                num_atoms_line = fin.readline()
                fout.write(num_atoms_line)
                num_atoms = int(num_atoms_line)
                fout.write(fin.readline())  # ITEM: BOX BOUNDS
                fout.write(fin.readline())
                fout.write(fin.readline())
                fout.write(fin.readline())

                atoms_header_line = fin.readline()
                fout.write(atoms_header_line)
                header_columns = atoms_header_line.strip().split()[2:]
                id_index = header_columns.index('c_id_track')
                type_index = header_columns.index('c_type_track')

                # Compute index range
                half_width = width_per_frame * frame
                index_low = max(0, center_index - half_width)
                index_high = min(max_index, center_index + half_width)

                # Read and process atoms
                for _ in range(num_atoms):
                    atom_line = fin.readline()
                    fields = atom_line.strip().split()
                    atom_id = int(fields[id_index])

                    if index_low <= atom_id <= index_high:
                        fields[type_index] = "13"
                    fout.write(" ".join(fields) + "\n")
                frame += 1
            else:
                fout.write(line)
            line = fin.readline()

# Usage
modify_lammpstrj(
    input_file="summerschool.lammpstrj",
    output_file="modified.lammpstrj",
    center_index=27169,
    width_per_frame=45
)
