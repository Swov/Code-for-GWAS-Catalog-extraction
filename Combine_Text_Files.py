import io

def append_file_excluding_first_line(file1_path, file2_path, output_path):
    # Read the contents of the first file
    with open(file1_path, 'r', encoding = "utf-8") as file1:
        content1 = file1.read()
    
    # Read the contents of the second file, excluding the first line
    with open(file2_path, 'r', encoding = "utf-8") as file2:
        lines = file2.readlines()
        # Skip the first line
        content2 = ''.join(lines[1:])
    
    # Combine the contents
    combined_content = content1 + content2
    
    # Write the combined content to the output file
    with open(output_path, 'w', encoding = "utf-8") as output_file:
        output_file.write(combined_content)

i = 1
append_file_excluding_first_line(f'Trial_File_Page{0}.txt',f'Trial_File_Page{i}.txt','Combined.txt')
while i < 141:
    append_file_excluding_first_line('Combined.txt',f'Trial_File_Page{i}.txt','Combined.txt')
    i += 1
print('Finished')

