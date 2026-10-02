The new vision for the workflow is as follows:
1. Add germline variant column to .full.mut.tsv with code in converting_from_old_format_to_new
2. Have one file that contains a function to create the tumor and normal from the .full.mut.tsv
3. Have a second file that runs on a single sample. It should offer accepting a .full.mut.tsv, in which case it runs the function from the first file before running, or just the normal and tumor, in which case it runs from file directly
4. Have a third file that runs over an entire directory by calling the second file
5. 
The idea is to create a pyramid of files calling each other (within the same directory) to maxmize code reuse
