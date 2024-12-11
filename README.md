# Blastn-C
## Description
This is the final PJ for the course of CLRS(H) at FDU.
The project has implemented a simplified version of the BLAST algorithm to **compare a query sequence with a database of sample sequences**.

## Usage

- Place your query sequences and sample sequences in the `seq` folder.
  - The query sequence should be named `query.fasta`
  - The sample sequence should be named `sample.fasta` 
  - the sample sequence file can contain multiple sequences.
- You can modify the macros in `main.c` to change the location where your sequences are stored.
- The alignment results will be output in the `output` folder, with the same name as the corresponding sample sequences.
- It is recommended to run the main program in the `bin` directory. 

```bash
$ make all
$ cd bin
$ ./main
```

