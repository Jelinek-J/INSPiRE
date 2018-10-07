# INSPiRE
knowledge-based protein-protein INteraction Sites PREdictor. For details about how the algorithm works please read https://doi.org/10.1186/s12859-017-1921-4. 

Documentation is under construction, stay tuned. But briefly:

# Use of the framework #
You can use the framework in three ways:
1. Header-only library;
2. Individual programs to run individual parts of the INSPiRE algorithm
3. Single program to run INSPiRE algorithm

## Header-only library ##
This use is ideal if you want to update INSPiRE algorithm or extend it for some new usecases.
For this use, you only need to add header files from 'src/backend' and 'src/elemental' directories to your project (or set corresponding paths). 'src/backend' containts the basic parts of the algorithm. 'src/elemental' contains some elemental methods and constants and also hides all usage of boost library, so if you have some problems with boost, you need to rewrite this files only.

## Fragmented INSPiRE tools ##
This possibility is focused on the situation, when you want to optimize the INSPiRE algorithm's configuration, because some temporary files can be reused, so it is useless to compute them again and again.
Using the files mentioned in the previous subsection, you need to compile (each separately) following files from 'src/frontend' directory:
1. For a constraction of knowledge-base and queries:  
    1. 'frontend/index.cpp': creates index file of proteins and their residues, chains and models. The file is used by other tools.
    2. 'frontend/aminoacids.cpp': creates transformation file to convert aminoacids' three-letter codes to one-letter codes.
    3. 'frontend/features.cpp': creates files with extracted features. Actually, it is possible to extract coordinates, amino acid type (here can be used the file created in the previous point), temperature, and interfaces.
    4. 'frontend/subgraphs.cpp': extracts chosen types of subgraphs and edges in them.
    5. 'frontend/fingerprints.cpp': create fingerprints from files generated in previous steps.
2. 'frontend/mine.cpp': searches in knowledge-base for fingerprints similar to query fingerprints.
3. 'frontend/classify.cpp': classifies fingerprints found by the tool from the previous step.
4. 'frontend/predict.cpp': pronounces a prediction based on the statistics computed using the tool from the previous step.
5. 'frontend/optimize.cpp': finds the best parameters for predictors in 'frontend/predict.cpp'.

## Single INSPiRE tool ##
You will probably prefer this approach if you want to use the INSPiRE for its original purpose, i.e. prediction of new protein-protein interaction interfaces. Instead of all the '\*.cpp' you need only a single file 'frontend/inspire.cpp' that is roughly equivalent to pipeline of tools from the previous subsection that corresponds to the INSPiRE algorithm. At most, you can optionally use 'frontend/aminoacids.cpp' to change the original mapping of three-letters codes to one-letter codes. 

There is also a premaked knowledge-base to make the use of INSPiRE easier. To make this knowledge-base:
1. All proteins without DNA and RNA were downloaded from Protein Data Bank in pdb format.
2. Only complexes with 'REMARK 350' were taken.
3. All complexes with monomer as 'BIOMOLECULE 1' were filtered out.
4. All complexes with more than 50 000 residues were filtered out.
5. All complexes with 10 or more residues without specified carbon alpha were filtered out.
6. All complexes containing something else than aminoacids (modified aminacids are allowed) were filtered out.
7. All complexes containing a chain shorter than 20 aminoacids were filtered out.
8. All complexes containing more than 1% of residues without specified carbon alpha were filtered out.
9. Reduction of redundancies:
    1. Fingerprints were grouped by aminoacid type, fingerprint and interface.
    2. For each group, only one fingerprint from each protein was preserved. (So e.g. if it is a symmetric tetramer, three fingerprints were preserved and only one was preserved.)
    3. For each group of size k, leave only ⌈log_6(k)+1⌉, except the situation where k_I != k_N but ⌈log_6(k_I)+1⌉ == ⌈log_6(k_n)+1⌉ - in such a case, leave one more item in the bigger group.
