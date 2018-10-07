# INSPiRE
knowledge-based protein-protein INteraction Sites PREdictor. For details about how the algorithm works please read https://doi.org/10.1186/s12859-017-1921-4. 

Documentation is under construction, stay tuned. But briefly:

# Use of the framework #
You can use the framework in three ways:
1. Header-only library;
2. Individual programs to run individual parts of the INSPiRE algorithm
3. Single program to run INSPiRE algorithm

## Header-only library ##
This use is ideal if you want to update INSPiRE algorithm or extend for some new usecases.
For this use, you only need to add header files from src/backend and src/elemental directories to your project (or set corresponding paths).

## Fragmented INSPiRE tools ##
This possibility is focused on the situation, if you want to optimize the INSPiRE algorithm configuration, because some temporary files can be reused and it is useless to compute them again and again.
In addition to files from the previous usecase, you need to compile (each separately) following files from src/frontend directory:
1. For a constraction of knowledge-base and queries:  
  1. 'frontend/index.cpp': creates index file of proteins and their residues, chains and models. The file is used by other tools.
  2. 'frontend/aminoacids.cpp': creates transformation file to convert aminoacids three-letter codes to one-letter codes.
  3. 'frontend/features.cpp': creates files with extracted features. Actually, it is able to extract coordinates, amino acid type (here can be used the file created in the previous point), temperature, and interfaces.
  4. 'frontend/subgraphs.cpp': generate chosen types of subgraphs and edges in them.
  5. 'frontend/fingerprints.cpp': create fingerprints from files generated in previous steps.
2. 'frontend/mine.cpp': searches in knowledge-base for fingerprints similar to query fingerprints.
3. 'frontend/classify.cpp': classifies fingerprints found using the tool from the previous step.
4. 'frontend/predict.cpp': pronounces a prediction based on the statistics computed using the tool from the previous step.
5. 'frontend/optimize.cpp': finds the best parameters for predictors in 'frontend/predict.cpp'.

## Single INSPiRE tool ##
You will probably prefer this approach if you want to use the INSPiRE for its original purpose, i.e. prediction of new protein-protein interaction interfaces. Instead of all the '\*.cpp' you need only a single file 'frontend/inspire.cpp' that is roughly equivalent to pipeline of tools from the previous subsection that corresponds to the INSPiRE algorithm. At most, you can optionally use 'frontend/aminoacids.cpp' to change the original mapping of three-letters codes to one-letter codes. 
