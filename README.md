# INSPiRE
knowledge-based protein-protein INteraction Sites PREdictor

Documentation is under construction, stay tuned.
Briefly:
1 For a knowledge-base:
1.1 create index file with 'frontend/index.cpp';
1.2 optionally, create transformation file for aminoa acids codes with 'frontend/aminoacids.cpp';
1.3 create chosen feature files with 'frontend/features.cpp'; requisites are coordinates (and interfaces);
1.4 generate chosen types of subgraphs and edges with 'frontend/subgraphs.cpp';
1.5 create fingerprints from files generated in previous steps with 'frontend/fingerprints.cpp'.
2 Repeat the same for each query protein, beware of a different switch in 'frontend/fingerprints.cpp'.
3 Search for similar fingerprints with 'frontend/mine.cpp'.
4 Classify results of the previous step with 'frontend/classify.cpp'
