FISHR2:
	g++ -O2 -w -o ./binaries/FISHR2 -I /work/KellerLab/opt/bin/include -I include   ./src/Ped.cpp ./src/ErrorFinderManager.cpp ./src/Consolidator.cpp ./src/ErrorCalculator.cpp ./src/Compute.cpp ./src/Gen.cpp ./src/Haps.cpp ./src/Sample.cpp ./src/GERMLINE_0001.cpp ./src/GERMLINE.cpp ./src/Share.cpp ./src/Chromosome.cpp ./src/ChromosomePair.cpp ./src/HMIndividualsExtractor.cpp ./src/MarkerSet.cpp ./src/Individual.cpp ./src/Individuals.cpp ./src/InputManager.cpp ./src/MatchFactory.cpp ./src/MatchesBuilder.cpp ./src/NucleotideMap.cpp ./src/PEDIndividualsExtractor.cpp ./src/Match.cpp ./src/PolymorphicIndividualsExtractor.cpp ./src/SNP.cpp ./src/SNPPositionMap.cpp ./src/SNPs.cpp

PARAMETER_FINDER:
	g++ -O2 -w  -I /work/KellerLab/opt/bin/include -o  ./utilities/parameter_finder_binaries/parameter_finder ./utilities/parameter_finder/ErrorFinderMain.cpp

IE_CALCULATOR:
	g++ -O2 -w -I /work/KellerLab/opt/bin/include -o ./utilities/ie_calculator_binaries/ie_calculator ./utilities/ie_calculator/src/IBG-ProjectFISHR-B.cpp ./utilities/ie_calculator/src/HandleFlags.cpp ./utilities/ie_calculator/src/ReadFiles.cpp ./utilities/ie_calculator/src/Compute.cpp ./utilities/ie_calculator/src/Ibd.cpp ./utilities/ie_calculator/src/Bmid.cpp ./utilities/ie_calculator/src/Ped.cpp 

GAP:
	g++ -O2 -w -std=c++0x ./utilities/gap/gap.cpp -o ./utilities/gap_binaries/gap

clean:
	rm -rf ./binaries/FISHR2 ./utilities/parameter_finder_binaries/parameter_finder ./utilities/ie_calculator_binaries/ie_calculator
