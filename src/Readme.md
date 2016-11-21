TO COMPILE:

g++ -O2 -o GERMLINE2 -I include GERMLINE_0001.cpp GERMLINE.cpp Share.cpp Chromosome.cpp ChromosomePair.cpp HMIndividualsExtractor.cpp MarkerSet.cpp Individual.cpp Individuals.cpp InputManager.cpp MatchFactory.cpp MatchesBuilder.cpp NucleotideMap.cpp PEDIndividualsExtractor.cpp Match.cpp PolymorphicIndividualsExtractor.cpp SNP.cpp SNPPositionMap.cpp SNPs.cpp






TO RUN:
./GERMLINE2 -mapfile ./Beagle.Phased.Group2.1k.map -pedfile ./Beagle.Phased.Group2.1k.ped -outfile GL_OUT -bin_out -bits 20 -err_hom 0 -err_het 0 -min_m 3 -homoz  -w_extend -h_extend 

./GERMLINE2 -mapfile ./22.test.match.map -pedfile ./22.test.ped -outfile GL_OUT -bin_out -bits 20 -err_hom 0 -err_het 0 -min_m 3 -homoz  -w_extend -h_extend 

./GERMLINE2 -mapfile ./18.q.phased.cM.map -pedfile ./18.q.phased.ped -outfile GL_OUT -bin_out -bits 20 -err_hom 0 -err_het 0 -min_m 3 -homoz -w_extend -h_extend

./GERMLINE2 -mapfile ./22.q.phased.cM.map -pedfile ./22.q.phased.ped -outfile GL_OUT -bin_out -bits 20 -err_hom 0 -err_het 0 -min_m 3 -homoz -w_extend -h_extend
