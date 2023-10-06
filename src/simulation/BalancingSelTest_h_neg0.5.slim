// set up a simple neutral simulation
initialize() {
	//Set constants
	defineConstant("outDir","[Specify your output directory]");
	defineConstant("selCoef", -0.01);
	defineConstant("balCoef", -0.5);
	defineConstant("N", 10000); //sample size
	defineConstant("simID", getSeed());
	defineConstant("seedID", getSeed());
	initializeMutationRate(2.5e-8);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1", m1, 1.0);
	initializeMutationType("m2", balCoef, "f", selCoef);    // balanced
	initializeGenomicElementType("g2", m2, 1.0);
	initializeGenomicElement(g1, 0, 9999);
	initializeRecombinationRate(2.5e-8);
}
1 {
	// save this run's identifier, used to save and restore
	sim.addSubpop("p1", N);
}
100000 late() {
	sim.addSubpopSplit("p2", N, p1);
	sim.addSubpopSplit("p3", N, p1);
	sim.outputFull(outDir+"temp_" + simID+"_" + seedID+ ".txt");
	target = sample(p2.genomes, 1);
	target.addNewDrawnMutation(m2, 4999);
}
///////////////////////////////////////////////////////////
101000 late() {
	mut = sim.mutationsOfType(m2);
	if (size(mut) == 1){cat(seedID + "101000: Established\n");}
	else
	{
		cat(seedID + "101000: LOST – RESTARTING\n");
		// go back to generation of split
		sim.readFromPopulationFile(outDir+"temp_" + simID+"_" + seedID+ ".txt");
		// start a newly seeded run
		setSeed(rdunif(1, 0, asInteger(2^62) - 1));
		// re-introduce the sweep mutation
		target = sample(p2.genomes, 1);
		target.addNewDrawnMutation(m2, 4999);
	}
}
200000 late() {
	mut = sim.mutationsOfType(m2);
	if (size(mut) == 1){cat(seedID + "200000: Established\n");}
	else
	{
		cat(seedID + "200000: LOST – RESTARTING\n");
		// go back to generation of split
		sim.readFromPopulationFile(outDir+"temp_" + simID+"_" + seedID+ ".txt");
		// start a newly seeded run
		setSeed(rdunif(1, 0, asInteger(2^62) - 1));
		// re-introduce the sweep mutation
		target = sample(p2.genomes, 1);
		target.addNewDrawnMutation(m2, 4999);
	}
}
///////////////////////////////////////////////////////////
350000 late() {
	mut = sim.mutationsOfType(m2);
	if (size(mut) == 1)
	{
		if (sim.mutationFrequencies(NULL, mut) > 0.001)
		{
			cat(seedID + "350000: ESTABLISHED\n");
			cat(sim.mutationFrequencies(NULL, mut));
			sim.deregisterScriptBlock(self);
			genome2 = p2.sampleIndividuals(100).genomes;
			genome2.outputVCF(filePath= outDir + "Human_slim_" +simID + ".txt",simplifyNucleotides=T);
			genome3 = p3.sampleIndividuals(1).genomes;
			genome3.outputVCF(filePath=outDir + "Chimp_slim_" + simID + ".txt",simplifyNucleotides=T);
			sim.simulationFinished();
		}
		else{
			cat(sim.mutationFrequencies(NULL, mut));
			cat(seedID + "350000: Not Established – RESTARTING\n");
			// go back to generation of split
			sim.readFromPopulationFile(outDir+"temp_" + simID+"_" + seedID+ ".txt");
			// start a newly seeded run
			setSeed(rdunif(1, 0, asInteger(2^62) - 1));
			// re-introduce the sweep mutation
			target = sample(p2.genomes, 1);
			target.addNewDrawnMutation(m2, 4999);
		}
	}
	else
	{
		cat(seedID + "350000: LOST – RESTARTING\n");
		// go back to generation of split
		sim.readFromPopulationFile(outDir+"temp_" + simID+"_" + seedID+ ".txt");
		// start a newly seeded run
		setSeed(rdunif(1, 0, asInteger(2^62) - 1));
		// re-introduce the sweep mutation
		target = sample(p2.genomes, 1);
		target.addNewDrawnMutation(m2, 4999);
	}
}
