/**
 * 
 */
package valuepostprocessing;

import java.io.IOException;

import alignmentblock.AlignmentBlock;

/**
 * @author hk3
 * this is to  map a given gff files into   alignment blocks (files) and printout the results as  a gff file.
 *
 */
public class MapGFFToAlignmentBlockMain {

	/**
	 * @param args
	 */
	public static void main(String[] args)  throws IOException, Exception{
		// TODO Auto-generated method stub
		int    blockNumber                = Integer.parseInt(args[0]);
		String  mapFile                   = args[1];
		String  enhancerInputGffFile         = args[2];
		String dnaseIAllStagesInputGffFile   = args[3];
		String genesInputGffFile             = args[4];
		String enhancerOutputFile         = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/Composure_Output_Model_20051212_With_SegmentWindow_350/Enhancers.gff";
		String dnaseIAllStagesOutputFile  = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/Composure_Output_Model_20051212_With_SegmentWindow_350/DNaseI.gff";
		String genesOutputFile            = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/Composure_Output_Model_20051212_With_SegmentWindow_350/Genes.gff";
		String alignmentFile              =  "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/Alignments/fly_CAF1.1/alignments/"+args[0]+"/mavid.fa";

		MapGFFToAlignmentBlock alignmentblock = new MapGFFToAlignmentBlock();
		
		alignmentblock.setBlockNumber(blockNumber);
		alignmentblock.setMapFile(mapFile);
		alignmentblock.setEnhancersInputGffFile(enhancerInputGffFile);
		alignmentblock.setEnhancersOutputFile(enhancerOutputFile);
		alignmentblock.setGenesInputGfffile(genesInputGffFile);
		alignmentblock.setGenesOutputFile(genesOutputFile);
		alignmentblock.setDnaseIInputGffFile(dnaseIAllStagesInputGffFile);
		alignmentblock.setDnaseIOutputFile(dnaseIAllStagesOutputFile);
		alignmentblock.setAlignmentFile(alignmentFile);
				
		// to make sure we dont use old files or we dont append new files to the old ones.
		alignmentblock.deleteFile(enhancerOutputFile);
		alignmentblock.deleteFile(genesOutputFile);
		alignmentblock.deleteFile(dnaseIAllStagesOutputFile);
		
		//get this alignment block specifications:
		String [] specifications = alignmentblock.getSpecifications(mapFile, blockNumber);
		//run  
		alignmentblock.mapGffToAlignment(specifications, enhancerInputGffFile, alignmentFile, enhancerOutputFile);		
		alignmentblock.mapGffToAlignment(specifications, dnaseIAllStagesInputGffFile, alignmentFile, dnaseIAllStagesOutputFile);
		alignmentblock.mapGffToAlignment(specifications, genesInputGffFile, alignmentFile, genesOutputFile);

	}/*main*/

}/*MapGFFToAlignmentBlockMain*/
