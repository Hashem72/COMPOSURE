/**
 * 
 */
package alignmentblock;

import java.io.IOException;

/**
 * @author hk3
 *
 */
public class AlignmentBlockTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stub

		int    blockNumber                = Integer.parseInt(args[0]);
		String  mapFile                   = args[1];
		String  enhancerInputFile         = args[2];
		String dnaseIAllStagesInputFile   = args[3];
		String genesInputFile             = args[4];
		String enhancerOutputFile         = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/Composure_Output_Model_20051212/Enhancers_Overlapped_With_Blocks.gff";
		String dnaseIAllStagesOutputFile  = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/Composure_Output_Model_20051212/DNaseI_Overlapped_With_Blocks.gff";
		String genesOutputFile            = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/Composure_Output_Model_20051212/genes_Overlapped_With_Blocks.gff";

				
		AlignmentBlock block = new AlignmentBlock();
		
		block.setEnhancerInputFile(enhancerInputFile);
		block.setEnhancerOutputFile(enhancerOutputFile);
		block.setDnaseIInputFile(dnaseIAllStagesInputFile);
		block.setDnaseIOutputFile(dnaseIAllStagesOutputFile);
		block.setGenesInputFile(genesInputFile);
		block.setGenesOutputFile(genesOutputFile);
	
		
		//make sure we dont use old files or we dont append new files to the old ones.
		block.deleteFile(enhancerOutputFile);
		block.deleteFile(dnaseIAllStagesOutputFile);
		block.deleteFile(genesOutputFile);

		
		block.setBlockNumber(blockNumber);
		block.setMapFile(mapFile);
		String [] specifications = block.getSpecifications(mapFile, blockNumber);
		
		//print overlaps into the given files
		block.printBlockOverlaps(specifications, enhancerInputFile, enhancerOutputFile);
		block.printBlockOverlaps(specifications, dnaseIAllStagesInputFile, dnaseIAllStagesOutputFile);
		block.printBlockOverlaps(specifications, genesInputFile, genesOutputFile);
		
	}/*main*/

}/*AlignmentBlockTest*/


