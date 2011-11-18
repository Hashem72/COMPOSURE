package alignmentblock;
import java.io.*;





/**
 * @author hk3
 *17-03-2011
 */


public class AlignmentBlock {
	//instances
	private int blockNumber;
	private String mapFile;
	private String enhancersInputFile;
	private String dnaseIInputFile;
	private String genesInputfile;
	private String enhancersOutputFile;
	private String dnaseIOutputFile;
	private String genesOutputFile;
	
	
	//getters
	public int getBlockNumber(){
		return blockNumber;		
	}/*getBlockNumber*/
	
	public String getMapFile(){
		return mapFile;
	}/*getMapFile*/
	
	public String getEnhancerInputFile(){
		return enhancersInputFile;
	}/*getEnhancerInputFile*/
	
	public String getDnaseIInputFile(){
		return dnaseIInputFile;
	} /*getDnaseIInputFile*/
	
	public String getEnhancerOutputFile(){
		return enhancersOutputFile;
	}/*getEnhancerOutputFile*/
	
	public String getDnaseIOutputFile(){
		return dnaseIOutputFile;
	}/*getDnaseIOutputFile*/
	
	public String getGenesInputFile(){
		return genesInputfile;
	}/*getGenesInputFile*/
	
	public String getGenesOutputFile(){
		return genesOutputFile;
	}/*getGenesOutputFile*/
	
	//setters
	public void setBlockNumber(int i){
		this.blockNumber =i;
	}/*setBlockNumber*/
	
	public void setMapFile(String mF){
		this.mapFile = mF;
		
	}/*setMapFile*/
	
	public void setEnhancerInputFile(String eif){
		this.enhancersInputFile = eif;
	}/*setEnhancerInputFile*/
	
	public void setEnhancerOutputFile(String eof){
		this.enhancersOutputFile = eof;
	}/*setEnhancerOutputFile*/
	
	public void setDnaseIInputFile(String dif){
		this.dnaseIInputFile = dif;
	}/*setDnaseIInputFile*/
	
	public void setDnaseIOutputFile(String dof){
		this.dnaseIOutputFile = dof;
	}/*setDnaseIOutputFile*/
	
	public void setGenesInputFile(String gif){
		this.genesInputfile = gif;
	}/*setGenesInputFile*/
	public void setGenesOutputFile(String gof){
		this.genesOutputFile = gof;
	}/*setGenesOutputFile*/
	
	//delet file
	public void deleteFile(String file){
		File f = new File(file);
		boolean exists = f.exists();
		if(exists){
			boolean success = f.delete();
			if(!success){
				System.out.println("Deletion failed.");
				System.exit(0);
			}
		}
	}/*deleteFile*/
	
	//get block specifications from the map file. Note that the block number is the line number in the map file.
	public String [] getSpecifications(String file, int lineNumber) throws IOException{
		String [] specifications = null;
		try{
			FileInputStream fs = new FileInputStream(file);
			BufferedReader br = new BufferedReader(new InputStreamReader(fs));
			String lineOfInterest = null;
			for(int i=0; i<lineNumber; i++){
				lineOfInterest = br.readLine();
			}
			specifications = lineOfInterest.split("\t");
		}
		catch(FileNotFoundException fe){
			fe.fillInStackTrace();	
	}/*getSpecifications*/
		return specifications;
}/*AlignmentBlock*/
	
	//prints the block overlaps with any given features such as enchancers, DnaseI, genes and ect into the give file
	public void printBlockOverlaps(String [] blockSpecifications, String inputFile, String outputFile) throws IOException{
		String segmentChr    =  blockSpecifications[1];
		int segmentStart     = Integer.parseInt(blockSpecifications[2]);
		int segmentEnd       = Integer.parseInt(blockSpecifications[3]);
		try{
			FileReader      fr =  new  FileReader(inputFile);
			BufferedReader  br =  new  BufferedReader(fr);
			FileWriter      fw =  new  FileWriter(outputFile, true);
			BufferedWriter  bw =  new  BufferedWriter(fw);
			String  line;
			while ( (line = br.readLine()) != null){
				String [] Features = line.split("\t");
				String featureChr = Features[0];
				int    featureStart = Integer.parseInt(Features[3]);
				int    featureEnd   = Integer.parseInt(Features[4]);
				if ( (segmentChr.equals(featureChr) ) && (featureEnd > segmentStart )  && (segmentEnd  >  featureStart) ){
					String [] overlappedFeatures = line.split("\t");
					//make segmentStart as origin
					int featureStart_new = featureStart- segmentStart;
					int featureEnd_new   = featureEnd  - segmentStart;
					String modifiedLine  = overlappedFeatures[0] +"\t"+ overlappedFeatures[1] + "\t" + overlappedFeatures[2] + "\t" + featureStart_new + "\t" + featureEnd_new + "\t.\t+\t.\t.\n";
					bw.append(modifiedLine);
					}
				}
				bw.close();			
		}
		catch(FileNotFoundException fe){
			fe.fillInStackTrace();
		}
	}/*printBlockOverlaps*/
}/*AlignmentBlock*/


//featureEnd >= segmentStart)  && (featureStart <= segmentEnd )

