package hmmwithdirichletprior;

import java.io.File;

import net.derkholm.nmica.utils.CliTools;
import net.derkholm.nmica.utils.ConfigurationException;

public class SHORTSEQUENCERMOMOVAL {

	/**
	 * @param args 24-11-2011 Hashem Koohy	
	 * some of the alignment blocks in fly multiple alignment data are really short. This is to remove short alignment blocks from either comp or ser  or any other directory that contains block block related data.
	 * the form of files in this directory is assumed to be as chr_start_end.fileExtension.
	 * directory name and and a length threshold are passed as arguments. blocks with a length shorter than length threshold will be deleted.
	 */
	
	public String dirName;
	public int  lengthThreshold;
	
	public void setDirName(String dn){
		this.dirName = dn;
	}/*setDirName*/
	
	public void setLengthThreshold(int l){
		this.lengthThreshold = l;
	}/*setLengthThreshold*/
	
	public String getDirName(){
		return dirName;
	}/*getDirName*/
	
	public int getLengthThreshold(){
		return lengthThreshold;
	}/*getLengthThreshold*/
	
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		SHORTSEQUENCERMOMOVAL app = new SHORTSEQUENCERMOMOVAL();
		args = CliTools.configureBean(app, args);
		app.run(args);
	}/*main*/
	
	private void run(String [] args) throws Exception{
		deleteShortFilesFromShortBlocksInthisDirectory(lengthThreshold, dirName);
	}/*run*/
	
	public void deleteShortFilesFromShortBlocksInthisDirectory(int thresholdForLength, String fullNameOfDir){
		File dir = new File(fullNameOfDir);
		String [] children = dir.list();
		
		
		int  numberOfDeletedFiles = 0;
		int numberOfFiles = 0;
		if(children == null){
			System.out.println("Either directory is empty or is not a directory!");
		}
		else{
			
			for(String child:children){
				if(!child.startsWith(".")){
					numberOfFiles++;
					String blockSplited  [] = child.split("\\.");
					String blockId = blockSplited[0];
					String blockSpecifications [] = blockId.split("_");
					int blockStart = Integer.valueOf(blockSpecifications[1]);
					int blockEnd   = Integer.valueOf(blockSpecifications[2]);
					int blockLength = blockEnd - blockStart;
					if(blockLength < 0){
						System.out.println("Found block wiht negativ length!!!!  look at + " + child);
					}
					
					if(blockLength < lengthThreshold ){
						System.out.println("found some files to delete!");
						String fileFullName = dirName+"/"+child;
						deletFile(fileFullName);
						numberOfDeletedFiles++;
					}
				}
			}
		}
		System.out.println("There were " + numberOfFiles + " file(s) and from them " + numberOfDeletedFiles + " file(s) deleted." );
		
	}/*deleteShortFilesFromShortBlocksInthisDirectory*/
	
	public void deletFile(String file){
		File  f = new File(file);
		boolean fileExists = f.exists();
		if(fileExists){
			boolean success = f.delete();
			if(!success){
				System.out.println("Deletion failed!");
				System.exit(0);
			}
		}
	}/*deletFile*/

}/*SHORTSEQUENCERMOMOVAL*/
