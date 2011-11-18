/**
 * 
 */
package valuepostprocessing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;



/**
 * @author hk3 04-05-2011
 *Files in Alignment_Eisenlab are named based on their coordinates from r4 ie each file in the directory has a name like: chr_startr4_endr4. This code is to rename them based on
 *their coordinates in r5, ie chr_startr5_endr5. conversion tool is what that exists in flybase. This conversion tool takes chr:start..end format as its  input and output, therefore
 *we will read all files in a directory and then print them  with chr:start..end format into a file. This file will is used in flybase to convert into r5 coordinates. The output is 
 *a file in which each line is: "chr:startr4..endr4 chr:startr5..endr5" (note there is a tab between old and new coordinates). Then this file is read here in our code to rename the files.  
 *
 */
public class ChangeFileNamesInADirectory {

	/**
	 * @param args
	 */
	
	private String dirName;
	private String outputFileName;
	private String fileOfListOfConvertedNames ;
	private String chr;
	
	//GETTERS
	public String getDirName(){
		return dirName;
	}/*getDirName*/
	
	public String getOutputFileName(){
		return outputFileName;
	}/*getOutputFileName*/
	
	public String getFileOfListOfConvertedNames(){
		return fileOfListOfConvertedNames;
	}/*getFileOfListOfConvertedNames*/
	public String getChr(){
		return chr;
	}/*getChr*/
	
	//SETTERS
	public void setDirName(String dn){
		this.dirName = dn;
	}/*setDirName*/
	
	public void setOutputFileName(String ofn){
		this.outputFileName = ofn;
	}/*setOutputFileName*/
	
	public void setFileOfListOfConvertedNames(String flcn){
		this.fileOfListOfConvertedNames = flcn;
	}/*setFileOfListOfConvertedNames*/
	
	public void setChr(String chr){
		this.chr = chr;
	}/*setChr*/
	
	
	//OTHERS
	public void readNameOfFilesInADirectoryReformatThemAndWriteThemIntoAFile(String dirName, String ouputFileName, String chr) throws IOException{
		
		ArrayList<String> newlyFormattedFileNames = new ArrayList<String>();
		File dir = new File(dirName);
		String [] files = dir.list();
		if( files == null){
			System.out.println("No file found!");
		}
		else{
			for(int i= 0; i< files.length; i++){
				String oneFileName = files[i];
				
				
				if(oneFileName.startsWith(chr)){
					String [] fileNameAndExtensionSplitted = oneFileName.split("\\.");
					String fileNameWithoudExtension        = fileNameAndExtensionSplitted[0];
					//System.out.println(fileNameWithoudExtension);
					String [] fileNameSplitted = fileNameWithoudExtension.split("\\_");
					String fileNameNewFormat  = fileNameSplitted[0]+":"+fileNameSplitted[1]+".."+fileNameSplitted[2];
					newlyFormattedFileNames.add(fileNameNewFormat);
				
				}
			}
		}
		
		File outputFile = new File(ouputFileName);
		PrintWriter pw  = new PrintWriter( new FileWriter(outputFile));
		
		System.out.println("number of files is: " + newlyFormattedFileNames.size());
		for(int i= 0;i<newlyFormattedFileNames.size();i++){
			pw.println(newlyFormattedFileNames.get(i));
			//System.out.println(newlyFormattedFileNames.get(i));
		}
		pw.close(); 		
	}/*readNameOfFilesInADirectoryReformatThemAndWriteThemIntoAFile*/
	
	
 public void getListOfFilesInADirecotryWithExtensionRemoved(String dirName)throws IOException{
	 String outputFileName = dirName+"_List_of_Alignment_Files.txt";
	 ArrayList<String> filesList = new ArrayList<String>();
	 File dir  = new File(dirName);
	 String [] files = dir.list();
	 if(files == null){
		 System.err.println("No file found!");
	 }
	 else{
		 for(int i=0; i<files.length; i++){
			 String []oneFileSplitted = files[i].split("\\.");
			 filesList.add(oneFileSplitted[0]);
		 }
	 }
	// System.out.println(dirName);
	 //System.out.println(outputFileName);
	 File outputFile = new File(outputFileName);
	 PrintWriter pw = new PrintWriter(new FileWriter(outputFile));
	 for(int i= 0; i< filesList.size();i++){
		 pw.println(filesList.get(i));
	 }
	 pw.close();
	 
 }	/*getListOfFilesInADirecotryWithExtensionRemoved*/
	
	public void renameAFile(String oldFileName, String newFileName , String pathToTheDirectory){
		String oldFileFullName = pathToTheDirectory+"/"+oldFileName;
		String newFileFullName = pathToTheDirectory+"/"+newFileName;
		File oldFile = new File(oldFileFullName);
		File newFile = new File(newFileFullName);
		if(!oldFile.exists()){
			System.err.println("File " + oldFile + " not exists");
		}
		oldFile.renameTo(newFile);
	}/*renameAFile*/
	
	public static String makeNewFileName(String str){//Note: here each file is of form chr:start..end and this method has to reformat it as chr_start_end
		String string = str.replace(":", "_");
		string = string.replaceAll("\\..", "\\_");
		return string;
	}/*makeNewFileName*/
	
	
	public void readMapFileForEachChroAndAcordinglyChangeAlignmentFileNameFromOldReleasToNewRelease(String pathToFile, String chr) throws IOException{
		String mapFileName      = pathToFile+"/"+chr+"_Coordinates_in_r4_and_r5.txt";
		File mapFile            = new File(mapFileName);
		if(!mapFile.exists()){
			System.err.println("File "+ mapFileName + "not exists");
			return;
		}
		
		FileReader fr = new FileReader(mapFile);
		BufferedReader br = new BufferedReader(fr);
		
		String strLine;
	    while ((strLine = br.readLine()) != null)   {
	    	String [] coordinatesInTwoRelease = strLine.split("\\t");
	    	String coordinatesInReleaseFour  = makeNewFileName(coordinatesInTwoRelease[0]);
	    	String coordinatesInReleaseFive  = makeNewFileName(coordinatesInTwoRelease[1]);
	    	String oldFileName = coordinatesInReleaseFour+".fa";
	    	String newFileName = coordinatesInReleaseFive+".fa";
	    	if(!oldFileName.equals(newFileName)){
	    		System.out.println("old = " + oldFileName + "new = "  + newFileName);
	    	}
	    	String path = pathToFile+"/"+chr;
	    	String pathToTest = "/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB/TEST";
	    	renameAFile(oldFileName,newFileName,path);
	    	
	    }
	    fr.close();

		
	}
	
	public static void main(String[] args) throws IOException{
		
		// TODO Auto-generated method stub
		
		String dirName       = args[0];
		String ouputFileName = args[1];
		String chr           = args[2];
		
		ChangeFileNamesInADirectory directory = new ChangeFileNamesInADirectory();
		//ArrayList<String> newlyFormattedFileNames = new ArrayList<String>();
		
		directory.setDirName(dirName);
		//directory.setOutputFileName(ouputFileName);
		directory.setChr(chr);
		//directory.getListOfFilesInADirecotryWithExtensionRemoved(dirName);
		//directory.readNameOfFilesInADirectoryReformatThemAndWriteThemIntoAFile(dirName, ouputFileName, chr);
		//directory.readMapFileForEachChroAndAcordinglyChangeAlignmentFileNameFromOldReleasToNewRelease("/Users/hk3/Desktop/Main/Composure_Droshophila_Model/ALIGNMENTS_EISENLAB", chr);
		
	}/*main*/

}/*ChangeFileNamesInADirectory*/

