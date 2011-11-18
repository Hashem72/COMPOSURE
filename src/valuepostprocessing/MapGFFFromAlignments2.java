/**
 * 
 */
package valuepostprocessing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.Symbol;

import net.derkholm.nmica.utils.CliTools;


/**
 * @author 
 * Note: This is in fact MapGFFFromAlignments from Thomas Down biobits project. I have made only a slight changes in order to make it applicable to my data.
 *
 */
public class MapGFFFromAlignments2 {

	/**
	 * @param args
	 */
	   private File aliDir;
	    private File aliGffDir;
	    private String gffPattern = "(([^_]+)_([0-9]+)_([0-9]+))(\\.comp)?\\.gff";
	    private String aliPattern = "%s.fa";
	    private String aliGffOutputFile = null;
	    
	    
	    public void setAliPattern(String aliPattern) {
			this.aliPattern = aliPattern;
		}

		public void setGffPattern(String gffPattern) {
			this.gffPattern = gffPattern;
		}

		public void setAliDir(File aliDir) {
	        this.aliDir = aliDir;
	    }

	    public void setAliGffDir(File aliGffDir) {
	        this.aliGffDir = aliGffDir;
	    }
	    
	    public void setAliGffOutputFile(String afo){
	    	this.aliGffOutputFile = afo;
	    }/*setAliGffOutputFile*/

	
	public static void main(String[] args)  throws Exception {
		// TODO Auto-generated method stub
        MapGFFFromAlignments2 app = new MapGFFFromAlignments2();
        args = CliTools.configureBean(app, args);
        app.run(args);
	}/*main*/
	
	   private void run(String[] args)
       throws Exception
   {
       Pattern fnPattern = Pattern.compile(gffPattern);
              
       PrintWriter pw;
       if(aliGffOutputFile != null){//given the output file name, the data must be printed into it.
    	   File gffOutputFile  = new File(this.aliGffOutputFile);
    	   pw = new PrintWriter(new FileWriter(gffOutputFile));
       	}
       else{
    	   pw = new PrintWriter(new OutputStreamWriter(System.out));
       }
       final GFFWriter gffw = new GFFWriter(pw);
       for (File f : aliGffDir.listFiles()) {
    	   String fileName = f.getName();
    	   if(fileName.startsWith(".")){
    		   continue;
    	   }
           Matcher m = fnPattern.matcher(f.getName());
           if (m.matches()) {
               String baseName = m.group(1);
               final String chr = m.group(2);
               int blockMin = Integer.parseInt(m.group(3));
               int blockMax = Integer.parseInt(m.group(4));
               
               File aliFile = new File(aliDir, String.format(aliPattern, baseName));
               if (!aliFile.exists()) {
                   System.err.println("Couldn't find alignment for " + baseName);
                   continue;
               }
               
               BufferedReader seqStream;
               if (aliFile.getName().endsWith(".gz")) {
               	seqStream = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(aliFile))));
               } else {
               	seqStream = new BufferedReader(new FileReader(aliFile));
               }
               Sequence masterSeq = SeqIOTools.readFastaDNA(seqStream).nextSequence();
               seqStream.close();
               final int[] ungapMap = makeUngapMap(masterSeq, blockMin);
               
               GFFParser parser = new GFFParser();
               BufferedReader gffStream = new BufferedReader(new FileReader(f));
               parser.parse(gffStream, new GFFDocumentHandler() {

                   public void startDocument(String locator) {
                   }

                   public void endDocument() {
                   }

                   public void commentLine(String comment) {
                   }

                   public void recordLine(GFFRecord record) {
                       SimpleGFFRecord r2 = new SimpleGFFRecord(record);
                       r2.setSeqName(chr);
                       r2.setStart(ungapMap[record.getStart()]);
                       r2.setEnd(ungapMap[record.getEnd()]);
                       gffw.recordLine(r2);
                   }
                   
               });
               gffStream.close();
           } else {
               System.err.println("Bad name " + f.getName());
           }
       }
       pw.flush();
   }

   private int[] makeUngapMap(Sequence masterSeq, int blockMin) {
       int[] ugm = new int[masterSeq.length() + 1];
       int pos = blockMin - 1;
       Symbol gap = masterSeq.getAlphabet().getGapSymbol();
       for (int p = 1; p <= masterSeq.length(); ++p) {
           if (masterSeq.symbolAt(p) != gap) {
               ++pos;
           }
           ugm[p] = pos;
       }
       return ugm;
   }
}/*MapGFFFromAlignments2*/
