package hmmwithdirichletprior;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.DPFactory;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.SimpleEmissionState;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.State;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.dp.onehead.SingleDP;
import org.biojava.bio.dp.onehead.SingleDPMatrix;
import org.biojava.bio.symbol.AbstractAlphabet;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.DoubleAlphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SimpleAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

import java.util.*;

import org.biojava.bio.seq.*;
import org.biojava.bio.symbol.*;

public class CASINO {

	/**
	 * @param args
	 * 
	 * THIS IS A TEST, TO SEE HOW BIOJAVA HMM WORKS
	 * @throws Exception 
	 */
	
	public static void main(String[] args) throws Exception {
		
		
		CASINO app = new CASINO();
		MarkovModel casino = app.createCasino();
		DP dp=DPFactory.DEFAULT.createDP(casino);
	      StatePath obs_rolls = dp.generate(10);
	      
	      for(int i =1;i<obs_rolls.length();i++){
	    	  System.out.println(obs_rolls.symbolAt(i).getName());
	      }
	      System.exit(0);
	      
	      SymbolList roll_sequence = obs_rolls.symbolListForLabel(StatePath.SEQUENCE);
	      SymbolList x             = obs_rolls.symbolListForLabel(StatePath.STATES);
	      SymbolList[] res_array = {roll_sequence};
	      StatePath v = dp.viterbi(res_array, ScoreType.PROBABILITY);
	      
	      
	      ///
	      
	      DP dp2 = new SingleDP(casino);
	      SingleDPMatrix forwardMatrix = (SingleDPMatrix) dp2.forwardMatrix(res_array, ScoreType.PROBABILITY);
	      double forwardScore    = forwardMatrix.getScore();
	      System.out.printf("forward core: %g%n", forwardScore);
	      SingleDPMatrix backwordMatrix = (SingleDPMatrix) dp2.backwardMatrix(res_array, ScoreType.PROBABILITY);
	      
	      State[] states = forwardMatrix.states();
	      
	      ////
	      
	      //print out obs_sequence, output, state symbols.
	      for(int i = 1; i <= obs_rolls.length()/60; i++) {
	    	  for(int j= (i-1)*60+1; j<=i*60; j++){
	    	  
	    	  System.out.print(  roll_sequence.symbolAt(j).getName() );
	    	  
	    	 // System.out.print( x.symbolAt(j).getName().charAt(0));	    	 
	    	  //System.out.print( v.symbolAt(j).getName().charAt(1));
	    	  
	    	  }
	    	  System.out.println();
	    	  for(int j= (i-1)*60+1; j<=i*60; j++){
	    		  System.out.print( x.symbolAt(j).getName().charAt(0));
	    		  
	    	  }
	    	  System.out.println();
	    	  for(int j= (i-1)*60+1; j<=i*60; j++){
	    		  System.out.print( v.symbolAt(j).getName().charAt(1));
	    	  }
	    	  
	    	  //System.out.println();
	    	  System.out.println();
	    	  System.out.println();
	        //for(int j=i*60; 
	      } 

	}/*main*/

	public static MarkovModel createCasino() throws Exception {
		    Symbol[] rolls=new Symbol[6];
		    

		    //set up the dice alphabet
		    SimpleAlphabet diceAlphabet=new SimpleAlphabet();
		    diceAlphabet.setName("DiceAlphabet");
		    

		    for(int i=1;i<7;i++) {
		      try {
		        rolls[i-1]= AlphabetManager.createSymbol( Integer.toString(i),  Annotation.EMPTY_ANNOTATION);
		        diceAlphabet.addSymbol(rolls[i-1]);
		        //System.out.println(rolls[i-1].getName());
		      } catch (Exception e) {
		        throw new Exception("Can't create symbols to represent dice rolls"
		        		);
		      }
		    }
		   
		   // System.exit(0);
		    
		   //////////////////////// 
		    
		    int [] advance = { 1 };
		    Distribution fairD;
		    Distribution loadedD;
		    try {
		      fairD = DistributionFactory.DEFAULT.createDistribution(diceAlphabet);
		      loadedD = DistributionFactory.DEFAULT.createDistribution(diceAlphabet);
		    } catch (Exception e) {
		      throw new Exception("Can't create distributions");
		    }
		    EmissionState fairS = new SimpleEmissionState("fair", Annotation.EMPTY_ANNOTATION, advance, fairD);
		    EmissionState loadedS = new SimpleEmissionState("loaded", Annotation.EMPTY_ANNOTATION, advance, loadedD);
		    
		    ////////////////////////
		    
		    SimpleMarkovModel casino = new SimpleMarkovModel(1, diceAlphabet, "Casino");
		    try {
		      casino.addState(fairS);
		      casino.addState(loadedS);
		    } catch (Exception e) {
		      throw new Exception("Can't add states to model");
		    }
		    
		    //////
		    //define transitions
		    try {
		        casino.createTransition(casino.magicalState(),fairS);
		        casino.createTransition(casino.magicalState(),loadedS);
		        casino.createTransition(fairS,casino.magicalState());
		        casino.createTransition(loadedS,casino.magicalState());
		        casino.createTransition(fairS,loadedS);
		        casino.createTransition(loadedS,fairS);
		        casino.createTransition(fairS,fairS);
		        casino.createTransition(loadedS,loadedS);
		      } catch (Exception e) {
		        throw new Exception("Can't create transitions");
		      }
		      
		      ////////////////////
		      //set emission probabilities 
		      try {
		          for(int i=0;i<rolls.length;i++)   {
		            fairD.setWeight(rolls[i],1.0/6.0);
		            loadedD.setWeight(rolls[i], 0.1);
		          }
		          loadedD.setWeight(rolls[5],0.5);
		        } catch (Exception e) {
		          throw new Exception("Can't set emission probabilities");
		        }
		        
		        //set up transition scores.
		        try {
		          Distribution dist;
		          //from magical to other defined  transitions
		          dist = casino.getWeights(casino.magicalState());
		          dist.setWeight(fairS, 0.8);
		          dist.setWeight(loadedS, 0.2);
		          
		          //from fairS to other defined transitions
		          dist = casino.getWeights(fairS);
		          dist.setWeight(loadedS,               0.04);
		          dist.setWeight(fairS,                 0.95);
		          dist.setWeight(casino.magicalState(), 0.01);
		          
		          //from loadedS to other defined transitions
		          dist = casino.getWeights(loadedS);
		          dist.setWeight(fairS,                 0.09);
		          dist.setWeight(loadedS,               0.90);
		          dist.setWeight(casino.magicalState(), 0.01);
		        } catch (Exception e) {
		          throw new Exception("Can't set transition probabilities");
		        }
		        ///////////////////
		        
		        return casino;
	  }/*createCasino*/
	
	

}/*CASION*/
