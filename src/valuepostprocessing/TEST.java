/**
 * 
 */
package valuepostprocessing;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * @author hk3
 *
 */





public class TEST {
	public static double getMax(int i,double[][] X){
		double x_max = Double.MAX_VALUE;
		int numOfRows = X[1].length;
		ArrayList<Double> AL = new ArrayList<Double>();
		
			for(int k= 0; k< numOfRows; k++){
				AL.add(X[i-1][k]);
			}
		
		int maxIndex = -1;
		Object obj_max = Collections.max(AL);
		maxIndex =  AL.indexOf(obj_max); 
		System.out.println("for i = " + i + " max is " + obj_max + " and index is " + maxIndex);
		return x_max;
		
	}

	
	public  static class PAIR{
		public double first;
		public int second;
		PAIR(double x , int y){
			this.first = x;
			this.second = y;
		}
	}/*PAIR*/

	public static PAIR getMax_2(int i, double [][]x){
		
		double x_max = Double.MAX_VALUE;
		int numberOfRows = x[1].length;
		ArrayList<Double> AL = new ArrayList<Double>();
		for (int k =0; k<numberOfRows; k++){
			AL.add(x[i-1][k]);
		}
		int maxIndex = -1;
		Object obj_max = Collections.max(AL);
		maxIndex = AL.indexOf(obj_max);
		double max = (Double) obj_max;
		PAIR result = new PAIR(max,maxIndex);
		return result;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double aMatrix [][] = new double[5][5];
		Random generator = new Random();
		
		for(int i =0; i<5; i++){
			for(int j = 0; j<5; j++){
				aMatrix[j][i]= j*i + generator.nextInt( 10);
				System.out.print(aMatrix[j][i] + "   ");
			}
			System.out.println();
			if(i>0){
				PAIR  aPair = 	getMax_2(i,aMatrix );
				System.out.println("i = "+ i + " max = " + aPair.first + " index = " + aPair.second);

			}
		}
		
		System.out.println("*********************");
		
	/*	for(int i = 0;i<6;i++){
			if(i>0){
				PAIR  aPair = 	getMax_2(i,aMatrix );
				System.out.println("i = "+ i + " max = " + aPair.first + " index = " + aPair.second);
			}
		} */
		
	}/*main*/
	
}/*TEST*/
