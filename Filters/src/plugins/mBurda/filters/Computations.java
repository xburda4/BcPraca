package plugins.mBurda.filters;

import flanagan.complex.Complex;
import flanagan.complex.ComplexMatrix;
import flanagan.math.FourierTransform;

public class Computations {
	public double[][] getGaborKernel(){
		
		
		return null;
	}
	
	
	public static ComplexMatrix FourierTransform2D(double[][] input){//double[polohaVStlpci][polohaVRiadku]
		if(input == null) return null;
		ComplexMatrix output = new ComplexMatrix(input.length,input[0].length);
		/*Fourierova transformácia nad každým riadkom a následne priradenie hodnôt do matice output
		 * 
		 * Počítam s tým,že input je vo formáte [y][x]?
		 */
		ComplexMatrix rowMatrix;
		for(int y = 0;y<output.getNrow();y++){
			FourierTransform fourierHorizontal = new FourierTransform(input[y]);
			fourierHorizontal.transform();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal.getTransformedDataAsComplex());
			for(int x = 0;x<output.getNcol();x++){
				output.setElement(y, x, rowMatrix.getElementCopy(0, x));//malo by to byť dobre
			}
		}
		/*tmp je output matica otočená o 90 stupňov*/
		ComplexMatrix tmp = new ComplexMatrix(output.getNcol(), output.getNrow());
		for(int y = 0;y<output.getNrow();y++){
			for(int x=0;x<output.getNcol();x++){
				tmp.setElement(x, y, output.getElementCopy(y, x));//robí to to,čo chcem
			}
		}
		Complex[] tmpRow = new Complex[tmp.getNrow()];
		//výpočet FFT nad stĺpcami
		for(int y = 0;y<tmp.getNcol();y++){
			for(int x = 0;x<tmp.getNrow();x++){
				tmpRow[x] = tmp.getElementCopy(x,y);
			}
			FourierTransform fourierHorizontal = new FourierTransform(tmpRow);
			fourierHorizontal.transform();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal.getTransformedDataAsComplex());
			for(int x = 0;x<tmp.getNrow();x++){
				tmp.setElement(x, y, rowMatrix.getElementCopy(0, x));
			}
		}
		//priradenie do output
		for(int y=0;y<output.getNcol();y++){
			for(int x=0;x<output.getNrow();x++){
				output.setElement(x, y, tmp.getElementCopy(y, x));
			}
		}
		return output;
	}

	public static double[][] InverseFourierTransform2D(ComplexMatrix input){
		if(input == null) return null;
		/*Fourierova transformácia nad každým riadkom a následne priradenie hodnôt do matice output
		 */
		ComplexMatrix rowMatrix;
		Complex[] tmpRow = new Complex[input.getNrow()];
		for(int y = 0;y<input.getNcol();y++){
			for(int x = 0;x<input.getNrow();x++){
				tmpRow[x] = input.getElementCopy(x,y);
			}
			FourierTransform fourierHorizontal = new FourierTransform(tmpRow);
			fourierHorizontal.inverse();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal.getTransformedDataAsComplex());
			for(int x = 0;x<rowMatrix.getNrow();x++){
				input.setElement(x, y, rowMatrix.getElementCopy(0, x));
			}
		}
		/*tmp je output matica otočená o 90 stupňov*/
		ComplexMatrix tmp = new ComplexMatrix(input.getNcol(), input.getNrow());
		for(int y = 0;y<input.getNcol();y++){
			for(int x=0;x<input.getNrow();x++){
				tmp.setElement(y, x, input.getElementCopy(x, y));
			}
		}
		tmpRow = new Complex[tmp.getNrow()];
		//výpočet FFT nad stĺpcami
		for(int y = 0;y<tmp.getNcol();y++){
			for(int x = 0;x<tmp.getNrow();x++){
				tmpRow[x] = tmp.getElementCopy(x,y);
			}
			FourierTransform fourierHorizontal=new FourierTransform(tmpRow);
			fourierHorizontal.inverse();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal.getTransformedDataAsComplex());
			for(int x = 0;x<rowMatrix.getNrow();x++){
				tmp.setElement(x, y, rowMatrix.getElementCopy(0, x));
			}
		}
		//priradenie do output
		double[][] output = new double[input.getNcol()][input.getNrow()];
		for(int y=0;y<output.length;y++){
			for(int x=0;x<output[y].length;x++){
				output[y][x] = tmp.getElementCopy(y, x).getReal();
			}
		}
		return output;
	}
	
}
