package plugins.mBurda.filters;

//import org.jdesktop.swingx.image.GaussianBlurFilter;

import flanagan.complex.Complex;
import flanagan.complex.ComplexMatrix;
import flanagan.math.FourierTransform;

public class Computations {
	
	public static double[][] getGaborKernel(int radius,double wavelength,double angle){
		return getGaborKernel(radius, wavelength,angle,8.5,0.75);
	}
	
	public static double[][] getGaborKernel(int radius,double wavelength,double angle,double sigmaOnf,double thetaSigma){
		angle = angle*Math.PI/180;
		double centerFreq = 1.0/wavelength;
		double radialComp;
		double[][] kernel2D = new double[2*radius+1][2*radius+1];
		double theta,ds,dc,angulComp,dtheta;
		for(int x=-radius;x<radius;x++){
			for(int y=-radius;y<radius;y++){
				double position = Math.sqrt(x*x+y*y);
				radialComp = Math.exp((-(Math.log(position/centerFreq)*(Math.log(position/centerFreq)))/(2*((Math.log(sigmaOnf)*Math.log(sigmaOnf))))));
				kernel2D[y+radius][x+radius] = radialComp;
				
				theta = Math.atan2(-y, x);
				ds = Math.sin(theta)*Math.cos(angle)-Math.cos(theta)*Math.sin(angle);
				dc = Math.cos(theta)*Math.cos(angle)+Math.sin(theta)*Math.sin(angle);
				dtheta = Math.abs(Math.atan2(ds, dc));
				angulComp = Math.exp(-(dtheta*dtheta)/(2*(thetaSigma*thetaSigma)));
				kernel2D[y+radius][x+radius] = kernel2D[y+radius][x+radius] * angulComp;
			}
		}
		return kernel2D;
	}
	
	public static ComplexMatrix FourierTransform2D(double[][] input){
		if(input == null) return null;
		int height=1,width=1;
		while(height < input.length){
			height *= 2;
		}
		while(width < input[0].length){
			width *= 2;
		}
		ComplexMatrix output = new ComplexMatrix(height,width);
		/*Fourierova transformácia nad každým riadkom a následne priradenie hodnôt do matice output*/
		ComplexMatrix rowMatrix;
		Complex[] tmpRow;
		for(int y = 0;y<output.getNrow();y++){
			tmpRow = new Complex[width];
				for(int x = 0;x<tmpRow.length;x++){
					if(y<input.length && x<input[y].length){
						tmpRow[x] = new Complex(input[y][x]);
					} else {
						tmpRow[x] = new Complex(0);
					}
				}
				FourierTransform fourierHorizontal = new FourierTransform(tmpRow);
				fourierHorizontal.transform();
				rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal.getTransformedDataAsComplex());
				for(int x = 0;x<output.getNcol();x++){
					output.setElement(y, x, rowMatrix.getElementCopy(0, x));
					}
		}
		tmpRow = new Complex[output.getNrow()];
		//výpočet FFT nad stĺpcami
		for(int x = 0;x<output.getNcol();x++){
			for(int y = 0;y<output.getNrow();y++){
				tmpRow[y] = output.getElementCopy(y,x);
			}
			FourierTransform fourierHorizontal = new FourierTransform(tmpRow);
			fourierHorizontal.transform();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal.getTransformedDataAsComplex());
			for(int y = 0;y<output.getNrow();y++){
				output.setElement(y, x, rowMatrix.getElementCopy(0, y));
			}
		}
		return output;
	}

	public static double[][] InverseFourierTransform2D(ComplexMatrix input){
		if(input == null) return null;
		/*Fourierova transformácia nad každým riadkom a následne priradenie hodnôt do matice output
		 */
		ComplexMatrix rowMatrix;
		Complex[] tmpRow = new Complex[input.getNcol()];
		for(int y = 0;y<input.getNrow();y++){
			for(int x = 0;x<input.getNcol();x++){
				tmpRow[x] = input.getElementCopy(y,x);
			}
			FourierTransform fourierHorizontal = new FourierTransform(tmpRow);
			fourierHorizontal.inverse();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal.getTransformedDataAsComplex());
			for(int x = 0;x<rowMatrix.getNcol();x++){
				input.setElement(y, x, rowMatrix.getElementCopy(0, x));
			}
		}
		double[][] output = new double[input.getNrow()][input.getNcol()];
		tmpRow = new Complex[input.getNrow()];
		//výpočet FFT nad stĺpcami
		for(int x = 0;x<input.getNcol();x++){
			for(int y = 0;y<input.getNrow();y++){
				tmpRow[y] = input.getElementCopy(y, x);
			}
			FourierTransform fourierHorizontal=new FourierTransform(tmpRow);
			fourierHorizontal.inverse();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal.getTransformedDataAsComplex());
			for(int y = 0;y<input.getNrow();y++){
				input.setElement(y, x, rowMatrix.getElementCopy(0, y));
		
			}
		}
		//priradenie do output
		
		for(int y=0;y<input.getNrow();y++){
			for(int x=0;x<input.getNcol();x++){
				output[y][x] = input.getElementCopy(y, x).getReal();
			}
		}
		return output;
	}
	
}
