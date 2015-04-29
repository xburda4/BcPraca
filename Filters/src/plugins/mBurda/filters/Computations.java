package plugins.mBurda.filters;

//import org.jdesktop.swingx.image.GaussianBlurFilter;

import java.awt.image.BufferedImage;
import java.util.ArrayList;

import flanagan.complex.Complex;
import flanagan.complex.ComplexMatrix;
import flanagan.math.FourierTransform;

public class Computations {
	
	public static double[][] getGaborKernel(int width,int height,double wavelength,double angle){
		return getGaborKernel(width,height, wavelength,angle,8.5,0.75);
	}
	
	public static double[][] getGaborKernel(int width,int height,double wavelength,double angle,double sigmaOnf,double thetaSigma){
		angle = angle*Math.PI/180;
		double centerFreq = 1.0/wavelength;
		double radialComp;
		double[][] kernel2D = new double[height][width];
		double theta,ds,dc,angulComp,dtheta;
		for(int x=-width/2;x<width/2;x++){
			for(int y=-height/2;y<height/2;y++){
				double position = Math.sqrt(x*x+y*y);
				radialComp = Math.exp((-(Math.log(position/centerFreq)*(Math.log(position/centerFreq)))/(2*((Math.log(sigmaOnf)*Math.log(sigmaOnf))))));
				kernel2D[y+height/2][x+width/2] = radialComp;
				
				theta = Math.atan2(-y, x);
				ds = Math.sin(theta)*Math.cos(angle)-Math.cos(theta)*Math.sin(angle);
				dc = Math.cos(theta)*Math.cos(angle)+Math.sin(theta)*Math.sin(angle);
				dtheta = Math.abs(Math.atan2(ds, dc));
				angulComp = Math.exp(-(dtheta*dtheta)/(2*(thetaSigma*thetaSigma)));
				kernel2D[y+height/2][x+width/2] = kernel2D[y+height/2][x+width/2] * angulComp;
			}
		}
		return kernel2D;
	}
	
//	public static BufferedImage multFreqDomain2D(BufferedImage img,double angleStep,int numberOfAngleSteps,double scaleStep,int numberOfScaleSteps){
//		ComplexMatrix ft = FourierTransform2D(img);
//		ArrayList list = new ArrayList<double[][]>();
//		for(int i = 0;i<numberOfAngleSteps;i++){
//			for(int j=0;j<numberOfScaleSteps;j++){
//				double[][] gk = getGaborKernel(ft.getNcol(), ft.getNrow(),0.08,0,1.2,0.25);
//				for(int y= 0;y<ft.getNrow();y++){
//					for(int x=0;x<ft.getNcol();x++){
//						gk[y][x] = gk[y][x]*ft.get
//					}
//				}
//			}
//		}
//		
//		return null;
//	}
	
	public static ComplexMatrix FourierTransform2D(BufferedImage img){
		double[][] matrix = new double[img.getHeight()][img.getWidth()];
		for(int y=0;y<img.getHeight();y++){
			for(int x=0;x<img.getWidth();x++){
				matrix[y][x] = img.getRaster().getSample(x, y, 0);
			}
		}
		return FourierTransform2D(matrix);
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
