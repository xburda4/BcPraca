package plugins.mBurda.filters;

import java.awt.image.BufferedImage;

import Jama.Matrix;
import flanagan.complex.Complex;
import flanagan.complex.ComplexMatrix;
import flanagan.math.FourierTransform;
import icy.gui.dialog.MessageDialog;
import icy.image.IcyBufferedImage;

public class Computations {

	public static double[][] getGaborKernel2D(int width, int height,
			double wavelength, double angle) {
		return getGaborKernel2D(width, height, wavelength, angle, 8.5, 0.75);
	}

	public static double[][] getGaborKernel2D(int width, int height,
		double wavelength, double angle, double sigmaOnf, double thetaSigma) {
		angle = angle * Math.PI / 180;
		double centerFreq = 1.0 / wavelength;
		double radialComp;
		double[][] kernel2D = new double[height][width];
		double theta, ds, dc, angulComp, dtheta;
		long time = System.nanoTime();
		for (int x = -width / 2; x < width / 2; x++) {
			for (int y = -height / 2; y < height / 2; y++) {
				double position = Math.sqrt(x * x + y * y);
				radialComp = Math
						.exp((-(Math.log(position / centerFreq) * (Math
								.log(position / centerFreq))) / (2 * ((Math
								.log(sigmaOnf) * Math.log(sigmaOnf))))));
				kernel2D[y + height / 2][x + width / 2] = radialComp;

				theta = Math.atan2(-y, x);
				ds = Math.sin(theta) * Math.cos(angle) - Math.cos(theta)
						* Math.sin(angle);
				dc = Math.cos(theta) * Math.cos(angle) + Math.sin(theta)
						* Math.sin(angle);
				dtheta = Math.abs(Math.atan2(ds, dc));
				angulComp = Math.exp(-(dtheta * dtheta)
						/ (2 * (thetaSigma * thetaSigma)));
				kernel2D[y + height / 2][x + width / 2] = kernel2D[y + height
						/ 2][x + width / 2]
						* angulComp;
			}
		}
		time = System.nanoTime() - time;
		System.out.println("Kernel time is " + time + " ns");
		return kernel2D;
	}

	public static double getLogGaborKernelPoint(int x, int y, double scale,
			double orientation) {
		double pointValue = Math
				.exp((-(Math.log(Math.sqrt(x * x + y * y) / (1/scale)) * (Math
						.log(Math.sqrt(x * x + y * y) / (1/scale)))) / (2 * ((Math
						.log(8.5) * Math.log(8.5))))));

		double theta = Math.atan2(-y, x);
		//orientation = orientation;
		double ds = Math.sin(theta) * Math.cos(orientation) - Math.cos(theta)
				* Math.sin(orientation);
		double dc = Math.cos(theta) * Math.cos(orientation) + Math.sin(theta)
				* Math.sin(orientation);
		double dtheta = Math.abs(Math.atan2(ds, dc));
		pointValue *= Math.exp(-(dtheta * dtheta) / (2 * (0.75 * 0.75)));
		return pointValue;
	}

	public static double[][][] multiFTKernel(ComplexMatrix ftImg,double scale,double angle) {
		double[][][] output = new double[2][ftImg.getNrow()][ftImg.getNcol()];
		double logGabPoint;
		for(int y=0;y<ftImg.getNrow();y++){
			for(int x=0;x<ftImg.getNcol();x++){
				logGabPoint = getLogGaborKernelPoint(x - ftImg.getNcol()/2,y - ftImg.getNrow()/2, scale, angle);
				//možné použitie fcie getLogGaborKernel() pred cyklom
				output[0][y][x] = logGabPoint * ftImg.getElementCopy(y, x).getReal();
				output[1][y][x] = /*logGabPoint * */ ftImg.getElementCopy(y, x).getImag();
			}
		}
		return output;
	}
	
	public static Matrix[][] getPhaseCong(ComplexMatrix ftImg, double[] scales,
			double[] orientations, double threshold, double cutoffValue,
			double gainFactor) {
		if(ftImg == null || scales == null || orientations == null) {
			MessageDialog.showDialog("Incorrect argument");
			return null;
		}

		double[][][][][] componentSO = new double[orientations.length][scales.length][2][ftImg
				.getNrow()][ftImg.getNcol()];
		double temp0,temp1;
		for (int orients = 0; orients < orientations.length; orients++) {
			for (int scs = 0; scs < scales.length; scs++) {
				componentSO[orients][scs] = multiFTKernel(ftImg, scales[scs], orientations[orients]);
				componentSO[orients][scs] = InverseFourierTransform2D(createComplexMatrix(componentSO[orients][scs]), false);
				for(int y=0;y<componentSO[orients][scs][0].length;y++){
					for(int x=0;x<componentSO[orients][scs][0][y].length;x++){
						temp0 = componentSO[orients][scs][0][y][x];
						temp1 = componentSO[orients][scs][1][y][x];
						componentSO[orients][scs][0][y][x] = Math.sqrt(temp0*temp0+temp1*temp1);
						componentSO[orients][scs][1][y][x] = Math.atan2(temp1, temp0);
						}
					}
			}
		}
		
		double denominator,numerator,weight,phaseDevMeasure,Amax=0,sumOfAmplis,meanPhase,tmp;
		Matrix[][] phaseCongMatrix = new Matrix[ftImg.getNrow()][ftImg.getNcol()];
		for (int y = 0; y < ftImg.getNrow(); y++) {
			for (int x = 0; x < ftImg.getNcol(); x++) {
				//phaseCongMatrix[y][x] = new Matrix(2,2);
				tmp = 0;
				for (int orients = 0; orients < orientations.length; orients++) {
					
					sumOfAmplis = 0; meanPhase = 0; numerator = 0; 					
					for (int scs = 0; scs < scales.length; scs++) {
						sumOfAmplis+=componentSO[orients][scs][0][y][x];
						if(scs == 0) Amax = componentSO[orients][scs][0][y][x];
						else Amax = Math.max(Amax,componentSO[orients][scs][0][y][x]);
						meanPhase += componentSO[orients][scs][1][y][x];
					}
					meanPhase /= scales.length;
					denominator = sumOfAmplis;

					for (int scs = 0; scs < scales.length; scs++) {
						phaseDevMeasure = Math.cos(componentSO[orients][scs][1][y][x]-meanPhase)-
								Math.abs(Math.sin(componentSO[orients][scs][1][y][x])-meanPhase);
						tmp = componentSO[orients][scs][0][y][x]*phaseDevMeasure - threshold;
						
						weight = 1 + Math.exp(gainFactor*(cutoffValue-(1/scales.length)*(sumOfAmplis/(Amax+0.000000000001))));
						tmp *= weight;
						numerator += tmp;
					}
					denominator += 0.0000000000000000001;
					
					Matrix mat = new Matrix(2,2);
					mat.set(0, 0, (Math.cos(orients)*Math.cos(orients)+1/2));
					mat.set(0, 1, Math.cos(orients)*Math.sin(orients));
					mat.set(1, 0, mat.get(0, 1));
					mat.set(1, 1, (Math.sin(orients)*Math.sin(orients))+1/2);
					
					if(phaseCongMatrix[y][x] == null) phaseCongMatrix[y][x] = mat.times(numerator/denominator);
					phaseCongMatrix[y][x].plus(mat.times(numerator/denominator));
					//phaseCongMatrix[y][x] += (numerator/denominator);
				}
			}
		}
		return phaseCongMatrix;
	}
	
	
	
	
	
//	public static double[][] getPhaseCong(ComplexMatrix ftImg, double[] scales,
//			double[] orientations, double threshold, double cutoffValue,
//			double gainFactor) {
//		if(ftImg == null || scales == null || orientations == null) {
//			MessageDialog.showDialog("Incorrect argument");
//			return null;
//		}
//		double[][][][][] componentSO = new double[orientations.length][scales.length][2][ftImg
//				.getNrow()][ftImg.getNcol()];
//		
//		for (int orients = 0; orients < orientations.length; orients++) {
//			for (int scs = 0; scs < scales.length; scs++) {
//				componentSO[orients][scs] = multiFTKernel(ftImg, scales[scs], orientations[orients]);
//			}
//		}
//		
//		double denominator,numerator,weight,phaseDevMeasure,Amax=0,sumOfAmplis,meanPhase;
//		double[][] phaseCongMatrix = new double[ftImg.getNrow()][ftImg.getNcol()];
//		for (int y = 0; y < ftImg.getNrow(); y++) {
//			for (int x = 0; x < ftImg.getNcol(); x++) {
//				phaseCongMatrix[y][x] = 0;
//				for (int orients = 0; orients < orientations.length; orients++) {
//					
//					sumOfAmplis = 0; meanPhase = 0; numerator = 0; denominator = 0;
//					
//					for (int scs = 0; scs < scales.length; scs++) {
//						sumOfAmplis+=componentSO[orients][scs][0][y][x];
//						if(scs == 0) Amax = componentSO[orients][scs][0][y][x];
//						else Amax = Math.max(Amax,componentSO[orients][scs][0][y][x]);
//						meanPhase += componentSO[orients][scs][1][y][x];
//					}
//					meanPhase /= scales.length;
//					
//					for (int scs = 0; scs < scales.length; scs++) {
//						if(componentSO[orients][scs][0][y][x] == 0){
//							numerator = 0;
//							continue;
//						}
//						phaseDevMeasure = Math.cos(componentSO[orients][scs][1][y][x]-meanPhase)-
//								Math.abs(Math.sin(componentSO[orients][scs][1][y][x])-meanPhase);
//						numerator = componentSO[orients][scs][0][y][x]*phaseDevMeasure - threshold;
//						if(numerator < 0){
//							numerator = 0;
//							continue;
//						}
//						
//						weight = 1 + Math.exp(gainFactor*(cutoffValue-(1/scales.length)*(sumOfAmplis/(Amax+0.000000000001))));
//						numerator *= weight;
//					}
//					denominator += 0.000000000001;
//					phaseCongMatrix[y][x] += (numerator/denominator);
//				}
//			}
//		}
//		return phaseCongMatrix;
//	}

	public static ComplexMatrix FourierTransform2D(IcyBufferedImage img,
			boolean isAlt){
		double[][] matrix = new double[img.getHeight()][img.getWidth()];
		for (int y = 0; y < img.getHeight(); y++) {
			for (int x = 0; x < img.getWidth(); x++) {
				matrix[y][x] = img.getRaster().getSample(x, y, 0);
			}
		}
		return FourierTransform2D(matrix, isAlt);
	}
	
	public static ComplexMatrix FourierTransform2D(BufferedImage img,
			boolean isAlt){
		double[][] matrix = new double[img.getHeight()][img.getWidth()];
		for (int y = 0; y < img.getHeight(); y++) {
			for (int x = 0; x < img.getWidth(); x++) {
				matrix[y][x] = img.getRaster().getSample(x, y, 0);
			}
		}
		return FourierTransform2D(matrix, isAlt);
	}

	public static ComplexMatrix FourierTransform2D(double[][] input,
			boolean isAlt) {
		if (input == null)
			return null;
		int height = 1, width = 1;
		double real, img;
		while (height < input.length) {
			height *= 2;
		}
		while (width < input[0].length) {
			width *= 2;
		}
		ComplexMatrix output = new ComplexMatrix(height, width);
		/*
		 * Fourierova transformácia nad každým riadkom a následne priradenie
		 * hodnôt do matice output
		 */
		long time = System.nanoTime();
		ComplexMatrix rowMatrix;
		Complex[] tmpRow;
		for (int y = 0; y < output.getNrow(); y++) {
			tmpRow = new Complex[width];
			for (int x = 0; x < tmpRow.length; x++) {
				if (y < input.length && x < input[y].length) {
					tmpRow[x] = new Complex(input[y][x]);
				} else {
					tmpRow[x] = new Complex(0);
				}
			}
			FourierTransform fourierHorizontal = new FourierTransform(tmpRow);
			fourierHorizontal.transform();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal
					.getTransformedDataAsComplex());
			for (int x = 0; x < output.getNcol(); x++) {
				output.setElement(y, (x+output.getNcol()/2)%output.getNcol(), rowMatrix.getElementCopy(0, x));
			}
		}
		tmpRow = new Complex[output.getNrow()];
		// výpočet FFT nad stĺpcami
		for (int x = 0; x < output.getNcol(); x++) {
			for (int y = 0; y < output.getNrow(); y++) {
				tmpRow[y] = output.getElementCopy(y, x);
			}
			FourierTransform fourierHorizontal = new FourierTransform(tmpRow);
			fourierHorizontal.transform();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal
					.getTransformedDataAsComplex());
			for (int y = 0; y < output.getNrow(); y++) {
				output.setElement(y, x, rowMatrix.getElementCopy(0, (y+rowMatrix.getNcol()/2)%rowMatrix.getNcol()));
				if (isAlt) {
					real = output.getElementCopy(y, x).getReal();
					img = output.getElementCopy(y, x).getImag();
					output.setElement(y, x, Math.sqrt(real * real + img * img),
							Math.atan(img / real));
				}
			}
		}
		time = System.nanoTime() - time;
		System.out.println("Fourier transform time is " + time + " ns");
		return output;
	}

	public static double[][] fttToDoubleArr2D(ComplexMatrix input,boolean isAlt){
		if (input == null)
			return null;
		double[][] output = new double[input.getNrow()][input.getNcol()];
		for (int y = 0; y < input.getNrow(); y++) {
			for (int x = 0; x < input.getNcol(); x++) {
				if(isAlt)
				output[y][x] = input.getElementCopy(y, x).getReal();
				else output[y][x] = input.getElementCopy(y, x).getImag();
			}
		}
		return output;
	}
	
	public static double[][][] InverseFourierTransform2D(ComplexMatrix input,boolean isAlt) {
		if (input == null)
			return null;
		
		/*
		 * Fourierova transformácia nad každým riadkom a následne priradenie
		 * hodnôt do matice output
		 */
		double[][][] output = new double[2][input.getNrow()][input.getNcol()];
		
		if(isAlt){
			for (int y = 0; y < input.getNrow(); y++) {
				for (int x = 0; x < input.getNcol(); x++) {
					output[0][(y+input.getNrow()/2)%input.getNrow()][(x+input.getNcol()/2)%input.getNcol()] = input.getElementCopy(y, x).getReal()*Math.exp(input.getElementCopy(y, x).getImag());
				}	
			}
		}
		
		
		ComplexMatrix rowMatrix;
		Complex[] tmpRow = new Complex[input.getNcol()];
		for (int y = 0; y < input.getNrow(); y++) {
			for (int x = 0; x < input.getNcol(); x++) {
				tmpRow[x] = input.getElementCopy(y, (x+input.getNcol()/2)%input.getNcol());
			}
			FourierTransform fourierHorizontal = new FourierTransform(tmpRow);
			fourierHorizontal.inverse();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal
					.getTransformedDataAsComplex());
			for (int x = 0; x < rowMatrix.getNcol(); x++) {
				input.setElement(y, x, rowMatrix.getElementCopy(0, x));
			}
		}
		
		tmpRow = new Complex[input.getNrow()];
		// výpočet FFT nad stĺpcami
		for (int x = 0; x < input.getNcol(); x++) {
			for (int y = 0; y < input.getNrow(); y++) {
				tmpRow[y] = input.getElementCopy((y+input.getNrow()/2)%input.getNrow(), x);
			}
			FourierTransform fourierHorizontal = new FourierTransform(tmpRow);
			fourierHorizontal.inverse();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal
					.getTransformedDataAsComplex());
			for (int y = 0; y < input.getNrow(); y++) {
				input.setElement(y, x, rowMatrix.getElementCopy(0, y));

			}
		}
		// priradenie do output

		for (int y = 0; y < input.getNrow(); y++) {
			for (int x = 0; x < input.getNcol(); x++) {
				output[0][y][x] = input.getElementCopy(y, x).getReal();
				output[1][y][x] = input.getElementCopy(y, x).getImag();
			}
		}
		return output;
	}

	
	public static ComplexMatrix createComplexMatrix(double[][][] input){
		ComplexMatrix output = new ComplexMatrix(input[0].length,input[0][0].length);
		for(int y=0;y<input[0].length;y++){
			for(int x=0;x<input[0][y].length;x++){
				output.setElement(y, x, input[0][y][x], input[1][y][x]);
			}
		}
		return output;
	}
}
