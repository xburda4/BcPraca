package plugins.mBurda.filters;

//import org.jdesktop.swingx.image.GaussianBlurFilter;

import icy.gui.dialog.MessageDialog;

import java.awt.image.BufferedImage;

import flanagan.complex.Complex;
import flanagan.complex.ComplexMatrix;
import flanagan.math.FourierTransform;

public class Computations {

	public static double[][] getGaborKernel(int width, int height,
			double wavelength, double angle) {
		return getGaborKernel(width, height, wavelength, angle, 8.5, 0.75);
	}

	public static double[][] getGaborKernel(int width, int height,
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
				.exp((-(Math.log(Math.sqrt(x * x + y * y) / scale) * (Math
						.log(Math.sqrt(x * x + y * y) / scale))) / (2 * ((Math
						.log(8.5) * Math.log(8.5))))));

		double theta = Math.atan2(-y, x);
		orientation = orientation * Math.PI / 180;
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
				logGabPoint = getLogGaborKernelPoint(x, y, scale, angle);
				//možné použitie fcie getLogGaborKernel() pred cyklom
				output[0][y][x] = logGabPoint * ftImg.getElementCopy(y, x).getReal();
				output[1][y][x] = logGabPoint * ftImg.getElementCopy(y, x).getImag();
			}
		}
		return output;
	}
	
	public static double getPhaseCong(ComplexMatrix ftImg, double[] scales,
			double[] orientations, double threshold, double cutoffValue,
			double gainFactor) {
		if(ftImg == null || scales == null || orientations == null) {
			MessageDialog.showDialog("Incorrect argument");
			return 0;
		}
		double[][][][][] componentSO = new double[orientations.length][scales.length][2][ftImg
				.getNrow()][ftImg.getNcol()];
		
		for (int orients = 0; orients < orientations.length; orients++) {
			for (int scs = 0; scs < scales.length; scs++) {
				componentSO[orients][scs] = multiFTKernel(ftImg, scales[scs], orientations[orients]);
			}
		}
		
		double denominator,numerator,weight,phaseDevMeasure,Amax=0,sumOfAmplis,meanPhase;
		double[][] phaseCongMatrix = new double[ftImg.getNrow()][ftImg.getNcol()];
		for (int y = 0; y < ftImg.getNrow(); y++) {
			for (int x = 0; x < ftImg.getNcol(); x++) {
				phaseCongMatrix[y][x] = 0;
				for (int orients = 0; orients < orientations.length; orients++) {
					
					sumOfAmplis = 0; meanPhase = 0; numerator = 0; denominator = 0;
					
					for (int scs = 0; scs < scales.length; scs++) {
						sumOfAmplis+=componentSO[orients][scs][0][y][x];
						if(scs == 0) Amax = componentSO[orients][scs][0][y][x];
						else Amax = Math.max(Amax,componentSO[orients][scs][0][y][x]);
					}
					
					for (int scs = 0; scs < scales.length; scs++) {
						if(componentSO[orients][scs][0][y][x] == 0){
							numerator = 0;
							continue;
						}
						phaseDevMeasure = Math.cos(componentSO[orients][scs][1][y][x]-meanPhase)-
								Math.abs(Math.sin(componentSO[orients][scs][1][y][x])-meanPhase);
						numerator = componentSO[orients][scs][0][y][x]*phaseDevMeasure - threshold;
						if(numerator < 0){
							numerator = 0;
							continue;
						}
						
						weight = 1 + Math.exp(gainFactor*(cutoffValue-(1/scales.length)*(sumOfAmplis/(Amax+0.000000000001))));
						numerator *= weight;
					}
					denominator += 0.000000000001;
					phaseCongMatrix[y][x] += (numerator/denominator);
				}
			}
		}
		return 0;
	}

	// public static double[][] getPhaseCongMatrix(ComplexMatrix ftImg,double[]
	// scales,double[] orientations,double threshold,double cutoffValue,double
	// gainFactor){
	// double[][] output = new double[ftImg.getNrow()][ftImg.getNcol()];
	// double phaseCong = 0,amplitude,deltaPhi,sumOfAmplis = 0;
	// for(int y=0;y<ftImg.getNrow();y++){
	// for(int x=0;x<ftImg.getNcol();x++){
	// for(int o=0;o<orientations.length;o++){
	//
	// for(int s=0;s<scales.length;s++){
	// amplitude = getLogGaborKernelPoint(x, y, scales[s], orientations[o]) *
	// ftImg.getElementCopy(y, x).getReal();
	// deltaPhi = Math.cos(ftImg.getElementCopy(y,
	// x).getImag()-TBD)-Math.sin(ftImg.getElementCopy(y, x).getImag()-TBD);
	// if(amplitude*deltaPhi-threshold<0) continue;
	// for(int i=0;i<scales.length;i++){
	// sumOfAmplis += ;
	// }
	// }
	// output[y][x] += phaseCong;
	// }
	// }
	// }
	//
	// return null;
	// }

	public static ComplexMatrix FourierTransform2D(BufferedImage img,
			boolean isAlt) {
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
				output.setElement(y, x, rowMatrix.getElementCopy(0, x));
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
				output.setElement(y, x, rowMatrix.getElementCopy(0, y));
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

	public static double[][] InverseFourierTransform2D(ComplexMatrix input) {
		if (input == null)
			return null;
		/*
		 * Fourierova transformácia nad každým riadkom a následne priradenie
		 * hodnôt do matice output
		 */
		ComplexMatrix rowMatrix;
		Complex[] tmpRow = new Complex[input.getNcol()];
		for (int y = 0; y < input.getNrow(); y++) {
			for (int x = 0; x < input.getNcol(); x++) {
				tmpRow[x] = input.getElementCopy(y, x);
			}
			FourierTransform fourierHorizontal = new FourierTransform(tmpRow);
			fourierHorizontal.inverse();
			rowMatrix = ComplexMatrix.rowMatrix(fourierHorizontal
					.getTransformedDataAsComplex());
			for (int x = 0; x < rowMatrix.getNcol(); x++) {
				input.setElement(y, x, rowMatrix.getElementCopy(0, x));
			}
		}
		double[][] output = new double[input.getNrow()][input.getNcol()];
		tmpRow = new Complex[input.getNrow()];
		// výpočet FFT nad stĺpcami
		for (int x = 0; x < input.getNcol(); x++) {
			for (int y = 0; y < input.getNrow(); y++) {
				tmpRow[y] = input.getElementCopy(y, x);
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
				output[y][x] = input.getElementCopy(y, x).getReal();
			}
		}
		return output;
	}

}
