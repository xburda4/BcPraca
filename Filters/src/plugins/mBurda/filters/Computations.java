package plugins.mBurda.filters;

import java.awt.image.BufferedImage;

import Jama.Matrix;
import flanagan.complex.Complex;
import flanagan.complex.ComplexMatrix;
import flanagan.math.FourierTransform;
import icy.gui.dialog.MessageDialog;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;

public class Computations {		
	/**
	 * Kontrola -> Výpočet fft
	 */
	public static Complex[] computeFFT(double[] input){
		FourierTransform four = new FourierTransform(input);
		four.transform();
		System.out.println("Hodnoty Fourierovej transformácie");
		for(int x = 0;x<four.getTransformedDataAsComplex().length;x++){
			System.out.print(four.getTransformedDataAsComplex()[x].getReal());
			System.out.println(" + "+four.getTransformedDataAsComplex()[x].getImag()+"i");
		}
		return four.getTransformedDataAsComplex();
	}
	public static void computeIFFT(Complex[] input){
		Complex tmp;
		for(int x=0;x<input.length/2;x++){
			tmp = input[x];
			input[x] = input[input.length-x-1];
			input[input.length-x-1] = tmp;
		}
		FourierTransform four = new FourierTransform(input);
		four.inverse();
		System.out.println("Hodnoty Inverznej FFT");
		for(int x = 0;x<four.getTransformedDataAsComplex().length;x++){
			System.out.print(four.getTransformedDataAsComplex()[x].getReal());
			System.out.println(" + "+four.getTransformedDataAsComplex()[x].getImag()+"i");
		}
	}
	
	public static double[][] multiFTKernel0(ComplexMatrix ftImg,double scale,double angle,int numberOfOrients) {
		double[][] output = new double[ftImg.getNrow()][ftImg.getNcol()];
//		double logGabPoint;
		
		final double[][] kernel = getGaborKernel2D(ftImg.getNcol(), ftImg.getNrow(), scale, angle,numberOfOrients);
		
		double tmp0=0,tmp1=0;
		final double[][] lp = lowPassFilter(ftImg.getNcol(),ftImg.getNrow());
		for(int y=0;y<ftImg.getNrow();y++){
			for(int x=0;x<ftImg.getNcol();x++){
				
//				logGabPoint = kernel[/*ftImg.getNrow()-1-*/y][/*ftImg.getNcol()-1-*/x];
				//logGabPoint = getLogGaborKernelPoint(x - ftImg.getNcol()/2,y - ftImg.getNrow()/2, scale, angle,0.75);
				//možné použitie fcie getLogGaborKernel() pred cyklom
				
				tmp0 = ftImg.getElementCopy(y, x).getReal();
				tmp1 = ftImg.getElementCopy(y, x).getImag();
				
				output[y][x] = kernel[ftImg.getNrow()-1-y][ftImg.getNcol()-1-x] * lp[y][x];
			}
		}
		return output;
	}
	
	public static double[][][] multiFTKernel1(ComplexMatrix ftImg,double scale,double angle,int numberOfOrients) {
		double[][][] output = new double[2][ftImg.getNrow()][ftImg.getNcol()];
//		double logGabPoint;
		
//		final double[][] kernel = getGaborKernel2D(ftImg.getNcol(), ftImg.getNrow(), scale, angle,numberOfOrients);
		
		
//		final double[][] lp = lowPassFilter(ftImg.getNcol(),ftImg.getNrow());
		for(int y=0;y<ftImg.getNrow();y++){
			for(int x=0;x<ftImg.getNcol();x++){
				
//				logGabPoint = kernel[/*ftImg.getNrow()-1-*/y][/*ftImg.getNcol()-1-*/x];
				//logGabPoint = getLogGaborKernelPoint(x - ftImg.getNcol()/2,y - ftImg.getNrow()/2, scale, angle,0.75);
				//možné použitie fcie getLogGaborKernel() pred cyklom
				
				
				output[0][y][x] = ftImg.getElementCopy(y, x).getReal();// * kernel[ftImg.getNrow()-1-y][ftImg.getNcol()-1-x] * lp[y][x];
				output[1][y][x] = ftImg.getElementCopy(y, x).getImag();// * kernel[ftImg.getNrow()-1-y][ftImg.getNcol()-1-x] * lp[y][x];
			}
		}
		return output;
	}
	
	public static double[][][] multiFTKernel2(ComplexMatrix ftImg,double scale,double angle,int numberOfOrients) {
		double[][][] output = new double[2][ftImg.getNrow()][ftImg.getNcol()];
//		double logGabPoint;
		
		final double[][] kernel = getGaborKernel2D(ftImg.getNcol(), ftImg.getNrow(), scale, angle,numberOfOrients);
		
		double tmp0=0,tmp1=0;
		final double[][] lp = lowPassFilter(ftImg.getNcol(),ftImg.getNrow());
		for(int y=0;y<ftImg.getNrow();y++){
			for(int x=0;x<ftImg.getNcol();x++){
				
//				logGabPoint = kernel[/*ftImg.getNrow()-1-*/y][/*ftImg.getNcol()-1-*/x];
				//logGabPoint = getLogGaborKernelPoint(x - ftImg.getNcol()/2,y - ftImg.getNrow()/2, scale, angle,0.75);
				//možné použitie fcie getLogGaborKernel() pred cyklom
				
				output[0][y][x] = ftImg.getElementCopy(y, x).getReal() * kernel[ftImg.getNrow()-1-y][ftImg.getNcol()-1-x] * lp[y][x];
				output[1][y][x] = ftImg.getElementCopy(y, x).getImag() * kernel[ftImg.getNrow()-1-y][ftImg.getNcol()-1-x] * lp[y][x];
			}
		}
		return output;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 * Returns real part of logGabor kernel in frequency domain
	 * @param width of image
	 * @param height of image
	 * @param wavelength 1/frequency
	 * @param angle of signal
	 * @return image of LogGabor kernel in frequency domain as 2D double array
	 * */
	public static double[][] getGaborKernel2D(int width, int height,
			double wavelength, double angle,int numOfOrients) {
		return getGaborKernel2D(width, height, wavelength, angle, 0.75, numOfOrients);
	}

	public static double[][] lowPassFilter(int width,int height){
		double[][] lowpass = new double[height][width];
		for(int y= -(height/2);y<height/2;y++){
			for(int x = -(width/2);x<width/2;x++){
				//lowpass[y + (height/2)][x + (width/2)] = Math.sqrt(x*x+y*y) < radius ? 1 : 0 ;
				lowpass[y + (height/2)][x + (width/2)] = 1.0/ (1+Math.pow((Math.sqrt(x*x+y*y)/(0.45*width)),30));
			}
		}
		return lowpass;
	}
	
	public static double[][] highPassFilter(int width,int height,int radius){
		double[][] highpass = new double[height][width];
		for(int y= -(height/2);y<height/2;y++){
			for(int x = -(width/2);x<width/2;x++){
				highpass[y + (height/2)][x + (width/2)] = Math.sqrt(x*x+y*y) > radius ? 1 : 0 ;
			}
		}
		return highpass;
	}
	
	public static double[][] bandwidthPassFilter(int width,int height,int min,int max){
		double[][] highpass = new double[height][width];
		for(int y= -(height/2);y<height/2;y++){
			for(int x = -(width/2);x<width/2;x++){
				highpass[y + (height/2)][x + (width/2)] = (Math.sqrt(x*x+y*y)< max)&&(Math.sqrt(x*x+y*y) > min) ? 1 : 0 ;
			}
		}
		return highpass;
	}
	
	/**
	 * Returns  real part of logGabor kernel in frequency domain.
	 * @param width of image
	 * @param height of image
	 * @param wavelength 1/frequency
	 * @param angle of signal,in radians
	 * @param sigmaOnf determines scale of centre frequency
	 * @param thetaSigma standard deviation of gauss with respect to angle
	 * @return image of LogGabor kernel in frequency domain as 2D double array
	 * */
	public static double[][] getGaborKernel2D(int width, int height,
		double wavelength, double angle, double sigmaOnf, int numberOfOrients) {
//		angle = angle * Math.PI / 180;
		double centerFreq = 1.0 / wavelength;
		double radialComp;
		double[][] kernel2D = new double[height][width];
		double theta, ds, dc, angulComp, dtheta;
//		long time = System.nanoTime();
		for (int x = -(width / 2); x < (width / 2); x++) {
			for (int y = -(height / 2); y < (height / 2); y++) {
				double position;
				if(x != 0 || y != 0)position = Math.sqrt((double)(x)/width * (double)(x)/width + (double)(y)/height * (double)(y)/height);
				//else position = 1;
				else {
					kernel2D[y][x] = 0;
					continue;
					}
				radialComp = Math
						.exp(-(Math.log(position / centerFreq) * (Math
								.log(position / centerFreq))) / (2 * ((Math
								.log(sigmaOnf) * Math.log(sigmaOnf)))));

				theta = Math.atan2(-(double)(y)/height, (double)(x)/width);
				ds = Math.sin(theta) * Math.cos(angle) - Math.cos(theta)
						* Math.sin(angle);
				dc = Math.cos(theta) * Math.cos(angle) + Math.sin(theta)
						* Math.sin(angle);
				dtheta = Math.abs(Math.atan2(ds, dc));
//				angulComp = Math.exp(-(dtheta * dtheta)
//						/ (2 * (thetaSigma * thetaSigma)));
				
				dtheta = Math.min(dtheta*numberOfOrients/2.0,Math.PI);
				angulComp = (Math.cos(dtheta)+1)/2.0;
				kernel2D[y + height / 2][x + width / 2] = radialComp *
						 angulComp;
			}
		}
		
		
//		time = System.nanoTime() - time;
//		System.out.println("Kernel time is " + time + " ns");
		return kernel2D;
	}

	/**
	 * Returns one pixel of logGabor kernel
	 * @param x coordinate
	 * @param y coordinate
	 * @param scale of frequency of logGabor filter
	 * @param angle of signal 
	 * @return value of logGabor kernel at given coordinates
	 * */
	public static double getLogGaborKernelPoint(int x, int y, double scale,
			double orientation,double sigmaOnf) {
		double pointValue = Math
				.exp((-(Math.log(Math.sqrt(x * x + y * y) / (1/scale)) * (Math
						.log(Math.sqrt(x * x + y * y) / (1/scale)))) / (2 * ((Math
						.log(sigmaOnf) * Math.log(sigmaOnf))))));

		double theta = Math.atan2(-y, x);
		double ds = Math.sin(theta) * Math.cos(orientation) - Math.cos(theta)
				* Math.sin(orientation);
		double dc = Math.cos(theta) * Math.cos(orientation) + Math.sin(theta)
				* Math.sin(orientation);
		double dtheta = Math.abs(Math.atan2(ds, dc));
		pointValue *= Math.exp(-(dtheta * dtheta) / (2 * (0.75 * 0.75)));
		return pointValue;
	}

	/**
	 * Multiplies image in frequency domain and logGabor kernel in frequency domain
	 * @param ComplexMatrix image in frequency domain
	 * @param double scale of wavelength of logGabor kernel
	 * @param double angle of wave of logGabor kernel
	 * @return double[][][] real and imaginary component of every pixel in image
	 * */
	public static double[][][] multiFTKernel(ComplexMatrix ftImg,double scale,double angle,int numberOfOrients) {
		double[][][] output = new double[2][ftImg.getNrow()][ftImg.getNcol()];
//		double logGabPoint;
		
		final double[][] kernel = getGaborKernel2D(ftImg.getNcol(), ftImg.getNrow(), scale, angle,numberOfOrients);
		
		final double[][] lp = lowPassFilter(ftImg.getNcol(),ftImg.getNrow());
		for(int y=0;y<ftImg.getNrow();y++){
			for(int x=0;x<ftImg.getNcol();x++){
				
//				logGabPoint = kernel[/*ftImg.getNrow()-1-*/y][/*ftImg.getNcol()-1-*/x];
				//logGabPoint = getLogGaborKernelPoint(x - ftImg.getNcol()/2,y - ftImg.getNrow()/2, scale, angle,0.75);
				//možné použitie fcie getLogGaborKernel() pred cyklom
				
			
				output[0][y][x] = ftImg.getElementCopy(y, x).getReal() * kernel[ftImg.getNrow()-1-y][ftImg.getNcol()-1-x] * lp[y][x];
				output[1][y][x] = ftImg.getElementCopy(y, x).getImag() * kernel[ftImg.getNrow()-1-y][ftImg.getNcol()-1-x] * lp[y][x];
			}
		}
		return output;
	}
	
	/**
	 * Computes phase congruency with specified variables.
	 * @param ftImg input image in frequency domain
	 * @param scales scales of wavelengths of logGabor kernels
	 * @param orientations angles of waves of logGabor kernels
	 * @param threshold
	 * @param cutoffValue
	 * @param gainFactor 
	 * @return Returns matrix array of tensors
	 */
	public static Matrix[][] getPhaseCong(ComplexMatrix ftImg, double[] scales,
			double[] orientations, double threshold, double cutoffValue,
			double gainFactor) {
		if(ftImg == null || scales == null || orientations == null) {
			MessageDialog.showDialog("Incorrect argument");
			return null;
		}

		/*double[rotácie][škály][real/imag][y][x]
		 * */
		double[][][][][] componentSO = new double[orientations.length][scales.length][2][ftImg
				.getNrow()][ftImg.getNcol()];
		double eps = 0.000000000001;
		Filter.phaseValues = new double[ftImg.getNrow()][ftImg.getNcol()];
		
		double temp0,temp1;
		/*
		 * Test
		 * */
		
		double[][][][] absSO= new double[orientations.length][scales.length][ftImg.getNrow()][ftImg.getNcol()];
		
		for (int orients = 0; orients < orientations.length; orients++) {
			for (int scs = 0; scs < scales.length; scs++) {
				componentSO[orients][scs] = multiFTKernel(ftImg, scales[scs], orientations[orients],orientations.length);
				componentSO[orients][scs] = InverseFourierTransform2D(createComplexMatrix(componentSO[orients][scs]));
//				componentSO[orients][scs] = newInv(createComplexMatrix(componentSO[orients][scs]),false);
				//componentSO[orients][scs][1] = shiftArray(componentSO[orients][scs][1]);
				//componentSO[orients][scs][0] = shiftArray(componentSO[orients][scs][0]);
				for(int y=0;y<componentSO[orients][scs][0].length;y++){
					for(int x=0;x<componentSO[orients][scs][0][y].length;x++){
						temp0 = componentSO[orients][scs][0][y][x];
						temp1 = componentSO[orients][scs][1][y][x];
						
						
						absSO[orients][scs][y][x] = Math.sqrt((temp0*temp0)+(temp1*temp1));
//						componentSO[orients][scs][0][y][x] = Math.sqrt((temp0*temp0)+(temp1*temp1));
//						if(orients ==0 && scs ==0)Filter.phaseValues[y][x] = componentSO[0][0][0][y][x];
//						componentSO[orients][scs][1][y][x] = Math.atan2(temp1,temp0);
						}
					}
			}
		}
		
		
//		Filter.phaseValues = new double[ftImg.getNrow()][ftImg.getNcol()];
		
		double sumE,sumO;
		double denominator,numerator,weight,phaseDevMeasure,Amax=0,sumOfAmplis,tmp;
		Matrix[][] phaseCongMatrix = new Matrix[ftImg.getNrow()][ftImg.getNcol()];
		for (int y = 0; y < ftImg.getNrow(); y++) {
			for (int x = 0; x < ftImg.getNcol(); x++) {
				//phaseCongMatrix[y][x] = new Matrix(2,2);
				tmp = 0;
				for (int orients = 0; orients < orientations.length; orients++) {
					sumE = 0; sumO = 0;
					sumOfAmplis = 0; numerator = 0; 					
					for (int scs = 0; scs < scales.length; scs++) {
						sumOfAmplis += absSO[orients][scs][y][x];
						sumE += componentSO[orients][scs][0][y][x];
						sumO += componentSO[orients][scs][1][y][x];
						if(scs == 0) Amax = componentSO[orients][scs][0][y][x];
						else Amax = Math.max(Amax,componentSO[orients][scs][0][y][x]);
//						meanPhase += componentSO[orients][scs][1][y][x];
						
					}
					weight = 1 + Math.exp(gainFactor*(cutoffValue-((1/scales.length)*(sumOfAmplis/(Amax+eps)-1))));
					//meanPhase /= scales.length;
					denominator = sumOfAmplis;

					double XE = Math.sqrt(sumE*sumE+sumO*sumO)+ eps;
					double mE = sumE / XE;
					double mO = sumO / XE;
					phaseDevMeasure = 0;
					for (int scs = 0; scs < scales.length; scs++) {
						phaseDevMeasure += mE*componentSO[orients][scs][0][y][x] + mO*componentSO[orients][scs][1][y][x] - 
								Math.abs(componentSO[orients][scs][0][y][x]*mO - componentSO[orients][scs][1][y][x]*mE);
						
//						phaseDevMeasure = Math.cos(componentSO[orients][scs][1][y][x]-meanPhase)-
//								Math.abs(Math.sin(componentSO[orients][scs][1][y][x]-meanPhase));
						tmp = /*componentSO[orients][scs][0][y][x]**/phaseDevMeasure - threshold;
						if(tmp < 0) tmp = 0;
						numerator += (tmp/weight);
					}
					denominator += eps;
					
					Matrix mat = new Matrix(2,2);
					mat.set(0, 0, (Math.cos(orientations[orients])*Math.cos(orientations[orients])-1));
					mat.set(0, 1, Math.cos(orientations[orients])*Math.sin(orientations[orients]));
					mat.set(1, 0, mat.get(0, 1));
					mat.set(1, 1, (Math.sin(orientations[orients])*Math.sin(orientations[orients]))-1);
					
					if(phaseCongMatrix[y][x] == null) phaseCongMatrix[y][x] = mat.times(numerator/denominator);
					else phaseCongMatrix[y][x].plus(mat.times(numerator/denominator));
					
					Filter.phaseValues[y][x] += numerator/denominator;
					
				}
			}
		}
		return phaseCongMatrix;
	}
	/**
	 * Transforms image to frequency domain.
	 * @param img input image
	 * @param isAlt if true, method returns magnitude and phase,otherwise real and imaginary part
	 * @return Returns image in frequency domain.
	 */
	public static ComplexMatrix FourierTransform2D(IcyBufferedImage img,boolean isAlt){
		double[][] matrix = new double[img.getHeight()][img.getWidth()];
		for (int y = 0; y < img.getHeight(); y++) {
			for (int x = 0; x < img.getWidth(); x++) {
				matrix[y][x] = img.getRaster().getSample(x, y, 0);
			}
		}
//		return newFour(matrix,isAlt);
		return FourierTransform2D(matrix,isAlt);
	}
	
	/**
	 * Transforms image to frequency domain
	 * @param img input image
	 * @param isAlt if true, method returns magnitude and phase,otherwise real and imaginary part 
	 * @return image in frequency domain
	 */
	public static ComplexMatrix FourierTransform2D(BufferedImage img,boolean isAlt){
		double[][] matrix = new double[img.getHeight()][img.getWidth()];
		for (int y = 0; y < img.getHeight(); y++) {
			for (int x = 0; x < img.getWidth(); x++) {
				matrix[y][x] = img.getRaster().getSample(x, y, 0);
			}
		}
		return FourierTransform2D(matrix,isAlt);
	}
	/**
	 * Transforms image to frequency domain
	 * @param input 2D double array
	 * @param isAlt if true, method returns magnitude and phase,otherwise real and imaginary part
	 * @return image in frequency domain
	 */
	public static ComplexMatrix FourierTransform2D(double[][] input,boolean isAlt) {
		if (input == null)
			return null;
		int height = 1, width = 1;
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
//				verzia spravená 26.3.
//				output.setElement(y, x, rowMatrix.getElementCopy(0, x));
				
//				verzia pred 26.3.
				output.setElement(y, /*x,rowMatrix.getElementCopy(0, x) );*/(x+output.getNcol()/2)%output.getNcol(), rowMatrix.getElementCopy(0, x));

				
//				output.setElement(y, x, output.getElementCopy(y, x).getReal(),output.getElementCopy(y, x).getImag());
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
//				verzia spravená 26.3.
//				output.setElement(y, x, rowMatrix.getElementCopy(0, y));
				
//				verzia pred 26.3.
				output.setElement(y, x, /*rowMatrix.getElementCopy(0, y));*/rowMatrix.getElementCopy(0, (y+rowMatrix.getNcol()/2)%rowMatrix.getNcol()));
				
				
				//output.setElement(y, x, output.getElementCopy(y, x).getReal(),output.getElementCopy(y, x).getImag()*(-1));
			
			
			}
		}
		double tmp0,tmp1;
		if(isAlt)
		for(int y=0;y<output.getNrow();y++){
			for(int x=0;x<output.getNcol();x++){
				tmp0 = output.getElementCopy(y, x).getReal();
				tmp1 = output.getElementCopy(y, x).getImag();
				output.setElement(y, x, Math.sqrt(tmp0*tmp0+tmp1*tmp1), Math.atan2(tmp1,tmp0));
			}
		}
		return output;
	}

	public static void getThatFourier(double[][] array){
		FourierTransform ft;
		ComplexMatrix complex = new ComplexMatrix(2, 4);
		for(int j = 0;j<2;j++){
			ft = new FourierTransform(array[j]);
			ft.transform();
			for(int i = 0;i<ft.getTransformedDataAsComplex().length;i++){
				complex.setElement(j, i, ft.getTransformedDataAsComplex()[i]);
//				System.out.println(j+":Real:"+ft.getTransformedDataAsComplex()[i].getReal() + " Imag:"+ft.getTransformedDataAsComplex()[i].getImag());
			}
		}
		Complex[] tmpRow = new Complex[2];
		for(int j = 0;j<4;j++){
			for(int i = 0;i<2;i++){
				tmpRow[i] = complex.getElementCopy(i, j);
			}
			ft = new FourierTransform(tmpRow);
			ft.transform();
			for(int i = 0;i<2;i++){
				complex.setElement(i, j, ft.getTransformedDataAsComplex()[i]);
//				System.out.println(j+":Real:"+ft.getTransformedDataAsComplex()[i].getReal() + " Imag:"+ft.getTransformedDataAsComplex()[i].getImag());
			}
		}
	}
	
	/**
	 * Creates double array from complex matrix
	 * @param input complex matrix
	 * @param isAlt if true returns real part,if false returns imaginary part
	 * @return 2D double array
	 */
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
	
	/**
	 * Transforms image from frequency domain.
	 * @param input image in frequency domain.
	 * @param isAlt if true,computes from magnitude and phase
	 * @return image in spatial domain
	 */
	public static double[][][] InverseFourierTransform2D(ComplexMatrix input) {
		if (input == null)
			return null;
		
		/*
		 * Inverzná Fourierova transformácia nad každým riadkom a následne priradenie
		 * hodnôt do matice output
		 */
		double[][][] output = new double[2][input.getNrow()][input.getNcol()];
		
		ComplexMatrix rowMatrix;
		Complex[] tmpRow = new Complex[input.getNcol()];
		Complex number;
		for (int y = 0; y < input.getNrow(); y++) {
			for (int x = 0; x < input.getNcol(); x++) {
				//pozor na zbytočný kód	
				//26.3.
//				number = input.getElementCopy(y, x);
				
//				verzia pred 26.3.
				number = input.getElementCopy(y, /*x);*/(x+input.getNcol()/2)%input.getNcol());
					
					
				//number.setImag(number.getImag()*(-1));
					
					
				tmpRow[x] = number;
				
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
//				verzia spravená 26.3.
//				number = input.getElementCopy(y, x);
				
//				verzia pred 26.3.
				number = input.getElementCopy((y+input.getNrow()/2)%input.getNrow(), x);

				
//				number.setImag(number.getImag()*(-1));
				tmpRow[y] = number;
				//tmpRow[y] = input.getElementCopy(y, x);
				
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
				output[1][y][x] = input.getElementCopy(y, x).getImag();//*(-1);
			}
		}
		return output;
	}

	/**
	 * Create complex matrix from double array.
	 * @param input 3D double array
	 * @return 
	 */
	public static ComplexMatrix createComplexMatrix(double[][][] input){
		ComplexMatrix output = new ComplexMatrix(input[0].length,input[0][0].length);
		for(int y=0;y<input[0].length;y++){
			for(int x=0;x<input[0][y].length;x++){
				output.setElement(y, x, input[0][y][x], input[1][y][x]);
			}
		}
		return output;
	}
	
	public static double[][] shiftArray(double[][] input){
		double[][] output = new double[input.length][input[0].length];
		for(int y = 0;y<input.length;y++){
			for(int x = 0;x<input[y].length;x++){
				output[y][x] = input[input.length-1-y][input[y].length-1-x];
			}
		}
		return output;
	}
}
