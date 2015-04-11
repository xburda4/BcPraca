package plugins.mBurda.filters;

//import org.jdesktop.swingx.image.GaussianBlurFilter;

import flanagan.complex.Complex;
import flanagan.complex.ComplexMatrix;
import flanagan.math.FourierTransform;

public class Computations {
	
	/**Method from swingx.image.GaussianBlurFilter
	 * */
	static double[] createGaussianKernel(int radius) {
        if (radius < 1) {
            throw new IllegalArgumentException("Radius must be >= 1");
        }

        double[] data = new double[radius * 2 + 1];

        float sigma = radius / 3.0f;
        float twoSigmaSquare = 2.0f * sigma * sigma;
        float sigmaRoot = (float) Math.sqrt(twoSigmaSquare * Math.PI);
        float total = 0.0f;

        for (int i = -radius; i <= radius; i++) {
            float distance = i * i;
            int index = i + radius;
            data[index] = (float) Math.exp(-distance / twoSigmaSquare) / sigmaRoot;
            total += data[index];
        }

        for (int i = 0; i < data.length; i++) {
            data[i] /= total;
        }

        return data;
    }

	
	
//	public double[][] getGaborKernel(int radiusOfGaussian){
//		double[] gausKernel = createGaussianKernel(radiusOfGaussian);
//		
//		return null;
//	}
	
	public static ComplexMatrix FourierTransform2D(double[][] input){//double[polohaVStlpci][polohaVRiadku]
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
			//	input.setElement(x, y, input.getElementCopy(x, y).getReal()/rowMatrix.getNrow(), input.getElementCopy(x, y).getImag()/rowMatrix.getNrow());
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
		//		input.setElement(x, y, input.getElementCopy(x, y).getReal()/rowMatrix.getNrow(), input.getElementCopy(x, y).getImag()/rowMatrix.getNrow());
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
