package plugins.mBurda.filters;

import org.jdesktop.swingx.image.GaussianBlurFilter;

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

	
	
	public double[][] getGaborKernel(int radiusOfGaussian){
		double[] gausKernel = createGaussianKernel(radiusOfGaussian);
		
		return null;
	}
	
	
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
		/*Fourierova transformácia nad každým riadkom a následne priradenie hodnôt do matice output
		 * 
		 * */
		ComplexMatrix rowMatrix;
		Complex[] tmpRow;
		for(int y = 0;y<output.getNrow();y++){
			tmpRow = new Complex[height];
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
		/*tmp je output matica otočená o 90 stupňov*/
		ComplexMatrix tmp = new ComplexMatrix(output.getNcol(), output.getNrow());
		for(int y = 0;y<output.getNrow();y++){
			for(int x=0;x<output.getNcol();x++){
				tmp.setElement(x, y, output.getElementCopy(y, x));
			}
		}
		tmpRow = new Complex[tmp.getNrow()];
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
				output.setElement(x, y, tmp.getElementCopy(x, y));
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
				tmp.setElement(y, x, input.getElementCopy(x,y));
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
				output[y][x] = tmp.getElementCopy(x, y).getReal();
			}
		}
		return output;
	}
	
}
