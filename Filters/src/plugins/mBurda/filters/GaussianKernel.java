package plugins.mBurda.filters;

public class GaussianKernel {
	
	public double[] computeGaussian(int radius){
		double[] gaussian = new double[radius*2 +1];
		double sigma = radius/3.0;
		for(int x=0;x<gaussian.length;x++){
				double r = radius - x;
				gaussian[x] = Math.exp(-0.5 * (r*r)/(sigma*sigma));
		}
		return gaussian;
	}
}
