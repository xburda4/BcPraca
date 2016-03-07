package plugins.mBurda.filters;

import java.awt.image.WritableRaster;

import Jama.Matrix;
import icy.image.IcyBufferedImage;

public abstract class Filter {
	protected static IcyBufferedImage source;
	protected static Matrix[][] phaseCong;
	public static double[][] phaseValues;
	/**
	 * Computes eigenvalues of image at point (x,y)
	 * @param x x-coordinate
	 * @param y y-coordinate
	 * @return array of eigenvalues
	 */
	public static double[] getEigenValues(int x,int y){
		Matrix hessian = new Matrix(2,2);
		WritableRaster raster = source.getRaster();
		
		if(x+1<source.getWidth() && x-1>=0)
		{
			//Horizontal approximation
			double a = (raster.getSample(x+1, y, 0)+raster.getSample(x-1, y, 0)-2*raster.getSample(x, y, 0));
			hessian.set(0, 0, a); 	
		} else hessian.set(0, 0, 0);
		if(y+1 < source.getHeight() && y-1>=0)
		{
			//Vertical approximation
			double a = (raster.getSample(x, y+1, 0)+raster.getSample(x, y-1, 0)-2*raster.getSample(x, y, 0));
			hessian.set(1, 1, a);
		} else hessian.set(1, 1, 0);
		if(x+1<source.getWidth() && x-1>=0 && y+1 < source.getHeight() && y-1>=0)
		{
			//Diagonal approximation
			double a = (raster.getSample(x+1, y+1, 0)+raster.getSample(x-1, y-1, 0)-raster.getSample(x+1, y-1, 0)-raster.getSample(x-1, y+1, 0))/4;
			hessian.set(0,1,a);
		} else hessian.set(0, 1, 0);
		hessian.set(1, 0, hessian.get(0,1));
		
		double[] eig = hessian.eig().getRealEigenvalues();
		double tmp = 0;
		if(Math.abs(eig[0]) > Math.abs(eig[1])){
			tmp = eig[1];
			eig[1] = eig[0];
			eig[0] = tmp;
		}
		return eig;
	}
	
	/**
	 * Computes eigen-values of tensor of phase congruencies
	 * @param x x-coordinate
	 * @param y y-coordinate
	 * @return array of eigen-values
	 */
	public double[] getEigenValuesWithPhase(int x,int y){
		double[] eig = phaseCong[y][x].eig().getRealEigenvalues();
		double tmp = 0;
		if(Math.abs(eig[0]) > Math.abs(eig[1])){
			tmp = eig[1];
			eig[1] = eig[0];
			eig[0] = tmp;
		}
		return eig;
	}
}
