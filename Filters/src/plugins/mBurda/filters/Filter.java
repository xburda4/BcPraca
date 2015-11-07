package plugins.mBurda.filters;

import java.awt.image.WritableRaster;

import Jama.Matrix;
import icy.image.IcyBufferedImage;

public abstract class Filter {
	protected static IcyBufferedImage source;
	protected static double[][] phaseCong;
	
	
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
	
//	public double[] getEigenValuesWithPhase(int x,int y){
//		Matrix mat = new Matrix(2,2);
//		
//	}
	
	public double[] getEigenValuesWithPhase(int x,int y){
		Matrix hessian = new Matrix(2,2);
		WritableRaster raster = source.getRaster();
		
		if(x+1<source.getWidth() && x-1>=0)
		{
			//Horizontal approximation
			double a = (raster.getSample(x+1, y, 0)*phaseCong[x+1][y]+raster.getSample(x-1, y, 0)*phaseCong[x-1][y]-2*raster.getSample(x, y, 0)*phaseCong[x][y]);
			hessian.set(0, 0, a); 	
		} else hessian.set(0, 0, 0);
		if(y+1 < source.getHeight() && y-1>=0)
		{
			//Vertical approximation
			double a = (raster.getSample(x, y+1, 0)*phaseCong[x][y+1]+raster.getSample(x, y-1, 0)*phaseCong[x][y-1]-2*raster.getSample(x, y, 0)*phaseCong[x][y]);
				hessian.set(1, 1, a);
		} else hessian.set(1, 1, 0);
		if(x+1<source.getWidth() && x-1>=0 && y+1 < source.getHeight() && y-1>=0)
		{
			//Diagonal approximation
			double a = (raster.getSample(x+1, y+1, 0)*phaseCong[x+1][y+1]+raster.getSample(x-1, y-1, 0)*phaseCong[x-1][y-1]-raster.getSample(x+1, y-1, 0)*phaseCong[x+1][y-1]-raster.getSample(x-1, y+1, 0)*phaseCong[x-1][y+1])/4;
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
}
