package plugins.mBurda.filters;


import java.awt.Graphics;
import java.awt.image.BufferedImage;

import org.jdesktop.swingx.image.GaussianBlurFilter;

import flanagan.complex.ComplexMatrix;

import icy.gui.dialog.MessageDialog;
import icy.image.IcyBufferedImage;
import icy.image.colormodel.IcyColorModel;
import icy.plugin.abstract_.PluginActionable;
import icy.sequence.Sequence;
import icy.type.DataType;

public class Filters extends PluginActionable {
	@Override
	public void run() {
//		if(getActiveSequence() != null) 
//		System.out.println("Depth is " + getActiveSequence().getSizeT() );
//		
//		Vesselness3D ves = new Vesselness3D(getActiveSequence(),3,0.5,30);
//		ves.makeImage3D();
//		addSequence(ves.ret);
		
		//rozmazanie obrazu s okolim 
		GaussianBlurFilter gauss = new GaussianBlurFilter(3);
		BufferedImage image = getGrayScale(gauss.filter(getActiveSequence().getImage(0, 0), null));
//
//		Vesselness ves = new Vesselness(blurred,1.5,0.08);
//		ves.makeImage2D();
//		addSequence(new Sequence("Vesselness",ves.ret));
//
//		Neuriteness neu = new Neuriteness(blurred,-5);
//		neu.makeImage2D();
//		addSequence(new Sequence("Neuriteness",neu.ret));
		
		
//		BufferedImage img = getGrayScale(getActiveSequence().getImage(0, 0));
//		
//		
//		double[][] image = new double[img.getRaster().getHeight()][img.getRaster().getWidth()];
//		for(int y = 0;y<img.getRaster().getHeight();y++){
//			for(int x=0;x<img.getRaster().getWidth();x++){
//				image[y][x] = img.getRaster().getSampleDouble(x, y, 0);
//			}
//		}
		
		ComplexMatrix matrix = Computations.FourierTransform2D(image,false);
		
		
		double[] scs = {0.0015,0.002,0.005,0.009,0.013};
		double[] ors = {0,90,180,270};
		
		Vesselness2D ves = new Vesselness2D(image,5.5,0.007);
		
		Neuriteness2D neu = new Neuriteness2D(image,-5);
		
		ves.phaseCong = Computations.getPhaseCong(matrix, scs, ors, -2, -55 ,8);
		neu.phaseCong = Computations.getPhaseCong(matrix, scs, ors, -2, -55 ,8);
//		ves.phaseCong = Computations.getPhaseCong(matrix, scs, ors, -2, 0.4 ,10);
		ves.makeImage2D();
		addSequence(new Sequence("Vesselness",ves.ret));
		ves.makeImageWithPhase2D();
		addSequence(new Sequence("Vesselness with phase",ves.ret));
		neu.makeImage2D();
		addSequence(new Sequence("Neuriteness",neu.ret));
		neu.makeImageWithPhase2D();
		addSequence(new Sequence("Neuriteness with phase",neu.ret));
//		addSequence(new Sequence("invFTT",makeImage2D(Computations.InverseFourierTransform2D(matrix, true))));
//		
//		addSequence(new Sequence("Magnitude",makeImage2D(Computations.fttToDoubleArr(matrix,true))));
//		addSequence(new Sequence("Phase",makeImage2D(Computations.fttToDoubleArr(matrix, false))));
		
//		ves.makeImage2D();
//		addSequence(new Sequence("Final product - Vesselness",ves.ret));

		
//		MessageDialog.showDialog("Filt is done !");
	}

	public IcyBufferedImage makeImage2D(double[][] source){
		long start = System.nanoTime();
		IcyBufferedImage ret = new IcyBufferedImage(source[0].length,source.length,IcyColorModel.createInstance(1, DataType.DOUBLE));
		
		//Without this,it would take forever
		ret.beginUpdate();
		
		for(int y=0;y< source.length;y++){
			for(int x=0;x<source[y].length;x++){
				ret.setData(x, y, 0, source[y][x]);
			}
		}
		ret.endUpdate();
		
		System.out.println("Time of making an image is " + (System.nanoTime() - start));
		
		return ret;
	 }
	
	private BufferedImage getGrayScale(BufferedImage original){
		BufferedImage image = new BufferedImage(original.getWidth(), original.getHeight(),
				BufferedImage.TYPE_BYTE_GRAY);  
		Graphics g = image.getGraphics();
		g.drawImage(original, 0, 0, null);
		g.dispose(); 
		return image;
	}
}


