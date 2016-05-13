package plugins.mBurda.filters;

//import java.awt.Graphics;
//import java.awt.Image;
//import java.awt.Toolkit;
import java.awt.image.BufferedImage;
//import java.awt.image.FilteredImageSource;
//import java.awt.image.ImageFilter;
//import java.awt.image.ImageProducer;
//import java.io.File;
//import java.io.IOException;

//import javax.imageio.ImageIO;
//import javax.swing.GrayFilter;

//import org.jdesktop.swingx.image.GaussianBlurFilter;

//import flanagan.complex.ComplexMatrix;
import icy.gui.dialog.MessageDialog;
import icy.gui.frame.progress.AnnounceFrame;
//import icy.image.IcyBufferedImage;
//import icy.image.colormodel.IcyColorModel;
import icy.sequence.Sequence;
//import icy.type.DataType;
import plugins.adufour.ezplug.*;

public class CurvilinearStructures extends EzPlug {

	Vesselness2D ves = null;
	Neuriteness2D neu = null;
	
	
//	EzVarInteger blur = new EzVarInteger("Radius of blur",1,0,10,1);
	EzVarBoolean vesselness = new EzVarBoolean("Compute Vesselness?",true);
	EzVarBoolean neuriteness = new EzVarBoolean("Compute neuriteness?",true);
	EzVarDouble betaThreshold =new EzVarDouble("Disparsity control",0,-50,50,0.05);
	EzVarDouble gammaThreshold = new EzVarDouble("Relative brightness control",0,-30,30,0.02);
//	EzVarDouble alphaSteer = new EzVarDouble("Steerable filter equivalent",0,-50,50,0.05);
	EzVarDouble cutoffValue = new EzVarDouble("Cutoff value",0,-100,100,0.05);
	EzVarDouble threshold = new EzVarDouble("Threshold",0,-100,100,1);
	EzVarDouble gainFactor = new EzVarDouble("Gain factor",0,-50,50,1);
	EzVarBoolean vesPhase = new EzVarBoolean("Compute vesselness with phase congruency?",true);
	EzVarBoolean neuPhase = new EzVarBoolean("Compute neuriteness with phase congruency?",true);
	EzVarInteger angleKernels = new EzVarInteger("Number of kernels for angle computation.",2,2,8,2);
	EzVarInteger scaleKernels = new EzVarInteger("Number of kernels for scale computation.",1,1,8,1);
	
	EzGroup vesGroup = new EzGroup("Vesselness",betaThreshold,gammaThreshold);
//	EzGroup neuGroup = new EzGroup("Neuriteness",alphaSteer);
	EzGroup phaseCongGroup = new EzGroup("Phase congruency parameters",threshold,cutoffValue,gainFactor);
	EzGroup outputGroup = new EzGroup("Show output",vesselness,vesPhase,neuriteness,neuPhase);
	EzGroup kernelGroup = new EzGroup("Kernel parametrization",angleKernels,scaleKernels);
	
	@Override
	protected void initialize() {
//		addEzComponent(blur);
		addEzComponent(outputGroup);
		addEzComponent(vesGroup);
//		addEzComponent(neuGroup);
		addEzComponent(phaseCongGroup);
		addEzComponent(kernelGroup);
	}

	@Override
	protected void execute() {
		ves = null;
		neu = null;
		triggerStart();
//		GaussianBlurFilter gauss = new GaussianBlurFilter(blur.getValue());
//		BufferedImage img = (BufferedImage)getGrayScale(gauss.filter(getActiveImage(), null));
		BufferedImage img = /*gauss.filter(*/getActiveImage();//, null);
		if(img == null) {
			MessageDialog.showDialog("This plugin needs opened image.");
			return;
		} else {
			new AnnounceFrame("Plugin runs on image "+ getActiveSequence().getName());
		}
		double[] scs = {scaleKernels.getValue()};
		{
			int tmp = 1;
			for(int i=2;i<=scaleKernels.getValue();i*=2){
				tmp=i;
			}
			scs = new double[tmp];
		}
		for(int i=0;i<scs.length;i++){
			scs[i] = 3* Math.pow(2.1,i);
		}
		
		double[] ors = {angleKernels.getValue()};
		{
			int tmp = 2;
			for(int i=2;i<=angleKernels.getValue();i*=2){
				tmp=i;
			}
			ors = new double[tmp];
		}
		for(int i=0;i<ors.length;i++){
			ors[i] = ((double)(i)/ors.length*Math.PI);
		}
//		long time;
		if(vesselness.getValue()){
//			time = System.nanoTime();
			
			ves = new Vesselness2D(img, betaThreshold.getValue(), gammaThreshold.getValue());
			addSequence(new Sequence("Vesselness",ves.makeImage2D()));
			
//			time = System.nanoTime() - time;
//			System.out.println("Time of computing Vesselness is " + time/(10.0*100000) + "ns");
		}
		if(vesPhase.getValue()){
//			time = System.nanoTime();
			if(ves == null) ves = new Vesselness2D(img, betaThreshold.getValue(), gammaThreshold.getValue());
			Filter.phaseCong = Computations.getPhaseCong(Computations.FourierTransform2D(Filter.source,false), 
					scs, ors, threshold.getValue(), cutoffValue.getValue(), gainFactor.getValue(),img.getWidth(),img.getHeight());
			addSequence(new Sequence("Vesselness with phase congurency",ves.makeImageWithPhase2D()));
			
//			time = System.nanoTime() - time;
//			System.out.println("Time of computing Vesselness with PCT is " + time/(10.0*100000) + "ns");
		}
		if(neuriteness.getValue()){
//			time = System.nanoTime();
			
			neu = new Neuriteness2D(img, -1/3);
			addSequence(new Sequence("Neuriteness",neu.makeImage2D()));
			
//			time = System.nanoTime() - time;
//			System.out.println("Time of computing Neuriteness is " + time/(10.0*100000) + "ns");
		}
		if(neuPhase.getValue()){
//			time = System.nanoTime();
			
			if(neu == null) neu = new Neuriteness2D(img, -1/3);
			if(!vesPhase.getValue()) Filter.phaseCong = Computations.getPhaseCong(Computations.FourierTransform2D(Filter.source,true), 
					scs, ors, threshold.getValue(), cutoffValue.getValue(), gainFactor.getValue(),img.getWidth(),img.getHeight());
			addSequence(new Sequence("Neuriteness with phase",neu.makeImageWithPhase2D()));
			
//			time = System.nanoTime() - time;
//			System.out.println("Time of computing Neuriteness with PCT is " + time/(10.0*100000) + " ns");
		}
//		System.out.println("Memory usage - max : "+Runtime.getRuntime().maxMemory() + " bytes");
		triggerStop();
		/*
		 * Toto všetko je tu len na testovanie a výpis medzivýpočtov
		 * 
		 * */
//		double[][] inp = {{7,13,25,69},{854,1,444,32}};
//		Computations.getThatFourier(inp);
		
//		IcyBufferedImage in = getActiveImage();
//		addSequence(new Sequence("FTImag",Filters.makeImage2D(Computations.fttToDoubleArr2D(Computations.FourierTransform2D(in,false),false))));
//		addSequence(new Sequence("FTReal",Filters.makeImage2D(Computations.fttToDoubleArr2D(Computations.FourierTransform2D(in,false),true))));
//		double[][] kernel = new double[256][256];// = (Computations.getGaborKernel2D(256, 256, 6.3, 0,4));
//		double[][] lp = Computations.lowPassFilter(256, 256);
//		
//		for(int y=0;y<kernel.length;y++)
//			for(int x=0;x<kernel[y].length;x++)
//			{
//				kernel[y][x] = Computations.getLogGaborKernelPoint(x, y, 256, 256, 6.3, 0, 4, 0.75) * lp[y][x];
//			}
//		addSequence(new Sequence("Kernel",Filters.makeImage2D(kernel)));
		
//		addSequence(new Sequence("InverseF",Filters.makeImage2D(Computations.InverseFourierTransform2D(Computations.createComplexMatrix(Computations.multiFTKernel(Computations.FourierTransform2D(in,false), 0.07, 45)),false)[0])));
//		addSequence(new Sequence("MultiReal",Filters.makeImage2D(Computations.shiftArray(Computations.fttToDoubleArr2D(Computations.createComplexMatrix(Computations.multiFTKernel(Computations.FourierTransform2D(in,false), 6.3, 0,4)),true)))));
//		addSequence(new Sequence("MultiImag",Filters.makeImage2D(Computations.shiftArray(Computations.fttToDoubleArr2D(Computations.createComplexMatrix(Computations.multiFTKernel(Computations.FourierTransform2D(in,false), 6.3, 0,4)),false)))));
//
//		addSequence(new Sequence("InversMultiReal",Filters.makeImage2D(/*Computations.shiftArray*/(Computations.InverseFourierTransform2D((Computations.createComplexMatrix(Computations.multiFTKernel(Computations.FourierTransform2D(in,false), 6.3, 0, 4))))[0]))));         
//		addSequence(new Sequence("InversMultiImag",Filters.makeImage2D(/*Computations.shiftArray*/(Computations.InverseFourierTransform2D((Computations.createComplexMatrix(Computations.multiFTKernel(Computations.FourierTransform2D(in,false), 6.3, 0, 4))))[1]))));
//		//		
//		addSequence(new Sequence("InverseF",Filters.makeImage2D(Computations.InverseFourierTransform2D(Computations.FourierTransform2D(img,false))[0])));
		
		
//		addSequence(new Sequence("gab*lp",Filters.makeImage2D(Computations.shiftArray(Computations.multiFTKernel0(Computations.FourierTransform2D(in, false), 6.3, 0, 4)))));
//		addSequence(new Sequence("ftR",Filters.makeImage2D(/*Computations.shiftArray*/(Computations.multiFTKernel1(Computations.FourierTransform2D(in, false), 6.3, 0, 4)[0]))));
//		addSequence(new Sequence("ftI",Filters.makeImage2D(/*Computations.shiftArray*/(Computations.multiFTKernel1(Computations.FourierTransform2D(in, false), 6.3, 0, 4)[1]))));
//		addSequence(new Sequence("multiR",Filters.makeImage2D(Computations.shiftArray(Computations.multiFTKernel2(Computations.FourierTransform2D(in, false), 6.3, 0, 4)[0]))));
//		addSequence(new Sequence("multiI",Filters.makeImage2D(Computations.shiftArray(Computations.multiFTKernel2(Computations.FourierTransform2D(in, false), 6.3, 0, 4)[1]))));

		
		
		
//		getPhase();
		}
	/*
	 * Metóda bude zmazaná,slúži na výpis fázových kongruencií
	 * */
//	private void getPhase(){
//		IcyBufferedImage ret = new IcyBufferedImage(Filter.phaseValues.length, Filter.phaseValues[0].length, IcyColorModel.createInstance(1, DataType.DOUBLE));
//		ret.beginUpdate();
//		for(int y=0;y<Filter.phaseValues.length;y++){
//			for(int x=0;x<Filter.phaseValues[y].length;x++){
//				ret.setData(x, y, 0, Filter.phaseValues[y][x]);
//			}
//		}
//		ret.endUpdate();
//		addSequence(new Sequence("Phase congruency",ret));
//	}
//	private IcyBufferedImage getImage(){
////		IcyBufferedImage ret = new IcyBufferedImage(Filter.phaseCong[0].length, Filter.phaseCong.length, IcyColorModel.createInstance(1, DataType.DOUBLE));
//		IcyBufferedImage ret = new IcyBufferedImage(200, 200, IcyColorModel.createInstance(1, DataType.DOUBLE));
//		ret.beginUpdate();
//		double[][] dob = Computations.getGaborKernel2D(200, 200, 0.13, 90, 0.75, 3);
//		for(int y=0;y<200;y++){
//			for(int x=0;x<200;x++){
//				ret.setDataAsDouble(x, y, 0, dob[x][y]);
//			}
//		}
//		ret.endUpdate();
//		return ret;
//	}
	private void triggerStart(){
		System.out.println("Trigger start");
	}
	
	private void triggerStop(){
		System.out.println("Trigger Stop");
	}
	
	@Override
	public void clean() {
		// TODO Auto-generated by Icy4Eclipse
	}
}
