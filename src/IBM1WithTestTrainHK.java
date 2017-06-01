import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * 
 * @author anahita
 *
 */
public class IBM1WithTestTrainHK extends IBM_Model1{
	
	private static final long serialVersionUID = 1L;
	//t parameter of the 20K model
	HashMap<Pair, Double> t_ibm1_20k;
	
	public IBM1WithTestTrainHK(HashMap<Pair, Double> t_ibm1_20k){
		this.t_ibm1_20k = t_ibm1_20k;
	}
	
	/**
	 * Smoothes t(f|e) parameter by backing-off to a uniform probability 1/V
	 * @param f  source word
	 * @param e  target word
	 * @param t_table  t parameter
	 * @return		   the somoothed t(f|e)
	 */
	@Override
	public double tProbability(String f, String e, HashMap<Pair, Double> t_table) {
		// TODO Auto-generated method stub
		Pair fe = new Pair(f,e);
		
		int V = 163303;
		if (t_ibm1_20k.containsKey(fe))
			return t_ibm1_20k.get(fe);
		else if (t_table.containsKey(fe))
			return t_table.get(fe);
		else
			return 1.0/V; 
	}
	public static void mainIBM1(int trainingSize, int testSize, String trainPrefix
			, String sourceLang, String targetLang, String testPrefix, String referenceFile, String augmentedTrainPrefix, int augmentedTrainSize, String alignmentFile) throws IOException{
		IBM_Model1 ibm1 = new IBM_Model1();
		//int trainingSize = 20357;//2373245;//18000;2373245
        //int testSize = 1956;//20357;//2357;20357
		
        
        //String[][] trainBitext = ibm1.readBitext("data-te/train.20k.seg.cln.cn","data-te/train.20k.seg.cln.en",trainingSize);
        //String[][] testBitext = ibm1.readBitext("data-te/test.seg.cln.cn","data-te/test.seg.cln.en",testSize);
		
		String[][] trainBitext = IBM_Model1.readBitext(trainPrefix+"."+sourceLang, trainPrefix+"."+targetLang, trainingSize);
		String[][] testBitext = IBM_Model1.readBitext(testPrefix+"."+sourceLang, testPrefix+"."+targetLang, testSize);
        //trainingSize = 1020357;
        //String[][] train20KHKBitext = ibm1.readBitext("data-te/train.20k+hk.cn","data-te/train.20k+hk.en",trainingSize);
		String[][] train20KHKBitext = IBM_Model1.readBitext(augmentedTrainPrefix+"."+sourceLang, augmentedTrainPrefix+"."+targetLang, augmentedTrainSize);
		
		ibm1.initializeCountsWithoutSets(trainBitext);
		System.out.println("length " + trainBitext.length);
		
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainBitext);
		
	/*	ArrayList<String> ibmModelAlignment = ibm1.print_alignment_SD_ibm1(testBitext, t_fe, alignmentFile);
		//ArrayList<String> reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		ArrayList<String> reference = ibm1.convertFileToArrayList(referenceFile);
		ibm1.gradeAlignmentWithType(testSize, testBitext, reference, ibmModelAlignment);
	*/	
		IBM1WithTestTrainHK ibm120KHk = new IBM1WithTestTrainHK(t_fe);
		ibm120KHk.initializeCountsWithoutSets(train20KHKBitext);
		System.out.println("length " + train20KHKBitext.length);
		HashMap<Pair, Double> t_fe_20KHK = ibm120KHk.EM_IBM1(ibm120KHk.f_count, ibm120KHk.e_count, ibm120KHk.fe_count, train20KHKBitext);
		
		ArrayList<String> ibmModelAlignment = ibm120KHk.print_alignment_SD_ibm1(testBitext, t_fe_20KHK, alignmentFile);
		if (referenceFile != null){
			ArrayList<String> reference = ibm120KHk.convertFileToArrayList(referenceFile);
			ibm120KHk.gradeAlignmentWithType(testSize, testBitext, reference, ibmModelAlignment);
		}
	}
	public static void main(String[] args) throws IOException {
		IBM_Model1 ibm1 = new IBM_Model1();
		int trainingSize = 20357;//2373245;//18000;2373245
        int testSize = 1956;//20357;//2357;20357
		
        
        String[][] trainBitext = IBM_Model1.readBitext("data-te/train.20k.seg.cln.cn","data-te/train.20k.seg.cln.en",trainingSize);
        String[][] testBitext = IBM_Model1.readBitext("data-te/test.seg.cln.cn","data-te/test.seg.cln.en",testSize);
		
        trainingSize = 1020357;
        String[][] train20KHKBitext = IBM_Model1.readBitext("data-te/train.20k+hk.cn","data-te/train.20k+hk.en",trainingSize);
		
		
		ibm1.initializeCountsWithoutSets(trainBitext);
		System.out.println("length " + trainBitext.length);
		
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainBitext);
		
		ArrayList<String> ibmModelAlignment = ibm1.print_alignment_SD_ibm1(testBitext, t_fe, "alignment");
		ArrayList<String> reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		ibm1.gradeAlignmentWithType(testSize, testBitext, reference, ibmModelAlignment);
	
		IBM1WithTestTrainHK ibm120KHk = new IBM1WithTestTrainHK(t_fe);
		ibm120KHk.initializeCountsWithoutSets(train20KHKBitext);
		System.out.println("length " + train20KHKBitext.length);
		HashMap<Pair, Double> t_fe_20KHK = ibm120KHk.EM_IBM1(ibm120KHk.f_count, ibm120KHk.e_count, ibm120KHk.fe_count, train20KHKBitext);
		
		ibmModelAlignment = ibm120KHk.print_alignment_SD_ibm1(testBitext, t_fe_20KHK,"alignment");
		reference = ibm120KHk.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		ibm120KHk.gradeAlignmentWithType(testSize, testBitext, reference, ibmModelAlignment);
		
	}
	
}
