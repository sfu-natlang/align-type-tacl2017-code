import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class HMMWithNullTestTrainHKTe extends HMMWithNullForTuningTe{
	private static final long serialVersionUID = 1L;
	
	HashMap<String, Integer> f_count = new HashMap<String, Integer>();
    HashMap<String, Integer> e_count = new HashMap<String, Integer>();
    HashMap<Pair, Integer> fe_count = new HashMap<Pair, Integer>();
	
    HMMWithNullForTuningTe hmm20K;
    /**
     * 
     * @param t_table
     * @param e_count
     * @param fe_count
     * @param hmm20K initialize this model with parameters of the HMM model created from the 20K data. 
     */
	public HMMWithNullTestTrainHKTe(HashMap<Pair, Double> t_table, HashMap<String, Integer> e_count, HashMap<Pair, Integer> fe_count, HMMWithNullForTuningTe hmm20K) {
		super(t_table);
		this.t_table = t_table;
		this.e_count = e_count;
		this.fe_count = fe_count;
		this.hmm20K = hmm20K;
	}
	/**
	 * Smoothes t(f|e) parameter by using the t parameters of 20K model.
	 * If the pair was not seen, uses the t parameter of the model created from 1M sentences.
	 * Otherwise, back-off to a uniform probability 1/V.
	 * @param f  source word
	 * @param e  target word
	 * @return   the somoothed t(f|e)
	 */
	public double tProbability(String f, String e){
		Pair fe = new Pair(f,e);
		
		int V = 163303;
		 if (hmm20K.t_table.containsKey(fe))
				return hmm20K.t_table.get(fe);
		 else if (t_table.containsKey(fe)){
			
			return t_table.get(fe);
		}else{
			if (e.equals("null"))
				return nullEmissionProb;
			else
				return 1.0/V;
		}
			 
	}
	public static void mainHMM(int trainingSize, int testSize, String trainPrefix
			, String sourceLang, String targetLang, String testPrefix, String referenceFile, String alignmentFile, String augmentedTrainPrefix, int augmentedTrainSize) throws IOException{
		IBM_Model1 ibm1 = new IBM_Model1();
		//int trainingSize = 19987;//20357;
		//int testSize = 1956;
		//String[][] trainBitext = ibm1.readBitext("data-te/train.20k.seg.cln.cn","data-te/train.20k.seg.cln.en",trainingSize);
		//String[][] testBitext = ibm1.readBitext("data-te/test.seg.cln.cn","data-te/test.seg.cln.en",testSize);
		
		//String[][] trainBitext = ibm1.readBitext("data-jasneet/train.20K.clean.cn","data-jasneet/train.20K.clean.en",trainingSize);
		
		String[][] trainBitext = IBM_Model1.readBitext(trainPrefix+"."+sourceLang, trainPrefix+"."+targetLang, trainingSize);
		String[][] testBitext = IBM_Model1.readBitext(testPrefix+"."+sourceLang, testPrefix+"."+targetLang, testSize);
		
		
		ibm1.initializeCountsWithoutSets(trainBitext);
		System.out.println("length " + trainBitext.length);
		
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainBitext);
		
		HMMWithNullForTuningTe hmm20K = new HMMWithNullForTuningTe(t_fe);
		hmm20K.baumWelch(trainBitext, ibm1.f_count, ibm1.fe_count);
		hmm20K.multiplyOneMinusP0H(trainBitext, hmm20K.p0H);
		System.out.println("p0H " + hmm20K.p0H + " nullEmissionProb" + hmm20K.nullEmissionProb);
		
		System.out.println("Baseline HMM for 20K train, 2k test:");
/*		ArrayList<String> hmmModelAlignment1 = hmm20K.findBestAlignmentsForAll_AER(testBitext, testSize, "alignment_baseline_20K_2K_Aug26");
		ArrayList<String> reference1 = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		double fScore1 = ibm1.gradeAlignmentWithType(testSize, testBitext, reference1, hmmModelAlignment1);
*/
		
		trainingSize = 1019987;//1020357;//2373245;//1020357; 
		IBM_Model1 ibm120KHK = new IBM_Model1();
		//String[][] train20KHKBitext = ibm120KHK.readBitext("data-te/train.20k+hk.cn","data-te/train.20k+hk.en",trainingSize);
		//String[][] train20KHKBitext = ibm120KHK.readBitext("data-jasneet/train.20K+1M.word.clean.cn","data-jasneet/train.20K+1M.word.clean.en",trainingSize);
		String[][] train20KHKBitext = IBM_Model1.readBitext(augmentedTrainPrefix+"."+sourceLang, augmentedTrainPrefix+"."+targetLang, augmentedTrainSize);
		//String[][] train20KHKBitext = ibm120KHK.readTritext("data-jasneet/train.20K+2M.word.cn","data-jasneet/train.20K+2M.word.en","data-te/train.20k.seg.cln.wa",trainingSize);
		ibm120KHK.initializeCountsWithoutSets(train20KHKBitext);
		System.out.println("length " + train20KHKBitext.length);
		HashMap<Pair, Double> t_fe_20KHK = ibm1.EM_IBM1(ibm120KHK.f_count, ibm120KHK.e_count, ibm120KHK.fe_count, train20KHKBitext);
		
		HMMWithNullTestTrainHKTe hmm = new HMMWithNullTestTrainHKTe(t_fe_20KHK, ibm120KHK.e_count, ibm120KHK.fe_count, hmm20K);
		hmm.baumWelch(train20KHKBitext, ibm120KHK.f_count, ibm120KHK.fe_count);
		hmm.multiplyOneMinusP0H(train20KHKBitext,hmm.p0H);
      	 	
/*		FileOutputStream fos = new FileOutputStream("hmmWithNull20KHK.ser");
                ObjectOutputStream oos = new ObjectOutputStream(fos);
                oos.writeObject(hmm);
                fos.close();
                oos.close(); 
  */      	
		
/*		FileInputStream fis = new FileInputStream("hmmWithNull20KHK.ser");
                ObjectInputStream ois = new ObjectInputStream(fis);
                HMMWithNullTestTrainHKTe hmm = (HMMWithNullTestTrainHKTe) ois.readObject();	
*/		hmm.p0H = 0.3;
		hmm.nullEmissionProb = 0.000005;
		System.out.println("p0H " + hmm.p0H + " nullEmissionProb" + hmm.nullEmissionProb);		
		
		System.out.println("Baseline HMM for Augmented model HK+20K train, 2k test:");

		ArrayList<String> hmmModelAlignment = hmm.findBestAlignmentsForAll_AER(testBitext, testSize, alignmentFile);
		//ArrayList<String> reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		if (referenceFile != null){
			ArrayList<String> reference = ibm1.convertFileToArrayList(referenceFile);
			ibm1.gradeAlignmentWithType(testSize, testBitext, reference, hmmModelAlignment);
		}
		
		//hmm.findBestAlignmentsForAll_AER(train20KHKBitext, trainingSize, "alignment_baseline_train_Jan25");
	}
	public static void testBaselineHMM() throws IOException{
		IBM_Model1 ibm1 = new IBM_Model1();
		int trainingSize = 19987;//20357;
		//int testSize = 1956;
		//String[][] trainBitext = ibm1.readBitext("data-te/train.20k.seg.cln.cn","data-te/train.20k.seg.cln.en",trainingSize);
		//String[][] testBitext = ibm1.readBitext("data-te/test.seg.cln.cn","data-te/test.seg.cln.en",testSize);
		
		String[][] trainBitext = IBM_Model1.readBitext("data-jasneet/train.20K.clean.cn","data-jasneet/train.20K.clean.en",trainingSize);
		ibm1.initializeCountsWithoutSets(trainBitext);
		System.out.println("length " + trainBitext.length);
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainBitext);
		HMMWithNullForTuningTe hmm20K = new HMMWithNullForTuningTe(t_fe);
		hmm20K.baumWelch(trainBitext, ibm1.f_count, ibm1.fe_count);
		hmm20K.multiplyOneMinusP0H(trainBitext, hmm20K.p0H);
		System.out.println("p0H " + hmm20K.p0H + " nullEmissionProb" + hmm20K.nullEmissionProb);
		
		System.out.println("Baseline HMM for 20K train, 2k test:");
/*		ArrayList<String> hmmModelAlignment1 = hmm20K.findBestAlignmentsForAll_AER(testBitext, testSize, "alignment_baseline_20K_2K_Aug26");
		ArrayList<String> reference1 = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		double fScore1 = ibm1.gradeAlignmentWithType(testSize, testBitext, reference1, hmmModelAlignment1);
*/
		
		trainingSize = 1019987;//1020357;//2373245;//1020357; 
		IBM_Model1 ibm120KHK = new IBM_Model1();
		//String[][] train20KHKBitext = ibm120KHK.readBitext("data-te/train.20k+hk.cn","data-te/train.20k+hk.en",trainingSize);
		String[][] train20KHKBitext = IBM_Model1.readBitext("data-jasneet/train.20K+1M.word.clean.cn","data-jasneet/train.20K+1M.word.clean.en",trainingSize);
		//String[][] train20KHKBitext = ibm120KHK.readTritext("data-jasneet/train.20K+2M.word.cn","data-jasneet/train.20K+2M.word.en","data-te/train.20k.seg.cln.wa",trainingSize);
		ibm120KHK.initializeCountsWithoutSets(train20KHKBitext);
		System.out.println("length " + train20KHKBitext.length);
		HashMap<Pair, Double> t_fe_20KHK = ibm1.EM_IBM1(ibm120KHK.f_count, ibm120KHK.e_count, ibm120KHK.fe_count, train20KHKBitext);
		
		HMMWithNullTestTrainHKTe hmm = new HMMWithNullTestTrainHKTe(t_fe_20KHK, ibm120KHK.e_count, ibm120KHK.fe_count, hmm20K);
		hmm.baumWelch(train20KHKBitext, ibm120KHK.f_count, ibm120KHK.fe_count);
		hmm.multiplyOneMinusP0H(train20KHKBitext,hmm.p0H);
      	 	
/*		FileOutputStream fos = new FileOutputStream("hmmWithNull20KHK.ser");
                ObjectOutputStream oos = new ObjectOutputStream(fos);
                oos.writeObject(hmm);
                fos.close();
                oos.close(); 
  */      	
		
/*		FileInputStream fis = new FileInputStream("hmmWithNull20KHK.ser");
                ObjectInputStream ois = new ObjectInputStream(fis);
                HMMWithNullTestTrainHKTe hmm = (HMMWithNullTestTrainHKTe) ois.readObject();	
*/		hmm.p0H = 0.3;
		hmm.nullEmissionProb = 0.000005;
		System.out.println("p0H " + hmm.p0H + " nullEmissionProb" + hmm.nullEmissionProb);		
		
		System.out.println("Baseline HMM for Augmented model HK+20K train, 2k test:");

		/*ArrayList<String> hmmModelAlignment = hmm.findBestAlignmentsForAll_AER(testBitext, testSize, "alignment_baseline_test_train_HK_tProb_all_Aug26");
		ArrayList<String> reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		double fScore = ibm1.gradeAlignmentWithType(testSize, testBitext, reference, hmmModelAlignment);
		hmm.findBestAlignmentsForAll_AER(train20KHKBitext, trainingSize, "alignment_baseline_train_Aug12");
		*/
		hmm.findBestAlignmentsForAll_AER(train20KHKBitext, trainingSize, "alignment_baseline_train_Jan25");
	}
	
	/**
	 * Computes F-Score of a given alignment file 
	 * @throws IOException
	 */
	public static double fscoreOfTest() throws IOException{
		IBM_Model1 ibm1 = new IBM_Model1();
		int testSize = 1956;
		String[][] testBitext = IBM_Model1.readBitext("data-te/test.seg.cln.cn","data-te/test.seg.cln.en",testSize);
		ArrayList<String> reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		ArrayList<String> hmmModelAlignment = ibm1.convertFileToArrayList("../../working-project/trainMosesWithGizaOnTestTrain/model/aligned.non-giza.srctotgt");///data-jasneet/fscoreWithIncGiza/20k-2k/test.non-giza.aligned.srctotgt");//alignment.non-giza.srctotgt");//20K+1M-2k/test.aligned.non-giza.step08.r.srctotgt");//test.non-giza.aligned.gdf");//test.aligned2.non-giza.srctotgt");//test.non-giza.aligned.srctotgt");
		double fScore = ibm1.gradeAlignmentWithType(testSize, testBitext, reference, hmmModelAlignment);
		return fScore;

	}
	public static void main(String[] args) throws IOException, ClassNotFoundException{
		System.out.println("F-score = " + fscoreOfTest());
		//System.out.println("F-score = " + fScore());
		//testTrainData();
		//testReverseBaselineHMM();
		//testBaselineHMM();
	}
}
