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
		
		
		String[][] trainBitext = IBM_Model1.readBitext(trainPrefix+"."+sourceLang, trainPrefix+"."+targetLang, trainingSize);
		String[][] testBitext = IBM_Model1.readBitext(testPrefix+"."+sourceLang, testPrefix+"."+targetLang, testSize);
		
		
		ibm1.initializeCountsWithoutSets(trainBitext);
		
		System.out.println("Training an IBM1 model on the bitext of size " + trainBitext.length + ".");
		
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainBitext);
		
		HMMWithNullForTuningTe hmm20K = new HMMWithNullForTuningTe(t_fe);
		System.out.println("Training an HMM model on the bitext of size " + trainBitext.length + ".");
		hmm20K.baumWelch(trainBitext, ibm1.f_count, ibm1.fe_count);
		hmm20K.multiplyOneMinusP0H(trainBitext, hmm20K.p0H);
		 
		IBM_Model1 ibm120KHK = new IBM_Model1();
		
		String[][] train20KHKBitext = IBM_Model1.readBitext(augmentedTrainPrefix+"."+sourceLang, augmentedTrainPrefix+"."+targetLang, augmentedTrainSize);
		
		ibm120KHK.initializeCountsWithoutSets(train20KHKBitext);
		System.out.println("Training an IBM1 model on the bitext of size " + train20KHKBitext.length + ".");
		HashMap<Pair, Double> t_fe_20KHK = ibm1.EM_IBM1(ibm120KHK.f_count, ibm120KHK.e_count, ibm120KHK.fe_count, train20KHKBitext);
		
		HMMWithNullTestTrainHKTe hmm = new HMMWithNullTestTrainHKTe(t_fe_20KHK, ibm120KHK.e_count, ibm120KHK.fe_count, hmm20K);
		System.out.println("Training an HMM model on the bitext of size " + train20KHKBitext.length + ".");
		hmm.baumWelch(train20KHKBitext, ibm120KHK.f_count, ibm120KHK.fe_count);
		hmm.multiplyOneMinusP0H(train20KHKBitext,hmm.p0H);
      	 	
		hmm.p0H = 0.3;
		hmm.nullEmissionProb = 0.000005;
		

		ArrayList<String> hmmModelAlignment = hmm.findBestAlignmentsForAll_AER(testBitext, testSize, alignmentFile);
		
		if (referenceFile != null){
			ArrayList<String> reference = ibm1.convertFileToArrayList(referenceFile);
			ibm1.gradeAlignmentWithType(testSize, testBitext, reference, hmmModelAlignment);
		}
		
	}	

}
