import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;


public class HMMWithNullForTuningTe implements Serializable{
	private static final long serialVersionUID = 1L;
	HashSet<Integer> targetLengthSet;
	
	
	protected double[][][] a;
	
	protected HashMap<Pair, Double> t_table;
	
	protected double[] pi;
	
	double p0H = 0.3;
	double nullEmissionProb = 0.000005;
	double smoothFactor = 0.1;	
	
	/**
	 * Constructor
	 * @param t_table   translation table which gives the probability p(f|e).
	 * 	                probabilities are stored in this HashMap of Source-Target Pairs (f, e).
	 */
	public HMMWithNullForTuningTe(HashMap<Pair, Double> t_table) {
		this.t_table = t_table;
	}

	/**
	 * Computes forward algorithm for HMM. Note that we re-scale the probabilities to avoid 0 probabilities for long sentences.  
	 * References:
	 * (1) wikipedia Baum Welch algorithm 
	 * (2) Mathematics of HMM and
	 * @param a  HMM transition parameters
	 * @param pi HMM start parameters
	 * @param y  source sentence
	 * @param N  length of the target sentence
	 * @param T  length of the source sentence
	 * @param d  target sentence
	 * @param t_table  translation table
	 * @return 		   a pair of alphaHat probabilities and c_hat (for alphaHat c_hat refer to reference (2))
	 */
	public GenericPair<double[][], double[]> forwardWithTScaled(double[][][] a, double[] pi, String[] y, int N, int T, String[] d, HashMap<Pair, Double> t_table){
		double[] c_scaled = new double[T+1];
		double[][] alphaHat = new double[N+1][];
		for (int i = 0; i < alphaHat.length ; i++){
			alphaHat[i] = new double[T+1];
		}
		double totalAlphaDoubleDot = 0;
		for (int i = 1; i < N+1 ; i++){
			
			alphaHat[i][1] = pi[i]*t_table.get(new Pair(y[0], d[i-1])).doubleValue();
			totalAlphaDoubleDot += alphaHat[i][1];
		}
		c_scaled[1] = 1.0/totalAlphaDoubleDot;
		for (int i = 1; i < N+1 ; i++){
			alphaHat[i][1] = c_scaled[1]*alphaHat[i][1];
		}
		for (int t = 1; t < T; t++){
			totalAlphaDoubleDot = 0;
			for (int j = 1; j < N + 1 ; j++){
				double total = 0;
				for (int i = 1; i < N + 1 ; i++)
					total += alphaHat[i][t]*a[i][j][N];
				alphaHat[j][t+1] = t_table.get(new Pair(y[t],d[j-1]) )*total;
				totalAlphaDoubleDot += alphaHat[j][t+1];
			}
			c_scaled[t+1] = 1.0/totalAlphaDoubleDot;
			for (int i = 1; i < N + 1 ; i++)
	            alphaHat[i][t+1] = c_scaled[t+1]*alphaHat[i][t+1];
		}
		return new GenericPair<double[][], double[]>(alphaHat, c_scaled);
	}
	/**
	 * Computes backward algorithm for HMM
	 * References:
	 * (1) wikipedia Baum Welch algorithm 
	 * (2) Mathematics of HMM and
	 * @param a  HMM transition parameters
	 * @param pi  HMM start parameters
	 * @param y  source sentence
	 * @param N  length of the target sentence
	 * @param T  length of the source sentence
	 * @param d  target sentence
	 * @param t_table  translation table
	 * @param c_scaled  c scaler (computed using forwardWithTScaled method) 
	 * @return betaHat probabilties
	 */
	public double[][] backwardWithTScaled(double[][][] a, double[] pi, String[] y, int N, int T, String[] d, HashMap<Pair, Double> t_table, double[] c_scaled){
		double[][] betaHat = new double[N+1][T+1]; 
		/*for (int i = 0; i < betaHat.length ; i++){
			betaHat[i] = new double[T+1];
		}*/
		for (int i = 1; i < N+1 ; i++){
			betaHat[i][T] = c_scaled[T];
		}
		for (int t = T-1; t > 0; t--){
			for (int i = 1; i < N+1 ; i++){
				double total = 0;
				for (int j = 1; j < N + 1 ; j++){
					if ( i > a.length || j > a[0].length){
						System.out.println("exiting " + i + " " + j + " " + a.length + " " + a[0].length);
						System.exit(1);
					}
					total += betaHat[j][t+1]*a[i][j][N]*t_table.get(new Pair(y[t], d[j-1]));
				}
				betaHat[i][t] = c_scaled[t]*total;
			}
		}
		return betaHat;
	}
	

	/**
	 * Initialize HMM parameters (a and pi). t_table is initialized from the t_table of IBM1. 
	 * @param N  length of the target sentence
	 */
	public void  initializeModel(int N){
		int twoN = 2*N;
		a = new double[twoN+1][twoN+1][N+1];
		
		pi = new double[twoN+1];
		for (int i = 1; i < twoN+1; i++)
			pi[i] = 1.0/twoN;
		for (int i = 1; i < N+1; i++)
			for (int j = 1; j < N + 1; j++)
				a[i][j][N] = 1.0/N;
		
	}
	/**
	 * Computes maximum length among target sentences
	 * @param bitext
	 * @return maximum length among target sentences
	 */
	public GenericPair<Integer, HashSet<Integer>> maxTargetSentenceLength(String[][] bitext){
		int maxLength = 0;
		
		HashSet<Integer> targetLengthSet = new HashSet<Integer>();
		for (String[] fe : bitext){
			int tempLength = fe[1].split(" ").length;
			if ( tempLength > maxLength)
				maxLength = tempLength;
			targetLengthSet.add(tempLength);
		}
		return new GenericPair<Integer, HashSet<Integer>>(maxLength, targetLengthSet);
	}
	/**
	 * This method reverse a bitext of <Pair, Integer> such that it becomes <Integer, Pair>.
	 * This is done for performance issues.  
	 * @param sd_count  source-target count map 
	 * @return the reversed bitext
	 */
	public GenericPair<HashMap<Pair, Integer>, HashMap<Integer, Pair>> mapBitextToInt(HashMap<Pair, Integer>sd_count){
		HashMap<Pair, Integer> index = new HashMap<Pair, Integer>();
		HashMap<Integer, Pair> biword = new HashMap<Integer, Pair>();
		int i = 0;
		for (Pair pair : sd_count.keySet()){
			index.put(pair, i);
			biword.put(i, pair);
			i++;
		}
		return new GenericPair<HashMap<Pair, Integer>, HashMap<Integer, Pair>>(index, biword);
	}
	public void addToMap(HashMap<Integer, Double> map, Integer key, double value){
		if (map.containsKey(key)){
			double currentValue = map.get(key);
			map.put(key, currentValue+value);
		}else{
			map.put(key, value);
		}
	}
	/**
	 * 
	 * @param tritext_sd tritext of train data 
	 * * 	
	 * @param s_count  source words count table
	 * @param sd_count  source-target words count table	 
	 */
	public void baumWelch(String[][] tritext_sd, HashMap<String, Integer>s_count,HashMap<Pair, Integer>sd_count) {
		
		
		GenericPair<Integer, HashSet<Integer>> NandtargetLengthSet = maxTargetSentenceLength(tritext_sd);
		
		int N = NandtargetLengthSet.a;
		targetLengthSet = NandtargetLengthSet.b;
		
		//System.out.println("N "+N);
		GenericPair<HashMap<Pair, Integer>, HashMap<Integer, Pair>> indexBiwordPair = mapBitextToInt(sd_count);
		HashMap<Pair, Integer> indexMap = indexBiwordPair.a;
		HashMap<Integer, Pair> biword = indexBiwordPair.b;
		
		int L = tritext_sd.length;
		int sd_size = indexMap.size();
		
		double[] totalGammaDeltaOverAllObservations_t_i = null;
		HashMap<String, Double>totalGammaDeltaOverAllObservations_t_overall_states_over_dest = null;
		
		for (int iteration = 0; iteration < 5; iteration++){
			
			double logLikelihood = 0;
			
			totalGammaDeltaOverAllObservations_t_i = new double[sd_size];
			
	        totalGammaDeltaOverAllObservations_t_overall_states_over_dest = new HashMap<String, Double>();

	        double[] totalGamma1OverAllObservations = new double[N+1];
	        double[][][] totalC_j_Minus_iOverAllObservations = new double[N+1][N+1][N+1];
	        
	        double[][] totalC_l_Minus_iOverAllObservations = new double[N+1][N+1];
	        
	        long start0 = System.currentTimeMillis();
	        
	        for (String[] f_e_wa : tritext_sd){
	        		
	        	String[] y = f_e_wa[0].split(" ");
	        	String[] x = f_e_wa[1].split(" "); 
	        	
	        	
	        	int T = y.length;
	        	N = x.length;
	        	HashMap<Integer, Double> c = new HashMap<Integer, Double>();
	        	
	        	if (iteration == 0)
	        		initializeModel(N);

	        	GenericPair<double[][], double[]> alphaHatC = forwardWithTScaled(a, pi, y, N, T, x, t_table);
	        	
	        	
	        	double[][] alpha_hat = alphaHatC.a;
	        	double[] c_scaled = alphaHatC.b;
	        	
	        	
	        	double[][] beta_hat = backwardWithTScaled(a, pi, y, N, T, x, t_table, c_scaled);
	        	
	        	
	        	
	        	double[][] gamma = new double[N+1][T+1];
	        	
	        	
	        	//Setting gamma
	        	for (int t = 1; t < T; t++){
	        		logLikelihood += -1*Math.log(c_scaled[t]);
	        		for (int i = 1; i < N+1; i++){
	        			gamma[i][t] = (alpha_hat[i][t]*beta_hat[i][t])/c_scaled[t];
	        			totalGammaDeltaOverAllObservations_t_i[indexMap.get(new Pair(y[t-1], x[i-1])).intValue()] += gamma[i][t];
	        		}
	        	}
	        	int t = T;
	        	logLikelihood += -1*Math.log(c_scaled[t]);
	        	for (int i = 1; i < N+1; i++){
        			gamma[i][t] = (alpha_hat[i][t]*beta_hat[i][t])/c_scaled[t];
        			
        			totalGammaDeltaOverAllObservations_t_i[indexMap.get(new Pair(y[t-1], x[i-1])).intValue()] += gamma[i][t];
        		}
	        	
	        	
	        	
	        	
	        	for (t = 1; t < T; t++)
	        		for (int i = 1; i < N+1; i++)
	        			for (int j = 1; j < N+1; j++){
	        				double xi = alpha_hat[i][t]*a[i][j][N]*t_table.get(new Pair(y[t],x[j-1]))*beta_hat[j][t+1];
	        				if (c.containsKey(j-i)){
        						double value = c.get(j-i);
        						c.put(j-i, value+xi);
        					}else{
        						c.put(j-i, xi);
        					}
	        			}
	        	
	        	for (int i = 1; i < N+1; i++){
	        		for (int j = 1; j < N+1; j++)
	        			totalC_j_Minus_iOverAllObservations[i][j][N] += (c.containsKey(j-i)? c.get(j-i) : 0);
	        		for (int l = 1; l < N+1; l++)
	        			totalC_l_Minus_iOverAllObservations[i][N] += (c.containsKey(l-i)? c.get(l-i) : 0);
	        	}
	        	
	        	
	        	for (int i = 1; i < N+1; i++)
	        		totalGamma1OverAllObservations[i] += gamma[i][1];

	        	
			
	        }//end of loop over bitext


	        
	        
	        System.out.println("log likelihood " + logLikelihood);
			N = totalGamma1OverAllObservations.length - 1;
			
				for (int k = 0; k < sd_size ; k++){
					Pair pair = biword.get(k);
					
					String e = pair.right;
					if (totalGammaDeltaOverAllObservations_t_overall_states_over_dest.containsKey(e)){
						double value = totalGammaDeltaOverAllObservations_t_overall_states_over_dest.get(e);
						totalGammaDeltaOverAllObservations_t_overall_states_over_dest.put(e, value + totalGammaDeltaOverAllObservations_t_i[k]);
					}else{
						totalGammaDeltaOverAllObservations_t_overall_states_over_dest.put(e, totalGammaDeltaOverAllObservations_t_i[k]);
					}
				}
			
			
			long end = System.currentTimeMillis();
			
			System.out.println("time spent in E-step: " + (end-start0)/1000.0 );
			
			int twoN = 2*N;
			
			//M-Step
			
			a = new double[twoN+1][twoN+1][N+1];
			pi = new double[twoN+1];
			t_table = new HashMap<Pair, Double>();
			
			
			//System.out.println("set " + targetLengthSet);
			for (int I : targetLengthSet){
				for (int i = 1; i < I+1; i++){
					for (int j = 1; j < I+1; j++)
						a[i][j][I] = totalC_j_Minus_iOverAllObservations[i][j][I]/totalC_l_Minus_iOverAllObservations[i][I];
				}
			}
			for (int i = 1; i < N+1; i++){
				//We can try uniform
				pi[i] = totalGamma1OverAllObservations[i]*(1.0/L);
			}
			
			for (int k = 0; k < sd_size ; k++){
				Pair pair = biword.get(k);
				String f = pair.left;
				String e = pair.right;
				t_table.put(new Pair(f, e), totalGammaDeltaOverAllObservations_t_i[k]/totalGammaDeltaOverAllObservations_t_overall_states_over_dest.get(e));
			}
			
			long end2 = System.currentTimeMillis();
			System.out.println("time spent in M-step: " + (end2-end)/1000.0 );
			System.out.println("iteration " + iteration);
			
		}	
		
		
	}
	/**
	 * This method checks the validity of the transition parameters (probabilities sum to 1)
	 * @param N
	 * @return True if a parameter is a probability and False otherwise
	 */
	public boolean checkProbability(int N){
		boolean isProb = true;
		for (int i = 1; i < a.length; i++){
			double sum = 0.0;
			for (int j = 1; j < a[0].length; j++){
				sum += a[i][j][N];
			}
			System.out.println("sum " + sum);
			if (sum != 1){
				System.out.println("not a probability");
				isProb = false;
			}else{
				System.out.println("Yessss");
			}
		}
		return isProb;
	}
	/**
	 * Smoothes t(f|e) parameter by backing-off to a uniform probability 1/V
	 * @param f  source word
	 * @param e  target word
	 * @return  the somoothed t(f|e)
	 */
	public double tProbability(String f, String e){
		Pair fe = new Pair(f,e);
		
		int V = 163303; //size of vocabulary
		if (t_table.containsKey(fe)){
			return t_table.get(fe);
			
		}
		else{
			if (e.equals("null"))
				return nullEmissionProb;
			else
				return 1.0/V;
		} 
	}
	/**
	 * Smoothes p(i|i',I) probabilities by backing-off to a uniform probability 1/I
	 * This is an alternative to Och and Ney, 2003 interpolation: \alpha*p(i|i',I) + (1-\alpha)*1/I.
	 * @param iPrime  previous target position
	 * @param i  current target position
	 * @param I  length of target sentence
	 * @return either p(i|i',I) or back-off to 1/I
	 */
	//p(i|i',I) is smoothed to uniform distribution for now --> p(i|i',I) = 1/I
	//we can make it interpolation form like what Och and Ney did
	public double aProbability(int iPrime, int i, int I){
		if (targetLengthSet.contains(I))
			return a[iPrime][i][I];
		else{
			return 1.0/I;
		}
	}
	/**
	 * Computes Viterbi decoding algorithm for HMM-based word alignment for a given sentence pair.
	 * For handling source words that have no alignment in the target sentence, we add N empty words (NULLs)
	 * to the target sentence and follow Och and Ney (2000a) approach. This is implemented in multiplyOneMinusP0H method and is called before
	 * the logViterbi method is called for each sentence pair. 
	 * @param N  target sentence length
	 * @param o  source sentence
	 * @param d  target sentence
	 * @return the best target positions (which could be NULL) for source positions 1 to N.
	 */
	public ArrayList<Integer> logViterbi(int N, String[] o, String[] d){
		//checkProbability(N);
		int twoN = 2*N;
		double[][] V = new double[twoN+1][o.length];
		int[][]ptr = new int[twoN+1][o.length];
		String[] newd = new String[2*d.length];
		int twoLend = 2*d.length;
		System.arraycopy(d, 0, newd, 0, d.length);
		for (int i = N; i < twoLend; i ++)
			newd[i] = "null";
		
		for (int q = 1 ; q < twoN + 1; q++){
				double t_o0_d_qMinus1 = tProbability(o[0],newd[q-1]);
				if (t_o0_d_qMinus1 == 0 || pi[q] == 0)
	                V[q][0] = Integer.MIN_VALUE;
	            else
	                V[q][0] = Math.log(pi[q])+Math.log(t_o0_d_qMinus1);
		}
		
		//printArray(V);
		for (int t = 1; t < o.length; t++){
			for (int q = 1 ; q < twoN + 1; q++){
				double maximum = Integer.MIN_VALUE;
				int max_q = Integer.MIN_VALUE;
				double t_o_d_qMinus1 = tProbability(o[t], newd[q-1]);
				for (int q_prime = 1 ; q_prime < twoN + 1; q_prime++){
					double a_q_prime_q_N = aProbability(q_prime, q, N);
					if (a_q_prime_q_N != 0 && t_o_d_qMinus1 != 0){
						double temp = V[q_prime][t-1]+Math.log(a_q_prime_q_N)+Math.log(t_o_d_qMinus1);
						if ( temp > maximum){
							maximum = temp;
							max_q = q_prime; 
						}
					}
				}
				V[q][t] = maximum;
			    ptr[q][t] = max_q;
			}
		}
		
		double max_of_V = Integer.MIN_VALUE;
		int q_of_max_of_V = 0;
		for (int q = 1 ; q < twoN + 1; q++){
			if (V[q][o.length-1] > max_of_V){
				max_of_V = V[q][o.length-1];
			    q_of_max_of_V = q;
			}
		}
		ArrayList<Integer> trace = new ArrayList<Integer>();
		trace.add(q_of_max_of_V);
		int q = q_of_max_of_V;
		int i = o.length-1;
		while (i > 0){
			q = ptr[q][i];
			trace.add(0, q);
			i = i - 1;
		}
		return trace;
	}
	/**
	 * Following Och and Ney (2000a) approach, for handling source words that have no alignment in the target sentence, we add N empty words (NULLs)
	 * to the target sentence and enforce a set of constraints which are implemented in this method. 
	 * @param bitext  source-target bitext
	 * @param p0H  probability of transition to NULL (refer to Och and Ney, (2000))
	 */
	public void multiplyOneMinusP0H(String[][] bitext, double p0H){
		
		for (int I : targetLengthSet)
		    for (int i = 1; i < I+1; i++)
				for (int j = 1; j < I+1; j++)
				    a[i][j][I] *= (1-p0H);
		for (int I : targetLengthSet){
		    for (int i = 1; i < I+1; i++){
				for (int j = 1; j < I+1; j++){ 
				    a[i][i + I][I] = p0H;
		            a[i + I][i + I][I] = p0H;
		            a[i+I][j][I] = a[i][j][I];
				}
		    }
		}
	}
	/**
	 * Finds the best alignments (source-target positions (1-1,2-3,3-5, ..., j-i, ..., N-x) for each sentence pair in the given bitext) 
	 * @param bitext  	 	test bitext for which decoding should be done.
	 * @param num_lines  	number of lines or sentences from the beginning of the bitext on which decoding should be done.
	 * @param alignmentFile output alignment file.
	 * @return 				a list of alignments. Each alignment is a string of j-i (source-target) alignment positions. 
	 * @throws IOException
	 */
	public ArrayList<String> findBestAlignmentsForAll_AER(String[][] bitext, int num_lines, String alignmentFile) throws IOException{
		BufferedWriter alignment = new BufferedWriter(new FileWriter(alignmentFile));
		ArrayList<String> alignmentList = new ArrayList<String>();
		int n = 0;
		for (String[] pair : bitext){
			String[] S = pair[0].trim().split(" ");
			String[] D = pair[1].trim().split(" ");
			
			int N = D.length;
			ArrayList<Integer> bestAlignment = logViterbi(N, S, D);
			String line = "";
			for(int i = 0; i < bestAlignment.size(); i++){
				
				if (bestAlignment.get(i) <= N)
					line += (i+1) +"-"+bestAlignment.get(i) +" ";
			}
			alignmentList.add(line);
			alignment.write(line+"\n");
			
			if (n == num_lines-1){
				alignment.close();
				return alignmentList;
			}
			n++;
		}
		alignment.close();
		return alignmentList;
	}
	public static void mainHMM(int trainingSize, int testSize, String trainPrefix
			, String sourceLang, String targetLang, String testPrefix, String referenceFile, String alignmentFile) throws IOException{
		IBM_Model1 ibm1 = new IBM_Model1();
	
		String[][] trainBitext = IBM_Model1.readBitext(trainPrefix+"."+sourceLang, trainPrefix+"."+targetLang, trainingSize);
		String[][] testBitext = IBM_Model1.readBitext(testPrefix+"."+sourceLang, testPrefix+"."+targetLang, testSize);
		
		
		ibm1.initializeCountsWithoutSets(trainBitext);
		System.out.println("Training an IBM1 model on the bitext of size " + trainBitext.length + ".");
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainBitext);
		HMMWithNullForTuningTe hmmL = new HMMWithNullForTuningTe(t_fe);
		System.out.println("Training an HMM model on the bitext of size " + trainBitext.length + ".");
		hmmL.baumWelch(trainBitext, ibm1.f_count, ibm1.fe_count);
		
		hmmL.multiplyOneMinusP0H(trainBitext, hmmL.p0H);
		
		ArrayList<String> hmmModelAlignment = hmmL.findBestAlignmentsForAll_AER(testBitext, testSize, alignmentFile);
		if (referenceFile != null){
			ArrayList<String> reference = ibm1.convertFileToArrayList(referenceFile);
			ibm1.gradeAlignmentWithType(testSize, testBitext, reference, hmmModelAlignment);
		}
	}
	
	
}
