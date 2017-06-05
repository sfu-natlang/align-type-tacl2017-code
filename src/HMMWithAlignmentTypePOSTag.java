import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

/**
 * Implements an HMM model with alignment type for POS-tagged data.
 * The POS tag model for the generative model in the following paper:
 * Joint prediction of word alignment with alignment types. 
 * @author anahita
 *
 */
public class HMMWithAlignmentTypePOSTag implements Serializable{
	private static final long serialVersionUID = 1L;
 
	String[][] bitext_fe;

	HashMap<String, Integer> e_count = new HashMap<String, Integer>();
	
	protected double[][][] a;
	
	protected HashMap<Pair, Double> t_table;
	
	protected double[] pi;

	protected int H = 11;  /**number of alignment types**/
	
	protected HashMap<Triple<String, String, Integer>, Double> s;
	
	protected int I;
	
	HashMap<String,Integer> tagMap = new HashMap<String, Integer>();
	
	
	HashSet<Integer> targetLengthSet;
	
	HashMap<String, Double> count_e;
	double p0H = 0.3;
	double nullEmissionProb = 0.000005;
	double smoothFactor = 0.1;	
	
	double lambda1;
	double lambda2;
	double lambda3;
	
	
	double lambda = 1-1e-20;
	
	public void setLambda(double lambda){
		this.lambda = lambda;
	}
	
	public HMMWithAlignmentTypePOSTag(HashMap<Pair, Double> t_table) {
		this.t_table = t_table;
	}
	/**
	 * 
	 * @param t_table  			translation table
	 * @param e_count  			target words count table
	 * @param total_f_e_h  		total triples source word, target word, alignment type count table.
	 * @param fe_count  		source-target word count table
	 * @param fe_count_20K		source-target count table for the 20K data
	 */
	public HMMWithAlignmentTypePOSTag(HashMap<Pair, Double> t_table, HashMap<String, Integer> e_count, HashMap<Triple<String, String, Integer>, Double> total_f_e_h, HashMap<Pair, Integer> fe_count, HashMap<Pair, Integer> fe_count_20K){
		this.t_table = t_table;
		this.e_count = e_count;
		
		s = new HashMap<Triple<String,String,Integer>, Double>();
		for (Entry<Triple<String, String, Integer>, Double> entry : total_f_e_h.entrySet()){
			Triple<String , String, Integer> f_e_h = entry.getKey();
			String f = f_e_h.left;
			String e = f_e_h.middle;
			double value = entry.getValue();
			
			s.put(f_e_h, value/fe_count_20K.get(new Pair(f,e)));
		}
	}
	
	public void setLambdas(double lambda1,double lambda2,double lambda3){
		this.lambda1 = lambda1;
		this.lambda2 = lambda2;
		this.lambda3 = lambda3;
	}
	
	/**
	 * Computes forward algorithm for HMM. Note that we re-scale the probabilities to avoid 0 probabilities for long sentences.  
	 * References:
	 * (1) wikipedia Baum Welch algorithm 
	 * (2) Mathematics of HMM and
	 * @param a   HMM transition parameters
	 * @param pi  HMM start parameters
	 * @param y   source sentence
	 * @param N   length of the target sentence
	 * @param T   length of the source sentence
	 * @param d   target sentence
	 * @param t_table  translation table
	 * @return 			a pair of alphaHat probabilities and c_hat (for alphaHat c_hat refer to reference (2))
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
	 * @param a HMM transition parameters
	 * @param pi HMM start parameters
	 * @param y source sentence
	 * @param N length of the target sentence
	 * @param T length of the source sentence
	 * @param d target sentence
	 * @param t_table translation table
	 * @param c_scaled c scaler (computed using forwardWithTScaled method) 
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
	 * Also, tagMap which is a map of alignment types to unique numbers is initialized. 
	 * @param N length of the target sentence
	 */
	public void  initializeModel(int N){
		int twoN = 2*N;
		
		for (int i = 1; i < twoN+1; i++)
			pi[i] = 1.0/twoN;
		for (int i = 1; i < N+1; i++)
			for (int j = 1; j < N + 1; j++)
				a[i][j][N] = 1.0/N;
		
		tagMap.put("SEM", 1);
		tagMap.put("FUN", 2);
		tagMap.put("PDE", 3);
		tagMap.put("CDE", 4);
		tagMap.put("MDE", 5);
		tagMap.put("GIS", 6);
		tagMap.put("GIF", 7);
		tagMap.put("COI", 8);
		tagMap.put("TIN", 9);
		tagMap.put("NTR", 10);
		tagMap.put("MTA", 11);
		
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
	 * @param sd_count source-target count map 
	 * @return 		   the reversed bitext
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
	
	public static String[] strip(String englishWord){
		String indices = "";
		int wordLength = englishWord.length();
		for (int i = 0; i < wordLength; i++){
			char e_i = englishWord.charAt(i);
			if (Character.isDigit(e_i) || e_i == ','){
				indices += e_i;
			}
		}
		return indices.split(",");
	}
	/**
	 * set a to 0
	 */
	 public void clearA(){
         for (int i = 0; i < I+1; i++){
             for (int j = 0; j < I+1; j++){
                 a[i][j][I] = 0.0;
             }
         }
	 }
	/**
	 * 
	 * @param tritext_sd tritext of train data 
	 * * 	
	 * @param s_count  source words count table
	 * @param sd_count source-target words count table
	 * @throws IOException
	 */
	public void baumWelch(String[][] tritext_sd, HashMap<String, Integer>s_count,HashMap<Pair, Integer>sd_count) {
		
		
		GenericPair<Integer, HashSet<Integer>> NandtargetLengthSet = maxTargetSentenceLength(tritext_sd);
		
		int N = NandtargetLengthSet.a;
		targetLengthSet = NandtargetLengthSet.b;
		
		I = N;
		
		//System.out.println("N "+N);
		GenericPair<HashMap<Pair, Integer>, HashMap<Integer, Pair>> indexBiwordPair = mapBitextToInt(sd_count);
		HashMap<Pair, Integer> indexMap = indexBiwordPair.a;
		HashMap<Integer, Pair> biword = indexBiwordPair.b;
		
		int L = tritext_sd.length;
		int sd_size = indexMap.size();
		
		double[] totalGammaDeltaOverAllObservations_t_i = null;
		HashMap<String, Double>totalGammaDeltaOverAllObservations_t_overall_states_over_dest = null;
		
		int twoN = 2*N;
		a = new double[twoN+1][twoN+1][N+1];
		
		pi = new double[twoN+1];

		
		for (int iteration = 0; iteration < 5; iteration++){
			
			double logLikelihood = 0;
		
			totalGammaDeltaOverAllObservations_t_i = new double[sd_size];
		
	        totalGammaDeltaOverAllObservations_t_overall_states_over_dest = new HashMap<String, Double>();

	        double[] totalGamma1OverAllObservations = new double[N+1];
	        double[][][] totalC_j_Minus_iOverAllObservations = new double[N+1][N+1][N+1];
	    
	        double[][] totalC_l_Minus_iOverAllObservations = new double[N+1][N+1];
	        HashMap<Triple<String, String, Integer>, Double> total_f_e_h = new HashMap<Triple<String,String,Integer>, Double>();
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
	        	
	        	double[][][] xi = new double[N+1][N+1][T+1];
	        	
	        	
	        	
	        	//Setting gamma
	        	for (int t = 1; t < T; t++){
	        		logLikelihood += -1*Math.log(c_scaled[t]);
	        		for (int i = 1; i < N+1; i++){
	        			gamma[i][t] = (alpha_hat[i][t]*beta_hat[i][t])/c_scaled[t];
	        			
	        			totalGammaDeltaOverAllObservations_t_i[indexMap.get(new Pair(y[t-1], x[i-1])).intValue()] += gamma[i][t];
	        			for (int h = 1; h < H+1 ; h++){
		        			Triple<String, String, Integer> f_e_h = new Triple<String, String, Integer>(y[t-1], x[i-1], h);
		        			if (total_f_e_h.containsKey(f_e_h)){
		        				Double value = total_f_e_h.get(f_e_h);
		        				
		        				total_f_e_h.put(f_e_h, value+gamma[i][t]*sProbability(y[t-1], x[i-1], h));
		        				
		        			}else{
		        				total_f_e_h.put(f_e_h, gamma[i][t]*sProbability(y[t-1], x[i-1],h));
		        			}
	        			}
	        		}
	        	}
	        	int t = T;
	        	logLikelihood += -1*Math.log(c_scaled[t]);
	        	for (int i = 1; i < N+1; i++){
        			gamma[i][t] = (alpha_hat[i][t]*beta_hat[i][t])/c_scaled[t];
        			
        			totalGammaDeltaOverAllObservations_t_i[indexMap.get(new Pair(y[t-1], x[i-1])).intValue()] += gamma[i][t];
        			for (int h = 1; h < H+1 ; h++){
	        			Triple<String, String, Integer> f_e_h = new Triple<String, String, Integer>(y[t-1], x[i-1], h);
	        			if (total_f_e_h.containsKey(f_e_h)){
	        				Double value = total_f_e_h.get(f_e_h);
	        		
	        				total_f_e_h.put(f_e_h, value+gamma[i][t]*sProbability(y[t-1], x[i-1], h));
	        				
	        			}else{
	        				total_f_e_h.put(f_e_h, gamma[i][t]*sProbability(y[t-1], x[i-1],h));
	        			}
        			}
        		}
	        	
	        	
	        	//Setting xi
	        	
	        	for (t = 1; t < T; t++)
	        		for (int i = 1; i < N+1; i++)
	        			for (int j = 1; j < N+1; j++)
	        				xi[i][j][t] = alpha_hat[i][t]*a[i][j][N]*t_table.get(new Pair(y[t],x[j-1]))*beta_hat[j][t+1];
	        	
	        	
	        	
	        	
	        	for (int i = 1; i < N+1; i++)
	        		totalGamma1OverAllObservations[i] += gamma[i][1];
	        	for (int d = -N-1; d < N + 1; d++)
	        		for (t = 1; t < T+1; t++) //T+1 nabood alan taghir dadam
	        			for (int i = 1; i < N+1; i++)
	        				if (i+d <= N && i+d >= 1){
	        					if (c.containsKey(d)){
	        						double value = c.get(d);
	        						c.put(d, value+xi[i][i+d][t]);
	        					}else{
	        						c.put(d, xi[i][i+d][t]);
	        					}
	        				}
	        	for (int i = 1; i < N+1; i++){
	        		for (int j = 1; j < N+1; j++)
	        			totalC_j_Minus_iOverAllObservations[i][j][N] += c.get(j-i);
	        		for (int l = 1; l < N+1; l++)
	        			totalC_l_Minus_iOverAllObservations[i][N] += c.get(l-i);
	        	}
	        	
	        	
	        		
	        	//printArray(totalC_j_Minus_iOverAllObservations);
	        }//end of loop over bitext
	        
	        
	        
	        System.out.println("log likelihood " + logLikelihood);
			N = totalGamma1OverAllObservations.length - 1;
			
				for (int k = 0; k < sd_size ; k++){
					Pair pair = biword.get(k);
					//String f = pair.left;
					String e = pair.right;
					if (totalGammaDeltaOverAllObservations_t_overall_states_over_dest.containsKey(e)){
						double value = totalGammaDeltaOverAllObservations_t_overall_states_over_dest.get(e);
						totalGammaDeltaOverAllObservations_t_overall_states_over_dest.put(e, value + totalGammaDeltaOverAllObservations_t_i[k]);
					}else{
						totalGammaDeltaOverAllObservations_t_overall_states_over_dest.put(e, totalGammaDeltaOverAllObservations_t_i[k]);
					}
				}
		
			long end = System.currentTimeMillis();
			//System.out.println("time spent in the end of E-step: " + (end-start)/1000.0 );
			System.out.println("time spent in E-step: " + (end-start0)/1000.0 );
			
			twoN = 2*N;
			
			//M-Step
			
			a = new double[twoN+1][twoN+1][N+1];
			pi = new double[twoN+1];
			t_table = new HashMap<Pair, Double>();
			s = new HashMap<Triple<String,String,Integer>, Double>();
			
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
			for (Entry<Triple<String , String, Integer>, Double> entry : total_f_e_h.entrySet()){
				Triple<String , String, Integer> f_e_h = entry.getKey();
				String f = f_e_h.left;
				String e = f_e_h.middle;
				double value = entry.getValue();
				s.put(f_e_h, value/totalGammaDeltaOverAllObservations_t_i[indexMap.get(new Pair(f,e))]);
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
	 * @param f source word
	 * @param e target word
	 * @return  the somoothed t(f|e)
	 */
	public double tProbability(String f, String e){
		Pair fe = new Pair(f,e);
		
		int V = 163303;
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
	 * @param iPrime previous target position
	 * @param i 	 current target position
	 * @param I		 length of target sentence
	 * @return  	 either p(i|i',I) or back-off to 1/I
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
	 * This method is used in EM for smoothing p(h|t_f,t_e) parameters. We used to call these s parameters. 
	 * @param t_f	 POS tag of the source word f
	 * @param t_e	 POS tag of the target word e
	 * @param h		 alignment type
	 * @return 		 smoothed p(i|i',I)
	 */
	
	public double sProbability(String t_f, String t_e, Integer h){
		Triple<String, String, Integer> f_e_h = new Triple<String, String, Integer>(t_f,t_e,h);
		double first = (s.containsKey(f_e_h) ? s.get(f_e_h) : 0 );
		double[] tagDist = {0, 0.401, 0.264, 0.004, 0.004, 0.012, 0.205, 0.031, 0.008, 0.003, 0.086, 0.002};
		
		return lambda*first + (1-lambda)*tagDist[h];
			
	}
	/**
	 * Following Och and Ney (2000a) approach, for handling source words that have no alignment in the target sentence, we add N empty words (NULLs)
	 * to the target sentence and enforce a set of constraints which are implemented in this method. 
	 */
	public void multiplyOneMinusP0H(){
		
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
		
}
