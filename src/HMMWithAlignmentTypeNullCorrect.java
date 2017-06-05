import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;


public class HMMWithAlignmentTypeNullCorrect extends HMMWithAlignmentTypePOSTag implements Serializable{
	
	private static final long serialVersionUID = -2641697232475893704L;
	
	HashMap<Pair, Double> t_table_20K;
	
	//POS-tag model p(h|t_f,t_e) which comes from the POS-tagged data
	HashMap<Triple<String, String, Integer>, Double> tagModel_s;
	
	/**
	 * Constructor
	 * @param t_fe t 		parameters of the IBM1
	 * @param e_count 		target words count table (also comes from IBM1)
	 * @param total_f_e_h 	(source word, target word, alignment type) count table which is initialized in IBM1
	 * @param fe_count 		(source word, target word) count table 
	 * @param fe_count_20K 	(source word, target word) count table for a smaller model (eg when we want to first train on a small model such as 20k in the paper and then on a larger one)
	 * @param tagModel_s	p(h|t_f,t_e) initialization
	 */
	public HMMWithAlignmentTypeNullCorrect(HashMap<Pair, Double> t_fe,
			HashMap<String, Integer> e_count,
			HashMap<Triple<String, String, Integer>, Double> total_f_e_h,
			HashMap<Pair, Integer> fe_count, HashMap<Pair, Integer> fe_count_20K,
			HashMap<Triple<String, String, Integer>, Double> tagModel_s) {
		super(t_fe, e_count, total_f_e_h, fe_count, fe_count_20K);
		this.tagModel_s = tagModel_s;
		
	}
	/**
	 * sets the t_table_20K to the t parameters of the 20K data
	 * @param t_table_20K	t parameters of the 20K data
	 */
	public void setTTable20K(HashMap<Pair, Double> t_table_20K){
		this.t_table_20K = t_table_20K;
	}
	/**
	 * 
	 * @param tritext_sd tritext of train data 
	 * * 	
	 * @param s_count 	source words count table
	 * @param sd_count  source-target words count table	 
	 */
	public void baumWelch(String[][] tritext_sd, String[][] tritext_sd_tag, HashMap<String, Integer>s_count,HashMap<Pair, Integer>sd_count) {
		
		
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
		
		double[] totalGammaDeltaOverAllObservations_t_i = new double[sd_size];
		HashMap<String, Double>totalGammaDeltaOverAllObservations_t_overall_states_over_dest = null;
		Pair currentPair = new Pair("", "");
		
		int twoN = 2*N;
		a = new double[twoN+1][twoN+1][N+1];
		
		pi = new double[twoN+1];
		
		HashMap<Integer, Double> c = new HashMap<Integer, Double>();
		HashMap<Triple<String, String, Integer>, Double> total_f_e_h = new HashMap<Triple<String,String,Integer>, Double>();
		for (int iteration = 0; iteration < 5; iteration++){
			
			double logLikelihood = 0;
			
			Arrays.fill(totalGammaDeltaOverAllObservations_t_i, 0, sd_size, 0.0);
			
	        totalGammaDeltaOverAllObservations_t_overall_states_over_dest = new HashMap<String, Double>();

	        double[] totalGamma1OverAllObservations = new double[N+1];
	        double[][][] totalC_j_Minus_iOverAllObservations = new double[N+1][N+1][N+1];
	        
	        double[][] totalC_l_Minus_iOverAllObservations = new double[N+1][N+1];
	        
	        total_f_e_h.clear();
	        
	        long start0 = System.currentTimeMillis();
	        int sentenceNumber = 0;
	        for (String[] f_e_wa : tritext_sd){
	        	
	        	String[] y = f_e_wa[0].trim().split(" ");
	        	String[] x = f_e_wa[1].trim().split(" "); 
	        	
	        	
	        	String[] y_t = tritext_sd_tag[sentenceNumber][0].trim().split(" ");
	        	String[] x_t = tritext_sd_tag[sentenceNumber][1].trim().split(" ");
	        	int T = y.length;
	        	N = x.length;
	        	
	        	c.clear();
	        	
	        	if (iteration == 0)
	        		initializeModel(N);
	        	
	        	GenericPair<double[][], double[]> alphaHatC = forwardWithTScaled(a, pi, y, N, T, x, t_table);
	        	
	        	
	        	//long start = end;
	        	
	        	double[][] alpha_hat = alphaHatC.a;
	        	double[] c_scaled = alphaHatC.b;
	        	
	        	
	        	double[][] beta_hat = backwardWithTScaled(a, pi, y, N, T, x, t_table, c_scaled);
	        	
	        	
	        	double[][] gamma = new double[N+1][T+1];
	        	
	        	
	        	//Setting gamma
	        	for (int t = 1; t < T; t++){
	        		logLikelihood += -1*Math.log(c_scaled[t]);
	        		for (int i = 1; i < N+1; i++){
	        			gamma[i][t] = (alpha_hat[i][t]*beta_hat[i][t])/c_scaled[t];
	        			
	        			currentPair.setPair(y[t-1], x[i-1]);
	        			totalGammaDeltaOverAllObservations_t_i[indexMap.get(currentPair).intValue()] += gamma[i][t];
	        			for (int h = 1; h < H+1 ; h++){
		        			Triple<String, String, Integer> f_e_h = new Triple<String, String, Integer>(y[t-1], x[i-1], h);
		        			if (total_f_e_h.containsKey(f_e_h)){
		        				Double value = total_f_e_h.get(f_e_h);
		        				//we can precompute sProbability
		        				total_f_e_h.put(f_e_h, value+gamma[i][t]*sProbability(y[t-1], x[i-1], h, y_t[t-1], x_t[i-1]));
		        				
		        			}else{
		        				total_f_e_h.put(f_e_h, gamma[i][t]*sProbability(y[t-1], x[i-1], h, y_t[t-1], x_t[i-1]));
		        			}
	        			}
	        		}
	        	}
	        	int t = T;
	        	logLikelihood += -1*Math.log(c_scaled[t]);
	        	for (int i = 1; i < N+1; i++){
        			gamma[i][t] = (alpha_hat[i][t]*beta_hat[i][t])/c_scaled[t];
        			currentPair.setPair(y[t-1], x[i-1]);
        			
        			totalGammaDeltaOverAllObservations_t_i[indexMap.get(currentPair).intValue()] += gamma[i][t];
        			for (int h = 1; h < H+1 ; h++){
	        			Triple<String, String, Integer> f_e_h = new Triple<String, String, Integer>(y[t-1], x[i-1], h);
	        			if (total_f_e_h.containsKey(f_e_h)){
	        				Double value = total_f_e_h.get(f_e_h);
	        				//we can precompute sProbability
	        				total_f_e_h.put(f_e_h, value+gamma[i][t]*sProbability(y[t-1], x[i-1], h, y_t[t-1], x_t[i-1]));
	        				
	        			}else{
	        				total_f_e_h.put(f_e_h, gamma[i][t]*sProbability(y[t-1], x[i-1],h, y_t[t-1], x_t[i-1]));
	        			}
        			}
        		}
	        	
	        	
	        	//Setting xi
	        	
	        	for (t = 1; t < T; t++)
	        		for (int i = 1; i < N+1; i++)
	        			for (int j = 1; j < N+1; j++){
	        				currentPair.setPair(y[t],x[j-1]);
						//xi[i][j][t] = alpha_hat[i][t]*a[i][j][N]*t_table.get(new Pair(y[t],x[j-1]))*beta_hat[j][t+1];
	        				double xi = alpha_hat[i][t]*a[i][j][N]*t_table.get(currentPair)*beta_hat[j][t+1];
	        				if (c.containsKey(j-i)){
        						double value = c.get(j-i);
        						c.put(j-i, value+xi);
        					}else{
        						c.put(j-i, xi);
        					}
	        			}
	        	
	        	//System.out.println("time spent to set xi: " + (end-start)/1000.0 );
	        	
	        	
	        	for (int i = 1; i < N+1; i++)
	        		totalGamma1OverAllObservations[i] += gamma[i][1];
	        	
	        	for (int i = 1; i < N+1; i++){
	        		for (int j = 1; j < N+1; j++){
					
	        			totalC_j_Minus_iOverAllObservations[i][j][N] += (c.containsKey(j-i)? c.get(j-i) : 0);
				}
	        		for (int l = 1; l < N+1; l++)
	        			totalC_l_Minus_iOverAllObservations[i][N] += (c.containsKey(l-i)? c.get(l-i) : 0);
	        	}
	        	
	        	//System.out.println("time spent to set xi: " + (end-start)/1000.0 );
	      
	        		        	
				
	        	sentenceNumber++;
	        	
	        }//end of loop over bitext
	        long end = System.currentTimeMillis();
	        
	        

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
			
			
			end = System.currentTimeMillis();
			
			System.out.println("time spent in E-step: " + (end-start0)/1000.0 );
			
			twoN = 2*N;
			
			//System.out.println("lambdas " + lambda1 + " " + lambda2 + " " + lambda3);
			
			//M-Step
			clearA();
			Arrays.fill(pi, 0, twoN+1, 0.0);
			t_table.clear();
			s.clear();
			
			
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
				String e = pair.right;
				t_table.put(pair, totalGammaDeltaOverAllObservations_t_i[k]/totalGammaDeltaOverAllObservations_t_overall_states_over_dest.get(e));
			}
			for (Entry<Triple<String , String, Integer>, Double> entry : total_f_e_h.entrySet()){
				Triple<String , String, Integer> f_e_h = entry.getKey();
				String f = f_e_h.left;
				String e = f_e_h.middle;
				currentPair.setPair(f, e);
				double value = entry.getValue();
				s.put(f_e_h, value/totalGammaDeltaOverAllObservations_t_i[indexMap.get(currentPair)]);
			}
			
			long end2 = System.currentTimeMillis();
			//System.out.println("s size " + s.size());
			System.out.println("time spent in M-step: " + (end2-end)/1000.0 );
			System.out.println("iteration " + iteration);
		}	
		//This is added to write count(f,e) object in the file (for further easy read)
		
	}
	/**
	 * Smoothes t(f|e) parameter by backing-off to a uniform probability 1/V
	 * @param f 	source word
	 * @param e	 	target word
	 * @return  	the somoothed t(f|e)
	 */
	public double tProbability(String f, String e){
		Pair fe = new Pair(f,e);
		
		int V = 163303;
		
		if (t_table_20K!= null && t_table_20K.containsKey(fe))
			return t_table_20K.get(fe);
 		else if (t_table.containsKey(fe)){
			
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
	 * computes p(h|f,e) using the following linear interpolation:
	 * p(h|f,e) = lambda_1*p(h|f,e) + lambda_2*p(h|t_f,t_e) + lambda3*p(h)
	 * @param f	 	source word f
	 * @param e	 	target word e
	 * @param h	 	alignment type h
	 * @param t_f	POS tag of f
	 * @param t_e   POS tag of e
	 * @return 		p(h|f,e) 
	 */
	public double sProbability(String f, String e, Integer h, String t_f, String t_e){
		double[] tagDist = {0, 0.401, 0.264, 0.004, 0.004, 0.012, 0.205, 0.031, 0.008, 0.003, 0.086, 0.002};
		Triple<String, String, Integer> f_e_h = new Triple<String, String, Integer>(f,e,h);

		double first = (s.containsKey(f_e_h) ? s.get(f_e_h) : 0 );
		
		Triple<String, String, Integer> taggedFeh = new Triple<String, String, Integer>(t_f,t_e,h);
		
		double second = (tagModel_s.containsKey(taggedFeh) ? tagModel_s.get(taggedFeh) : 0 );
		
		
		first = lambda*first + (1-lambda)*tagDist[h];
		second = lambda*second + (1-lambda)*tagDist[h];
		
		return lambda1*first + lambda2*second + lambda3*tagDist[h];	
	}
	
	/**
	 * 
	 * Computes Viterbi decoding algorithm for HMM-based word alignment with alignment types for a given sentence pair.
	 * For handling source words that have no alignment in the target sentence, we add N empty words (NULLs)
	 * to the target sentence and follow Och and Ney (2000a) approach. This is implemented in multiplyOneMinusP0H method and is called before.
	 * the logViterbi method is called for each sentence pair. 
	 * For more information on how to find the best alignment and alignment type, see the paper.
	 * @param N	 	target sentence length
	 * @param o		source sentence
	 * @param d		target sentence
	 * @param t_o	POS-tagged source sentence
	 * @param t_d  	POS-tagged target sentence
	 * @return the best alignment, as well as best alignment types for source positions 1 to N.
	 */
	public GenericPair<ArrayList<Integer>, ArrayList<Integer>> logViterbi(int N, String[] o, String[] d,String[] t_o, String[] t_d){

		int twoN = 2*N;

		double[][][] V = new double[twoN+1][o.length][H+1];
		int[][][]ptr = new int[twoN+1][o.length][H+1];
		int[][][]ptr_h = new int[twoN+1][o.length][H+1];
		double[] tagDist = {0, 0.401, 0.264, 0.004, 0.004, 0.012, 0.205, 0.031, 0.008, 0.003, 0.086, 0.002};
		
		String[] newd = new String[2*d.length];
		String[] new_t_d = new String[2*t_d.length];
		
		int twoLend = 2*d.length;
		if (t_d.length != d.length){
			System.out.println("length of array " + t_d.length + " " + d.length);
		}
		System.arraycopy(d, 0, newd, 0, d.length);
		System.arraycopy(t_d, 0, new_t_d, 0, t_d.length);
		
		for (int i = N; i < twoLend; i ++){
			newd[i] = "null";
			new_t_d[i] = "null";
		}
		
		
		for (int q = 1 ; q < twoN + 1; q++){
			double t_o0_d_qMinus1 = tProbability(o[0],newd[q-1]);
			for (int h = 1; h < H + 1 ; h++){
							
            	Triple<String, String, Integer> feh = new Triple<String, String, Integer>(o[0],newd[q-1],h);
            	//write the tag of the data (add input parameters ...)
            	Triple<String, String, Integer> taggedFeh = new Triple<String, String, Integer>(t_o[0],new_t_d[q-1],h);
            	double first = (s.containsKey(feh) ? s.get(feh) : 0 );
            	
            	double second = (tagModel_s.containsKey(taggedFeh) ? tagModel_s.get(taggedFeh) : 0 );
            	
            	double s; 
            	first = first*(1-1e-20) + 1e-20*tagDist[h];
				second = second*(1-1e-20)+1e-20*tagDist[h];

				//July 10th
				s = lambda1*first + lambda2*second + lambda3*tagDist[h];
	            double temp_s = Math.log(s);
                V[q][0][h] = Math.log(pi[q])+ Math.log(t_o0_d_qMinus1) +temp_s;

			}

		}
		
	
		for (int t = 1; t < o.length; t++){
			for (int q = 1 ; q < twoN + 1; q++){
				double maximum = Integer.MIN_VALUE;
				int max_q = Integer.MIN_VALUE;
				int best_h = 0;
				double t_o_d_qMinus1 = tProbability(o[t], newd[q-1]);
				for (int q_prime = 1 ; q_prime < twoN + 1; q_prime++){
					double a_q_prime_q_N = aProbability(q_prime, q, N);
					if (a_q_prime_q_N != 0 && t_o_d_qMinus1 != 0){
						for (int hPrime = 1; hPrime < H + 1 ; hPrime++){
							double temp = V[q_prime][t-1][hPrime]+Math.log(a_q_prime_q_N)+Math.log(t_o_d_qMinus1);
							if ( temp > maximum){
								maximum = temp;
								max_q = q_prime; 
								best_h = hPrime;
							}
						}
					}
				}
				
            	for (int h = 1; h < H + 1 ; h++){
	            	Triple<String, String, Integer> feh = new Triple<String, String, Integer>(o[t],newd[q-1],h);
	            	//write the tag of the data (add input parameters ...)
	            	Triple<String, String, Integer> taggedFeh = new Triple<String, String, Integer>(t_o[t],new_t_d[q-1],h);
	            	double first = (s.containsKey(feh) ? s.get(feh) : 0 );
	            	double second = (tagModel_s.containsKey(taggedFeh) ? tagModel_s.get(taggedFeh) : 0 );
	            	
	            	
	            	double s; 
	            	first = first*(1-1e-20) + 1e-20*tagDist[h];
					second = second*(1-1e-20)+1e-20*tagDist[h];
					
					s = lambda1*first + lambda2*second + lambda3*tagDist[h];
					
	            	double temp_s = Math.log(s);
				    V[q][t][h] = maximum + temp_s;
				    ptr[q][t][h] = max_q;
				    ptr_h[q][t][h] = best_h;
			    
				}
			}
		}
		
		double max_of_V = Integer.MIN_VALUE;
		int q_of_max_of_V = 0;
		int h_of_max_of_V = 0;
		for (int q = 1 ; q < twoN + 1; q++){
			for (int h = 1; h < H + 1; h++){
				if (V[q][o.length-1][h] > max_of_V){
					max_of_V = V[q][o.length-1][h];
				    q_of_max_of_V = q;
				    h_of_max_of_V = h;
				}
			}
		}
		ArrayList<Integer> trace = new ArrayList<Integer>();
		ArrayList<Integer> bestLinkTrace = new ArrayList<Integer>();
		
		trace.add(q_of_max_of_V);
		bestLinkTrace.add(h_of_max_of_V);
		
		int q = q_of_max_of_V;
		int i = o.length-1;
		int h = h_of_max_of_V;

		while (i > 0){
			int temp_q = q;
			int temp_h = h;
			
			h = ptr_h[temp_q][i][temp_h];
			q = ptr[temp_q][i][temp_h];
			trace.add(0, q);
			bestLinkTrace.add(0,h);
			i = i - 1;
		}
		return new GenericPair<ArrayList<Integer>, ArrayList<Integer>>(trace, bestLinkTrace);
	}
	/**
	 * * Jointly finds the best alignment and alignment types (source-target positions and types (1-1(SEM),2-3(FUN),3-5(SEM), ..., j-i(h), ..., N-x(GIS)) for each sentence pair in the given bitext) 
	 * @param bitext		 	test bitext for which decoding should be done.
	 * @param tagModelBitext 	POS-tagged bitext
	 * @param num_lines		 	number of lines or sentences from the beginning of the bitext on which decoding should be done.
	 * @param alignmentFile	 	output alignment file.
	 * @param useAlignmentType 	set to True for the joint task and False for the word alignment task
	 * @return 					a list of alignments (if useAlignmentType is True) or a list of alignment and alignment types (if useAlignmentType is False). 
	 * For the former, each alignment is a string of j-i (source-target) alignment positions. For the latter, source-target positions and their corresponding types j-i(h) are computed for each source position in each sentence pair. 
	 * @throws IOException
	 */
	public ArrayList<String> findBestAlignmentsForAll_AER(String[][] bitext,String[][] tagModelBitext, int num_lines, String alignmentFile, boolean useAlignmentType) throws IOException{
		
		String[] linkMap = {"SEM", "FUN", "PDE", "CDE", "MDE", "GIS", "GIF", "COI", "TIN", "NTR", "MTA"};
		BufferedWriter alignment = new BufferedWriter(new FileWriter(alignmentFile));
		ArrayList<String> alignmentList = new ArrayList<String>();
		int n = 0;
		
		while (n <= num_lines){
			
			//System.out.println("VITERBI: line " + (n+1));
			String[] S = bitext[n][0].trim().split(" ");
			String[] D = bitext[n][1].trim().split(" ");
			
			//Source and target sentence in tagModel
			String[] taggedS = tagModelBitext[n][0].trim().split(" ");
			String[] taggedD = tagModelBitext[n][1].trim().split(" ");
			int N = D.length;
			
			GenericPair<ArrayList<Integer>,ArrayList<Integer>> bestAlignmentAndLink = logViterbi(N, S, D, taggedS, taggedD );
			ArrayList<Integer> bestAlignment = bestAlignmentAndLink.a;
			ArrayList<Integer> bestLink = bestAlignmentAndLink.b;
			
			String line = "";
			if (!useAlignmentType){
				for(int i = 0; i < bestAlignment.size(); i++){
					if (bestAlignment.get(i) <= N)
						line += (i+1) +"-"+bestAlignment.get(i) +" ";
				}
			}else{
				for(int i = 0; i < bestAlignment.size(); i++){
					if (bestAlignment.get(i) <= N){
						if (bestLink.get(i) == 0)
							line += (i+1) +"-"+bestAlignment.get(i) + "(" + linkMap[bestLink.get(i)] + ")" + " ";
						else
							line += (i+1) +"-"+bestAlignment.get(i) + "(" + linkMap[bestLink.get(i)-1] + ")" + " ";
					}
				}
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
	
	public static void augmentedModelExperimentMain() throws IOException, ClassNotFoundException{
		
		int trainingSize = 19987;//20357;//20357;//18000;
		
/**
 * For computing BLEU Score testTritexts are not important. We can even use nulls for them		
 */
		String[][] trainTritextTaggedModel = IBM_Model1.readTritext("data-jasneet/train.20K.tag.clean.cn","data-jasneet/train.20K.tag.clean.en","data-jasneet/train.20K.seg.cln.clean.wa",trainingSize);
		//String[][] testTritextTaggedModel = null;//IBM_Model1.readTritext("data-te/test.tags.cn","data-te/test.tags.en","data-te/test.seg.cln.wa",testSize);
		
		String[][] trainTritext = IBM_Model1.readTritext("data-jasneet/train.20K.clean.cn","data-jasneet/train.20K.clean.en","data-jasneet/train.20K.seg.cln.clean.wa",trainingSize);
		//String[][] testTritext = null;//IBM_Model1.readTritext("data-te/test.seg.cln.cn","data-te/test.seg.cln.en","data-te/test.seg.cln.wa",testSize);

		double lambda1 = 0.9999999999;
		double lambda2 = 9.999900827395436E-11;
		double lambda3 = 1.000000082740371E-15;
		
		IBM_Model1 ibm1TaggedModel = new IBM_Model1();
		IBM_Model1 ibm1 = new IBM_Model1();
		ibm1TaggedModel.initializeCounts(trainTritextTaggedModel, 0);
        HashMap<Pair, Double> t_fe_tagged = ibm1TaggedModel.EM_IBM1(ibm1TaggedModel.f_count, ibm1TaggedModel.e_count, ibm1TaggedModel.fe_count, trainTritextTaggedModel);
  
        HMMWithAlignmentTypePOSTag hmmTaggedModel = new HMMWithAlignmentTypePOSTag(t_fe_tagged, ibm1TaggedModel.e_count, ibm1TaggedModel.total_f_e_h,ibm1TaggedModel.fe_count,ibm1TaggedModel.fe_count);
        hmmTaggedModel.baumWelch(trainTritextTaggedModel, ibm1TaggedModel.f_count, ibm1TaggedModel.fe_count);
        hmmTaggedModel.multiplyOneMinusP0H();
        
        
        ibm1.initializeCounts(trainTritext, 0);
		
		System.out.println("length " + trainTritext.length);
		
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainTritext);
		
		HMMWithAlignmentTypeNullCorrect hmm20K = new HMMWithAlignmentTypeNullCorrect(t_fe, ibm1.e_count, ibm1.total_f_e_h,ibm1.fe_count, ibm1.fe_count,hmmTaggedModel.s);
		hmm20K.setLambdas(0.9999999999,9.999900827395436E-11,1.000000082740371E-15);
		hmm20K.baumWelch(trainTritext,trainTritextTaggedModel, ibm1.f_count, ibm1.fe_count);
		 
		hmm20K.multiplyOneMinusP0H();
		
		hmm20K.setLambdas(lambda1, lambda2, lambda3);
		
		HashMap<Pair,Double> hmm20k_t_table = new HashMap<Pair, Double>(hmm20K.t_table);
		
		
		//Augmented Model
		trainingSize = 1019987;//1020357; 
		String[][] train20KHKTritextTaggedModel = IBM_Model1.readTritext("data-jasneet/train.20K+1M.tag.clean.cn","data-jasneet/train.20K+1M.tag.clean.en","data-jasneet/train.20K.seg.cln.clean.wa",trainingSize);
		
		HashMap<Triple<String, String, Integer>, Double> s_hmmTaggedModelHK20K = runIBMHMMGAugmentedTaggedModel(ibm1TaggedModel, ibm1, train20KHKTritextTaggedModel);
		String[][] train20KHKTritext = IBM_Model1.readTritext("data-jasneet/train.20K+1M.word.clean.cn","data-jasneet/train.20K+1M.word.clean.en","data-jasneet/train.20K.seg.cln.clean.wa",trainingSize);

        runIBMHMMAugmentedModel(lambda1,lambda2, lambda3, ibm1, hmm20k_t_table, s_hmmTaggedModelHK20K,train20KHKTritext,train20KHKTritextTaggedModel); 
	
	}
	
	private static void testAndEvaluateAugmentedModel(int trainingSize,
			int testSize, String[][] testTritextTaggedModel,
			String[][] testTritext, 
			String[][] train20KHKTritextTaggedModel,
			String[][] train20KHKTritext, HMMWithAlignmentTypeNullCorrect hmm,
			String referenceFile, String referenceAlignmentAndType, String alignmentFile)
			throws IOException {
		
		IBM_Model1 ibm1 = new IBM_Model1();
		//Jan 18th/2017 2017/01/18
		
		ArrayList<String> hmmModelAlignment = hmm.findBestAlignmentsForAll_AER(testTritext, testTritextTaggedModel, testSize, alignmentFile, false);
		//reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		if (referenceFile != null){
			System.out.println("WA task");
			ArrayList<String> reference = ibm1.convertFileToArrayList(referenceFile);
			ibm1.gradeAlignmentWithType(testSize, testTritext, reference, hmmModelAlignment);
		}
		
		hmmModelAlignment = hmm.findBestAlignmentsForAll_AER(testTritext, testTritextTaggedModel, testSize, alignmentFile+".type", true);

		//ArrayList<String> reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.wa");
		if (referenceAlignmentAndType != null){
			System.out.println("WA+Type task");
			ArrayList<String> reference = ibm1.convertFileToArrayList(referenceAlignmentAndType);

			ibm1.gradeAlignmentWithTypeWAPlusTag(testSize, testTritext, reference, hmmModelAlignment);
			
		}
		
		

	//	hmm.findBestAlignmentsForAll_AER(train20KHKTritext, train20KHKTritextTaggedModel, trainingSize, "alignment_hmm_train_20k+2M_gen_WA_Aug5", false);

//		hmm.findBestAlignmentsForAll_AER(train20KHKTritext, train20KHKTritextTaggedModel, trainingSize, "alignment_generative.cn-en.clean.Jan28", false);
	}
	private static HMMWithAlignmentTypeNullCorrect runIBMHMMAugmentedModel(
			double lambda1,
			double lambda2,
			double lambda3,
			IBM_Model1 ibm1,
			HashMap<Pair, Double> hmm20K_t_table,
			HashMap<Triple<String, String, Integer>, Double> s_hmmTaggedModelHK20K,
			String[][] train20KHKTritext, String[][] train20khkTritextTaggedModel) {
		IBM_Model1 ibm120KHK = new IBM_Model1();
        
	    ibm120KHK.initializeCountsWithoutSets(train20KHKTritext);
	    //ibm120KHK.initializeCounts(train20KHKTritext,0);
		
	    System.out.println("Training an IBM1 model on the bitext of size " + train20KHKTritext.length + ".");
		//IBM Model1
		HashMap<Pair, Double> t_fe_20KHK = ibm1.EM_IBM1(ibm120KHK.f_count, ibm120KHK.e_count, ibm120KHK.fe_count, train20KHKTritext);
	        
        HMMWithAlignmentTypeNullCorrect hmm = new HMMWithAlignmentTypeNullCorrect(t_fe_20KHK, ibm120KHK.e_count, ibm1.total_f_e_h ,
        		ibm120KHK.fe_count, ibm1.fe_count, s_hmmTaggedModelHK20K);
        
        //This should not be IBM t_fe but hmm t_fe
        hmm.setTTable20K(hmm20K_t_table); 
        
        hmm.setLambdas(lambda1, lambda2, lambda3);
        System.out.println("Training an HMM model on the bitext of size " + train20KHKTritext.length + ".");
        hmm.baumWelch(train20KHKTritext,train20khkTritextTaggedModel, ibm120KHK.f_count, ibm120KHK.fe_count);
       	hmm.multiplyOneMinusP0H();
       	
		//System.out.println("p0H " + hmm.p0H + " nullEmissionProb" + hmm.nullEmissionProb + "lambda1 " + hmm.lambda1);
		return hmm;
	}
	private static HashMap<Triple<String, String, Integer>, Double> runIBMHMMGAugmentedTaggedModel(
			IBM_Model1 ibm1TaggedModel, IBM_Model1 ibm1,
			String[][] train20KHKTritextTaggedModel) {
		IBM_Model1 ibm120KHKTaggedModel = new IBM_Model1();
        
        //ibm120KHKTaggedModel.initializeCounts(train20KHKTritextTaggedModel, 0);
        ibm120KHKTaggedModel.initializeCountsWithoutSets(train20KHKTritextTaggedModel);
       
        System.out.println("Training an IBM1 model on the tagged bitext of size " + train20KHKTritextTaggedModel.length + ".");
        HashMap<Pair, Double> t_fe_20KHK_tagged = ibm1.EM_IBM1(ibm120KHKTaggedModel.f_count, ibm120KHKTaggedModel.e_count, ibm120KHKTaggedModel.fe_count, train20KHKTritextTaggedModel);
        
        //HMMWithAlignmentTypePOSTag hmmTaggedModelHK20K = new HMMWithAlignmentTypePOSTag(t_fe_20KHK_tagged, ibm120KHKTaggedModel.e_count, ibm120KHKTaggedModel.total_f_e_h,ibm120KHKTaggedModel.fe_count,);
        System.out.println("Training an HMM model on the tagged bitext of size " + train20KHKTritextTaggedModel.length + ".");
        HMMWithAlignmentTypePOSTag hmmTaggedModelHK20K = new HMMWithAlignmentTypePOSTag(t_fe_20KHK_tagged, ibm120KHKTaggedModel.e_count, ibm1TaggedModel.total_f_e_h, ibm120KHKTaggedModel.fe_count, ibm1TaggedModel.fe_count);
        hmmTaggedModelHK20K.baumWelch(train20KHKTritextTaggedModel, ibm120KHKTaggedModel.f_count, ibm120KHKTaggedModel.fe_count);
        hmmTaggedModelHK20K.multiplyOneMinusP0H();
		return hmmTaggedModelHK20K.s;
	}
	

	public static void FScoreExperiment1AlignmentTypeMain(int trainingSize, int testSize, String trainPrefix
			, String sourceLang, String targetLang, String testPrefix, String referenceFile, String output, 
			String goldTrainPrefix, String POSTaggedTrainPrefix, String POSTaggedTestPrefix, String referenceAlignmentAndType) throws IOException{
		IBM_Model1 ibm1 = new IBM_Model1();
        
        IBM_Model1 ibm1TaggedModel = new IBM_Model1();

		String[][] trainTritextTaggedModel = IBM_Model1.readTritext(POSTaggedTrainPrefix+"."+sourceLang, POSTaggedTrainPrefix+"."+targetLang, goldTrainPrefix, trainingSize);
		String[][] testBitextTaggedModel = IBM_Model1.readBitext(POSTaggedTestPrefix+"."+sourceLang, POSTaggedTestPrefix+"."+targetLang, testSize);
	
		ibm1TaggedModel.initializeCounts(trainTritextTaggedModel, 0);
		System.out.println("Training an IBM1 model on the tagged bitext of size " + trainTritextTaggedModel.length + ".");
        HashMap<Pair, Double> t_fe_tagged = ibm1TaggedModel.EM_IBM1(ibm1TaggedModel.f_count, ibm1TaggedModel.e_count, ibm1TaggedModel.fe_count, trainTritextTaggedModel);
  
        HMMWithAlignmentTypePOSTag hmmTaggedModel = new HMMWithAlignmentTypePOSTag(t_fe_tagged, ibm1TaggedModel.e_count, ibm1TaggedModel.total_f_e_h,ibm1TaggedModel.fe_count,ibm1TaggedModel.fe_count);
        System.out.println("Training an HMM model on the tagged bitext of size " + trainTritextTaggedModel.length + ".");
        hmmTaggedModel.baumWelch(trainTritextTaggedModel, ibm1TaggedModel.f_count, ibm1TaggedModel.fe_count);
        hmmTaggedModel.multiplyOneMinusP0H();
		
		String[][] trainTritext = IBM_Model1.readTritext(trainPrefix+"."+sourceLang, trainPrefix+"."+targetLang, goldTrainPrefix, trainingSize);
		String[][] testBitext = IBM_Model1.readBitext(testPrefix+"."+sourceLang, testPrefix+"."+targetLang, testSize);
				
		ibm1.initializeCounts(trainTritext, 0);
		System.out.println("Training an IBM1 model on the bitext of size " + trainTritext.length + ".");
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainTritext);
		HMMWithAlignmentTypeNullCorrect hmmL = new HMMWithAlignmentTypeNullCorrect(t_fe, ibm1.e_count, ibm1.total_f_e_h,ibm1.fe_count, ibm1.fe_count,hmmTaggedModel.s);
		hmmL.setLambdas(0.9999999999,9.999900827395436E-11,1.000000082740371E-15);
		System.out.println("Training an HMM model on the bitext of size " + trainTritext.length + ".");
		hmmL.baumWelch(trainTritext,trainTritextTaggedModel, ibm1.f_count, ibm1.fe_count);
		
 
		hmmL.multiplyOneMinusP0H();
		


		double lambda1 = 0.9999999999;
		double lambda2 = 9.999900827395436E-11;
		double lambda3 = 1.000000082740371E-15;
		
		hmmL.setLambdas(lambda1, lambda2, lambda3);
		
		ArrayList<String> hmmModelAlignment = hmmL.findBestAlignmentsForAll_AER(testBitext, testBitextTaggedModel, testSize, output, false);
		
		if (referenceFile != null){
			System.out.println("\nWA:");
			ArrayList<String> reference = ibm1.convertFileToArrayList(referenceFile);

			ibm1.gradeAlignmentWithType(testSize, testBitext, reference, hmmModelAlignment);
		}
		hmmModelAlignment = hmmL.findBestAlignmentsForAll_AER(testBitext, testBitextTaggedModel, testSize, output+".type", true);
		
		if (referenceAlignmentAndType != null){
			System.out.println("\nWA+Type:");
			ArrayList<String> reference = ibm1.convertFileToArrayList(referenceAlignmentAndType);
		
			ibm1.gradeAlignmentWithTypeWAPlusTag(testSize, testBitext, reference, hmmModelAlignment);
		}
						
	
	}
	
	public static void augmentedModelExperiment(int trainingSize, int testSize, String trainPrefix
			, String sourceLang, String targetLang, String testPrefix, String referenceFile, String alignmentFile, 
			String goldTrain, String POSTaggedTrainPrefix, String POSTaggedTestPrefix, String referenceAlignmentAndType
			, String augmentedTrainPrefix, String POSTaggedAugmentedTrainPrefix, int augmentedTrainSize) throws IOException, ClassNotFoundException{
/**
 * For computing BLEU Score testTritexts are not important. We can even use nulls for them		
 */
	
		String[][] trainTritextTaggedModel = IBM_Model1.readTritext(POSTaggedTrainPrefix+"."+sourceLang, POSTaggedTrainPrefix+"."+targetLang, goldTrain, trainingSize);
		String[][] testBitextTaggedModel = IBM_Model1.readBitext(POSTaggedTestPrefix+"."+sourceLang, POSTaggedTestPrefix+"."+targetLang, testSize);
		
		String[][] trainTritext = IBM_Model1.readTritext(trainPrefix+"."+sourceLang, trainPrefix+"."+targetLang, goldTrain, trainingSize);
		String[][] testBitext = IBM_Model1.readBitext(testPrefix+"."+sourceLang, testPrefix+"."+targetLang, testSize);
		

        double lambda1 = 0.9999999999;
        double lambda2 = 9.999900827395436E-11;
        double lambda3 = 1.000000082740371E-15;
        
        IBM_Model1 ibm1TaggedModel = new IBM_Model1();
        IBM_Model1 ibm1 = new IBM_Model1();
        ibm1TaggedModel.initializeCounts(trainTritextTaggedModel, 0);
        System.out.println("Training an IBM1 model on the tagged bitext of size " + trainTritextTaggedModel.length + ".");
        HashMap<Pair, Double> t_fe_tagged = ibm1TaggedModel.EM_IBM1(ibm1TaggedModel.f_count, ibm1TaggedModel.e_count, ibm1TaggedModel.fe_count, trainTritextTaggedModel);
  
        System.out.println("Training an HMM model on the tagged bitext of size " + trainTritextTaggedModel.length + ".");
        HMMWithAlignmentTypePOSTag hmmTaggedModel = new HMMWithAlignmentTypePOSTag(t_fe_tagged, ibm1TaggedModel.e_count, ibm1TaggedModel.total_f_e_h,ibm1TaggedModel.fe_count,ibm1TaggedModel.fe_count);
        hmmTaggedModel.baumWelch(trainTritextTaggedModel, ibm1TaggedModel.f_count, ibm1TaggedModel.fe_count);
        hmmTaggedModel.multiplyOneMinusP0H();
        
   
        ibm1.initializeCounts(trainTritext, 0);
                
        
        System.out.println("Training an IBM1 model on the bitext of size " + trainTritext.length + ".");
        HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainTritext);
        
        System.out.println("Training an HMM model on the bitext of size " + trainTritext.length + ".");
        HMMWithAlignmentTypeNullCorrect hmm20K = new HMMWithAlignmentTypeNullCorrect(t_fe, ibm1.e_count, ibm1.total_f_e_h,ibm1.fe_count, ibm1.fe_count,hmmTaggedModel.s);
        hmm20K.setLambdas(0.9999999999,9.999900827395436E-11,1.000000082740371E-15);
        hmm20K.baumWelch(trainTritext,trainTritextTaggedModel, ibm1.f_count, ibm1.fe_count);
         
        hmm20K.multiplyOneMinusP0H();
        
        
        hmm20K.setLambdas(lambda1, lambda2, lambda3);
        
                
        HashMap<Pair,Double> hmm20k_t_table = new HashMap<Pair, Double>(hmm20K.t_table);
                

        //Augmented Model
        String[][] train20KHKTritextTaggedModel = IBM_Model1.readTritext(POSTaggedAugmentedTrainPrefix+"."+sourceLang,POSTaggedAugmentedTrainPrefix+"."+targetLang,goldTrain,augmentedTrainSize);
        
        HashMap<Triple<String, String, Integer>, Double> s_hmmTaggedModelHK20K = runIBMHMMGAugmentedTaggedModel(ibm1TaggedModel, ibm1, train20KHKTritextTaggedModel);
       
        String[][] train20KHKTritext = IBM_Model1.readTritext(augmentedTrainPrefix+"."+sourceLang, augmentedTrainPrefix+"."+targetLang,goldTrain, augmentedTrainSize);

        HMMWithAlignmentTypeNullCorrect hmm = runIBMHMMAugmentedModel(lambda1,lambda2, lambda3, ibm1, hmm20k_t_table, s_hmmTaggedModelHK20K,train20KHKTritext,train20KHKTritextTaggedModel);
        testAndEvaluateAugmentedModel(trainingSize, testSize,testBitextTaggedModel, testBitext,train20KHKTritextTaggedModel, train20KHKTritext, hmm, referenceFile, referenceAlignmentAndType, alignmentFile);
		
	}
public static void augmentedReverseModelExperiment() throws IOException, ClassNotFoundException{
	int trainingSize = 19987;//20357;//18000;
        int testSize = 1956;//2357;

/**
 * For computing BLEU Score testTritexts are not important. We can even use nulls for them		
 */
		String[][] trainTritextTaggedModel = IBM_Model1.readTritext("data-jasneet/train.20K.tag.clean.en","data-jasneet/train.20K.tag.clean.cn","data-jasneet/train.20K.en-cn.seg.cln.clean.wa",trainingSize);
		String[][] testTritextTaggedModel = null;//IBM_Model1.readTritext("data-te/test.tags.cn","data-te/test.tags.en","data-te/test.seg.cln.wa",testSize);
		
		String[][] trainTritext = IBM_Model1.readTritext("data-jasneet/train.20K.clean.en","data-jasneet/train.20K.clean.cn","data-jasneet/train.20K.en-cn.seg.cln.clean.wa",trainingSize);
		String[][] testTritext = null;//IBM_Model1.readTritext("data-te/test.seg.cln.cn","data-te/test.seg.cln.en","data-te/test.seg.cln.wa",testSize);


        double lambda1 = 0.9999999999;
        double lambda2 = 9.999900827395436E-11;
        double lambda3 = 1.000000082740371E-15;
        
        
        IBM_Model1 ibm1TaggedModel = new IBM_Model1();
        IBM_Model1 ibm1 = new IBM_Model1();
        ibm1TaggedModel.initializeCounts(trainTritextTaggedModel, 0);
        HashMap<Pair, Double> t_fe_tagged = ibm1TaggedModel.EM_IBM1(ibm1TaggedModel.f_count, ibm1TaggedModel.e_count, ibm1TaggedModel.fe_count, trainTritextTaggedModel);
  
        HMMWithAlignmentTypePOSTag hmmTaggedModel = new HMMWithAlignmentTypePOSTag(t_fe_tagged, ibm1TaggedModel.e_count, ibm1TaggedModel.total_f_e_h,ibm1TaggedModel.fe_count,ibm1TaggedModel.fe_count);
        hmmTaggedModel.baumWelch(trainTritextTaggedModel, ibm1TaggedModel.f_count, ibm1TaggedModel.fe_count);
        hmmTaggedModel.multiplyOneMinusP0H();
        
        
        ibm1.initializeCounts(trainTritext, 0);
                
        System.out.println("length " + trainTritext.length);
        
        HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainTritext);
        
        HMMWithAlignmentTypeNullCorrect hmm20K = new HMMWithAlignmentTypeNullCorrect(t_fe, ibm1.e_count, ibm1.total_f_e_h,ibm1.fe_count, ibm1.fe_count,hmmTaggedModel.s);
        hmm20K.setLambdas(0.9999999999,9.999900827395436E-11,1.000000082740371E-15);
        hmm20K.baumWelch(trainTritext,trainTritextTaggedModel, ibm1.f_count, ibm1.fe_count);
         
        hmm20K.multiplyOneMinusP0H();
        
        System.out.println("p0H " + hmm20K.p0H + " nullEmissionProb" + hmm20K.nullEmissionProb);
        
        hmm20K.setLambdas(lambda1, lambda2, lambda3);
            
                
        HashMap<Pair,Double> hmm20k_t_table = new HashMap<Pair, Double>(hmm20K.t_table);
                
        //Augmented Model
        trainingSize = 1019987;

        String[][] train20KHKTritextTaggedModel = IBM_Model1.readTritext("data-jasneet/train.20K+1M.tag.clean.en","data-jasneet/train.20K+1M.tag.clean.cn","data-jasneet/train.20K.en-cn.seg.cln.clean.wa",trainingSize);
        
        HashMap<Triple<String, String, Integer>, Double> s_hmmTaggedModelHK20K = runIBMHMMGAugmentedTaggedModel(ibm1TaggedModel, ibm1, train20KHKTritextTaggedModel);
        String[][] train20KHKTritext = IBM_Model1.readTritext("data-jasneet/train.20K+1M.word.clean.en","data-jasneet/train.20K+1M.word.clean.cn","data-jasneet/train.20K.en-cn.seg.cln.clean.wa",trainingSize);

        HMMWithAlignmentTypeNullCorrect hmm = runIBMHMMAugmentedModel(lambda1,lambda2, lambda3, ibm1, hmm20k_t_table, s_hmmTaggedModelHK20K,train20KHKTritext,train20KHKTritextTaggedModel);

        testAndEvaluateAugmentedModel(trainingSize, testSize,testTritextTaggedModel, testTritext,train20KHKTritextTaggedModel, train20KHKTritext, hmm, "", "", "alignment");
	}
	
}
