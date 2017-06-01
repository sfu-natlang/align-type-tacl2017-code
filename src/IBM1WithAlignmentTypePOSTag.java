import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;


public class IBM1WithAlignmentTypePOSTag implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -2641697232475893704L;
	HashMap<String, Integer> f_count = new HashMap<String, Integer>();
	HashMap<String, Integer> e_count = new HashMap<String, Integer>();
	HashMap<Pair, Integer> fe_count = new HashMap<Pair, Integer>();
	HashMap<String,Integer> tagMap = new HashMap<String, Integer>();
	HashMap<Triple<String, String, Integer>, Double> total_f_e_h = new HashMap<Triple<String,String,Integer>, Double>();
	HashMap<Triple<String, String, Integer>, Double> s = new HashMap<Triple<String,String,Integer>, Double>();
	int H = 11;
	double lambda = 1-1e-20;
	
	public void setSProbabilities(HashMap<Pair, Integer> fe_count, HashMap<Triple<String, String, Integer>, Double> total_f_e_h) {
		for (Entry<Triple<String, String, Integer>, Double> entry : total_f_e_h.entrySet()){
			Triple<String , String, Integer> f_e_h = entry.getKey();
			String f = f_e_h.left;
			String e = f_e_h.middle;
			double value = entry.getValue();
			s.put(f_e_h, value/fe_count.get(new Pair(f,e)));
		}
	}
	
	public void initializeTagMap(){
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
	 * This method computes source and target counts as well as (source, target, alignment type) counts
	 * (f,e,h) counts are stored in total_f_e_h
	 * HMMWithAlignmentType initializes its s parameter from total_f_e_h 
	 * @param tritext
	 */
	public void initializeCounts(String[][] tritext){
		
		initializeTagMap();
		
		for (String[] fewa : tritext){
			String[] F = fewa[0].trim().split(" ");
			String[] E = fewa[1].trim().split(" ");
			String[] wa = fewa[2].trim().split(" ");
			
		
			//Setting f_count
			for (String f : F){ //I replaced F with setF ***********
				if (f_count.containsKey(f))
					f_count.put(f, new Integer(f_count.get(f).intValue()+1));
				else
					f_count.put(f, new Integer(1));
				
				//Setting fe_count
				for (String e : E){ 
					Pair fePair = new Pair(f,e);

					if (fe_count.containsKey(fePair))
						fe_count.put(fePair, new Integer(fe_count.get(fePair).intValue()+1));
					else
						fe_count.put(fePair, new Integer(1));
				}
			}
			//setting e_count
			for (String e : E)
				if (e_count.containsKey(e))
					e_count.put(e, new Integer(e_count.get(e).intValue()+1));
				else
					e_count.put(e, new Integer(1));
			
			//setting total_f_e_h count
			
				for (String alm: wa){
					String[] splittedAlignment = alm.split("-");
	        		String left = splittedAlignment[0];
	        		String[] leftPositions = strip(left);
	        		if (leftPositions.length == 1 && !leftPositions[0].equals("")){
		        		int chineseWordPos = Integer.parseInt(leftPositions[0]);
		        		String chineseWord = F[chineseWordPos-1];
		        		String right = splittedAlignment[1];
		        		int rightLength = right.length();
		        		
		        		String linkLabel = right.substring(rightLength-4, rightLength-1); 
		        		String[] engIndices = strip(right.substring(0, rightLength-5));
		        		if (!engIndices[0].equals("")){
			        		for (String wordIndex : engIndices){
			        			int engWordPos = Integer.parseInt(wordIndex);
			        			String engWord = E[engWordPos-1];
			        			int tagId = tagMap.get(linkLabel);
			        			Triple<String, String, Integer> f_e_h = new Triple<String, String, Integer>(chineseWord, engWord, tagId);
			        			
			        			if (total_f_e_h.containsKey(f_e_h)){
			        				Double value = total_f_e_h.get(f_e_h);
			        				total_f_e_h.put(f_e_h, value+1);
			        			}else{
			        				total_f_e_h.put(f_e_h, 1.0);
			        			}
			        		}
		        		}
		        		
		        	}
				}
			
			
			
		}
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
	
	public void initializeCountsOfAugmentedModel(String[][] bitext){
			
			initializeTagMap();
			
			for (String[] fe : bitext){
				String[] F = fe[0].trim().split(" ");
				String[] E = fe[1].trim().split(" ");
				
				
				//Setting f_count
				for (String f : F){ 
					if (f_count.containsKey(f))
						f_count.put(f, new Integer(f_count.get(f).intValue()+1));
					else
						f_count.put(f, new Integer(1));
					
					//Setting fe_count
					for (String e : E){ //I replaced E with setE ***********
						Pair fePair = new Pair(f,e);
						//System.out.println("pair "+f + " "  + e);
						if (fe_count.containsKey(fePair))
							fe_count.put(fePair, new Integer(fe_count.get(fePair).intValue()+1));
						else
							fe_count.put(fePair, new Integer(1));
					}
				}
				//setting e_count
				for (String e : E)
					if (e_count.containsKey(e))
						e_count.put(e, new Integer(e_count.get(e).intValue()+1));
					else
						e_count.put(e, new Integer(1));
				
				
			}
		}
	
	public HashMap<Pair, Double> EM_IBM1(HashMap<String, Integer> sourceCount, HashMap<String, Integer> targetCount, HashMap<Pair, Integer> st_count, String[][] bitext){
		HashMap<Pair, Double> t = new HashMap<Pair, Double>();
		
		double initialValue = new Double(1.0/sourceCount.size());
		//initialize t probabilities
		for (Pair p: st_count.keySet()){
			t.put(p, initialValue);
		}
		//initialize s probabilities
		//setSProbabilities();
		HashMap<Pair, Double> c = new HashMap<Pair, Double>();
		HashMap<String, Double> total = new HashMap<String, Double>();
		HashMap<Triple<String, String, Integer>, Double> c_feh = new HashMap<Triple<String,String,Integer>, Double>();
		
		Pair sdPair = new Pair("","");
		Pair newPair1 = new Pair("","");
		
		for (int i = 0; i < 5 ; i++){
			System.out.println("IBM iteration " + i);
			
			c.clear();
			total.clear();
			c_feh.clear();
			
			for (String[] pair : bitext){
				String[] S = pair[0].split(" ");
				String[] D = pair[1].split(" ");
				
				for (String s_i : S){
					double Z = 0;
					for (String d_j : D){
						
						sdPair.setPair(s_i, d_j);
						Z += t.get(sdPair).doubleValue();
					}
					for (String d_j : D){
						
						newPair1.setPair(s_i,d_j);
						
							
						if (c.containsKey(newPair1))
							c.put(newPair1, new Double(c.get(newPair1).doubleValue() + t.get(newPair1).doubleValue()/Z));
						else{
							Pair newPair = new Pair(s_i,d_j);
							c.put(newPair, new Double(t.get(newPair).doubleValue()/Z));
						}
						if (total.containsKey(d_j))
							total.put(d_j, new Double(total.get(d_j).doubleValue() + t.get(newPair1).doubleValue()/Z));
						else
							total.put(d_j, new Double( t.get(newPair1).doubleValue()/Z));
						for (int h = 1; h < H+1; h++){
							Triple<String, String, Integer> f_e_h = new Triple<String, String, Integer>(s_i, d_j, h);
							if (c_feh.containsKey(f_e_h))
								c_feh.put(f_e_h, c_feh.get(f_e_h) + (t.get(newPair1) * sProbability(s_i, d_j, h))/Z);
							else
								c_feh.put(f_e_h, (t.get(newPair1) * sProbability(s_i, d_j, h))/Z);
						}
					}
				}	
			}
			for (Pair sd : st_count.keySet()){
				t.put(sd, c.get(sd).doubleValue()/total.get(sd.right).doubleValue());
			}
			for (Entry<Triple<String , String, Integer>, Double> entry : c_feh.entrySet()){
				Triple<String , String, Integer> f_e_h = entry.getKey();
				String f = f_e_h.left;
				String e = f_e_h.middle;
				double c_fehValue = entry.getValue();
				s.put(f_e_h, c_fehValue/c.get(new Pair(f,e)) );
			}
		}
		
		return t;
	}
	
	public double tProbability(String f, String e, HashMap<Pair, Double> t_table){
		Pair fe = new Pair(f,e);
		
		int V = 163303;
		if (t_table.containsKey(fe))
			return t_table.get(fe);
		else
			return 1.0/V; 
	}
	public double sProbability(String t_f, String t_e, Integer h){
		Triple<String, String, Integer> f_e_h = new Triple<String, String, Integer>(t_f,t_e,h);
		double first = (s.containsKey(f_e_h) ? s.get(f_e_h) : 0 );
		double[] tagDist = {0, 0.401, 0.264, 0.004, 0.004, 0.012, 0.205, 0.031, 0.008, 0.003, 0.086, 0.002};
		
		return lambda*first + (1-lambda)*tagDist[h];
			
	}
	//returns an ArrayList of alignments results of aligning through IBM1 source to target
	public ArrayList<String> print_alignment_SD_ibm1(String[][] bitext, HashMap<Pair, Double> t){
		ArrayList<String> alignmentList = new ArrayList<String>();
		for (String[] pair : bitext){
			String[] S = pair[0].trim().split(" ");
			String[] D = pair[1].trim().split(" ");
			
			HashSet<IntPair> set = new HashSet<IntPair>();
			for (int i = 0; i < S.length; i++){
				double max_ts = 0;
				int argmax = -1;
				for (int j = 0; j < D.length; j++){
					//Pair p = new Pair(S[i],D[j]);
					double t_SiDj = tProbability(S[i], D[j], t);
					for (int h = 1; h < H+1; h++){
						double s_SiDj = sProbability(S[i], D[j], h);
						if (t_SiDj*s_SiDj > max_ts){
							max_ts = t_SiDj*s_SiDj;
							argmax = j;
						}
					}
				}
				
				set.add(new IntPair(i,argmax));
			}
			String line = "";
			for (IntPair ip : set){
				System.out.print((ip.left+1) + "-" + (ip.right+1) + " ");
				line += (ip.left+1) + "-" + (ip.right+1) + " ";
			}
			alignmentList.add(line);
			System.out.println();
		}
		return alignmentList;
	}
	public double gradeAlignmentWithType(int alignementTestSize, String[][] bitext, ArrayList<String> reference, ArrayList<String> systemAlignment){
		int size_a = 0;
		int size_s = 0;
		int size_a_and_s = 0;
		int size_a_and_p = 0;
		
		int min = (alignementTestSize < reference.size() ? alignementTestSize : reference.size()); 
		for (int i = 0; i < min; i++){
			//alignment of sentence i
			String a = systemAlignment.get(i);
			String g = reference.get(i);
			
			int size_f = bitext[i][0].trim().split(" ").length;
			int size_e = bitext[i][1].trim().split(" ").length;
			
			
			HashSet<IntPair> alignment = new HashSet<IntPair>();
			String[] pairsWithDash = a.trim().split(" ");
			for (String pwd : pairsWithDash){
				int index = pwd.indexOf('-');
				if (index != -1){					
						alignment.add(new IntPair(Integer.parseInt(pwd.substring(0,index)), Integer.parseInt(pwd.substring(index+1)) ));
				}
			}
			
			for (IntPair ip : alignment){
				if (ip.left > size_f || ip.right > size_e){
					System.err.print("NOT A VALID LINK");
					System.out.println(i + " " + ip.left +" " + size_f + " " + ip.right + " " + size_e);
				}
			}
			//grade
			HashSet<IntPair> sure = new HashSet<IntPair>();
			HashSet<IntPair> possible = new HashSet<IntPair>();
			String[] surePairsWithDash = g.trim().split(" ");
			for (String spwd : surePairsWithDash){
				int index = spwd.indexOf('-');
				if (index != -1){
					String[] engPositions = spwd.substring(index+1).split(",");
					for (String engPos : engPositions)
						sure.add(new IntPair(Integer.parseInt(spwd.substring(0,index)), Integer.parseInt(engPos) ));
				}
				index = spwd.indexOf('?');
				if (index != -1)
					possible.add(new IntPair(Integer.parseInt(spwd.substring(0,index)), Integer.parseInt(spwd.substring(index+1)) ));
			}
			
			size_a += alignment.size();
			size_s += sure.size();
			HashSet<IntPair> aAnds = new HashSet<IntPair>(sure);
			aAnds.retainAll(alignment);
			size_a_and_s += aAnds.size();
			
			HashSet<IntPair> aAndp = new HashSet<IntPair>(possible);
			aAndp.retainAll(alignment);
			size_a_and_p += (aAndp.size() + aAnds.size());
			
			
			
		}
		double precision = (double)size_a_and_p / size_a;
		double recall = (double)size_a_and_s / size_s;
		double aer = 1 - ((double)(size_a_and_s + size_a_and_p) / (size_a + size_s));
		double fScore = 2*precision*recall/(precision+recall);
		System.out.println("Precision = " + precision + "\n Recall = " + recall + "\n AER = " + aer + "\n F-score = " + fScore);
		return fScore;
	}
	
}
