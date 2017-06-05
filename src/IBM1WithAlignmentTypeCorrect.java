import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;


public class IBM1WithAlignmentTypeCorrect extends IBM1WithAlignmentTypePOSTag implements Serializable{
	
	private static final long serialVersionUID = -2641697232475893704L;
	HashMap<Triple<String, String, Integer>, Double> tagModel_s;
	double lambda1;
	double lambda2;
	double lambda3;
	String[][] bitext_tag_fe;
	
	HashMap<Pair, Double> t_table_20K;
	
	public IBM1WithAlignmentTypeCorrect(HashMap<Triple<String, String, Integer>, Double> tagModel_s, String[][] bitext_tag_fe){
		this.tagModel_s = tagModel_s;
		this.bitext_tag_fe = bitext_tag_fe;
	}
	public void setLambdas(double lambda1,double lambda2,double lambda3){
		this.lambda1 = lambda1;
		this.lambda2 = lambda2;
		this.lambda3 = lambda3;
	}
	public void setTTable20K(HashMap<Pair, Double> t_table_20K){
		this.t_table_20K = t_table_20K;
	}
	
	public HashMap<Pair, Double> EM_IBM1(HashMap<String, Integer> sourceCount, HashMap<String, Integer> targetCount, HashMap<Pair, Integer> st_count, String[][] bitext){
		HashMap<Pair, Double> t = new HashMap<Pair, Double>();
		
		double initialValue = new Double(1.0/sourceCount.size());
		//initialize t probabilities
		for (Pair p: st_count.keySet()){
			t.put(p, initialValue);
		}
		
		HashMap<Pair, Double> c = new HashMap<Pair, Double>();
		
		HashMap<String, Double> total = new HashMap<String, Double>();
		HashMap<Triple<String, String, Integer>, Double> c_feh = new HashMap<Triple<String,String,Integer>, Double>();


		Pair sdPair = new Pair("","");
		Pair newPair1 = new Pair("","");
		
		for (int iteration = 0; iteration < 5 ; iteration++){
			System.out.println("IBM iteration " + iteration);
			
			c.clear();
			total.clear();
			c_feh.clear();
			
			int sentenceNumber = 0;
			for (String[] pair : bitext){
				
				String[] S = pair[0].trim().split(" ");
				String[] D = pair[1].trim().split(" ");
				
				String[] tagS = bitext_tag_fe[sentenceNumber][0].trim().split(" ");
				String[] tagD = bitext_tag_fe[sentenceNumber][1].trim().split(" ");
				
				int i = 0;
				for (String s_i : S){
					double Z = 0;
					for (String d_j : D){
						
						sdPair.setPair(s_i, d_j);
						Z += t.get(sdPair).doubleValue();
					}
					int j = 0;
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
								c_feh.put(f_e_h, c_feh.get(f_e_h) + (t.get(newPair1) * sProbability(s_i, d_j, h, tagS[i], tagD[j]))/Z);
							else
								c_feh.put(f_e_h, (t.get(newPair1) * sProbability(s_i, d_j, h, tagS[i], tagD[j]))/Z);
						}
						j++;
					}
					i++;
				}	
				sentenceNumber++;
			}
			for (Pair sd : st_count.keySet()){
				if (!total.containsKey(sd.right))
					System.out.println(sd + " and " + sd.right);
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
		
		if (t_table_20K != null && t_table_20K.containsKey(fe))
			return t_table_20K.get(fe);
 		else if (t_table.containsKey(fe)){
			return t_table.get(fe);
		}
		else{
				return 1.0/V;
		} 
	}
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
	//returns an ArrayList of alignments results of aligning through IBM1 source to target
		public ArrayList<String> print_alignment_SD_ibm1(String[][] bitext, String[][] bitextTaggedModel, HashMap<Pair, Double> t, String alignmentFile, boolean useAlignmentType) throws IOException{
			String[] linkMap = {"SEM", "FUN", "PDE", "CDE", "MDE", "GIS", "GIF", "COI", "TIN", "NTR", "MTA"};
			ArrayList<String> alignmentList = new ArrayList<String>();
			BufferedWriter alignment = new BufferedWriter(new FileWriter(alignmentFile));
			int sentenceNumber = 0;
			for (String[] pair : bitext){
				
				String[] S = pair[0].trim().split(" ");
				String[] D = pair[1].trim().split(" ");
				
				String[] tagS = bitextTaggedModel[sentenceNumber][0].trim().split(" ");
				String[] tagD = bitextTaggedModel[sentenceNumber][1].trim().split(" ");
				
				HashSet<IntPair> set = new HashSet<IntPair>();
				String lineWATag = "";
				for (int i = 0; i < S.length; i++){
					double max_ts = 0;
					int argmax = -1;
					int bestH = 0;
					for (int j = 0; j < D.length; j++){
						
						double t_SiDj = tProbability(S[i], D[j], t);
						for (int h = 1; h < H+1; h++){
							double s_SiDj = sProbability(S[i], D[j], h, tagS[i], tagD[j]);
							if (t_SiDj*s_SiDj > max_ts){
								max_ts = t_SiDj*s_SiDj;
								argmax = j;
								bestH = h;
							}
						}
					}
					
					set.add(new IntPair(i,argmax));
					
					lineWATag += (i+1) + "-" + (argmax+1) + "(" + linkMap[bestH-1] + ") ";	
				}
				String line = "";
				if (!useAlignmentType)
					for (IntPair ip : set){
						line += (ip.left+1) + "-" + (ip.right+1) + " ";
					}
				if (!useAlignmentType){
					alignmentList.add(line);
					alignment.write(line+"\n");
				}
				else{
					alignmentList.add(lineWATag);
					alignment.write(lineWATag+"\n");
				}
				sentenceNumber++;
			}
			alignment.close();
			return alignmentList;
		}
		public ArrayList<String> convertFileToArrayList(String fileName){
			ArrayList<String> list = new ArrayList<String>();
			BufferedReader reader = null;
			try{
				 reader = new BufferedReader(new FileReader(fileName));
				String line = null;
				while ((line = reader.readLine()) != null){
					list.add(line);
				}
				reader.close();
			}catch(IOException ioe){
				ioe.printStackTrace();
			}finally{
				try{
					reader.close();
				}catch(IOException ioe){
					ioe.printStackTrace();
				}
				
			}
			return list;
		}
		public double gradeAlignmentWithTypeWAPlusTag(int alignementTestSize, String[][] bitext, ArrayList<String> reference, ArrayList<String> systemAlignment){
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
				
				
				HashSet<Triple<Integer, Integer, String>> alignment = new HashSet<Triple<Integer, Integer, String>>();
				String[] pairsWithDash = a.trim().split(" ");
				
				for (String pwd : pairsWithDash){
					int index = pwd.indexOf('-');
					
					if (index != -1){
						String right = pwd.substring(index+1);
						int rightLength = right.length();
						int engPos = Integer.parseInt(right.substring(0,rightLength-5));
		        		String linkLabel = right.substring(rightLength-4, rightLength-1); 
						alignment.add(new Triple<Integer, Integer, String>(Integer.parseInt(pwd.substring(0,index)), engPos, linkLabel ) );
					}
					
				}
				for (Triple<Integer, Integer, String> ip : alignment){
					if (ip.left > size_f || ip.middle > size_e){
						System.err.print("NOT A VALID LINK");
						System.out.println(i + " " + ip.left +" " + size_f + " " + ip.right + " " + size_e);
					}
				}
				//grade
				HashSet<Triple<Integer, Integer, String>> sure = new HashSet<Triple<Integer, Integer, String>>();
				//HashSet<Triple<Integer, Integer, String>> possible = new HashSet<Triple<Integer, Integer, String>>();
				String[] surePairsWithDash = g.trim().split(" ");
				for (String spwd : surePairsWithDash){

					int index = spwd.indexOf('-');
					
					if (index != -1){
						String left = spwd.substring(0,index);
						String[] leftPositions = strip(left);
		        		if (leftPositions.length == 1 && !leftPositions[0].equals("")){
							String right = spwd.substring(index+1);
			        		int rightLength = right.length();
			        		String linkLabel = right.substring(rightLength-4, rightLength-1); 
			        		String[] engPositions = strip(right);
			        		
							
							for (String engPos : engPositions){
								if (!engPos.equals(""))
									sure.add(new Triple<Integer, Integer, String>(Integer.parseInt(leftPositions[0]), Integer.parseInt(engPos), linkLabel ) );
							}
						}
					}
					
				}
				
				size_a += alignment.size();
				size_s += sure.size();
				HashSet<Triple<Integer, Integer, String>> aAnds = new HashSet<Triple<Integer, Integer, String>>(sure);
				aAnds.retainAll(alignment);
				size_a_and_s += aAnds.size();
				
				
				size_a_and_p += (aAnds.size());
			
			}
			double precision = (double)size_a_and_p / size_a;
			double recall = (double)size_a_and_s / size_s;
			double aer = 1 - ((double)(size_a_and_s + size_a_and_p) / (size_a + size_s));
			double fScore = 2*precision*recall/(precision+recall);
			System.out.println("Precision = " + precision + "\n Recall = " + recall + "\n AER = " + aer + "\n F-score = " + fScore);
			return fScore;
		}
		
		
		public static void augmentedModelExperimentMain(int trainingSize, int testSize, String trainPrefix
		, String sourceLang, String targetLang, String testPrefix, String referenceFile, String output, 
		String goldTrain, String POSTaggedTrainPrefix, String POSTaggedTestPrefix, String referenceAlignmentAndType
		, String augmentedTrainPrefix, String POSTaggedAugmentedTrainPrefix, int augmentedTrainSize) throws IOException, ClassNotFoundException{

			IBM1WithAlignmentTypePOSTag ibm1TaggedModel = new IBM1WithAlignmentTypePOSTag();
			
			String[][] trainTritextTaggedModel = IBM_Model1.readTritext(POSTaggedTrainPrefix+"."+sourceLang, POSTaggedTrainPrefix+"."+targetLang, goldTrain, trainingSize);
			String[][] testBitextTaggedModel = IBM_Model1.readBitext(POSTaggedTestPrefix+"."+sourceLang, POSTaggedTestPrefix+"."+targetLang, testSize);
			
			String[][] trainTritext = IBM_Model1.readTritext(trainPrefix+"."+sourceLang, trainPrefix+"."+targetLang, goldTrain, trainingSize);
			String[][] testBitext = IBM_Model1.readBitext(testPrefix+"."+sourceLang, testPrefix+"."+targetLang, testSize);
			
			
			ibm1TaggedModel.initializeCounts(trainTritextTaggedModel);
			ibm1TaggedModel.setSProbabilities(ibm1TaggedModel.fe_count, ibm1TaggedModel.total_f_e_h);
			
			System.out.println("Training an IBM1 model on the tagged bitext of size " + trainTritextTaggedModel.length + ".");
			ibm1TaggedModel.EM_IBM1(ibm1TaggedModel.f_count, ibm1TaggedModel.e_count, ibm1TaggedModel.fe_count, trainTritextTaggedModel);
			
			
			IBM1WithAlignmentTypeCorrect ibm1 = new IBM1WithAlignmentTypeCorrect(ibm1TaggedModel.s, trainTritextTaggedModel);
			
			ibm1.initializeCounts(trainTritext);
			
			//set lambdas
			ibm1.setLambdas(0.9999999999,9.999900827395436E-11,1.000000082740371E-15);
			ibm1.setSProbabilities(ibm1.fe_count, ibm1.total_f_e_h);
			System.out.println("Training an IBM1 model on the bitext of size " + trainTritextTaggedModel.length + ".");
			HashMap<Pair, Double> t_table_20K = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainTritext);
			
			//Augmented Model
			
			String[][] train20KHKTritextTaggedModel = IBM_Model1.readTritext(POSTaggedAugmentedTrainPrefix+"."+sourceLang,POSTaggedAugmentedTrainPrefix+"."+targetLang,goldTrain,augmentedTrainSize);
			
			IBM1WithAlignmentTypePOSTag ibm120KHKTagged = new IBM1WithAlignmentTypePOSTag();
			ibm120KHKTagged.initializeCountsOfAugmentedModel(train20KHKTritextTaggedModel);
			ibm120KHKTagged.setSProbabilities(ibm1TaggedModel.fe_count, ibm1TaggedModel.total_f_e_h);
			System.out.println("Training an IBM1 model on the tagged bitext of size " + train20KHKTritextTaggedModel.length + ".");
			ibm120KHKTagged.EM_IBM1(ibm120KHKTagged.f_count, ibm120KHKTagged.e_count, ibm120KHKTagged.fe_count, train20KHKTritextTaggedModel);
			
			
			String[][] train20KHKTritext = IBM_Model1.readTritext(augmentedTrainPrefix+"."+sourceLang, augmentedTrainPrefix+"."+targetLang,goldTrain, augmentedTrainSize);
			
			IBM1WithAlignmentTypeCorrect ibm120KHK = new IBM1WithAlignmentTypeCorrect(ibm120KHKTagged.s, train20KHKTritextTaggedModel);
			ibm120KHK.initializeCountsOfAugmentedModel(train20KHKTritext);
			ibm120KHK.setTTable20K(t_table_20K);
			ibm120KHK.setLambdas(0.9999999999,9.999900827395436E-11,1.000000082740371E-15);
			ibm120KHK.setSProbabilities(ibm1.fe_count, ibm1.total_f_e_h);
			System.out.println("Training an IBM1 model on the bitext of size " + train20KHKTritext.length + ".");
			HashMap<Pair, Double> t_fe = ibm120KHK.EM_IBM1(ibm120KHK.f_count, ibm120KHK.e_count, ibm120KHK.fe_count, train20KHKTritext);
			
			System.out.println("\nWA:");			
			ArrayList<String> ibmModelAlignment = ibm120KHK.print_alignment_SD_ibm1(testBitext, testBitextTaggedModel, t_fe, output, false);
			
			if (referenceFile != null){
				ArrayList<String> reference = ibm120KHK.convertFileToArrayList(referenceFile);
				ibm120KHK.gradeAlignmentWithType(testSize, testBitext, reference, ibmModelAlignment);
			}
			System.out.println("\nWA+Type:");
			ibmModelAlignment = ibm120KHK.print_alignment_SD_ibm1(testBitext, testBitextTaggedModel, t_fe, output+".type", true);
			
			if (referenceAlignmentAndType != null){
				ArrayList<String> reference = ibm120KHK.convertFileToArrayList(referenceAlignmentAndType);
				ibm120KHK.gradeAlignmentWithTypeWAPlusTag(testSize, testBitext, reference, ibmModelAlignment);
			}
			
		}
		
		public static void fScoreExperiment(int trainingSize, int testSize, String trainPrefix
				, String sourceLang, String targetLang, String testPrefix, String referenceFile, String output, 
				String goldTrainPrefix, String POSTaggedTrainPrefix, String POSTaggedTestPrefix, String referenceAlignmentAndType) throws IOException{
			IBM1WithAlignmentTypePOSTag ibm1TaggedModel = new IBM1WithAlignmentTypePOSTag();
			
			String[][] trainTritextTaggedModel = IBM_Model1.readTritext(POSTaggedTrainPrefix+"."+sourceLang, POSTaggedTrainPrefix+"."+targetLang, goldTrainPrefix, trainingSize);
			String[][] testBitextTaggedModel = IBM_Model1.readBitext(POSTaggedTestPrefix+"."+sourceLang, POSTaggedTestPrefix+"."+targetLang, testSize);
			
			ibm1TaggedModel.initializeCounts(trainTritextTaggedModel);
			ibm1TaggedModel.setSProbabilities(ibm1TaggedModel.fe_count, ibm1TaggedModel.total_f_e_h);
			System.out.println("Training an IBM1 model on the tagged bitext of size " + trainTritextTaggedModel.length + ".");
			ibm1TaggedModel.EM_IBM1(ibm1TaggedModel.f_count, ibm1TaggedModel.e_count, ibm1TaggedModel.fe_count, trainTritextTaggedModel);
			
			String[][] trainTritext = IBM_Model1.readTritext(trainPrefix+"."+sourceLang, trainPrefix+"."+targetLang, goldTrainPrefix, trainingSize);
			String[][] testBitext = IBM_Model1.readBitext(testPrefix+"."+sourceLang, testPrefix+"."+targetLang, testSize);
			
			IBM1WithAlignmentTypeCorrect ibm1 = new IBM1WithAlignmentTypeCorrect(ibm1TaggedModel.s, trainTritextTaggedModel);
			
			ibm1.initializeCounts(trainTritext);
			
			//set lambdas
			ibm1.setLambdas(0.9999999999,9.999900827395436E-11,1.000000082740371E-15);
			
			ibm1.setSProbabilities(ibm1.fe_count, ibm1.total_f_e_h);
			System.out.println("Training an IBM1 model on the bitext of size " + trainTritext.length + ".");
			HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainTritext);
			
			
			System.out.println("\nWA:");			
			ArrayList<String> ibmModelAlignment = ibm1.print_alignment_SD_ibm1(testBitext, testBitextTaggedModel, t_fe, output , false);
			
			if (referenceFile != null){
				ArrayList<String> reference = ibm1.convertFileToArrayList(referenceFile);
				ibm1.gradeAlignmentWithType(testSize, testBitext, reference, ibmModelAlignment);
			}
			
			System.out.println("\nWA+Type:");
			ibmModelAlignment = ibm1.print_alignment_SD_ibm1(testBitext, testBitextTaggedModel, t_fe, output+".type", true);
			
			if (referenceAlignmentAndType != null){
				ArrayList<String> reference = ibm1.convertFileToArrayList(referenceAlignmentAndType);
				ibm1.gradeAlignmentWithTypeWAPlusTag(testSize, testBitext, reference, ibmModelAlignment);
			}

		}
		
}
