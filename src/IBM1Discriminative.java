import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;


public class IBM1Discriminative extends IBM_Model1{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	HashMap<Triple<Integer, Integer, Integer>, String> alignmentTypeMap = new HashMap<Triple<Integer, Integer, Integer>, String>();
	HashMap<Triple<Integer, Integer, Integer>, Double> probabilityMap = new HashMap<Triple<Integer, Integer, Integer>, Double>();
	
	double lambda = 1e-20;
	
	HashMap<Pair, Double> t_ibm1_20k;
	
	public void setTTable20k(HashMap<Pair, Double> t_ibm1_20k){
		this.t_ibm1_20k = t_ibm1_20k;
	}
	
	public void readLogisticRegFile(String fileName) throws IOException{
		BufferedReader br = null;
		try {
	    	File fileDir = new File(fileName);	    	 
			br = new BufferedReader(new FileReader(fileDir));
			String line = br.readLine();
			
			while (line != null){
				String[] info = line.split(" ");
				
				Triple<Integer, Integer, Integer> lineEnChPair = new Triple<Integer, Integer, Integer>(Integer.parseInt(info[0]), Integer.parseInt(info[1]), Integer.parseInt(info[2]));
				alignmentTypeMap.put(lineEnChPair, info[3]);
				probabilityMap.put(lineEnChPair, Double.parseDouble(info[4]));
				
				line = br.readLine();
			}
		}finally {
	        br.close();
	    }
	}
	
	public int findLink(int fPos, int ePos, int lineNumber, int N){
		
		Triple<Integer, Integer, Integer> triple = new  Triple<Integer, Integer, Integer>(lineNumber, ePos, fPos);
    	
		int h = 1;
		if (ePos < N) {
			String type = alignmentTypeMap.get(triple);

			if (type == null)
				System.out.println("nukkkkkkllll" + triple + " N " + N );
			h = tagMap.get(type);
			return h;
		}
		else return -1;
	}
	public double tProbability(String f, String e, HashMap<Pair, Double> t_table){
		Pair fe = new Pair(f,e);
		
		int V = 163303;
		if (t_ibm1_20k.containsKey(fe))
		return t_ibm1_20k.get(fe);
		else if (t_table.containsKey(fe))
			return t_table.get(fe);
		else
			return 1.0/V; 
	}
	//returns an ArrayList of alignments results of aligning through IBM1 source to target
	public ArrayList<String> print_alignment_SD_ibm1(String[][] bitext, HashMap<Pair, Double> t, boolean useAlignmentType){
		double[] tagDist = {0, 0.401, 0.264, 0.004, 0.004, 0.012, 0.205, 0.031, 0.008, 0.003, 0.086, 0.002};
		String[] linkMap = {"SEM", "FUN", "PDE", "CDE", "MDE", "GIS", "GIF", "COI", "TIN", "NTR", "MTA"};
		ArrayList<String> alignmentList = new ArrayList<String>();
		int lineNumber = 0;
		for (String[] pair : bitext){
			String[] S = pair[0].trim().split(" ");
			String[] D = pair[1].trim().split(" ");
			
			HashSet<IntPair> set = new HashSet<IntPair>();
			String lineWATag = "";
			for (int i = 0; i < S.length; i++){
				double max_ts = 0;
				int argmax = -1;
				int bestH = 0;
				for (int j = 0; j < D.length; j++){
					
					Triple<Integer, Integer, Integer> triple = new  Triple<Integer, Integer, Integer>(lineNumber, j, i);
					String type = alignmentTypeMap.get(triple);
					double first = 0;
					
					double prob = 0;
					
					int h = tagMap.get(type);
            		prob = probabilityMap.get(triple);
					first = prob;
					double s_SiDj = (1-lambda)*(first);
            		
					s_SiDj += lambda*tagDist[h];
	            	
					double t_SiDj = tProbability(S[i], D[j], t);
					
						if (t_SiDj*s_SiDj > max_ts){
							max_ts = t_SiDj*s_SiDj;
							argmax = j;
							bestH = h;
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
			
			
			if (!useAlignmentType)
				alignmentList.add(line);
			else
				alignmentList.add(lineWATag);
			
			lineNumber++;
		}
		return alignmentList;
	}
	
	public static void fScoreExperiment() throws IOException{
		
		int trainingSize = 203;//2373245;//18000;2373245
        int testSize = 19;//20357;//2357;20357
		
		String[][] trainTritext = IBM_Model1.readTritext("data-te/train.20k.seg.cln.cn","data-te/train.20k.seg.cln.en","data-te/train.20k.seg.cln.wa",trainingSize);
		String[][] testTritext = IBM_Model1.readTritext("data-te/test.seg.cln.cn","data-te/test.seg.cln.en","data-te/test.seg.cln.wa",testSize);
		
		IBM1Discriminative ibm1 = new IBM1Discriminative();
		
		ibm1.initializeCounts(trainTritext,0);
		
		
		//ibm1.setSProbabilities(ibm1.fe_count, ibm1.total_f_e_h);
		ibm1.readLogisticRegFile("logisticRegFile");
		
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainTritext);
		System.out.println("length " + trainTritext.length);
		
		System.out.println("\n WA:");			
		ArrayList<String> ibmModelAlignment = ibm1.print_alignment_SD_ibm1(testTritext, t_fe, false);
		ArrayList<String> reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		ibm1.gradeAlignmentWithType(testSize, testTritext, reference, ibmModelAlignment);
		
		System.out.println("\n WA+Tag:");
		ibmModelAlignment = ibm1.print_alignment_SD_ibm1(testTritext, t_fe, true);
		reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.wa");
		ibm1.gradeAlignmentWithTypeWAPlusTag(testSize, testTritext, reference, ibmModelAlignment);
		
	}
	
public static void augmentedModelExperimentMain() throws IOException{
		
		int trainingSize = 20357;//2373245;//18000;2373245
        int testSize = 1956;//20357;//2357;20357
		
		String[][] trainTritext = IBM_Model1.readTritext("data-te/train.20k.seg.cln.cn","data-te/train.20k.seg.cln.en","data-te/train.20k.seg.cln.wa",trainingSize);
		String[][] testTritext = IBM_Model1.readTritext("data-te/test.seg.cln.cn","data-te/test.seg.cln.en","data-te/test.seg.cln.wa",testSize);
		
		IBM1Discriminative ibm1 = new IBM1Discriminative();
		
		ibm1.initializeCountsWithoutSets(trainTritext);
		
		//ibm1.setSProbabilities(ibm1.fe_count, ibm1.total_f_e_h);
		ibm1.readLogisticRegFile("logisticRegFile");
		
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainTritext);
		System.out.println("length " + trainTritext.length);
		
		System.out.println("\n WA:");			
	/*	ArrayList<String> ibmModelAlignment = ibm1.print_alignment_SD_ibm1(testTritext, t_fe, false);
		ArrayList<String> reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		ibm1.gradeAlignmentWithType(testSize, testTritext, reference, ibmModelAlignment);
		
		System.out.println("\n WA+Tag:");
		ibmModelAlignment = ibm1.print_alignment_SD_ibm1(testTritext, t_fe, true);
		reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.wa");
		ibm1.gradeAlignmentWithTypeWAPlusTag(testSize, testTritext, reference, ibmModelAlignment);
	*/	
		trainingSize = 1020357;
		
		IBM1Discriminative ibm120KHK = new IBM1Discriminative();
		ibm120KHK.setTTable20k(t_fe);
		String[][] train20KHKTritext = IBM_Model1.readTritext("data-te/train.20k+hk.cn","data-te/train.20k+hk.en","data-te/train.20k.seg.cln.wa",trainingSize);
        ibm120KHK.initializeCountsWithoutSets(train20KHKTritext);
		System.out.println("length " + train20KHKTritext.length);
		ibm120KHK.readLogisticRegFile("logisticRegFile");
		HashMap<Pair, Double> t_fe_20KHK = ibm120KHK.EM_IBM1(ibm120KHK.f_count, ibm120KHK.e_count, ibm120KHK.fe_count, train20KHKTritext);
		
		System.out.println("\n WA:");			
		ArrayList<String> ibmModelAlignment = ibm120KHK.print_alignment_SD_ibm1(testTritext, t_fe_20KHK, false);
		ArrayList<String> reference = ibm120KHK.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		ibm120KHK.gradeAlignmentWithType(testSize, testTritext, reference, ibmModelAlignment);
		
		System.out.println("\n WA+Tag:");
		ibmModelAlignment = ibm120KHK.print_alignment_SD_ibm1(testTritext, t_fe_20KHK, true);
		reference = ibm120KHK.convertFileToArrayList("data-te/test.seg.cln.wa");
		ibm120KHK.gradeAlignmentWithTypeWAPlusTag(testSize, testTritext, reference, ibmModelAlignment);
		
	}
	public static void main(String[] args) throws IOException {
		//fScoreExperiment();
		augmentedModelExperimentMain();
	}
	
}
