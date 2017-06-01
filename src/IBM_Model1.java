import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.io.*;

public class IBM_Model1 implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	/**
	 * This class has a problem. It depends on Triple (Actually IBM1 should be independent of alignment types)
	 * change this ------------
	 */
	
	HashMap<String, Integer> f_count = new HashMap<String, Integer>();
	HashMap<String, Integer> e_count = new HashMap<String, Integer>();
	HashMap<Pair, Integer> fe_count = new HashMap<Pair, Integer>();
	
	HashMap<String,Integer> tagMap = new HashMap<String, Integer>();
	HashMap<Triple<String, String, Integer>, Double> total_f_e_h = new HashMap<Triple<String,String,Integer>, Double>();
	
		
	/**
	 * Initializes a map of tags to integers
	 */
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
	 * 
	 * @param fileName     source file name
	 * @param fileName2    target file name
	 * @param linesToRead  number of lines to read from both files to create the bitext
	 * @return bitext of source-target
	 * @throws IOException
	 */
	public static String[][] readBitext(String fileName, String fileName2, int linesToRead) throws IOException{
		BufferedReader br = null;
		BufferedReader br2 = null;
		String[][] text = new String[linesToRead][2];
		System.out.print("Reading the training bitext");
	    try {
	    	File fileDir = new File(fileName);
	    	 
			br = new BufferedReader(new InputStreamReader(new FileInputStream(fileDir), "UTF8"));
			
			File fileDir2 = new File(fileName2);
	    	 
			br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fileDir2), "UTF8"));
	        
	        String line = br.readLine();
	        String line2 = br2.readLine();
	        
	        int cnt = 0;
	        while (line != null) {
	        	text[cnt][0] = line;
	        	text[cnt][1] = line2;
	        	cnt++;
	        	//System.out.println(cnt);
	           
	            line = br.readLine();
	            line2 = br2.readLine();
	            if (cnt == linesToRead){
	            	return text;
	            }
	        }
	        
	    } finally {
	        br.close();
	        br2.close();
	    }
		System.out.println("Done");
		return text;
	}
	/**
	 * 
	 * @param fileName    source file name
	 * @param fileName2   target file name
	 * @param fileName3   gold word alignment and alignment type file name 
	 * @param linesToRead number of lines to read from both files to create the bitext
	 * @return bitext of source-target
	 * @throws IOException
	 */
	public static String[][] readTritext(String fileName, String fileName2, String fileName3, int linesToRead) throws IOException{
		BufferedReader br = null;
		BufferedReader br2 = null;
		BufferedReader br3 = null;
		
		String[][] text = new String[linesToRead][3];
		
	    try {
	    	File fileDir = new File(fileName);
	    	 
			br = new BufferedReader(new InputStreamReader(new FileInputStream(fileDir), "UTF8"));
			
			File fileDir2 = new File(fileName2);
	    	 
			br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fileDir2), "UTF8"));
			
			File fileDir3 = new File(fileName3);
	    	 
			br3 = new BufferedReader(new InputStreamReader(new FileInputStream(fileDir3), "UTF8"));
	        
	        String line = br.readLine();
	        String line2 = br2.readLine();
	        String line3 = br3.readLine();
	        
	        int cnt = 0;
	        while (line != null) {
	        	text[cnt][0] = line;
	        	text[cnt][1] = line2;
	        	text[cnt][2] = line3;
	        	cnt++;
	        	System.out.println(cnt);
	            
	            line = br.readLine();
	            line2 = br2.readLine();
	            line3 = br3.readLine();
	            if (cnt == linesToRead){
	            	return text;
	            }
	        }
	       
	    } finally {
	        br.close();
	        br2.close();
	        br3.close();
	    }
		return text;
	}
	
	
	/**
	 * initializes source count, target count and source-target count tables (maps)
	 * @param bitext  bitext of source-target
	 */
public void initializeCountsWithoutSets(String[][] bitext){
		 initializeTagMap();	
		for (String[] fewa : bitext){
			String[] F = fewa[0].split(" ");
			String[] E = fewa[1].split(" ");
			
			//Setting f_count
			for (String f : F){ 
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
			
		}
		
	}
	
	/**
	 * This method computes source and target counts as well as (source, target, alignment type) counts
	 * (f,e,h) counts are stored in total_f_e_h
	 * HMMWithAlignmentType initializes its s parameter from total_f_e_h 
	 * @param tritext
	 * @param testSize
	 */
	public void initializeCounts(String[][] tritext, int testSize){
		
		initializeTagMap();
		int sentenceNumber = 1;
		for (String[] fewa : tritext){
			String[] F = fewa[0].split(" ");
			String[] E = fewa[1].split(" ");
			String[] wa = fewa[2].split(" ");
			
			
			//Setting f_count
			for (String f : F){ 
				if (f_count.containsKey(f))
					f_count.put(f, new Integer(f_count.get(f).intValue()+1));
				else
					f_count.put(f, new Integer(1));
				
				//Setting fe_count
				for (String e : E){ 
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
			
			//setting total_f_e_h count
			if (sentenceNumber > testSize){
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
			sentenceNumber++;
			
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
	/**
	 * EM for IBM1
	 * @param sourceCount
	 * @param targetCount
	 * @param st_count
	 * @param bitext
	 * @return t parameters
	 */
	public HashMap<Pair, Double> EM_IBM1(HashMap<String, Integer> sourceCount, HashMap<String, Integer> targetCount, HashMap<Pair, Integer> st_count, String[][] bitext){
		//Pair p = new Pair("a", "b");
		HashMap<Pair, Double> t = new HashMap<Pair, Double>();
		
		double initialValue = new Double(1.0/sourceCount.size());
		
		for (Pair p: st_count.keySet()){
			t.put(p, initialValue);
		}

		Pair sdPair = new Pair("","");
		Pair newPair1 = new Pair("","");

		for (int i = 0; i < 5 ; i++){
			System.out.println("IBM iteration " + i);
			HashMap<Pair, Double> c = new HashMap<Pair, Double>();
			HashMap<String, Double> total = new HashMap<String, Double>();
			
			for (String[] pair : bitext){
				String[] S = pair[0].split(" ");
				String[] D = pair[1].split(" ");
				
				for (String s_i : S){
					double Z = 0;
					for (String d_j : D){
						//Pair sdPair = new Pair(s_i,d_j);	
						//System.out.println("t is null"+ s_i + d_j);
						sdPair.setPair(s_i, d_j);
						Z += t.get(sdPair).doubleValue();
					}
					for (String d_j : D){
						newPair1.setPair(s_i,d_j);
						//Pair newPair = new Pair(s_i,d_j);
						if (c.containsKey(newPair1))
							c.put(newPair1, c.get(newPair1).doubleValue() + t.get(newPair1).doubleValue()/Z);
						else{
							Pair newPair = new Pair(s_i,d_j);
							c.put(newPair, t.get(newPair).doubleValue()/Z);
						}
						if (total.containsKey(d_j))
							total.put(d_j, total.get(d_j).doubleValue() + t.get(newPair1).doubleValue()/Z);
						else
							total.put(d_j, t.get(newPair1).doubleValue()/Z);
					}
				}	
			}
			for (Pair sd : st_count.keySet()){
				t.put(sd, c.get(sd).doubleValue()/total.get(sd.right).doubleValue());
			}
		}
		
		return t;
	}
	/**
	 * Smoothes t(f|e) parameter by backing-off to a uniform probability 1/V
	 * @param f  source word
	 * @param e  target word
	 * @return   the somoothed t(f|e)
	 */
	public double tProbability(String f, String e, HashMap<Pair, Double> t_table){
		Pair fe = new Pair(f,e);
		
		int V = 163303;
		if (t_table.containsKey(fe))
			return t_table.get(fe);
		else
			return 1.0/V; 
	}
	//returns an ArrayList of alignments results of aligning through IBM1 source to target
	public ArrayList<String> print_alignment_SD_ibm1(String[][] bitext, HashMap<Pair, Double> t, String alignmentFile) throws IOException{
		ArrayList<String> alignmentList = new ArrayList<String>();
		BufferedWriter alignment = new BufferedWriter(new FileWriter(alignmentFile));
		for (String[] pair : bitext){
			String[] S = pair[0].split(" ");
			String[] D = pair[1].split(" ");
			
			HashSet<IntPair> set = new HashSet<IntPair>();
			for (int i = 0; i < S.length; i++){
				double max_t = 0;
				int argmax = -1;
				for (int j = 0; j < D.length; j++){
					double t_SiDj = tProbability(S[i], D[j], t);
					if (t_SiDj > max_t ){
						max_t = t_SiDj;
						argmax = j;
					}
				}
				
				set.add(new IntPair(i,argmax));
			}
			String line = "";
			for (IntPair ip : set){
				//System.out.print((ip.left+1) + "-" + (ip.right+1) + " ");
				line += (ip.left+1) + "-" + (ip.right+1) + " ";
			}
			alignmentList.add(line);
			//System.out.println();
			alignment.write(line+"\n");
		}
		alignment.close();
		return alignmentList;
	}
	/**
	 * 
	 * @param alignementTestSize  number of lines to do test on
	 * @param bitext  	 		test bitext
	 * @param reference	 		gold reference alignment. Including 1-to-many alignments.
	 * @param systemAlignment	 my output alignment
	 * @return F1 score
	 */
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
	/**
	 * 
	 * @param alignementTestSize number of lines to do test on
	 * @param bitext			 test bitext
	 * @param reference			 gold reference alignment with alignment types
	 * @param systemAlignment    my output alignment
	 * @return F1 score
	 */
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
	/**
	 * 
	 * @param alignementTestSize	 number of lines to do test on
	 * @param bitext				 test bitext
	 * @param reference				 gold reference with sure and possible links
	 * @param systemAlignment		 my output alignment
	 */
	public void gradeAlign(int alignementTestSize, String[][] bitext, ArrayList<String> reference, ArrayList<String> systemAlignment){
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
				if (index != -1)
					alignment.add(new IntPair(Integer.parseInt(pwd.substring(0,index)), Integer.parseInt(pwd.substring(index+1)) ));
				
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
				if (index != -1)
					sure.add(new IntPair(Integer.parseInt(spwd.substring(0,index)), Integer.parseInt(spwd.substring(index+1)) ));
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
		System.out.println("Precision = " + precision + "\n Recall = " + recall + "\n AER = " + aer);
	}
	/**
	 * Converts a file to an ArrayList of String (line by line)
	 * @param fileName
	 * @return  an ArrayList of String where each element of the list is a line of fileName
	 */
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
	
	public static void mainIBM1(int trainingSize, int testSize, String trainPrefix
			, String sourceLang, String targetLang, String testPrefix, String referenceFile, String alignmentFile) throws IOException{
		IBM_Model1 ibm1 = new IBM_Model1();
		//int trainingSize = 20357;//2373245;//18000;2373245
        //int testSize = 1956;//20357;//2357;20357
		
        
        //String[][] trainBitext = ibm1.readBitext("data-te/train.20k.seg.cln.cn","data-te/train.20k.seg.cln.en",trainingSize);
        //String[][] testBitext = ibm1.readBitext("data-te/test.seg.cln.cn","data-te/test.seg.cln.en",testSize);
		
		String[][] trainBitext = IBM_Model1.readBitext(trainPrefix+"."+sourceLang, trainPrefix+"."+targetLang, trainingSize);
		String[][] testBitext = IBM_Model1.readBitext(testPrefix+"."+sourceLang, testPrefix+"."+targetLang, testSize);
        //trainingSize = 1020357;
        //String[][] train20KHKBitext = ibm1.readBitext("data-te/train.20k+hk.cn","data-te/train.20k+hk.en",trainingSize);
		
		
		ibm1.initializeCountsWithoutSets(trainBitext);
		System.out.println("length " + trainBitext.length);
		
		HashMap<Pair, Double> t_fe = ibm1.EM_IBM1(ibm1.f_count, ibm1.e_count, ibm1.fe_count, trainBitext);
		
		ArrayList<String> ibmModelAlignment = ibm1.print_alignment_SD_ibm1(testBitext, t_fe, alignmentFile);
		//ArrayList<String> reference = ibm1.convertFileToArrayList("data-te/test.seg.cln.gold.wa");
		if (referenceFile != null){
			ArrayList<String> reference = ibm1.convertFileToArrayList(referenceFile);
			ibm1.gradeAlignmentWithType(testSize, testBitext, reference, ibmModelAlignment);
		}
	
	}
}

class Pair implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	String left;
	String right;
	
	public Pair(String left, String right){
		this.left = left;
		this.right = right;
	}
	public void setPair(String left, String right){
		this.left = left;
		this.right = right;
	}
	public boolean equals(Object obj){
		if (!(obj instanceof Pair))
            return false;
		if (obj == this)
            return true;
		Pair p = (Pair) obj;
		if (p.left.equals(this.left) && p.right.equals(this.right))
			return true;
		return false;
	}
    public int hashCode() {
		 final int prime = 31;
         int result = 1;
         result = prime * result
                      + ((left == null) ? 0 : left.hashCode());
         result = prime * result + ((right == null) ? 0 : right.hashCode());
         return result;
	 }
    public String toString(){
    	return "(" + left + "," + right + ")";
    }
}

class IntPair{
	int left;
	int right;
	
	public IntPair(int left, int right){
		this.left = left;
		this.right = right;
	}
	public boolean equals(Object obj){
		if (!(obj instanceof IntPair))
            return false;
		if (obj == this)
            return true;
		IntPair p = (IntPair) obj;
		if (p.left == this.left && p.right == this.right)
			return true;
		return false;
	}
	 public int hashCode() {
		 final int prime = 31;
         int result = 1;
         result = prime * result + left;
         result = prime * result + right;
         return result;
	 }
	 public String toString(){
		 return "(" + this.left + ", " + this.right + ")";
	 }
    
}
