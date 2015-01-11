import gnu.getopt.Getopt;
import gnu.getopt.LongOpt;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Random;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class IsoformExpressionLearning {

	private static final double EPSILON = 2.2204e-16;
	private static final double NOISE_STD = 0.0001;
	private static int VERBOSE_LEVEL = 10;
	private static final boolean DEBUG = true;
	private static int TOTAL_READS_MAPPED = 0;

	private static PrintStream ERR_STREAM;


	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		String filename = "",transcripts_filename="";
		boolean printUsage = false;
		LongOpt[] longopts = new LongOpt[1];
		longopts[0] = new LongOpt("help", LongOpt.NO_ARGUMENT, null, 'h');

		Getopt g = new Getopt("IsoformExpressionLearning", args, "T:F:V:N:h",longopts);
		int c;

		while ((c = g.getopt()) != -1)
		{
			switch(c)
			{
			case 'h':
				printUsage = true;
				break;
			case 'F':
				filename = g.getOptarg();
				break;
			case 'T':
				transcripts_filename = g.getOptarg();
				break;
			case 'V':
				VERBOSE_LEVEL = Integer.parseInt(g.getOptarg());
				break;
			case 'N':
				TOTAL_READS_MAPPED = Integer.parseInt(g.getOptarg());
				break;

			case '?':
				printUsage = true;
				break; 
				//
			default:
				printUsage = true;
			}
		}
		printUsage = printUsage || filename.equals("") || transcripts_filename.equals("") || TOTAL_READS_MAPPED==0; 

		if (printUsage)
		{
			System.err.println("");
			System.err.println("########################################################################################");
			System.err.println("#");
			System.err.println("# Required:");
			System.err.println("#  -T  <string>  transcript fasta file");
			System.err.println("#  -F  <string>  prefix for reads2Transcripts file");
			System.err.println("#  -N  <int>     total number of reads mapped");
			System.err.println("#");
			System.err.println("# Optional:");
			System.err.println("#  -V <int>                           verbosity level ");
			System.err.println("#                                        (default: 10 - progress of method + some stats)");
			System.err.println("#                                        (20 - maximum verbosity)");
			System.err.println("#  ");
			System.err.println("#");
			System.err.println("########################################################################################");
			System.err.println("");
			System.exit(1);

		}
		ERR_STREAM = new PrintStream(new FileOutputStream(filename + ".ExpEM.err"));
		debugMes("Started",10);
		DecimalFormat df = new DecimalFormat("#.###");

		// initialize
		Vector<String> TranscriptNames = readTranscripts(transcripts_filename);
		debugMes("Read "+TranscriptNames.size()+" transcripts",10);

		HashMap<String,Integer> TranscriptIndices = new HashMap<String,Integer>();
		for (Integer i=0; i<TranscriptNames.size(); i++){
			TranscriptIndices.put(TranscriptNames.get(i), i);
		}

		Vector<Double> totalReadsUsedInThisComp = new Vector<Double>();
		initVecToZero(totalReadsUsedInThisComp, 1);

		Vector<Double> theta = initTheta(filename+".reads2Transcripts.txt",TranscriptIndices,totalReadsUsedInThisComp);
		debugMes("theta = "+vectorToString(theta),20);
		//		Vector<Double> theta = initThetaUniformly(TranscriptIndices);
		Vector<Double> like_term2 = new Vector<Double>();
		initVecToZero(like_term2, 1);
		Vector<Double> E_Nt = calcE_Nt(filename+".reads2Transcripts.txt",theta,TranscriptIndices,like_term2);

//		debugMes("E_Nt = "+E_Nt+";",10);
		debugMes("sum E_Nt = "+sumVector(E_Nt)+" and N="+totalReadsUsedInThisComp.get(0),20);
		Vector<Double> lengthsToTheMinus1 = calcLengths(TranscriptNames);

		Vector<Double> initial_EFPKM = calcEFPKM(E_Nt,lengthsToTheMinus1);

		Vector<Double> likeHist = new Vector<Double>();
		likeHist.add(calcLikelihood(E_Nt,theta,like_term2));
		// run iteratively
		for (int iter=0; iter<=40 ; iter++)
		{
			debugMes("iteration = "+iter+", log likelihood = "+df.format(likeHist.lastElement()),10);
			// E-step
			if (iter>0)
			{
				initVecToZero(like_term2, 1);
				E_Nt = calcE_Nt(filename+".reads2Transcripts.txt",theta,TranscriptIndices,like_term2);
			}
			debugMes("E_t = "+vectorToString(E_Nt),20);
			// M-step
			theta = calcTheta(E_Nt);
			debugMes("theta = "+vectorToString(theta),20);
			//// re-calc the likelihood
			likeHist.add(calcLikelihood(E_Nt,theta,like_term2));
		}

		Vector<Double> final_EFPKM = calcEFPKM(E_Nt,lengthsToTheMinus1);

		printOutput(filename,likeHist,initial_EFPKM,final_EFPKM,TranscriptNames);

		debugMes("Done",10);
	}





	private static void printOutput(String filename, Vector<Double> likeHist,
			Vector<Double> initial_EFPKM, Vector<Double> final_EFPKM,
			Vector<String> TranscriptNames) throws IOException {

		DecimalFormat df = new DecimalFormat("#.###");
		PrintStream p;
//		p= new PrintStream(new FileOutputStream(filename + "_like_hist.txt"));
//		printVector(p,likeHist);
//		p.close();

//		p= new PrintStream(new FileOutputStream(filename + "_initial_EFPKM.txt"));
//		printVector(p,initial_EFPKM);
//		p.close();
//
//		p= new PrintStream(new FileOutputStream(filename + "_final_EFPKM.txt"));
//		printVector(p,final_EFPKM);
//		p.close();

		p= new PrintStream(new FileOutputStream(filename + "_out.txt"));
		p.println("TranscriptName\tinitialEFPKM\tfinalEFPKM");
		for (int i=0; i<TranscriptNames.size(); i++)
		{
			String tName = TranscriptNames.get(i);
			String init_EFPKM = df.format(initial_EFPKM.get(i));
			String fin_EFPKM = df.format(final_EFPKM.get(i));
			p.println(tName+"\t"+init_EFPKM+"\t"+fin_EFPKM);
			debugMes("transcript "+ tName + " had "+init_EFPKM+" and now has "+fin_EFPKM,10);
		}
		p.close();

	}

//	/**
//	 * print the given vector as a column vector
//	 * @param p
//	 * @param vec
//	 */
//	private static void printVector(PrintStream p, Vector<Double> vec) {
//		DecimalFormat df = new DecimalFormat("#.####");
//
//		for (int i=0; i<vec.size(); i++)
//			p.println(df.format(vec.get(i)));
//	}

	/**
	 * return the formatted toString
	 * @param vec
	 * @return
	 */
	private static String vectorToString(Vector<Double> vec) {
		DecimalFormat df = new DecimalFormat("#.####");
		String res = "[";
		for (int i=0; i<vec.size(); i++)
			res = res + df.format(vec.get(i)) + ", ";
		res = res.substring(0, res.length()-2);
		res = res + "]";
		
		return res;
	}
	
	
	
	/**
	 * Read the transcript names from the given filename
	 * @param string
	 * @return
	 * @throws IOException 
	 */
	private static Vector<String> readTranscripts(String file) throws IOException {
		Vector<String> transNames = new Vector<String>();
		BufferedReader fileB = 	new BufferedReader(new FileReader(file)); 
		String l;

		while (fileB.ready())
		{
			l = fileB.readLine();
			if (l.startsWith(">")) //read only transcript names
				transNames.add(l.substring(1));
		}		

		return transNames;
	}


	/**
	 * Go over the reads to transcript mapping, and initialize theta according to how many reads were mapped to each transcript
	 * @param file
	 * @param transcriptIndices 
	 * @param totalReadsUsedInThisComp - enter here the total number of reads mapped (N)
	 * @return
	 * @throws IOException 
	 */
	private static Vector<Double> initTheta(String file, HashMap<String,Integer> transcriptIndices, Vector<Double> totalReadsUsedInThisComp) 
	throws IOException {
		int numT = transcriptIndices.size();
		Vector<Double> theta = new Vector<Double>();
		initVecToZero(theta,numT);
		BufferedReader fileB = 	new BufferedReader(new FileReader(file)); 
		String l,read,trans;
		Vector<Double> tmpPerReadName = new Vector<Double>();
		initVecToZero(tmpPerReadName,numT);

		String prevRead = "";
		while (fileB.ready())

		{
			l = fileB.readLine();
			//			61G9EAAXX100520:5:100:10000:5699        comp3190_c716_seq2_FPKM_all:11.742_FPKM_rel:2.677_len:2401_path:[40631,40349,40358,2726]        376-444 558-444 -0.0305
			//			61G9EAAXX100520:5:100:10000:5699        comp3190_c716_seq4_FPKM_all:13.841_FPKM_rel:6.304_len:2888_path:[0,40349,40358,2726]    863-931 1045-931        -0.0305
			//			61G9EAAXX100520:5:100:10000:7578        comp754_c304_seq1_FPKM_all:27.333_FPKM_rel:2.886_len:901_path:[28182,65285,65293,65295,517,535,538,1140,1147,1289]      417-485 535-485 -0.0111
			//			readName transName locRead2 locRead2 prob			
			String[] fields = l.split("\t");
			read = fields[0];
			trans = fields[1];

			int t_id = transcriptIndices.get(trans);

			if (read.equals(prevRead))
			{
				// deal with this read
				tmpPerReadName.set(t_id,tmpPerReadName.get(t_id)+1);

			}else{
				// new read, close the prevRead
				if (!prevRead.equals("")){

					totalReadsUsedInThisComp.set(0, totalReadsUsedInThisComp.get(0)+1);

					normVector(tmpPerReadName);
					addToVector(theta,tmpPerReadName);

					// prepare for the next read
					initVecToZero(tmpPerReadName,numT);
				}
				// deal with this read
				tmpPerReadName.set(t_id,tmpPerReadName.get(t_id)+1);

				prevRead = read;

			}

		}
		// handle the last read
		totalReadsUsedInThisComp.set(0, totalReadsUsedInThisComp.get(0)+1);

		normVector(tmpPerReadName);
		addToVector(theta,tmpPerReadName);

		normVector(theta);

		//add noise
		addNormalNoise(theta,NOISE_STD);

		return theta;
	}

	/**
	 * given the reads to transcript mapping, and theta, calc E_Nt: (which basically calcs the posterior)
	 * $$E(N_t) = \sum_i \frac{ P(r_i | t)\cdot \theta_t} {\sum_{t'} P(r_i | t')\cdot \theta_{t'}}$$
	 * @param string
	 * @param theta
	 * @return
	 */
	private static Vector<Double> calcE_Nt(String file, Vector<Double> theta,
			HashMap<String,Integer> transcriptIndices, Vector<Double> likeTerm2) 
			throws IOException {
		int numT = theta.size();
		Vector<Double> E_Nt = new Vector<Double>();
		initVecToZero(E_Nt,numT);
		BufferedReader fileB = 	new BufferedReader(new FileReader(file)); 
		String l,read,trans;
		Double prob; 

		Vector<Double> tmpPerReadName = new Vector<Double>();
		initVecToZero(tmpPerReadName,numT);
		String prevRead = "";
		while (fileB.ready())

		{
			l = fileB.readLine();
			//			61G9EAAXX100520:5:100:10000:5699        comp3190_c716_seq2_FPKM_all:11.742_FPKM_rel:2.677_len:2401_path:[40631,40349,40358,2726]        376-444 558-444 -0.0305
			//			61G9EAAXX100520:5:100:10000:5699        comp3190_c716_seq4_FPKM_all:13.841_FPKM_rel:6.304_len:2888_path:[0,40349,40358,2726]    863-931 1045-931        -0.0305
			//			61G9EAAXX100520:5:100:10000:7578        comp754_c304_seq1_FPKM_all:27.333_FPKM_rel:2.886_len:901_path:[28182,65285,65293,65295,517,535,538,1140,1147,1289]      417-485 535-485 -0.0111
			//			readName transName locRead2 locRead2 prob			
			String[] fields = l.split("\t");
			read = fields[0];
			trans = fields[1];
			prob = Math.pow(10, Double.parseDouble(fields[4]));

			int t_id = transcriptIndices.get(trans);

			if (read.equals(prevRead))
			{
				// deal with this read
				tmpPerReadName.set(t_id,tmpPerReadName.get(t_id)+prob);

			}else{
				// new read, close the prevRead
				if (!prevRead.equals("")){
					//tmpPerReadName holds now P(r_i | T)
					likeTerm2.set(0, likeTerm2.get(0)+sumLogVector(tmpPerReadName));

					dotMultiply(tmpPerReadName,theta);
					normVector(tmpPerReadName);

					addToVector(E_Nt,tmpPerReadName);

					// prepare for the next read
					initVecToZero(tmpPerReadName,numT);
				}
				// deal with this read
				tmpPerReadName.set(t_id,tmpPerReadName.get(t_id)+prob);

				prevRead = read;

			}

		}

		//handle the last read
		//tmpPerReadName holds now P(r_i | T)
		likeTerm2.set(0, likeTerm2.get(0)+sumLogVector(tmpPerReadName));

		dotMultiply(tmpPerReadName,theta);
		normVector(tmpPerReadName);

		addToVector(E_Nt,tmpPerReadName);

		return E_Nt;
	}

	/**
	 * Return 1/len
	 * @param transcriptNames
	 * @return
	 */
	private static Vector<Double> calcLengths(Vector<String> transcriptNames) {
		Vector<Double> lengthsToTheMinus1 = new Vector<Double>();
		initVecToZero(lengthsToTheMinus1, transcriptNames.size());
		//extract lengths, then hold 1/|t|
		Pattern p = Pattern.compile("\\s*_len:([0-9]+)\\s*");
		Matcher m;
		String name="";
		Double len = new Double(0);
		for (int i=0 ; i<transcriptNames.size(); i++)
		{
			name = transcriptNames.get(i);
			m = p.matcher(name);
			if (m.find())
			{
				len = Double.parseDouble(m.group(1));
				lengthsToTheMinus1.add(i,1/len);
			}
			else
				debugMes("We have a problem, no length was found in "+name,10);
		}

		return lengthsToTheMinus1;
	}


	/**
	 * Given E_Nt, calc the posterior (10^9 /(|t|*N) * E(N_t) )
	 * @param E_Nt
	 * @param transcriptNames
	 * @param totalReadsMapped - N
	 * @return
	 */
	private static Vector<Double> calcEFPKM(Vector<Double> E_Nt,
			Vector<Double> lengthsToTheMinus1) {
		// post holds the normalized E(N_t)
		Vector<Double> post = new Vector<Double>(E_Nt);
		//FIXME !!!
//		normVector(post); 

		// factor holds 10^9 /(|t| * N)
		Vector<Double> factor = new Vector<Double>(lengthsToTheMinus1);

		Vector<Double> tmp = new Vector<Double>();
		initVecToZero(tmp, factor.size());

		Double tenToTheNinthByN = Math.pow(10, 9)/TOTAL_READS_MAPPED;

		addScalarToVector(tmp, tenToTheNinthByN);
		dotMultiply(factor, tmp);
		dotMultiply(post, factor);

		return post;
	}


	/**
	 * Given E_Nt, calc theta:
	 * $$\hat{\theta_t} = \frac{E(N_t)} {\sum_{t'} E(N_{t'})}$$
	 * @param E_Nt
	 * @return
	 */
	private static Vector<Double> calcTheta(Vector<Double> E_Nt) {
		Vector<Double> theta = new Vector<Double>();
		initVecToZero(theta, E_Nt.size());
		addToVector(theta, E_Nt);
		normVector(theta);

		return theta;
	}

	/**
	 * Given E_Nt and theta, calc the likelihood
	 * like = like_term1 + like_term2;
	 * like_term1 = sum(E_t.*log(theta)) 
	 * like_term2 = sum(log(P(R|T)+eps));
	 * @param E_Nt
	 * @param theta
	 * @param likeTerm2 
	 * @return
	 */
	private static double calcLikelihood(Vector<Double> E_Nt,
			Vector<Double> theta, Vector<Double> likeTerm2) {
		Vector<Double> tmpV = logVector(theta);
		dotMultiply(tmpV, E_Nt);
		double likeTerm1 = sumVector(tmpV);
		return likeTerm1+likeTerm2.get(0);
	}

	/**
	 * return log(vec)
	 * @param vec
	 * @return
	 */
	private static Vector<Double> logVector(Vector<Double> vec) {
		Vector<Double> res = new Vector<Double>();
		res.setSize(vec.size());
		for (int i=0; i<vec.size(); i++){
			double t = vec.get(i);
			if (t==0)
				res.set(i, Math.log10(EPSILON));
			else
				res.set(i, Math.log10(t));
		}
		return res;

	}

	/**
	 * sum(vec)
	 * @param vec
	 * @return
	 */
	private static double sumVector(Vector<Double> vec) {
		double sum = 0;
		for (Double t : vec){
			sum+=t;
		}
		return sum;
	}

	/**
	 * return sum(log10(vec))
	 * @param vec
	 * @return
	 */
	private static Double sumLogVector(Vector<Double> vec) {
		Double sum = new Double(0);
		for (Double t : vec){
			if (t==0)
				sum+=Math.log10(EPSILON);
			else
				sum+=Math.log10(t);
		}
		return sum;

	}


	/**
	 * vec = vec./sum(vec)
	 * @param vec
	 */
	private static void normVector(Vector<Double> vec) {
		double sum = sumVector(vec);
		if (sum>0)
			for (int i=0; i<vec.size(); i++)
				vec.set(i, vec.get(i)/sum);


	}


	/**
	 * vec += vecToBeAdded
	 * @param vec
	 * @param vecToBeAdded
	 */
	private static void addToVector(Vector<Double> vec, Vector<Double> vecToBeAdded) {
		if (vec.size()!=vecToBeAdded.size())
			assert(true);

		for (int i=0; i<vec.size(); i++){
			vec.set(i, vec.get(i)+vecToBeAdded.get(i));
		}
	}

	/**
	 * vec = vec+scalar
	 * @param vec
	 * @param scalar
	 */
	private static void addScalarToVector(Vector<Double> vec, double scalar) {
		for (int i=0; i<vec.size(); i++){
			vec.set(i, vec.get(i)+scalar);
		}

	}

	/**
	 * set size of vec to size, and init with zeros
	 * @param vec
	 * @param size
	 */
	private static void initVecToZero(Vector<Double> vec, int size) {
		vec.clear();
		vec.setSize(size);
		for (int i=0; i<vec.size(); i++)
			vec.set(i, new Double(0));

	}


	/**
	 * add a normally distributed noise, with the given std, to the given vec.
	 * @param vec
	 * @param std
	 */
	private static void addNormalNoise(Vector<Double> vec,double std) {
		Random generator = new Random();
		double noise;
		for (int i=0; i<vec.size(); i++)
		{
			noise = generator.nextGaussian()*std;
			vec.set(i, Math.max(0,vec.get(i)+noise));
		}
	}

	/**
	 * vec1 = vec1.*vec2
	 * @param vec1
	 * @param vec2
	 */
	private static void dotMultiply(Vector<Double> vec1,
			Vector<Double> vec2) {
		for (int i=0; i<vec1.size(); i++)
			vec1.set(i, vec1.get(i)*vec2.get(i));

	}

	//	/**
	//	 * Init theta uniformly
	//	 * @param transcriptIndices
	//	 * @return
	//	 */
	//	private static Vector<Double> initThetaUniformly(HashMap<String,Integer> transcriptIndices){
	//		int numT = transcriptIndices.size();
	//		Vector<Double> theta = new Vector<Double>();
	//		initVecToZero(theta,numT);
	//		addScalarToVector(theta, 1);
	//		normVector(theta);
	//		return theta;
	//	}


	/**
	 * print out the given error message, only if DEBUG=true
	 * @param mes Message
	 */
	private static void debugMes(String mes, int verbosityLevel)
	{
		if (DEBUG && verbosityLevel<=VERBOSE_LEVEL)
			ERR_STREAM.println(mes);
	}
}
