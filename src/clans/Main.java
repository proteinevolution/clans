package clans;

import java.io.*;
import java.text.ParseException;
import java.util.*;

public class Main {

	public static void main(String[] args) {

		if (args.length == 0) {
			docalc = false;

		} else if (args.length == 1) {
			if (args[0].charAt(0) != '-') {
				input_filename = args[0];
				docalc = false;
				args = new String[0];
			}
		}
		String[] inargs = args;

		// check if a config file is set
		if (!get_configuration_file_argument(args)) {
			System.err.println("unable to check args.");
			return;
		}

		// read configuration file if it exists
		File testfile = new File(conffilename);
		if (testfile.exists()) {
			if (!parse_configuration_file(conffilename)) {
				System.err.println("unable to read conffile " + conffilename + "; using defaults");
			}
		} else {
			System.err.println("Warning: " + conffilename
					+ " does not exists, will be using default and command-line options only");
		}

		if (!parse_arguments(inargs)) {
			System.err.println("unable to read args; exiting program.");
			print_usage_help();
			return;
		}

		if (!start_computations()) {
			System.err.println("Error in start_computations().");
		}
		return;
	}

	static StringBuffer errbuff = new StringBuffer();

	static String conffilename = "clans.conf"; // file name of the config file
	static String infilename = "stdin"; // default input
	static String cmd = "";// command to prepend to anything the system executes (i.e. nice -19)

	static String blastpath = "blastp ";// command necessary to start blast (if path location is needed don't forget it)
	static int blastblocks = 50;// the number of sequences to do at the same time in one blast
	static boolean addblastvbparam = true;// check to see whether I have more sequences than blast would normally return
											// hits for
	static String formatdbpath = "makeblastdb -dbtype prot";// command needed to execute formatdb (since blast+ 2.2.26
															// -dbtype is no longer an optional entry)
	static String[] referencedb = new String[0];// holds the databases to blast against to generate psiblast profiles

	static boolean skipcheckdone = true; // check for a DONE in tmpblasthsp and then skip all further checks (if false)
	static double eval = 10;// default maximum evalue to accept for hsp
	static double pval = 0.1;// default maximum pvalue to accept for hsp
	static float coverage = (float) 0;// necessary minimal coverage for blast hsp
	static float scval = (float) -1;// necessary minimal blast score/collumn for blast hsp
	static float ident = (float) 0;// necessary minimal identity to query for blast hsp
	static float minconf = 6;// what is the minimal confidence value to take as significant
	static int verbose = 1;// verbose mode
	static int cpu = 1; // how many cpu's should I use
	static boolean useallrounds = false;// do you want to use all or oly the last psiblast round for computing
										// confidences?
	static String saveblastname = "tmpblasthsp.txt";// name of blast savefile
	static boolean readblast = true;// do I want to read blasts from a savefile?
	static boolean lowmem = false;// do I want to temporarily save results to hdd?
	static boolean save_intermediate_results = false;// do I want to save my sequence positions in 3d?
	static boolean docalc = true;// do I want to do calculations or just open a clustertest window?
	static boolean nographics = false;// do I want to NOT start the graphical interface?
	static String savetoname = null;// define a savefile to save results to (only if used in conjunction with -dorounds
									// and -load)
	static boolean initialize = false; // if true will initialize the clustermap upon loading
	static int dorounds = -1; // how many rounds to cluster by (only if used in conjunction with -load)

	// variables used for adding new sequences to an already present dataset
	static String olddata = "";
	static String newseqs = "";
	static String input_filename = null;
	static boolean enrichseqs = false;
	static double gatherseqseval = 1e-10;
	static double rmseqseval = 1e-25;
	static int maxenrichseqsnum = -1;
	static int exhaustive = 1;// if I add new sequences do the fast version and only look for one way blast hits!

	static void print_usage_help() {
		System.out.println("USAGE: java -jar clans.jar [options]");
		System.out.println("If a outOfMemoryError occurs, try running it via java -Xmx###m -jar programname options");
		System.out.println("where ### is the number of megabytes of ram you are willing to allocate to this process");
		System.out.println("-------------------OPTIONS--------------------");
		System.out.println("-? print this information");
		System.out.println("-conf name of configuration file (def: clans.conf)");
		System.out.println("-infile name of input file");
		System.out.println("-load name-of-savefile");
		System.out.println("-rounds (int) (def:" + dorounds
				+ ") how many rounds to cluster for (only used in conjunction with -load)");
		System.out
				.println("-saveto String where to save the results to (only used in conjunction with -load and -rounds)");
		System.out.println("-initialize boolean (t/F) if true randomly initializes the graph in runs without GUI");
		System.out.println("-loadalt name-of-alternate-format-savefile");
		System.out.println("-lowmem t/F (doesn't do much at the moment)");
		System.out.println("-cmd String to prepend to all commands (i.e nice (unix) cmd (windows)) (def: \"\")");
		System.out.println("-blastpath \"path to blast [options]\" (def: " + blastpath + ")");
		System.out.println("-blastblocks \"number of sequences to search in one blast\" (def: " + blastblocks + ")");
		System.out
				.println("-addblastvb \"check to see whether you have more than the default 250/500 sequences blast returns hits for\" (def:"
						+ addblastvbparam + ")");
		System.out.println("-formatdbpath \"path to formatdb executable\" (def: " + formatdbpath + ")");
		System.out.println("-referencedb \"databases to psiblast against to generate the profile\"");
		System.out
				.println("-skipcheckdone \"check for a DONE in the tmpblasthsp file and skip further checks if found\" (def:"
						+ skipcheckdone + ")");
		System.out.println("-eval maximum e-value to collect hits for (def: 10)");
		System.out.println("-pval maximum p-value to collect hits for (def: 0.1)");
		System.out
				.println("-scval minimum score/query_length to collect hits for (def: -1);(if >=0 it disables the P-value and E-value filters and returns the scores normalized from 0-1 as attraction values!)");
		System.out.println("-verbose verbosity of the program (def:1)");
		System.out.println("-cpu number of CPU to use (def: 1)");
		System.out.println("-readblast T/f read former blast results");
		System.out.println("-savepos t/F save position data during each round (in file <INPUT_FILE>.savepos)");
		System.out.println("-docalc T/f do calculation or just load interface");
		System.out
				.println("-nographics t/F only do the blast runs; no graphical interface, results are saved to tmpblasthsp.txt");
		System.out.println("-olddata name of savefile (def: \"\")");
		System.out
				.println("-newseqs filename with new sequences to add to cluster (need to specify olddata) (def: \"\")");
		System.out.println("-enrichseqs t/F take newseqs as such or enrich with close relatives (blast/psiblast)");
		System.out
				.println("-gatherseqseval gather sequences for blast hits up to this evalue (in enrichment) (def: 1e-10)");
		System.out.println("-rmseqseval remove sequences more similar than rmseqseval (in enrichment) (def: 1e-20)");
		System.out
				.println("-maxenrichseqsnum get at most this number of sequences for each query (in enrichment) (def: unlimited)");
		System.out.println("-exhaustive (number) 0=one way search; 1=backvalidation; 2=redo all blast runs (def: 1)");
		System.out.println("            when adding sequences, calculate the pairwise blast values by(see above)");
		System.out.println("-------------------OPTIONS--------------------");
	}

	static void print_settings() {
		System.out.println("------------------SETTINGS-------------------");
		System.out.println("conffile=" + conffilename);
		System.out.println("infilename=" + infilename);
		System.out.println("loadname=" + input_filename);

		if (input_filename != null && dorounds >= 0) {
			System.out.println("rounds=" + dorounds);
			System.out.println("saveto=" + savetoname);
		}

		System.out.println("cmd=" + cmd);

		System.out.println("blastpath=" + blastpath);
		System.out.println("blastblocks=" + blastblocks);
		System.out.println("addblastvb=" + addblastvbparam);
		System.out.println("formatdbpath=" + formatdbpath);

		System.out.print("referencedb: ");
		for (int i = 0; i < referencedb.length; i++) {
			System.out.print(referencedb[i] + "; ");
		}

		System.out.println();
		System.out.println("skipcheckdone=" + skipcheckdone);
		System.out.println("eval=" + String.valueOf(eval));
		System.out.println("pval=" + String.valueOf(pval));
		System.out.println("scval=" + String.valueOf(scval));
		System.out.println("verbose=" + String.valueOf(verbose));
		System.out.println("cpu=" + String.valueOf(cpu));
		System.out.println("lowmem=" + lowmem);

		System.out.println("savepos=" + save_intermediate_results);
		System.out.println("docalc=" + docalc);

		System.out.println("nographics=" + nographics);

		System.out.println("readblast=" + readblast);
		System.out.println("olddata=" + olddata);
		System.out.println("newseqs=" + newseqs);
		System.out.println("enrichseqs=" + enrichseqs);
		if (enrichseqs) {
			System.out.println("gatherseqseval=" + gatherseqseval);
			System.out.println("rmseqseval=" + rmseqseval);
			System.out.println("maxenrichseqsnum=" + maxenrichseqsnum);
		}

		System.out.println("exhaustive=" + exhaustive + "; 0=one way search; 1=backvalidation; 2=redo all blast runs");
		System.out.println("------------------SETTINGS-------------------");
	}

	/**
	 * coordinates data preparation and starts the GUI
	 * 
	 * @return true if successful
	 */
	static boolean start_computations() {
		if (verbose > 0) {
			print_settings();
		}

		if (docalc) {

			if (olddata.length() == 0) { // generate similarities by running
											// BLAST on the input fasta file
				AminoAcidSequence[] sequences = AlignmentHandling.parse_fasta_format(infilename);

				if (verbose > 3) {
					System.out.println("sequences read:");
					for (int i = 0; i < sequences.length; i++) {
						System.out.println(i + " " + sequences[i].name);
						System.out.println(i + " " + sequences[i].seq);
					}
				}

				int sequence_number = sequences.length;

				if (sequence_number <= 1) {
					System.err.println("One or less sequences read, nothing to do.");
					return true;
				}

				// now set up a vector array that will hold the blast hsp's
				HashMap<String, Integer> sequence_name_internal_mapping = new HashMap<String, Integer>(
						(int) (sequence_number / 0.84), (float) 0.85);// mapping of names to array indices
				String[] sequence_names = new String[sequence_number];
				for (int i = 0; i < sequence_number; i++) {
					sequence_names[i] = sequences[i].name;
					sequences[i].name = new String("sequence" + i); // assign internally used names
					sequence_name_internal_mapping.put(sequences[i].name, new Integer(i));
					sequences[i].seq = sequences[i].seq.toUpperCase();
				}

				// run blast
				boolean isblastplus = true;
				if (blastpath.contains("blastall") || blastpath.contains("blastpgp")) {
					isblastplus = false;
				}
				searchblastv2 mysearchblast = new searchblastv2(errbuff, addblastvbparam, isblastplus, cpu,
						blastblocks, cmd, blastpath, formatdbpath, eval, pval, coverage, scval, ident, verbose,
						saveblastname, readblast, sequence_name_internal_mapping);
				MinimalHsp[] blasthits = mysearchblast.gethits(sequences);
				mysearchblast = null;

				// start the GUI
				if (nographics == false) {
					System.out.println("...reading data");
					ClusterData myclusterdata = new ClusterData(blasthits, sequences, sequence_names,
							sequence_name_internal_mapping, eval, pval, scval, verbose, cpu, save_intermediate_results,
							cmd, blastpath, addblastvbparam, formatdbpath, referencedb, errbuff, input_filename);
					myclusterdata.roundslimit = dorounds;// set the limit of how
															// often to run this
					ClusteringWithGui myclusterer = new ClusteringWithGui(myclusterdata);
					myclusterer.setVisible(true);

				} else {
					System.out.println("DONE. To visualize results restart program with the -nographics F option.");
				}

			} else if (newseqs.length() > 0) { // load data from olddata and add data from newseqs
				System.out.println("reading old data");
				
				
				saverunobject readdata;
				try {
					readdata = ClusterData.load_run_from_file(new java.io.File(olddata));
					
				} catch (FileNotFoundException e) {
					System.err.println("file not found: " + olddata);
					return false;
					
				} catch (ParseException e) {
					System.err.println("line " + e.getErrorOffset() + ": " + e.getMessage());
					return false;
					
				} catch (IOException e) {
					System.err.println(e.getMessage());
					return false;
				}

				System.out.println("reading new sequences");
				AminoAcidSequence[] newaln = AlignmentHandling.read(newseqs);

				if (enrichseqs) {
					if (referencedb.length == 0) {
						System.err.println("ERROR, no referencedb specified, unable to enrich dataset, skipping.");
					} else {
						System.out.println("starting sequences=" + newaln.length);
						newaln = new enrichutils().enrich(newaln, cmd, blastpath, formatdbpath, referencedb, cpu,
								gatherseqseval, rmseqseval, maxenrichseqsnum);
						System.out.println("enriched sequences=" + newaln.length);
					}
				}

				// now copy the position of the points and assign random positions to the new sequences
				int newelements = newaln.length;
				int readelements = readdata.inaln.length;
				int allelements = newelements + readelements;

				String[] sequence_names = new String[allelements];
				AminoAcidSequence[] sequences = new AminoAcidSequence[allelements];
				HashMap<String, Integer> sequence_name_internal_mapping = new HashMap<String, Integer>();

				float[][] allposarr = new float[allelements][3];
				for (int i = 0; i < readelements; i++) {
					sequence_names[i] = readdata.inaln[i].name;
					readdata.inaln[i].name = new String("sequence" + i);
					sequence_name_internal_mapping.put(readdata.inaln[i].name, new Integer(i));
					sequences[i] = readdata.inaln[i];
					allposarr[i][0] = readdata.posarr[i][0];
					allposarr[i][1] = readdata.posarr[i][1];
					allposarr[i][2] = readdata.posarr[i][2];
				}

				java.util.Random rand = new java.util.Random(System.currentTimeMillis());
				int[] newnumarr = new int[newelements];
				for (int i = 0; i < newelements; i++) {
					sequence_names[readelements + i] = newaln[i].name;
					newaln[i].name = new String("sequence" + (readelements + i));
					sequence_name_internal_mapping.put(newaln[i].name, new Integer(readelements + i));
					newnumarr[i] = readelements + i;
					sequences[readelements + i] = newaln[i];
					allposarr[readelements + i][0] = rand.nextFloat();
					allposarr[readelements + i][1] = rand.nextFloat();
					allposarr[readelements + i][2] = rand.nextFloat();
				}

				// run BLAST
				searchblast mysearchblast = new searchblast(errbuff, addblastvbparam);
				double mypval = readdata.pval;
				float mymaxmove = readdata.maxmove;
				MinimalHsp[] newblasthits = mysearchblast.gethits(readdata.inaln, readdata.blasthits, newaln, cmd,
						formatdbpath, blastpath, cpu, eval, pval, coverage, scval, ident, verbose,
						sequence_name_internal_mapping, useallrounds, lowmem, referencedb, exhaustive, readblast, true);

				// add only the new matches
				ArrayList<MinimalHsp> addblasthits = new ArrayList<MinimalHsp>();
                for (int i = newblasthits.length; --i >= 0;) {
					if (newblasthits[i].query >= readelements || newblasthits[i].hit >= readelements) {
						addblasthits.add(newblasthits[i]);
					}
				}

				// now I know which of the "new" blast hits to add
				int oldnum = readdata.blasthits.length;
				MinimalHsp[] blasthits = new MinimalHsp[oldnum + addblasthits.size()];
				System.arraycopy(readdata.blasthits, 0, blasthits, 0, oldnum);
				for (int i = addblasthits.size(); --i >= 0;) {
					blasthits[oldnum + i] = (MinimalHsp) addblasthits.get(i);
				}

				newaln = null;
				readdata = null;
				if (nographics == false) { // start the GUI
					System.out.println("...reading data");
					ClusterData myclusterdata = new ClusterData(new MinimalHsp[0], new AminoAcidSequence[0],
							new String[0], new HashMap<String, Integer>(), eval, pval, scval, verbose, cpu,
							save_intermediate_results, cmd, blastpath, addblastvbparam, formatdbpath, referencedb,
							errbuff, input_filename);
					myclusterdata.roundslimit = dorounds;// set the limit of how often to run this
					ClusteringWithGui myclusterer = new ClusteringWithGui(myclusterdata);
					myclusterer.initaddedseqs(blasthits, sequences, sequence_names, sequence_name_internal_mapping,
							newnumarr, allposarr, mymaxmove, mypval, true);
					readdata = null;
					myclusterer.setVisible(true);
				} else {
					System.out.println("DONE. To visualize results restart program with the -nographics F option");
				}

			} else { // load data from olddata
				System.out.println("Reading old data from " + olddata);
				saverunobject readdata;
				
				try {
					readdata = ClusterData.load_run_from_file(new java.io.File(olddata));
				} catch (FileNotFoundException e) {
					System.err.println("file not found: " + olddata);
					return false;
					
				} catch (ParseException e) {
					System.err.println("line " + e.getErrorOffset() + ": " + e.getMessage());
					return false;
					
				} catch (IOException e) {
					System.err.println(e.getMessage());
					return false;
				}
				
				int seqnum = readdata.inaln.length;
				HashMap<String, Integer> sequence_name_internal_mapping = new HashMap<String, Integer>(
						(int) (seqnum / 0.74), (float) 0.75); // mapping of names to array indices
				String[] sequence_names = new String[seqnum];

				for (int i = 0; i < seqnum; i++) {
					sequence_names[i] = readdata.inaln[i].name;
					readdata.inaln[i].name = new String("sequence" + i);// assign internally used names
					sequence_name_internal_mapping.put(readdata.inaln[i].name, new Integer(i));
					readdata.inaln[i].seq = readdata.inaln[i].seq.toUpperCase();
				}

				if (nographics == false) { // start the GUI
					System.out.println("...reading data");
					ClusterData myclusterdata = new ClusterData(new MinimalHsp[0], new AminoAcidSequence[0],
							new String[0], new HashMap<String, Integer>(), eval, pval, scval, verbose, cpu,
							save_intermediate_results, cmd, blastpath, addblastvbparam, formatdbpath, referencedb,
							errbuff, input_filename);

					myclusterdata.roundslimit = dorounds;

					ClusteringWithGui myclusterer = new ClusteringWithGui(myclusterdata);
					myclusterer.initaddedseqs(readdata.blasthits, readdata.inaln, sequence_names,
							sequence_name_internal_mapping, new int[0], readdata.posarr, readdata.maxmove,
							readdata.pval, false);
					readdata = null;
					myclusterer.setVisible(true);
				} else {
					System.out.println("DONE. To visualize results restart program with the -nographics F option");
				}
			}

		} else { // no need to prepare anything, just load an existing CLANS file

			if (dorounds >= 0 && savetoname != null) {
				// run in non-GUI mode if dorounds and savetoname are set
				if (!run_clans_without_gui()) {

				}

			} else { // if the reclustering is NOT the case

				if (nographics == false) { // load the input file and start the GUI

					ClusterData myclusterdata = new ClusterData(new MinimalHsp[0], new AminoAcidSequence[0],
							new String[0], new HashMap<String, Integer>(), eval, pval, scval, verbose, cpu,
							save_intermediate_results, cmd, blastpath, addblastvbparam, formatdbpath, referencedb,
							errbuff, input_filename);

					myclusterdata.roundslimit = dorounds;

					ClusteringWithGui myclusterer = new ClusteringWithGui(myclusterdata);
					myclusterer.setVisible(true);

				} else { // pretty useless case that's only educational to the user
					System.out
							.println("Nothing to do!, try starting the program with the -nographics option set to false");
				}
			}
		}
		return true;
	}

	/**
	 * run CLANS in command line mode. No GUI will be started. Results will be saved to a file.
	 * 
	 * @return true if run was successful
	 */
	private static boolean run_clans_without_gui() {
		System.out.println("non-graphical mode");

		if (input_filename == null) {
			System.err.println("input file is not set");
			return false;
		}

		ClusterData myclusterdata = new ClusterData(new MinimalHsp[0], new AminoAcidSequence[0], new String[0],
				new HashMap<String, Integer>(), eval, pval, scval, verbose, cpu, save_intermediate_results, cmd,
				blastpath, addblastvbparam, formatdbpath, referencedb, errbuff, input_filename);

		myclusterdata.roundslimit = dorounds;

		ClusteringWithoutGui myclusterer = new ClusteringWithoutGui(myclusterdata);

		try {
			myclusterer.data.load_clans_file(myclusterer.data.getAbsoluteInputfileName());
		} catch (FileNotFoundException e) {
			System.err.println("file not found: " + myclusterer.data.getAbsoluteInputfileName());
			return false;
			
		} catch (ParseException e) {
			System.err.println("line " + e.getErrorOffset() + ": " + e.getMessage());
			return false;
			
		} catch (IOException e) {
			System.err.println(e.getMessage());
			return false;
		}
		
		if (pval != -1) {
			myclusterer.data.pvalue_threshold = myclusterer.data.pval;
		}

		if (initialize) {
			myclusterer.initialize();
		} else {
			myclusterer.data.compute_attraction_values();
		}

		myclusterer.startstopthread(); // start the thread
		int waittime = 15000;// 15 seconds
		synchronized (myclusterer) {
			while (myclusterer.mythread.stop == false) {
				System.out.println("done clustering round " + myclusterer.data.roundsdone);
				try {
					myclusterer.wait(waittime);
				} catch (InterruptedException ie) {
					System.err.println("ERROR, interrupted wait in Main\n");
					System.exit(-5);
				}
			}
		}

		File savefile = new File(savetoname);
		
		try{
			myclusterer.data.save_to_file(savefile);
		} catch (IOException e) {
			System.err.println(e.getMessage() + "\n\n" + "YOUR DATA WAS NOT SAVED!!!\nTry saving to another location.");
			return false;
		} catch (IllegalStateException e) {
			System.err.println(e.getMessage() + "\n\n" + "YOUR DATA WAS NOT SAVED!!!\nTry saving to another location.");
			return false;
		}

		System.out.println("done clustering, saving results to file '" + savefile.getAbsolutePath() + "'");
		return true;
	}

	/**
	 * parse command line arguments
	 * 
	 * @param args
	 *            the command line arguments
	 * @return true if successful
	 */
	static boolean parse_arguments(String[] args) {
		int i = 0;
		while (i < args.length) {
			if (args[i].equals("?") || args[i].equals("-?")) {
				print_usage_help();
				System.exit(0);
			}
			if (args[i].equalsIgnoreCase("-conf") || args[i].equalsIgnoreCase("-c")) {
				// this should have been read in checkargs, so here just skip it
				i++;
				// now, i is the index of the config file
				if (i >= args.length) {
					System.err.println("Error reading -conf, missing argument.");
					return false;
				}
				i++;
				continue;
			}
			if ((args[i].equalsIgnoreCase("-infile")) || (args[i].equalsIgnoreCase("-i"))
					|| (args[i].equalsIgnoreCase("-in"))) {
				i++;
				if (i < args.length) {
					infilename = args[i];
				} else {
					System.err.println("Error reading -infile, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-load")) || (args[i].equalsIgnoreCase("-l"))) {
				i++;
				if (i < args.length) {
					input_filename = args[i];
					docalc = false;
				} else {
					System.err.println("Error reading -load, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-initialize"))) {
				i++;
				if (i < args.length) {
					initialize = (args[i].equalsIgnoreCase("TRUE") || args[i].equalsIgnoreCase("T"));
				} else {
					System.err.println("Error reading -initialize, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-saveto"))) {
				i++;
				if (i < args.length) {
					savetoname = args[i];
				} else {
					System.err.println("Error reading -saveto, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-referencedb")) || (args[i].equalsIgnoreCase("-refdb"))) {
				i++;
				if (i < args.length) {
					referencedb = args[i].split("\\s+", 0);
				} else {
					System.err.println("Error reading -referencedb, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if (args[i].equalsIgnoreCase("-cmd")) {
				int quotesfound = 0;
				cmd = "";
				if ((i + 1) < args.length) {
					if (args[i + 1].indexOf("\"") == -1) {// if the next elem is
															// not in quotes
						i++;
						cmd = args[i];
						i++;
						continue;
					} else {
						while (quotesfound < 2) {
							i++;
							if (i >= args.length) {
								System.err.println("Error reading -cmd, missing argument.");
								return false;
							}
							String curr = args[i];
							int qindex;
							if ((qindex = curr.indexOf("\"")) > -1) {
								quotesfound += 1;
								if (quotesfound == 1) {
									cmd += curr.substring(qindex + 1);
								}
								if (quotesfound == 2) {
									cmd += " " + curr.substring(0, qindex);
								}
								// now check if the second quote is somewhere in
								// command
								if ((qindex = cmd.indexOf("\"")) > -1) {
									quotesfound += 1;
									cmd = cmd.substring(0, cmd.indexOf("\""));
									continue;
								}
								continue;
							}
							if (quotesfound == 1) {
								cmd = cmd + " " + curr;
							}
						}
					}
				} else {
					System.err.println("Error reading -cmd, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-blastpath")) || (args[i].equalsIgnoreCase("-blast"))) {
				int quotesfound = 0;
				blastpath = "";
				if (args[i + 1].indexOf("\"") == -1) {// if the next elem is not in quotes
					i++;
					blastpath = args[i];
					i++;
					continue;
				} else {
					while (quotesfound < 2) {
						i++;
						if (i >= args.length) {
							System.err.println("Error reading -blastpath, missing argument.");
							return false;
						}
						String curr = args[i];
						int qindex;
						if ((qindex = curr.indexOf("\"")) > -1) {
							quotesfound++;
							if (quotesfound == 1) {// if this is the first quote
								blastpath = curr.substring(qindex + 1);
							}
							if (quotesfound == 2) {// if I found the second quote
								blastpath += " " + curr.substring(0, qindex);// add the rest
							}
							// now check if the second quote is somewhere in command
							if ((qindex = blastpath.indexOf("\"")) > -1) {
								quotesfound += 1;
								blastpath = blastpath.substring(0, blastpath.indexOf("\""));
								continue;
							}
						} else {// if this element does not contain a quote
							if (quotesfound == 1) {// but I already did find one
													// before
								blastpath += " " + curr;
							}
						}
					}
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-addblastvb"))) {
				i++;
				if (i < args.length) {
					if (args[i].equalsIgnoreCase("FALSE") || args[i].equalsIgnoreCase("F")) {
						addblastvbparam = false;
					} else {
						addblastvbparam = true;
					}
				} else {
					System.err.println("Error reading -addblastvb, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-formatdbpath")) || (args[i].equalsIgnoreCase("-fdb"))) {
				int quotesfound = 0;
				formatdbpath = "";
				if (args[i + 1].indexOf("\"") == -1) {// if the next elem is not in quotes
					i++;
					formatdbpath = args[i];
					i++;
					continue;
				} else {
					while (quotesfound < 2) {
						i++;
						if (i >= args.length) {
							System.err.println("Error reading -formatdbpath, missing argument.");
							return false;
						}
						String curr = args[i];
						int qindex;
						if ((qindex = curr.indexOf("\"")) > -1) {
							quotesfound += 1;
							if (quotesfound == 1) {
								formatdbpath += curr.substring(qindex + 1);
							}
							if (quotesfound == 2) {
								formatdbpath += " " + curr.substring(0, qindex);
							}
							// now check if the second quote is somewhere in formatdbpath
							if ((qindex = formatdbpath.indexOf("\"")) > -1) {
								quotesfound += 1;
								formatdbpath = formatdbpath.substring(0, formatdbpath.indexOf("\""));
								continue;
							}
							continue;
						}
						if (quotesfound == 1) {
							formatdbpath = formatdbpath + " " + curr;
						}
					}
					i++;
					continue;
				}
			}

			if ((args[i].equalsIgnoreCase("-blastblocks")) || (args[i].equalsIgnoreCase("-bb"))) {
				i++;
				if ((i) < args.length) {
					try {
						blastblocks = Integer.parseInt(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("unable to parse int from '" + args[i] + "' in -blastblocks");
						return false;
					}
				} else {
					System.err.println("Error reading -blastblocks, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-eval")) || (args[i].equalsIgnoreCase("-e"))) {
				i++;
				if ((i) < args.length) {
					try {
						eval = Double.parseDouble(args[i]);
						if (eval < pval) {
							pval = eval;
						}
					} catch (NumberFormatException e) {
						System.err.println("unable to parse double from '" + args[i] + "' in -eval");
						return false;
					}
				} else {
					System.err.println("Error reading -eval, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-pval")) || (args[i].equalsIgnoreCase("-p"))) {
				i++;
				if ((i) < args.length) {
					try {
						pval = Double.parseDouble(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("unable to parse double from '" + args[i] + "' in -pval");
						return false;
					}
				} else {
					System.err.println("Error reading -pval, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-scval")) || (args[i].equalsIgnoreCase("-sc"))) {
				i++;
				if ((i) < args.length) {
					try {
						scval = Float.parseFloat(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("unable to parse float from '" + args[i] + "' in -pval");
						return false;
					}
				} else {
					System.err.println("Error reading -pval, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-verbose")) || (args[i].equalsIgnoreCase("-v"))) {
				i++;
				if ((i) < args.length) {
					try {
						verbose = Integer.parseInt(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("unable to parse int from " + args[i] + " in -verbose.");
						return false;
					}
				} else {
					System.err.println("Error reading -verbose, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-cpu"))) {
				i++;
				if ((i) < args.length) {
					try {
						cpu = Integer.parseInt(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("unable to parse int from " + args[i] + " in -cpu.");
						return false;
					}
				} else {
					System.err.println("Error reading -cpu, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-dorounds"))) {
				i++;
				if ((i) < args.length) {
					try {
						dorounds = Integer.parseInt(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("unable to parse int from " + args[i] + " in -dorounds.");
						return false;
					}
				} else {
					System.err.println("Error reading -dorounds, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-readblast"))) {
				i++;
				if (i < args.length) {
					if (args[i].equalsIgnoreCase("FALSE") || args[i].equalsIgnoreCase("F")) {
						readblast = false;
					} else {
						readblast = true;// default
					}
				} else {
					System.err.println("Error reading -readblast, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-lowmem"))) {
				i++;
				if (i < args.length) {
					if (args[i].equalsIgnoreCase("TRUE") || args[i].equalsIgnoreCase("T")) {
						lowmem = true;
					} else {
						lowmem = false;// default
					}
				} else {
					System.err.println("Error reading -lowmem, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-savepos"))) {
				i++;
				if (i < args.length) {
					if (args[i].equalsIgnoreCase("TRUE") || args[i].equalsIgnoreCase("T")) {
						save_intermediate_results = true;
					} else {
						save_intermediate_results = false;// default
					}
				} else {
					System.err.println("Error reading -savepos, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-docalc"))) {
				i++;
				if (i < args.length) {
					if (args[i].equalsIgnoreCase("TRUE") || args[i].equalsIgnoreCase("T")) {
						docalc = true;
					} else {
						docalc = false;
					}
				} else {
					System.err.println("Error reading -docalc, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-nographics"))) {
				i++;
				if (i < args.length) {
					if (args[i].equalsIgnoreCase("TRUE") || args[i].equalsIgnoreCase("T")) {
						nographics = true;
					} else {
						nographics = false;
					}
				} else {
					System.err.println("Error reading -nographics, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-exhaustive"))) {
				i++;
				if ((i) < args.length) {
					try {
						exhaustive = Integer.parseInt(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("unable to parse int from " + args[i] + " in -exhaustive.");
						return false;
					}
					if (exhaustive > 2) {
						exhaustive = 2;
					} else if (exhaustive < 0) {
						exhaustive = 0;
					}
				} else {
					System.err.println("Error reading -exhaustive, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-olddata"))) {
				i++;
				if ((i) < args.length) {
					olddata = args[i];
				} else {
					System.err.println("Error reading -olddata, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-newseqs"))) {
				i++;
				if ((i) < args.length) {
					newseqs = args[i];
				} else {
					System.err.println("Error reading -newseqs, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if ((args[i].equalsIgnoreCase("-enrichseqs"))) {
				i++;
				if (i < args.length) {
					if (args[i].equalsIgnoreCase("true") || args[i].equalsIgnoreCase("t")) {
						enrichseqs = true;
					} else {
						enrichseqs = false;
					}
				} else {
					System.err.println("Error reading -enrichseqs, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if (args[i].equalsIgnoreCase("-gatherseqseval")) {
				i++;
				if (i < args.length) {
					try {
						gatherseqseval = Double.parseDouble(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("ERROR: unable to parse double from " + args[i] + " in -gatherseqseval");
						return false;
					}
				} else {
					System.err.println("Error reading -gatherseqseval, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if (args[i].equalsIgnoreCase("-rmseqseval")) {
				i++;
				if (i < args.length) {
					try {
						rmseqseval = Double.parseDouble(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("ERROR: unable to parse double from " + args[i] + " in -rmseqseval");
						return false;
					}
				} else {
					System.err.println("Error reading -rmseqseval, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			if (args[i].equalsIgnoreCase("-maxenrichseqsnum")) {
				i++;
				if (i < args.length) {
					try {
						gatherseqseval = Integer.parseInt(args[i]);
					} catch (NumberFormatException e) {
						System.err.println("ERROR: unable to parse int from " + args[i] + " in -maxenrichseqsnum");
						return false;
					}
				} else {
					System.err.println("Error reading -maxenrichseqsnum, missing argument.");
					return false;
				}
				i++;
				continue;
			}

			System.err.println("unknown option " + args[i]);
			return false;
		}
		return true;
	}

	/**
	 * parse a configuration file
	 * 
	 * @param filename
	 *            the configuration filename
	 * @return true if successful
	 */
	static boolean parse_configuration_file(String filename) {
		try {
			BufferedReader infile = new BufferedReader(new FileReader(filename));
			String inline;
			int enddata = 0;

			while ((inline = infile.readLine()) != null) {
				inline = inline.trim();

				if ((enddata = inline.indexOf("#")) > -1) {// if this is a line with a comment on it

					if (enddata > 1) {// if I have some data on this line
						inline = inline.substring(0, enddata);
						if ((parse_arguments(inline.split("\\s", 0))) == false) {
							System.err.println("Error reading on line " + inline);
							infile.close();
							return false;
						}
					} else {
						continue;
					}

				} else {// if this line has no comment on it
					if (inline.length() > 0) {
						if ((parse_arguments(inline.split("\\s", 0))) == false) {
							System.err.println("Error reading on line " + inline);
							infile.close();
							return false;
						}
					}
				}
			}
			infile.close();

		} catch (IOException e) {
			System.err.println("IOError reading from " + filename);
			return false;
		}
		return true;
	}

	/**
	 * check whether there is a configuration file in the argument list
	 * 
	 * @param args
	 *            the command line parameters
	 * @return true if there is a config file parameter followed by a configuration file name
	 */
	static boolean get_configuration_file_argument(String[] args) {
		for (int i = 0; i < args.length; i++) {

			if ((args[i].equalsIgnoreCase("-conf")) || (args[i].equalsIgnoreCase("-c"))) {
				if ((i + 1) < args.length) {
					conffilename = args[i + 1];
				} else {
					return false;
				}
			}

		}
		return true;
	}
}
