package clans.model;

// TODO: refactor public/protected/private

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.Vector;
import java.util.concurrent.CancellationException;

import javax.swing.SwingWorker;

import clans.algorithms.fruchtermanreingold.ClusterMethods;
import clans.algorithms.fruchtermanreingold.MovementComputerThread;
import clans.gui.Shapes;
import clans.io.ClusterDataLoadHelper;
import clans.io.ComfortableBufferedWriter;
import clans.misc.MyMath;
import clans.misc.OsUtils;
import clans.model.microarray.Replicates;
import clans.model.proteins.AminoAcidSequence;
import clans.model.proteins.MinimalAttractionValue;
import clans.model.proteins.MinimalHsp;

/**
 * Handles the data used in the clustering.
 */
public class ClusterData {

    // variables initialized on creation
    public MinimalHsp[] blasthits = null;
    public AminoAcidSequence[] sequences = null;
    public String[] sequence_names = null;
    public HashMap<String, Integer> nameshash = null;
    public double eval = -1;
    public double pval = -1;
    public float scval = -1;
    public int verbose = -1;
    public int cpu = -1;
    public boolean save_intermediate_results = false;
    public String cmd = null;
    public String blastpath = null;
    public boolean addblastvbparam = false;
    public String formatdbpath = null;
    public String[] referencedb = null;
    public StringBuffer errbuff = null;

    private String input_filename = null;

    // variables initialized later on
    public boolean nographics = false;
    public boolean complexatt = true;
    public int seqnum = 0;
    public float[] seqlengths = null;
    public MovementComputerThread[] movethreads = null;
    public boolean usescval = false;
    public MinimalAttractionValue[] attractionValues = null;

    // variables I use as part of the clustering
    public int rounds = 0;
    public float[][] positions = null;
    public float[][] posarr = null;
    public boolean cluster2d = false;
    public float maxmove = 0.1f;
    public double pvalue_threshold = 1;
    public double mineval = 1;
    public double cooling = 1;
    public double currcool = 1;
    public ClusterDataLoadHelper saveddata = null;
    public float attfactor = 10f;
    public float repfactor = 0.5f;
    public int attvalpow = 1;
    public int repvalpow = 1;
    public float dampening = 0.2f;
    public double minattract = 1;
    public float[] weights = null;
    public ArrayList<File> mapfiles = null;
    public ArrayList<File> lookupfiles = null;
    public java.util.Vector<Replicates> affyfiles = null;
    public boolean usefoldchange = false;
    public boolean avgfoldchange = false;
    public String namesdmp_file = "not_spcified";
    public String nodesdmp_file = "not_specified";
    public boolean showinfo = false;
    public int[] selectedSequencesIndices = new int[0];
    public int[] selectedSequencesIndicesStableCopy = new int[0];
    public float[][] movementsLastIteration = null;
    public float[][] movements = null;
    public float[][] posarrtmp = null;
    public int[][] drawarrtmp = null;
    public ArrayList<ArrayList<int[]>> draworder = null;
    public static int dimensions = 3; // TODO: use this field in the GUI too, currently GUI uses the button state
    public int elements = -1;
    /**
     * A 3x3 rotation matrix.
     */
    public double[][] rotmtx = MyMath.create3x3IdentityMatrix();
    /**
     * A 3x3 rotation matrix.
     */
    public double[][] myrotmtx = MyMath.create3x3IdentityMatrix();
    public MinimalAttractionValue[] orgattvals = null;
    public boolean attvalsimple = false;
    public boolean rescalepvalues = false;
    public double maxvalfound = 0;
    public float p2attfactor = 1;
    public float p2attoffset = 0;
    public int ovalsize = 4;
    public int dotsize = 2;
    public int groupsize = 4;
    public java.util.Vector<SequenceGroup> seqgroupsvec = new java.util.Vector<SequenceGroup>();
    public ArrayList<int[][]> polygons = null;
    public boolean showSequenceGroups = false;
    public boolean changedvals = false;
    public java.awt.Color[] colorarr = null;
    public float[] colorcutoffs = null;
    /**
	 * The number of rounds completed since the last press of the start/stop/resume button. This is used for running a
	 * fixed number of rounds as chosen by the user in the options window.
	 */
    public int roundsCompleted = 0;
    /**
	 * The number of rounds to run before stopping automatically. If set to -1, automatic stopping is disabled.
	 */
    private int roundsLimitBeforeStopping = -1;

    private int autosaveIntervalMinutes; // the autosave interval in minutes; -1 if value not in input file, 0 = disabled

    public ClusterData(MinimalHsp[] blasthits, AminoAcidSequence[] sequences, String[] namearr,
            HashMap<String, Integer> nameshash, double eval, double pval, float scval, int verbose, int cpu,
            boolean save_intermediate_results, String cmd, String blastpath, boolean addblastvbparam,
            String formatdbpath, String[] referencedb, StringBuffer errbuff, String input_filename) {

        this.sequences = ClusterMethods.removeGapsFromSequences(sequences);

        this.movethreads = new MovementComputerThread[cpu];

        this.blasthits = blasthits;
        this.sequence_names = namearr;
        this.nameshash = nameshash;
        this.eval = eval;

        if (pval != -1) {
            this.pval = pval;
        }

        this.scval = scval;
        this.verbose = verbose;
        this.cpu = cpu;
        this.save_intermediate_results = save_intermediate_results;
        this.cmd = cmd;
        this.blastpath = blastpath;
        this.addblastvbparam = addblastvbparam;
        this.formatdbpath = formatdbpath;
        this.referencedb = referencedb;
        this.errbuff = errbuff;

        this.seqnum = namearr.length;

        this.setInputFilename(input_filename);
    }

    /**
     * Sets the filename as an absolute one independent of whether the input filename is relative or absolute.
     * 
     * @param new_filename
     *            a relative or absolute input filename or null
     */
    public void setInputFilename(String new_filename) {

        input_filename = new_filename;

        if (new_filename == null) {
            return;
        }

        if (!new File(new_filename).isAbsolute()) {
            // the directory java was executed from is used as basis for all relative filenames
            input_filename = System.getProperty("user.dir") + System.getProperty("file.separator") + new_filename;
        }
    }

    /**
     * 
     * @return the absolute filename of the file associated with this clustering or null if the filename is not set
     */
    public String getAbsoluteInputfileName() {
        return input_filename;
    }

    /**
     * 
     * @return the base name of the file associated with this clustering or null if the filename is not set
     */
    public String getBaseInputfileName() {
        
        if (input_filename == null) {
            return null;
        }
        
        return new File(input_filename).getName();
    }
    
    /**
     * 
     * @return the name of the file used for writing intermediate results after each round if that flag is set
     */
    public String getIntermediateResultfileName() {
        return getAbsoluteInputfileName() + ".savepos";
    }

    /**
     * Resets the draw order so that it will be recomputed the next time it is used. 
     */
    public void resetDrawOrder() {
        draworder = new ArrayList<ArrayList<int[]>>();
    }
    
	public int getTotalNumberOfConnections() {
		// BLAST HSP data
		if (blasthits != null) {
			return blasthits.length;
		}

		// attvals mode data
		if (attractionValues == null) {
			return attractionValues.length;

		}

		// no data
		return 0;
	}
	
	public boolean knowsAutosave() {
		return autosaveIntervalMinutes != -1;
	}
	
	public int getAutosaveIntervalMinutes() {
		return autosaveIntervalMinutes;
	}
	
	public void setAutosaveIntervalMinutes(int new_value) {
		autosaveIntervalMinutes = new_value;
	}

	public void makeAutosaveAware() {
		if (!knowsAutosave()) {
			guessInitialAutosaveInterval();
		}
	}
	
	/**
	 * Guesses a default autosave interval for these data based on the number of connections. The rationale is that
	 * connections contribute most to file size and therefore also to the time it takes to save. These values should be
	 * relatively long intervals to not annoy users with frequent wait times while saving.
	 */
	private void guessInitialAutosaveInterval() {
		int total_connections = getTotalNumberOfConnections();

		if (total_connections == 0) {
			autosaveIntervalMinutes = 0;
			return;
		}

		if (total_connections < 1000000) { // < 1 million connections -> small file
			autosaveIntervalMinutes = 10;

		} else if (total_connections < 5000000) { // < 5 million connections -> medium file
			autosaveIntervalMinutes = 20;

		} else { // large file
			autosaveIntervalMinutes = 30;
		}
	}

	/**
	 * Convenience method that calls {@code load_run_from_file(File, SwingWorker)} with SwingWorker {@code null}. This
	 * method should be used in no-GUI mode.
	 */
	public static ClusterDataLoadHelper load_run_from_file(File infile) throws ParseException, IOException, FileNotFoundException {
		return load_run_from_file(infile, null);
	}
	
	/**
	 * load stuff from a savefile !!! NOTE!!!: this was edited so as to combine hsp's with the same query-hit
	 * combination irrespective of which sequence is the query and which the hit this is a valid approach in this case,
	 * as I later on anyways symmetrize the sequence interactions
	 * 
	 * @param infile
	 *            The file to load data from
	 * @param worker
	 *            If not null, this worker is used on a regular basis to check whether saving should be canceled. This
	 *            is used by the GUI to act on user aborts.
	 * @return A {@code saverunobject} instance containing the loaded data
	 * @throws ParseException
	 *             If the file is corrupted or has the wrong format
	 * @throws IOException
	 *             If an I/O error occurs
	 * @throws FileNotFoundException
	 *             If the file does not exist, cannot be opened, or is a directory
	 * @throws CancellationException
	 *             If worker is cancelled and loading should therefore be cancelled. In headless (no-GUI) mode this
	 *             Exception will not occur.
	 */
	public static ClusterDataLoadHelper load_run_from_file(File infile, SwingWorker<Void, Integer> worker) throws ParseException,
			IOException, FileNotFoundException, CancellationException {
         
		System.out.println("LOADING data from '" + infile.getAbsolutePath() + "'");

		ClusterDataLoadHelper myrun = new ClusterDataLoadHelper();
        myrun.file = infile;

        boolean is_biolayout = false;
		boolean errors_occurred = false;
		String error_message = "";
		
		boolean hasBlockPos = false;

		String inline;
		int expected_sequences = -1;
		int current_line = 0;
		int hspCount = 0;

		BufferedReader buffered_file_handle = new BufferedReader(new FileReader(infile));

		try {
			System.out.println("Count HSPs in '" + infile.getAbsolutePath() + "'");
			while ((inline = buffered_file_handle.readLine()) != null) {
				if (inline.equalsIgnoreCase("<hsp>")) {
					while (((inline = buffered_file_handle.readLine()) != null) && (inline.equalsIgnoreCase("</hsp>") == false)) {
						//skip empty lines
						if(inline.length()==0){
							continue;
						}
						hspCount++;
					}//end while count HSPs
				}
			}//end while for going through file for counting HSPs
			System.out.println("Counted " + hspCount + " HSPs");
		} catch (IOException e) {
			error_message = "IOError unable to read from " + infile.getAbsolutePath() + ":\n\t" + e.getMessage();
			System.err.println(error_message);
			throw new IOException(error_message);
		} finally {
			buffered_file_handle.close();
		}

		buffered_file_handle = new BufferedReader(new FileReader(infile));

		try {
            while ((inline = buffered_file_handle.readLine()) != null) {
            	current_line += 1;
                inline = inline.trim();
               
                if (inline.length() == 0) {
                    continue;
                
                } else if (inline.startsWith("#")) { // skip comment lines
                    continue;
                }
                
                if (inline.startsWith("sequences=")) {
                    try {
                        expected_sequences = Integer.parseInt(inline.substring(10));
                        System.out.println("sequences=" + expected_sequences);
                    } catch (NumberFormatException ne) {
                        System.err.println("Error parsing int from '" + inline.substring(10) + "'");
                    }
                    myrun.inaln = new AminoAcidSequence[expected_sequences];
                    myrun.posarr = new float[expected_sequences][3];
                    myrun.blasthits = new MinimalHsp[0];
                    continue;
                
                } else if (expected_sequences != -1) {
                	
                	checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)

                    if (inline.equalsIgnoreCase("<param>")) {
                        if (!myrun.parse_params_block(buffered_file_handle)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <param> block";
                        }

                    } else if (inline.equalsIgnoreCase("<function>")) {
                        if (!myrun.parse_function_block(buffered_file_handle)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <function> block";
                        }

                    } else if (inline.equalsIgnoreCase("<affyfiles>")) {
                        if (!myrun.parse_affyfiles_block(buffered_file_handle)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <affyfiles> block";
                        }

                    } else if (inline.equalsIgnoreCase("<rotmtx>")) {
                        if (!myrun.parse_rotmtx_block(buffered_file_handle)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <rotmtx> block";
                        }

                    } else if (inline.equalsIgnoreCase("<seq>")) {
                        if (!myrun.parse_seq_block(buffered_file_handle, expected_sequences)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <seq> block";
                        }

                    } else if (inline.equalsIgnoreCase("<weight>")) {
                        if (!myrun.parse_weight_block(buffered_file_handle, expected_sequences)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <weight> block";
                        }

                    } else if (inline.equalsIgnoreCase("<pos>")) {
                    	hasBlockPos = true;
                        if (!myrun.parse_pos_block(buffered_file_handle, expected_sequences)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <pos> block";
                        }

                    } else if (inline.equalsIgnoreCase("<hsp>")) {
                        if (!myrun.parse_hsp_block(buffered_file_handle, worker, hspCount)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <hsp> block";
                        }

                    } else if (inline.equalsIgnoreCase("<att>")) {
                        if (!myrun.parse_att_block(buffered_file_handle)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <att> block";
                        }

                    } else if (inline.equalsIgnoreCase("<seqgroups>")) {
                        if (!myrun.parse_seqgroups_block(buffered_file_handle)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <seqgroups> block";
                        }

                    } else if (inline.equalsIgnoreCase("<mtx>")) {
                        if (!myrun.parse_mtx_block(buffered_file_handle, expected_sequences)) {
                        	errors_occurred = true;
                        	error_message = "could not parse <mtx> block";
                        }

                    } else {
                    	errors_occurred = true;
                    	error_message = "unknown format in line \"" + inline + "\"";
                    }

                    if (errors_occurred) {
                    	System.err.println(error_message);
                    	throw new ParseException(error_message, current_line);
                    }

                } else {
                	is_biolayout = true; // let's try and parse it as biolayout file
                	break;
                }
            }
            
            checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)

        } catch (IOException e) {
        	error_message = "IOError unable to read from " + infile.getAbsolutePath() + ":\n\t" + e.getMessage();
        	System.err.println(error_message);
        	throw new IOException(error_message);

        } finally {
        	buffered_file_handle.close();
        }
		
		if (is_biolayout) {
            System.out.println("assuming BioLayout format");
            return load_biolayout_file(infile, worker);
		}

        if (!hasBlockPos) {
            // give it random values
            Random rand = new Random(System.currentTimeMillis());
            myrun.posarr = new float[myrun.inaln.length][3];
            for (int i = 0; i < myrun.inaln.length; i++) {
                myrun.posarr[i][0] = rand.nextFloat();
                myrun.posarr[i][1] = rand.nextFloat();
                myrun.posarr[i][2] = rand.nextFloat();
            }
        }

        return myrun;
    }

	/**
	 * Convenience method that calls {@code load_clans_file(String, SwingWorker)} with SwingWorker {@code null}. This
	 * method should be used in no-GUI mode.
	 */
	public void load_clans_file(String input_filename) throws ParseException,
	IOException, FileNotFoundException {
		load_clans_file(input_filename, null);
	}
	
	/**
	 * Read data from a CLANS format file. 
	 * 
	 * @param input_filename
	 *            The input filename.
	 * @param worker
	 *            If not null, this worker is used on a regular basis to check whether saving should be canceled. This
	 *            is used by the GUI to act on user aborts.
	 * @throws ParseException
	 *             If {@code ClusterData.load_run_from_file} throws it
	 * @throws IOException
	 *             If {@code ClusterData.load_run_from_file} throws it
	 * @throws FileNotFoundException
	 *             If {@code ClusterData.load_run_from_file} throws it
	 * @throws CancellationException
	 *             If {@code ClusterData.load_run_from_file} throws it
	 */
	public void load_clans_file(String input_filename, SwingWorker<Void, Integer> worker) throws ParseException,
			IOException, FileNotFoundException, CancellationException {
    	
        ClusterDataLoadHelper loaded_data = ClusterData.load_run_from_file(new java.io.File(input_filename), worker);
        
        checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
        
        injectLoadedDataIntoExisting(input_filename, loaded_data);
        
        System.out.println("File loaded:" + input_filename);
    }
    
    public void injectLoadedDataIntoExisting(String input_filename, ClusterDataLoadHelper loaded_data) {
        
    	this.input_filename = input_filename;
        
        this.autosaveIntervalMinutes = loaded_data.autosaveIntervalMinutes;

        sequences = ClusterMethods.removeGapsFromSequences(loaded_data.inaln);

        positions = loaded_data.posarr;
        blasthits = loaded_data.blasthits;
        usescval = loaded_data.usescval;

        if (blasthits == null) {
            // first time I load myattvals; cannot be anything else; don't need to sync
            attractionValues = loaded_data.attvals;
        }

        complexatt = loaded_data.complexatt;
        maxmove = loaded_data.maxmove;
        pvalue_threshold = loaded_data.pval;
        cooling = loaded_data.cooling;
        currcool = loaded_data.currcool;
        attfactor = loaded_data.attfactor;
        repfactor = loaded_data.repfactor;
        attvalpow = loaded_data.attvalpow;
        repvalpow = loaded_data.repvalpow;
        dampening = loaded_data.dampening;
        minattract = loaded_data.minattract;
        weights = loaded_data.weights;
        mapfiles = loaded_data.mapfiles;
        lookupfiles = loaded_data.lookupfiles;
        affyfiles = loaded_data.affyfiles;
        usefoldchange = loaded_data.usefoldchange;
        avgfoldchange = loaded_data.avgfoldchange;
        namesdmp_file = loaded_data.namesdmp_file;
        nodesdmp_file = loaded_data.nodesdmp_file;

        // be careful not to overwrite any blastpath and formatdbpath setting passed via command line
        if (blastpath.equals("blastall -p blastp")) {// if it was not changed via command line
            blastpath = loaded_data.blastpath;
        }
        if (formatdbpath.equals("formatdb -p T")) {// if it was not changed via command line
            formatdbpath = loaded_data.formatdbpath;
        }

        cluster2d = loaded_data.cluster2d;
        showinfo = loaded_data.showinfo;
        int number_of_sequences = sequences.length;
        nameshash = new HashMap<String, Integer>((int) (number_of_sequences / 0.75) + 1, (float) 0.75);// holds info
                                                                                                       // about which
                                                                                                       // name is which
                                                                                                       // array number
        sequence_names = new String[number_of_sequences];
        for (int i = 0; i < number_of_sequences; i++) {
            sequence_names[i] = sequences[i].name.trim();
            sequences[i].name = "sequence" + i;
            nameshash.put(sequences[i].name, new Integer(i));
        }
        elements = sequence_names.length;
        selectedSequencesIndices = new int[0];
        posarr = positions;
        movementsLastIteration = new float[elements][ClusterData.dimensions];
        movements = new float[elements][ClusterData.dimensions];
        posarrtmp = new float[elements][ClusterData.dimensions];
        drawarrtmp = new int[elements][ClusterData.dimensions];
        resetDrawOrder();

        myrotmtx = loaded_data.rotmtx;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                rotmtx[i][j] = myrotmtx[i][j];
            }
        }

        orgattvals = null;
        // first time I load myattvals; don't need to sync as nothing else can be using this yet

        compute_attraction_values();

        dotsize = loaded_data.dotsize;
        ovalsize = loaded_data.ovalsize;
        groupsize = loaded_data.groupsize;
        polygons = clans.gui.Shapes.get(groupsize);
        seqgroupsvec = loaded_data.seqgroupsvec;
        if (seqgroupsvec.size() > 0) {
            showSequenceGroups = true;
        }
        changedvals = true;

        if (loaded_data.colorarr != null) {
            colorarr = loaded_data.colorarr;
        }
        if (loaded_data.colorcutoffs != null) {
            colorcutoffs = loaded_data.colorcutoffs;
        }

        rounds = loaded_data.rounds;

        seqlengths = new float[number_of_sequences];
        float maxlength = 0;
        for (int i = 0; i < number_of_sequences; i++) {
            seqlengths[i] = sequences[i].seq.length();
            if (seqlengths[i] > maxlength) {
                maxlength = seqlengths[i];
            }
        }
        for (int i = 0; i < number_of_sequences; i++) {
            seqlengths[i] /= maxlength;
        }
    }

	/**
	 * 
	 * Read data from a biolayout format file.
	 * 
	 * @param input_file
	 *            The input file.
	 * @return The loaded data if loading successful.
	 * @throws ParseException
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	static ClusterDataLoadHelper load_biolayout_file(File input_file, SwingWorker<Void, Integer> worker) throws ParseException,
			IOException, FileNotFoundException, CancellationException {
    	ClusterDataLoadHelper myrun = new ClusterDataLoadHelper();
        float attval;
        String name1;
        String name2;
        String error_message;
        ArrayList<String> namelist = new ArrayList<String>();
        ArrayList<MinimalAttractionValue> datlist = new ArrayList<MinimalAttractionValue>();
        
        String inline;
        String[] tmparr;
        int vals = 0;
        int seq1num, seq2num, seqnum = 0;
        int line_number = 0;
        HashMap<String, Integer> nameshash = new HashMap<String, Integer>();
       
        BufferedReader inread = new BufferedReader(new FileReader(input_file));

        try {
            while ((inline = inread.readLine()) != null) {
            	line_number += 1;
            	
                if (inline.length() == 0) {
                    continue;
              
                } else if (inline.startsWith("//")) {
                    continue;
                
                } else {// if this is contact info
                	checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
                	
                    // split the line on spaces either in two or three
                    tmparr = inline.split("\\s+");
                    
                    if (tmparr.length == 2) {
                        if ((vals != 2) && (vals != 0)) {
                            System.out.println("ERROR found 2 elements (expected 3) on line '" + inline + "'");
                            System.err.println("Mixed 2 and 3 element entries in file " + input_file.getName()
                                    + " the graph layout WILL misrepresent distances");
                            System.err.println("NEVER mix 2 and 3 element entries!");
                        }
                        
                        vals = 2;
                        // then assign the connection a strength of 1
                        name1 = tmparr[0];
                        name2 = tmparr[1];
                        attval = 1;
                        // now see if the names are already known, else, add as new names
                        if (nameshash.containsKey(name1) == false) {
                            namelist.add(name1);
                            seq1num = seqnum;
                            nameshash.put(name1, new Integer(seq1num));
                            seqnum++;
                        } else {
                            seq1num = ((Integer) nameshash.get(name1)).intValue();
                        }

                        if (nameshash.containsKey(name2) == false) {
                            // System.out.println("new name "+name2);
                            namelist.add(name2);
                            seq2num = seqnum;
                            nameshash.put(name2, new Integer(seq2num));
                            seqnum++;
                        } else {
                            seq2num = ((Integer) nameshash.get(name2)).intValue();
                        }
                        datlist.add(new MinimalAttractionValue(seq1num, seq2num, attval));
                        
                    } else if (tmparr.length == 3) {
                        if ((vals != 3) && (vals != 0)) {
                            System.out.println("ERROR found 3 elements (expected 2) on line '" + inline + "'");
                            System.err.println("Mixed 2 and 3 element entries in file " + input_file.getName()
                                    + " the graph layout WILL misrepresent distances");
                            System.err.println("NEVER mix 2 and 3 element entries!");
                        }
                        
                        vals = 3;
                        name1 = tmparr[0];
                        name2 = tmparr[1];
                        
                        try {
                            attval = Float.parseFloat(tmparr[2]);
                        
                        } catch (NumberFormatException ne) {
                        	inread.close();
                        	
                        	error_message = "unable to parse float from '" + tmparr[2] + "' in '" + inline + "'";
                        	System.err.println(error_message);
                            throw new ParseException(error_message, line_number);
                        }
                        
                        // now see if the names are already known, else, add as new names
                        if (nameshash.containsKey(name1) == false) {
                            namelist.add(name1);
                            seq1num = seqnum;
                            nameshash.put(name1, new Integer(seq1num));
                            seqnum++;
                        } else {
                            seq1num = ((Integer) nameshash.get(name1)).intValue();
                        }

                        if (nameshash.containsKey(name2) == false) {
                            namelist.add(name2);
                            seq2num = seqnum;
                            nameshash.put(name2, new Integer(seq2num));
                            seqnum++;
                        } else {
                            seq2num = ((Integer) nameshash.get(name2)).intValue();
                        }
                        datlist.add(new MinimalAttractionValue(seq1num, seq2num, attval));

                    } else {
                    	inread.close();
                    	
                    	error_message = "unkonwn format: \"" + inline + "\"";
                    	System.err.println(error_message);
                        throw new ParseException(error_message, line_number);
                    }
                }
            }
        
        } catch (IOException e) {
        	error_message = "unable to read file \"" + input_file.getName() + "\"";
        	System.err.println(error_message);
            throw new IOException(error_message);

        } finally {
        	try{
        		inread.close();
        	} catch (IOException e) {
        		
        	}
        }
        
        // now make an array out of the arrayList
        myrun.attvals = (MinimalAttractionValue[]) datlist.toArray(new MinimalAttractionValue[0]);
        java.util.Random rand = new java.util.Random(System.currentTimeMillis());
        myrun.inaln = new AminoAcidSequence[namelist.size()];
        myrun.posarr = new float[namelist.size()][3];
        for (int i = namelist.size(); --i >= 0;) {
            myrun.inaln[i] = new AminoAcidSequence(((String) namelist.get(i)), "");
            myrun.posarr[i][0] = rand.nextFloat();
            myrun.posarr[i][1] = rand.nextFloat();
            myrun.posarr[i][2] = rand.nextFloat();
        }

        // now normalize the attraction values and symmetrize the array of attvals
        myrun.file = input_file; // marker for successful read
        return myrun;
    }

    /**
     * Read a file in matrix format
     * load data from a file with fasta format sequence input (sequence length can be null) and a matrix with pairwise
     * similarity values. Format:
     *     sequences=number_of_sequences
     *     <seqs> // FASTA format sequence records [headers are mandatory, sequences are optional]
     *     >header 1
     *     SEQUENCE1
     *     </seqs>
     *     <mtx>
     *     matrix_with_attracion-values
     *     </mtx>
     * @param input_filename
     * @return
     */
    public static ClusterDataLoadHelper load_matrix_file(File input_filename) {
        System.out.println("reading matrixdata");
      
        ClusterDataLoadHelper myrun = new ClusterDataLoadHelper();
        myrun.file = null;// if myrun has a filename all was read ok
        try {
            BufferedReader inread = new BufferedReader(new FileReader(input_filename));
            String inline;
            int seqs = -1;
            while ((inline = inread.readLine()) != null) {
                inline = inline.trim();
                if (inline.length() == 0) {
                    continue;
                } else if (inline.startsWith("#")) {// skip comment lines
                    continue;
                }
                if (inline.startsWith("sequences=")) {
                    try {
                        seqs = Integer.parseInt(inline.substring(10));
                    } catch (NumberFormatException ne) {
                        System.err.println("Error parsing int from '" + inline.substring(10) + "'");
                    }
                    myrun.inaln = new AminoAcidSequence[seqs];
                    myrun.posarr = new float[seqs][3];
                    myrun.blasthits = null;
                    continue;
                } else if (seqs != -1) {
                    inline = inline.trim();
                    if (inline.equalsIgnoreCase("<seqs>")) {
                        System.out.println("reading sequences");
                        int counter = 0;
                        String currname = "";
                        String currseq = "";
                        while (((inline = inread.readLine().trim()) != null)
                                && (inline.equalsIgnoreCase("</seqs>") == false)) {

                            if (inline.length() == 0) {
                                continue;
                            }
                            
                            // read the sequence data
                            if (inline.startsWith(">")) {
                                if (currname.length() > 0) {
                                    myrun.inaln[counter] = new AminoAcidSequence(currname, currseq);
                                    counter++;
                                }
                                currname = inline.substring(1);
                                currseq = "";
                            } else {
                                currseq += inline;
                            }
                        }
                        myrun.inaln[counter] = new AminoAcidSequence(currname, currseq);
                        counter++;
                        if (counter != seqs) {
                            System.err.println("ERROR, not found the number of specified sequences");
                            inread.close();
                            return myrun;
                        }
                        System.out.println("done reading sequences:" + seqs);

                    } else if (inline.equalsIgnoreCase("<mtx>")) {
                        System.out.println("reading matrix");
                        int counter = 0;
                        String[] tmparr;
                        HashMap<MinimalAttractionValue, MinimalAttractionValue> tmphash = new HashMap<MinimalAttractionValue, MinimalAttractionValue>();
                        while (((inline = inread.readLine()) != null) && (inline.equalsIgnoreCase("</mtx>") == false)) {
                            // skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            tmparr = inline.trim().split("\\s+");
                            if (tmparr.length != seqs) {
                                System.err.println("ERROR reading positions from " + inline + "; expecting " + seqs
                                        + " values");
                                inread.close();
                                return myrun;
                            }
                            try {
                                for (int i = 0; i < seqs; i++) {
                                    float tmpval = Float.parseFloat(tmparr[i]);
                                    if (tmpval != 0) {
                                        MinimalAttractionValue curratt = new MinimalAttractionValue(counter, i);
                                        if (tmphash.containsKey(curratt)) {
                                            tmphash.get(curratt).att += tmpval / 2;
                                        } else {
                                            curratt.att = tmpval / 2;
                                            tmphash.put(curratt, curratt);
                                        }
                                    }
                                }
                            } catch (NumberFormatException ne) {
                                System.err.println("ERROR, unable to parse float array from " + inline);
                                inread.close();
                                return myrun;
                            }
                            counter++;
                        }
                        myrun.attvals = (MinimalAttractionValue[]) (tmphash.values().toArray(new MinimalAttractionValue[0]));
                        System.out.println("done reading matrix:" + counter);
                        
                        if (counter != seqs) {
                            System.err.println("ERROR, not found the necessary number of positions");
                            inread.close();
                            return myrun;
                        }
                    
                    } else if (inline.equalsIgnoreCase("<seqgroups>")) {
                        // while I am reading the sequence groups
                        String[] tmparr;
                        SequenceGroup currgroup = null;
                        while (((inline = inread.readLine()) != null)
                                && (inline.equalsIgnoreCase("</seqgroups>") == false)) {
                            // skip empty lines
                            if (inline.length() == 0) {
                                continue;
                            }
                            tmparr = inline.split("=", 2);
                            
                            if (tmparr.length != 2) {
                                System.err.println("ERROR reading from savefile on line '" + inline + "'");
                                inread.close();
                                return myrun;
                            }
                            
                            if (tmparr[0].equalsIgnoreCase("name")) {
                                if (currgroup != null) {
                                    myrun.seqgroupsvec.addElement(currgroup);
                                }
                                currgroup = new SequenceGroup();
                                currgroup.name = tmparr[1];
                            } else if (tmparr[0].equalsIgnoreCase("color")) {
                                tmparr = tmparr[1].split(";");
                                try {
                                    int red = Integer.parseInt(tmparr[0]);
                                    int green = Integer.parseInt(tmparr[1]);
                                    int blue = Integer.parseInt(tmparr[2]);
                                    currgroup.color = new java.awt.Color(red, green, blue);
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing numbers from '" + inline + "'");
                                    inread.close();
                                    return myrun;
                                }
                            } else if (tmparr[0].equalsIgnoreCase("numbers")) {
                                tmparr = tmparr[1].split(";");
                                int[] retarr = new int[tmparr.length];
                                try {
                                    for (int i = 0; i < tmparr.length; i++) {
                                        retarr[i] = Integer.parseInt(tmparr[i]);
                                    }
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR parsing numbers from '" + inline + "'");
                                    inread.close();
                                    return myrun;
                                }
                                currgroup.sequences = retarr;
                            } else {
                                System.err.println("Error reading savefile in line" + inline);
                                inread.close();
                                return myrun;
                            }
                        }
                        if (currgroup != null) {
                            myrun.seqgroupsvec.addElement(currgroup);
                        }
                    } else {
                        System.err.println("ERROR, wrong format! unknown keyword on line " + inline);
                        inread.close();
                        return myrun;
                    }
                } else {
                    System.err
                            .println("ERROR, wrong format! Fist line has to start with: sequences=number_of_sequences");
                    inread.close();
                    return myrun;
                }
            }
            inread.close();
        } catch (IOException e) {
            System.err.println("IOError unable to read from " + input_filename.getName());
            return myrun;
        }
        myrun.file = input_filename;
        return myrun;
    }

	/**
	 * Convenience method that calls {@code save_to_file(String, SwingWorker)} with SwingWorker {@code null}. This
	 * method should be used in no-GUI mode.
	 */
	public void save_to_file(java.io.File output_file) throws IllegalStateException, IOException, CancellationException {
		save_to_file(output_file, null);
	}
    
	/**
	 * Writes a CLANS format output file to disk.
	 * 
	 * @param output_file
	 * @param worker
	 *            If not null, this worker is used on a regular basis to check whether saving should be canceled. This
	 *            is used by the GUI to act on user aborts.
	 * @throws IllegalStateException
	 *             If {@code ClusterData.saverun} throws it
	 * @throws IOException
	 *             If {@code ClusterData.saverun} throws it
	 * @throws CancellationException
	 *             If {@code ClusterData.saverun} throws it
	 */
	public void save_to_file(java.io.File output_file, SwingWorker<Void, Integer> worker) throws IllegalStateException,
			IOException, CancellationException {
        ClusterDataLoadHelper myrun = new ClusterDataLoadHelper();
        myrun.autosaveIntervalMinutes = autosaveIntervalMinutes;
        myrun.file = output_file;
        myrun.inaln = sequences;
        myrun.blasthits = blasthits;
        myrun.attvals = attractionValues;
        myrun.posarr = positions;
        myrun.maxmove = maxmove;
        myrun.pval = pvalue_threshold;
        myrun.usescval = usescval;
        if (attvalsimple) {
            myrun.complexatt = false;
        } else {
            myrun.complexatt = true;
        }
        myrun.rotmtx = rotmtx;
        myrun.seqgroupsvec = seqgroupsvec;
        myrun.cooling = cooling;
        myrun.currcool = currcool;
        myrun.attfactor = attfactor;
        myrun.attvalpow = attvalpow;
        myrun.repfactor = repfactor;
        myrun.repvalpow = repvalpow;
        myrun.dampening = dampening;
        myrun.minattract = minattract;
        myrun.blastpath = blastpath;
        myrun.formatdbpath = formatdbpath;
        myrun.dotsize = dotsize;
        myrun.ovalsize = ovalsize;
        myrun.groupsize = groupsize;
        myrun.mapfiles = mapfiles;
        myrun.lookupfiles = lookupfiles;
        myrun.usefoldchange = usefoldchange;
        myrun.avgfoldchange = avgfoldchange;
        myrun.affyfiles = affyfiles;
        myrun.namesdmp_file = namesdmp_file;
        myrun.nodesdmp_file = nodesdmp_file;
        if (cluster2d) {
            myrun.cluster2d = true;
        } else {
            myrun.cluster2d = false;
        }
        myrun.colorarr = colorarr;
        myrun.colorcutoffs = colorcutoffs;

        myrun.rounds = rounds;

        ClusterData.saverun(myrun, sequence_names, nographics, worker);
        myrun = null;
    }
    
	/**
	 * Convenience method that calls {@code safer_save_to_file(String, SwingWorker)} with SwingWorker {@code null}. This
	 * method should be used in no-GUI mode.
	 */
	public void safer_save_to_file(String output_filename) throws IllegalStateException, IOException {
		safer_save_to_file(output_filename, null);
	}
	
	/**
	 * Safer method to save a run that avoids source file corruption at the cost of temporary doubled disk space.
	 * <p>
	 * Internally, this writes to a temporary file in the parent folder of {@code output_file} and finally moves the
	 * temporary file to its final location, which is an instant operation. The almost instant move operation is very
	 * likely to prevent the source file from getting corrupted if CLANS crashes or is force-closed during the save
	 * operation.
	 * 
	 * @param output_file
	 *            the file to which the data should be written. The temporary file will be created in the parent folder
	 *            of output_file.
	 * @param worker
	 *            If not null, this worker is used on a regular basis to check whether saving should be canceled. This
	 *            is used by the GUI to act on user aborts.
	 * @throws IllegalStateException
	 *             if {@code save_to_file} throws it
	 * @throws IOException
	 *             if the temporary file cannot be created or cannot be moved to its final destination, or if
	 *             {@code save_to_file} throws it. All lead to unfinished saves.
	 */
	public void safer_save_to_file(String output_filename, SwingWorker<Void, Integer> worker)
			throws IllegalStateException, IOException {
		String error_message;
		
		File output_file = new File(output_filename);

		// generate temporary filename and save to it
		File temporary_file = null;
		try {
			temporary_file = File.createTempFile("###" + output_file.getName() + "_", ".save_in_progress", new File(
					output_file.getParent()));
		
		} catch (IOException e) {
			error_message = "unable to create temporary file. " + e.getMessage();
			throw new IOException(error_message);
		}

		try {
			save_to_file(temporary_file, worker);

		} catch (CancellationException e) {
			// cleanup temporary file
			if (temporary_file.exists() && !temporary_file.delete()) {
				System.err.println("Failed to delete the temporary save file. "
						+ "This is not critical but leaves a file in your CLANS folder.");
			}
			throw e;
		}

		// move original file to a backup location
		File backup_file = null;
		boolean create_backup = output_file.exists();
		if (create_backup) {
			
			try {
				// create empty temporary file to use its unique filename as backup name
				backup_file = File.createTempFile("###" + output_file.getName() + "_", ".temporary_backup", new File(
						output_file.getParent()));

			} catch (IOException e) {
				error_message = "unable to create file: " + e.getMessage();
				throw new IOException(error_message);

			} finally {
				/**
				 * Remove the empty file created by createTempFile if running on a Windows OS as the call to
				 * File.renameTo below fails on Windows if an empty file exists at the target location. On Linux and
				 * MacOS, this is not necessary.
				 */
				if (OsUtils.isWindows() && backup_file.exists() && !backup_file.delete()) {
					throw new IOException("Failed to delete the empty backup file, which is necessary on Windows.");
				}
			}
			
			try {
				checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
			
			} catch (CancellationException e) {
				throw e;
			}
			
			if (!output_file.renameTo(backup_file)) { // move original file to create backup
				error_message = "unable to move file\n\t" + output_filename + "\nto backup location\n\t"
						+ backup_file.getPath();
				throw new IOException(error_message);
			}
			
			try {
				checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
			
			} catch (CancellationException e) {
				if (!backup_file.renameTo(output_file)) { // try to restore original file by moving backup back
					error_message = "unable to move backup file\n\t" + backup_file.getPath() + "\nback to original location\n\t"
							+ output_filename;
					System.err.println(error_message);
				}
				throw e;
			}
		}
		
		// move temporary file to destination
		if (!temporary_file.renameTo(new File(output_filename))) {
			error_message = "unable to alter " + output_filename;
			throw new IOException(error_message);
		}

		// set source filename to newly saved file
		setInputFilename(output_filename);
		
		// remove backup file
		if (create_backup) {
			if (backup_file.exists() && !backup_file.delete()) {
				System.err.println("Failed to delete the backup save file. "
						+ "This is not critical but leaves a file in your CLANS folder.");
			}
		}
	}

    /**
     * read the seqgroup information from a separate file. the file has to either contain sequence names or numbers. if
     * numbers, the names of the sequences have to be specified as well. read the groups, assign sequence names where
     * necessary and find those names in namearr then assign groups based on the position of the sequence in namearr
     * 
     * @param input_filename
     */
    public void append_groups_or_clusters_from_file(File input_filename) {         
        try {
            BufferedReader inread = new BufferedReader(new FileReader(input_filename));
            String firstline = inread.readLine();
            inread.close();
            if (firstline.startsWith("<ids>")) {
                // I want to read cluster data not group data
                System.out.println("reading cluster data");
                append_clusters_from_file(input_filename);
            }else{
            	// else do the usual read
            	System.out.println("reading sequence group data");
            }
            inread = new BufferedReader(new FileReader(input_filename));
            Vector<String> tmpvec_strings = new Vector<String>();
            Vector<Integer> tmpvec_integers = new Vector<Integer>();
            HashMap<String, Integer> mynamehash = new HashMap<String, Integer>();
            String[] tmparr;
            String tmpstr = "";
            //System.out.println("current dataset seqnum:"+sequence_names.length);
            for (int i = 0; i < sequence_names.length; i++) {
                tmparr = sequence_names[i].split("\\s+");
                mynamehash.put(tmparr[0], new Integer(i));
            }            
            String[] mynamearr = null;
            int maxnum = 0;
            String inline;
            while ((inline = inread.readLine()) != null) {
                if (inline.length() == 0) {
                    continue;
                }
                if (inline.equalsIgnoreCase("<seq>")) {
                    // read the sequence names used in this file
                    while ((inline = inread.readLine()) != null) {
                        if (inline.length() == 0) {
                            continue;
                        }
                        if (inline.equalsIgnoreCase("</seq>")) {
                        	break;
                        } else if (inline.startsWith(">")) {
                            tmpstr = inline.substring(1).trim();
                            tmparr = tmpstr.split("\\s+");
                            tmpvec_strings.addElement(tmparr[0]);
                        }
                    }
                    maxnum = tmpvec_strings.size();
                    mynamearr = new String[maxnum];
                    tmpvec_strings.copyInto(mynamearr);
                    //System.out.println("Append dataset seqnum:"+tmpvec_strings.size());
                    tmpvec_strings.clear();
                } else if (inline.equalsIgnoreCase("<seqgroups>")) {
                    // read the groups defined in this file; format is CLANS
                	//System.out.println("in seqgroups");
                    if (maxnum == 0) {
                        System.err.println("ERROR; missing sequence names");
                        inread.close();
                        return;
                    }
                    int count = 0;
                    java.awt.Color color = java.awt.Color.red;
                    String name = null;
                    int type = 0;
                    int size = 0, tmpval;
                    boolean hide = false;
                    while ((inline = inread.readLine()) != null) {
                        if (inline.length() == 0) {
                            continue;
                        }
                        if (inline.equalsIgnoreCase("</seqgroups>")) {
                        	//System.out.println("END seqgroups");
                            break;
                        }
                        if (inline.startsWith("name=")) {
                        	name = inline.substring(5).trim();
                            type = 0;
                            //System.out.println("\tname="+name);
                        } else if (inline.startsWith("type=")) {
                            try {
                                tmpstr = inline.substring(5).trim();
                                type = Integer.parseInt(tmpstr);
                            } catch (NumberFormatException ne) {
                                System.err.println("ERROR, unable to parse int from '" + tmpstr + "'");
                            }
                        } else if (inline.startsWith("size=")) {
                            try {
                                tmpstr = inline.substring(5).trim();
                                size = Integer.parseInt(tmpstr);
                            } catch (NumberFormatException ne) {
                                System.err.println("ERROR, unable to parse int from '" + tmpstr + "'");
                            }
                        } else if (inline.startsWith("hide=")) {
                            try {
                                tmpstr = inline.substring(5).trim();
                                tmpval = Integer.parseInt(tmpstr);
                                if (tmpval == 1) {
                                    hide = true;
                                } else {
                                    hide = false;
                                }
                            } catch (NumberFormatException ne) {
                                System.err.println("ERROR, unable to parse int from '" + tmpstr + "'");
                            }
                        } else if (inline.startsWith("numbers=")) {
                            count++;
                            SequenceGroup mygroup = new SequenceGroup();
                            if (name != null) {
                                mygroup.name = name;
                                name = null;
                            } else {
                                mygroup.name = String.valueOf(count);
                            }
                            mygroup.type = type;
                            mygroup.color = color;
                            mygroup.hide = hide;
                            mygroup.size = size;
                            mygroup.polygon = Shapes.get(mygroup.type, mygroup.size);
                            color = java.awt.Color.red;
                            // now convert the numbers you read to the numbers used in namearr
                            tmparr = inline.substring(8).trim().split(";");
                            int[] tmp = new int[tmparr.length];
                            try {
                                for (int i = 0; i < tmparr.length; i++) {
                                    tmpstr = tmparr[i];
                                    tmp[i] = Integer.parseInt(tmpstr);
                                }
                            } catch (NumberFormatException ne) {
                                System.err.println("Unable to parse int from '" + tmpstr + "' in line '" + inline + "'");
                                continue;
                            }
                            //tmpvec_strings.clear();
                            tmpvec_integers.clear();
                            for (int i = 0; i < tmparr.length; i++) {
                                tmpstr = mynamearr[tmp[i]];
                                if (mynamehash.containsKey(tmpstr)) {
                                    tmpvec_integers.addElement(mynamehash.get(tmpstr));
                                } else {
                                    System.err.println("no name '" + tmpstr + "' found in current graph; skipping entry");
                                }
                            }
                            //elements = tmpvec_strings.size();
                            elements = tmpvec_integers.size();
                            mygroup.sequences = new int[elements];
                            for (int i = 0; i < elements; i++) {
                                mygroup.sequences[i] = ((Integer) tmpvec_integers.elementAt(i)).intValue();
                            }
                            if (elements > 0) {
                            	//System.out.println("Adding group '"+mygroup.name+"' ("+seqgroupsvec.size()+") with elements '"+elements);
                                seqgroupsvec.addElement(mygroup);
                            }else{
                            	//System.out.println("group size = 0");
                            }
                        } else if (inline.startsWith("color=")) {
                            tmpstr = inline.substring(6).trim();
                            tmparr = tmpstr.split(";");
                            if ((tmparr.length != 3) && (tmparr.length != 4)) {
                                System.err.println("ERROR reading color info from line '" + inline + "'");
                            }
                            int red, green, blue;
                            try {
                                tmpstr = tmparr[0];
                                red = Integer.parseInt(tmpstr);
                                tmpstr = tmparr[1];
                                green = Integer.parseInt(tmpstr);
                                tmpstr = tmparr[2];
                                blue = Integer.parseInt(tmpstr);
                                if (tmparr.length < 4) {
                                    color = new java.awt.Color(red, green, blue);
                                } else {
                                    tmpstr = tmparr[3];
                                    int alpha = Integer.parseInt(tmpstr);
                                    color = new java.awt.Color(red, green, blue, alpha);
                                }
                            } catch (NumberFormatException ne) {
                                System.err.println("unable to convert '" + tmpstr + "' to int in line '" + inline + "'");
                                color = java.awt.Color.red;
                            }
                        }
                    }
                }
            }
            inread.close();
        } catch (IOException ioe) {
            System.err.println("IOERROR; unable to read from file '" + input_filename.getAbsolutePath() + "'");
        }
    }//end append_groups_or_clusters_from_file(File input_filename)

    /**
     * read the cluster data generated by the iterative clustering approach (append groups to groupsvec)
     * 
     * @param input_filename
     */
    private void append_clusters_from_file(File input_filename) {
        HashMap<String, Integer> mynamehash = new HashMap<String, Integer>();
        String[] tmparr;
        String tmpstr = "";

        for (int i = 0; i < sequence_names.length; i++) {
            tmparr = sequence_names[i].split("\\s+");
            if (mynamehash.containsKey(tmparr[0])) {
                System.err.println("ERROR: non unique identifier in parent: name '" + tmparr[0]
                        + "' already defined as -->'" + mynamehash.get(tmparr[0]) + "'");
                return;
            } else {
                mynamehash.put(tmparr[0], new Integer(i));
            }
        }
        
        try {
            BufferedReader inread = new BufferedReader(new FileReader(input_filename));
            String inline;
            HashMap<String, Integer> clusternamehash = new HashMap<String, Integer>();
            String numberstring = "";
            String namestring = "";
            
            while ((inline = inread.readLine()) != null) {
            
                if (inline.equalsIgnoreCase("<ids>")) {
                    System.out.println("reading cluster IDS");
                    boolean stopids = false;
                    
                    while (stopids == false) {
                        inline = inread.readLine();
                        
                        if (inline == null) {
                            System.err.println("ERROR, reached EOF before end of <ids>");
                            inread.close();
                            return;
                        
                        } else if (inline.equalsIgnoreCase("</ids>")) {
                            stopids = true;
                        
                        } else {
                            tmparr = inline.split("\\s+", 2);
                            numberstring = tmparr[0];
                            namestring = tmparr[1];
                            int spaceidx = namestring.indexOf(" ", 1);
                            int quoteidx = namestring.indexOf("\"", 1);
                            if (spaceidx < quoteidx) {
                                namestring = namestring.substring(1, spaceidx);
                            } else {
                                namestring = namestring.substring(1, quoteidx);
                            }
                            if (mynamehash.containsKey(namestring)) {
                                clusternamehash.put(numberstring, mynamehash.get(namestring));
                            } else {
                                System.err.println("WARNING: name '" + namestring + "' is not known in parent");
                                clusternamehash.put(numberstring, null);
                            }
                        }
                    }
                    System.out.println("DONE reading cluster IDS; found:" + clusternamehash.size());
                
                } else if (inline.equalsIgnoreCase("<clusters>")) {
                    tmpstr = inread.readLine();// next line contains offset=number (substring position 7)
                    float offset = 0;
                    
                    try {
                        tmpstr = tmpstr.substring(7);
                        offset = Float.parseFloat(tmpstr);
                    
                    } catch (NumberFormatException ne) {
                        System.err.println("ERROR trying to parse float from '" + tmpstr + "'");
                        inread.close();
                        return;
                    }
                    
                    System.out.println("reading clusters; offset=" + offset);
                    boolean stopclusters = false;
              
                    while (stopclusters == false) {
                        inline = inread.readLine();

                        if (inline == null) {
                            System.err.println("ERROR, reached EOF before end of <clusters>");
                            inread.close();
                            return;

                        } else if (inline.equalsIgnoreCase("</clusters>")) {
                            stopclusters = true;

                        } else if (inline.equalsIgnoreCase("<cluster>")) {
                            SequenceGroup mycluster = new SequenceGroup();
                            int clustersize = -1;
                            String nameline = inread.readLine();
                            String sizeline = inread.readLine();
                            String elementline = inread.readLine();
                            tmpstr = inread.readLine();// this should be </cluster>

                            if ((tmpstr.equalsIgnoreCase("</cluster>") == false)
                                    || (nameline.startsWith("name=") == false)
                                    || (sizeline.startsWith("size=") == false)
                                    || (elementline.startsWith("elements=") == false)) {
                                System.err.println("ERROR parsing data from cluster:");
                                System.err.println("<cluster>");
                                System.err.println(nameline);
                                System.err.println(sizeline);
                                System.err.println(elementline);
                                System.err.println(tmpstr);
                            
                            } else {
                                mycluster.name = nameline.substring(5) + "_of:" + offset;
                                try {
                                    clustersize = Integer.parseInt(sizeline.substring(5));
                                } catch (NumberFormatException ne) {
                                    System.err.println("ERROR trying to parse int from line '" + sizeline + "'");
                                    inread.close();
                                    return;
                                }
                                tmparr = elementline.substring(9).split(";");
                                if (tmparr.length != clustersize) {
                                    System.err.println("WARNING: unequal cluster size: should=" + clustersize + " is="
                                            + tmparr.length);
                                    clustersize = tmparr.length;
                                }
                                mycluster.sequences = new int[clustersize];
                                Integer myint;
                                for (int i = 0; i < clustersize; i++) {
                                    if ((myint = clusternamehash.get(tmparr[i])) != null) {
                                        mycluster.sequences[i] = myint.intValue();
                                    } else {
                                        System.err.println("WARNING: undefined name for '" + tmparr[i] + "'");
                                    }
                                }

                                seqgroupsvec.add(mycluster);
                            }
                        } else {
                            System.out.println("read unknown line '" + inline + "'");
                        }
                    }
                }
            }
            inread.close();
        
        } catch (IOException ioe) {
            System.err.println("IOERROR; unable to read from file '" + input_filename.getAbsolutePath() + "'");
        }
        System.out.println("read " + seqgroupsvec.size() + " cluster entries");
    }

	/**
	 * Saves the data to the selected file
	 * 
	 * @param input
	 * @param sequence_names
	 * @param nographics
	 *            true if the GUI is used, false for commandline-only CLANS
	 * @param worker
	 *            If not null, this worker is used on a regular basis to check whether saving should be canceled. This
	 *            is used by the GUI to act on user aborts.
	 * @throws IllegalStateException
	 *             If neither HSPs nor attraction values are set in {@code input}
	 * @throws IOException
	 *             If FileWriter throws it
	 * @throws CancellationException
	 *             If worker is cancelled and saving should therefore be cancelled. In headless (no-GUI) mode this
	 *             Exception will not occur.
	 */
	private static void saverun(ClusterDataLoadHelper input, String[] sequence_names, boolean nographics,
			SwingWorker<Void, Integer> worker) throws IllegalStateException, IOException, CancellationException {

		checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
		
		// check if necessary data is available
		if (input.blasthits == null && input.attvals == null) {
			throw new IllegalStateException("Error while saving run: no attraction values or list of HSPs found");
		}

		ComfortableBufferedWriter outwrite = new ComfortableBufferedWriter(new FileWriter(input.file));
		
		try {
			
	        outwrite.println("sequences=" + input.inaln.length);
	
	        outwrite.println("<param>");
	
	        outwrite.println("attfactor=" + input.attfactor);
	        outwrite.println("attvalpow=" + input.attvalpow);
	        outwrite.println("autosaveinterval=" + input.autosaveIntervalMinutes);
	        outwrite.println("avgfoldchange=" + input.avgfoldchange);
	
	        outwrite.println("blastpath=" + input.blastpath);
	
	        outwrite.println("cluster2d=" + input.cluster2d);
	
	        outwrite.print("colorarr=");
	        for (int i = 0; i < input.colorarr.length; i++) {
	            java.awt.Color tmp = input.colorarr[i];
	            outwrite.print("(" + tmp.getRed() + ";" + tmp.getGreen() + ";" + tmp.getBlue() + "):");
	        }
	        outwrite.println();
	
	        outwrite.print("colorcutoffs=");
	        for (int i = 0; i < input.colorcutoffs.length; i++) {
	            outwrite.print(input.colorcutoffs[i] + ";");
	        }
	        outwrite.println();
	
	        outwrite.println("complexatt=" + input.complexatt);
	        outwrite.println("cooling=" + input.cooling);
	        outwrite.println("currcool=" + input.currcool);
	
	        outwrite.println("dampening=" + input.dampening);
	        outwrite.println("dotsize=" + input.dotsize);
	
	        outwrite.println("formatdbpath=" + input.formatdbpath);
	
	        outwrite.println("groupsize=" + input.groupsize);
	
	        outwrite.println("maxmove=" + input.maxmove);
	        outwrite.println("minattract=" + input.minattract);
	
	        if (input.namesdmp_file != null && input.nodesdmp_file != null) {
	            outwrite.println("namesdmp_file=" + input.namesdmp_file);
	            outwrite.println("nodesdmp_file=" + input.nodesdmp_file);
	        }
	
	        outwrite.println("ovalsize=" + input.ovalsize);
	
	        outwrite.println("pval=" + input.pval);
	
	        outwrite.println("repfactor=" + input.repfactor);
	        outwrite.println("repvalpow=" + input.repvalpow);
	        outwrite.println("rounds_done=" + input.rounds);
	
	        outwrite.println("showinfo=" + input.showinfo);
	
	        outwrite.println("usefoldchange=" + input.usefoldchange);
	        outwrite.println("usescval=" + input.usescval);
	
	        outwrite.println("zoom=" + input.zoom);
	
	        outwrite.println("</param>");
	        
	        checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
	
	        if ((input.mapfiles != null) && (input.mapfiles.size() > 0)) {
	            outwrite.println("<function>");
	            int num = input.mapfiles.size();
	            for (int i = 0; i < num; i++) {
	                if ((input.lookupfiles != null) && (input.lookupfiles.get(i) != null)) {
	                    outwrite.println(((File) input.mapfiles.get(i)).getAbsolutePath() + "';'"
	                            + ((File) input.lookupfiles.get(i)).getAbsolutePath());
	                } else {
	                    outwrite.println(((File) input.mapfiles.get(i)).getAbsolutePath() + "';'NONE");
	                }
	            }
	            outwrite.println("</function>");
	        }
	        
	        checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
	        
	        if (input.affyfiles != null) {
	            outwrite.println("<affyfiles>");
	            int repnum = input.affyfiles.size();
	            Replicates rep;
	            for (int i = 0; i < repnum; i++) {
	                rep = (Replicates) (input.affyfiles.get(i));
	                outwrite.println("<");
	                outwrite.println("abbreviation=" + rep.abbreviation);
	                outwrite.println("replicates=" + rep.replicates);
	                outwrite.println("wtreplicates=" + rep.wtreplicates);
	                outwrite.println("name=" + rep.name);
	                outwrite.println("wtname=" + rep.wtname);
	                outwrite.print("replicate=");
	                for (int j = rep.replicate.length; --j >= 0;) {
	                    outwrite.print(rep.replicate[j].getAbsolutePath() + "';'");
	                }
	                outwrite.println();
	                outwrite.print("wtreplicate=");
	                for (int j = rep.wtreplicate.length; --j >= 0;) {
	                    outwrite.print(rep.wtreplicate[j].getAbsolutePath() + "';'");
	                }
	                outwrite.println();
	                outwrite.println(">");
	            }
	            outwrite.println("</affyfiles>");
	        }
	        
	        checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
	        
	        outwrite.println("<rotmtx>");
	        for (int i = 0; i < 3; i++) {
	            for (int j = 0; j < 3; j++) {
	                outwrite.print(input.rotmtx[i][j] + ";");
	            }
	            outwrite.println();
	        }
	        outwrite.println("</rotmtx>");
	        
	        checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
	        
	        // first write the sequences to file
	        outwrite.println("<seq>");
	        for (int i = 0; i < input.inaln.length; i++) {
	            outwrite.println(">" + sequence_names[i]);
	            outwrite.println(input.inaln[i].seq);
	        }
	        outwrite.println("</seq>");
	
	        checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
	        
	        // write the sequence weights
	        if (input.weights != null) {
	            outwrite.println("<weight>");
	            for (int i = 0; i < input.weights.length; i++) {
	                outwrite.println(">" + sequence_names[i]);
	                outwrite.println(input.weights[i]);
	            }
	            outwrite.println("</weight>");
	        }
	        
	        checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
	        
	        // write the sequence groups
	        if ((input.seqgroupsvec != null) && (input.seqgroupsvec.size() > 0)) {
	            outwrite.println("<seqgroups>");
	            SequenceGroup mygroup;
	            for (int i = 0; i < input.seqgroupsvec.size(); i++) {
	                mygroup = input.seqgroupsvec.elementAt(i);
	                
	                outwrite.println("name=" + mygroup.name);
	                outwrite.println("type=" + mygroup.type);
	                outwrite.println("size=" + mygroup.size);
	                
	                if (mygroup.hide == true) {
	                    outwrite.println("hide=1");
	                } else {
	                    outwrite.println("hide=0");
	                }
	                
	                outwrite.println("color=" + mygroup.color.getRed() + ";" + mygroup.color.getGreen() + ";"
	                        + mygroup.color.getBlue() + ";" + mygroup.color.getAlpha());
	                outwrite.print("numbers=");
	
	                Arrays.sort(mygroup.sequences);
	                for (int j = 0; j < mygroup.sequences.length; j++) {
	                    outwrite.print(mygroup.sequences[j] + ";");
	                }
	                outwrite.println();
	            }
	            outwrite.println("</seqgroups>");
	        }
	        
	        checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
	        
	        // next write the sequence positions
	        outwrite.println("<pos>");
	        for (int i = 0; i < input.inaln.length; i++) {
	            outwrite.println(i + " " + input.posarr[i][0] + " " + input.posarr[i][1] + " " + input.posarr[i][2]);
	        }
	        outwrite.println("</pos>");
	        
	        checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
	
			if (input.blasthits != null) {
				// write blast HSPs
				int tmpsize = 0;
	
				outwrite.println("<hsp>");
	
				for (int i = 0; i < input.blasthits.length; i++) {
					checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
					
					outwrite.print(input.blasthits[i].query + " " + input.blasthits[i].hit + ":");
					tmpsize = input.blasthits[i].val.length;
	
					for (int j = 0; j < tmpsize; j++) {
						outwrite.print(input.blasthits[i].val[j] + " ");
					}
					outwrite.println();
				}
	
				outwrite.println("</hsp>");
	
			} else { // we check blasthits == null and attvals == null at the beginning, thus input.attvals is not null here
				// write attvals
				outwrite.println("<att>");
	
				for (int i = 0; i < input.attvals.length; i++) {
					outwrite.println(input.attvals[i].query + " " + input.attvals[i].hit + " " + input.attvals[i].att);
				}
	
				outwrite.println("</att>");
			}
			
			checkWorkerStatus(worker); // stop saving if told so by user (only in GUI mode)
			
		} finally {
			outwrite.close();
		}
	}

    /**
     * Saves the attraction values to a file.
     * @param output_file the file to which the output is written
     */
    public void save_attraction_values_to_file(File output_file) {
        try {
            PrintWriter outwrite = new PrintWriter(new BufferedWriter(new FileWriter(output_file)));
            MinimalAttractionValue myatt;
            for (int i = attractionValues.length; --i >= 0;) {
                myatt = attractionValues[i];
                outwrite.println(myatt.query + " " + myatt.hit + " " + myatt.att);
            }
            outwrite.close();
        } catch (IOException e) {
            System.err.println("IOError writing to " + output_file.getName());
            if (nographics == false) {
                javax.swing.JOptionPane.showMessageDialog(new javax.swing.JFrame(),
                        "IOERROR writing to '" + output_file.getName() + "'");
            }
        }
    }

    /**
     * Computes the attraction values.
     */
    public void compute_attraction_values() {

        if (blasthits == null) {// possible (if alternate data source was loaded)
            System.out.println("cannot compute attraction values; no BLAST HSPs present");
            return;
        }

        ArrayList<MinimalAttractionValue> tmpvec = new ArrayList<MinimalAttractionValue>();
        int number_of_blasthits = blasthits.length;
        HashMap<MinimalAttractionValue, MinimalAttractionValue> myhash = new HashMap<MinimalAttractionValue, MinimalAttractionValue>(number_of_blasthits);
        float maxattval = 0;
        maxvalfound = 0;// init to zero, is assigned value in getattvalsimple or mult
        // NOTE: this is not necessarily a symmetrical array. compute all values
        // and then symmetrize computing the average values

        if (rescalepvalues == false) {
            // make the attraction values
            if (attvalsimple) {
                for (int i = 0; i < number_of_blasthits; i++) {
                    MinimalAttractionValue curratt = new MinimalAttractionValue(blasthits[i].query, blasthits[i].hit);
                    if (myhash.containsKey(curratt)) {
                        curratt = myhash.get(curratt);
                        if (curratt.att == -1) {
                            // in this case keep the -1
                        } else {
                            float newatt = ClusterMethods.computeSimpleAttractionValue(blasthits[i].val, this);
                            if (newatt == -1) {
                                curratt.att = -1;
                            } else {
                                newatt /= 2;
                                curratt.att += newatt;
                            }
                        }
                    } else {
                        // if I've never encountered this query-hit pair before
                        curratt.att = ClusterMethods.computeSimpleAttractionValue(blasthits[i].val, this);
                        if (curratt.att != -1) {
                            curratt.att /= 2;
                        }
                        if (curratt.att != 0) {
                            myhash.put(curratt, curratt);
                            tmpvec.add(curratt);
                        }
                    }
                    if (curratt.att > maxattval) {
                        maxattval = curratt.att;
                    }
                }
            } else {
                for (int i = 0; i < number_of_blasthits; i++) {
                    MinimalAttractionValue curratt = new MinimalAttractionValue(blasthits[i].query, blasthits[i].hit);
                    if (myhash.containsKey(curratt)) {
                        curratt = myhash.get(curratt);
                        if (curratt.att == -1) {
                            // in this case keep the -1
                        } else {
                            float newatt = ClusterMethods.computeComplexAttractionValue(blasthits[i].val, this);
                            if (newatt == -1) {
                                curratt.att = -1;
                            } else {
                                newatt /= 2;
                                curratt.att += newatt;
                            }
                        }
                    } else {
                        // if I've never encountered this query-hit pair before
                        curratt.att = ClusterMethods.computeComplexAttractionValue(blasthits[i].val, this);
                        if (curratt.att != -1) {
                            curratt.att /= 2;
                        }
                        if (curratt.att != 0) {
                            myhash.put(curratt, curratt);
                            tmpvec.add(curratt);
                        }
                    }
                    if (curratt.att > maxattval) {
                        maxattval = curratt.att;
                    }
                }
            }
            // divide all vals by maxattval (-->range: 0-1)
            // standard, just divide all values by the maximum value
            // note, this does NOT symmetrize the attractions
            if (usescval == false) {
                for (int i = 0; i < tmpvec.size(); i++) {
                    if (tmpvec.get(i).att == -1) {
                        tmpvec.get(i).att = 1;

                    } else {
                        tmpvec.get(i).att /= maxattval;
                    }
                }
                
                p2attfactor = maxattval;
                p2attoffset = 0;
            } else {// if using scval
                p2attfactor = 1;
                p2attoffset = 0;
            }
        } else {// if rescalepvaluecheckbox==true
            float minattval = java.lang.Float.MAX_VALUE;
            // rescale the attraction values to range from 0 to 1 (with the smallest positive non-zero value as zero.
            if (attvalsimple) {
                for (int i = 0; i < number_of_blasthits; i++) {
                    MinimalAttractionValue curratt = new MinimalAttractionValue(blasthits[i].query, blasthits[i].hit);
                    if (myhash.containsKey(curratt)) {
                        curratt =  myhash.get(curratt);
                        if (curratt.att == -1) {
                            // in this case keep the -1
                        } else {
                            float newatt = ClusterMethods.computeSimpleAttractionValue(blasthits[i].val, this);
                            if (newatt == -1) {
                                curratt.att = -1;
                            } else {
                                newatt /= 2;
                                curratt.att += newatt;
                            }
                        }
                    } else {
                        // if I've never encountered this query-hit pair before
                        curratt.att = ClusterMethods.computeSimpleAttractionValue(blasthits[i].val, this);
                        if (curratt.att != -1) {
                            curratt.att /= 2;
                        }
                        if (curratt.att != 0) {
                            myhash.put(curratt, curratt);
                            tmpvec.add(curratt);
                        }
                    }
                    if (curratt.att > maxattval) {
                        maxattval = curratt.att;
                    }
                    if ((curratt.att > 0) && (curratt.att < minattval)) {
                        minattval = curratt.att;
                    }
                }
            } else {
                for (int i = 0; i < number_of_blasthits; i++) {
                    MinimalAttractionValue curratt = new MinimalAttractionValue(blasthits[i].query, blasthits[i].hit);
                    if (myhash.containsKey(curratt)) {
                        curratt = myhash.get(curratt);
                        if (curratt.att == -1) {
                            // in this case keep the -1
                        } else {
                            float newatt = ClusterMethods.computeComplexAttractionValue(blasthits[i].val, this);
                            if (newatt == -1) {
                                curratt.att = -1;
                            } else {
                                newatt /= 2;
                                curratt.att += newatt;
                            }
                        }

                    } else {
                        // if I've never encountered this query-hit pair before
                        curratt.att = ClusterMethods.computeComplexAttractionValue(blasthits[i].val, this);
                        if (curratt.att != -1) {
                            curratt.att /= 2;
                        }
                        if (curratt.att != 0) {
                            myhash.put(curratt, curratt);
                            tmpvec.add(curratt);
                        }

                    }
                    if (curratt.att > maxattval) {
                        maxattval = curratt.att;
                    }
                    if ((curratt.att > 0) && (curratt.att < minattval)) {
                        minattval = curratt.att;
                    }
                }
            }
            // and divide all vals by maxattval and offset by minattval(-->range: 0-1)
            float divval = maxattval - minattval;
            for (int i = 0; i < tmpvec.size(); i++) {
                if (tmpvec.get(i).att == -1) {
                    tmpvec.get(i).att = 1;
                } else {
                    tmpvec.get(i).att = (tmpvec.get(i).att - minattval) / divval;
                }
            }

            p2attfactor = divval;
            p2attoffset = minattval;
        }
        attractionValues = tmpvec.toArray(new MinimalAttractionValue[0]);

		String formatted_used_hsps = String.format("%" + Integer.toString(blasthits.length).length() + "s",
				attractionValues.length);
		String formatted_pvalue = String.format("%7s", new DecimalFormat("0.###E0").format(pvalue_threshold));
		String formatted_percentage = String.format("%7s",
				new DecimalFormat("0.00").format(100 * (float) (attractionValues.length) / blasthits.length) + "%");
		System.out.println("used/total HSPs at threshold " + formatted_pvalue + ": " + formatted_used_hsps + "/"
				+ blasthits.length + " (" + formatted_percentage + ")");
	}

    /**
     * Resets the graph with random sequence positions and recomputes "attraction" values from HSPs.
     */
    public void initialize() {
        rounds = 0;
        
        currcool = 1;

        reset_rotmtx();

        movements = null;
        movementsLastIteration = null;

        elements = sequence_names.length;

        if (elements == 0) { // no elements means nothing to do
            positions = new float[0][0];
            posarr = positions;
            return;
        }

        positions = new float[elements][dimensions];
        posarrtmp = new float[elements][dimensions];
        drawarrtmp = new int[elements][dimensions];
        
        movements = new float[elements][dimensions];
        movementsLastIteration = new float[elements][dimensions];
        for (int i = 0; i < elements; i++) {
            for (int j = 0; j < dimensions; j++) {
                movements[i][j] = 0;
                movementsLastIteration[i][j] = 0;
            }
        }

        // compute the "attraction values"
        if (blasthits != null) {
            if (attractionValues == null) {
                compute_attraction_values();
            }
        }

        for (int i = 0; i < positions.length; i++) {
            for (int j = 0; j < positions[j].length; j++) {
                
                if (cluster2d && j > 1) { // leave last entry 0
                    break;
                }
                
                positions[i][j] = ClusterMethods.rand.nextFloat() * 2 - 1;
            }
        }

        posarr = positions;
    }

    /**
     * Resets the rotation matrix to represent no rotation.
     */
    private void reset_rotmtx() {
    	MyMath.setTo3x3IdentityMatrix(rotmtx);
    	MyMath.setTo3x3IdentityMatrix(myrotmtx);
    }
    
    /**
     * adds a new group to the set of groups
     * 
     * @param name
     *            the new group's name
     * @param members
     *            indices of group members in the list of all entries of the clustering
     */
    public void add_group(String name, int[] members) {
        
        SequenceGroup group = new SequenceGroup();
        
        group.name = name;
        group.sequences = members;

        // initialize other values with defaults
        group.color = java.awt.Color.red;
        group.type = 0;
        group.size = groupsize;

        seqgroupsvec.addElement(group);
    }
    
    /**
     * removes the group at position index
     * 
     * @param index
     *            the index of the group to be removed
     * @throws IndexOutOfBoundsException
     *             if index > (number of existing SequenceGroups - 1)
     */
    public void remove_group(int index) throws IndexOutOfBoundsException {
        if (index > seqgroupsvec.size()) {
            throw new IndexOutOfBoundsException("SequenceGroup with index " + index + " cannot be removed as only "
                    + seqgroupsvec.size() + " groups exist.");
        }
        seqgroupsvec.remove(index);
    }
    
	/**
	 * @return true if the user set up a limit for the number of rounds to run before stopping automatically.
	 */
	public boolean hasRoundsLimit() {
		return roundsLimitBeforeStopping > 0;
	}
	
	/**
	 * @return the number of rounds to run before stopping automatically or -1 if no limit is currently in place.
	 */
	public int getRoundsLimit() {
		return roundsLimitBeforeStopping;
	}
	
	/**
	 * Sets the number of rounds to run before stopping automatically.
	 * 
	 * @param new_limit
	 *            The number of rounds to run before stopping.
	 */
	public void setRoundsLimit(int new_limit) {
		roundsLimitBeforeStopping = new_limit;
		roundsCompleted = 0;
	}
	
	/**
	 * Disables automatic stopping after a fixed number of rounds.
	 */
	public void disableRoundsLimit() {
		setRoundsLimit(-1);
	}
    
    /**
     * Used to indicate worker cancellation. 
     * @param worker null if no worker is present or the worker
     * @throws CancellationException If {@code worker} is not {@code null} and the worker is in cancelled state.
     */
    public static void checkWorkerStatus(SwingWorker<Void, Integer> worker) throws CancellationException {
		if (worker != null && worker.isCancelled()) {
			throw new CancellationException("worker has been canceled");
		}
    }
}
