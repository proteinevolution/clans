package clans;

import java.util.*;
import java.io.*;

public class searchblastv2 {

    //allow a few more parameters when searching for blast hits
    //1. allow block searches (i.e. run multiple queries per search)
    //2. always use a thread
    //3. allow use of blast+ as well as legacy blast
    public searchblastv2(StringBuffer errbuff, boolean addblastvbparam, boolean isblastplus, int cpu, int blastblocksize, String cmd, String blastpath, String formatdbpath, double eval, double pval, float coverage, float scval, float ident, int verbose, String saveblastname, boolean readblast, HashMap<String,Integer> nameshash) {
        this.errbuff = errbuff;
        this.addblastvbparam = addblastvbparam;
        this.isblastplus = isblastplus;
        this.cpu = cpu;
        this.blastblocksize = blastblocksize;
        this.cmd = cmd;
        this.blastpath = blastpath;
        this.formatdbpath = formatdbpath;
        this.eval = eval;
        this.pval = pval;
        this.coverage = coverage;
        this.scval = scval;
        this.ident = ident;
        this.verbose = verbose;
        this.saveblastname = saveblastname;
        this.nameshash = nameshash;
        this.readblast = readblast;
    }
    //the global variable that are the same for all blast runs
    StringBuffer errbuff;
    boolean addblastvbparam;
    boolean isblastplus;
    boolean readblast = true;//check for old blast results
    final Vector<AminoAcidSequence> queryvec = new Vector<AminoAcidSequence>();//will synchronize over this, hence final!
    final Integer syncdone = new Integer(1);//synchronize the "thread done" checks on this one
    final Integer filesync = new Integer(1);//synchronize writing to a file on this one
    int cpu = 1;
    int blastblocksize = 1;
    String cmd;
    String blastpath;
    String formatdbpath;
    String saveblastname;
    double eval;
    double pval;
    float coverage;
    float scval;
    float ident;
    int verbose;
    PrintWriter outwrite;
    HashMap<String, Integer> nameshash;

    //I need a method to get hits for a full set of sequences
    public MinimalHsp[] gethits(AminoAcidSequence[] queryseqs) {
        int allseqnum = queryseqs.length;
        //first set up the queryvec with the sequences to search for
        for (int i = queryseqs.length; --i >= 0;) {
            queryvec.add(queryseqs[i]);
        }//end for i
        //now read through former BLAST results and ONLY repeat those sequences that were not already completed
        //completed searches should be present with a line # done for sequence_name
        if (new File(saveblastname).canRead() && readblast) {
            System.out.println("checking old blast output found in '" + saveblastname + "'");
            HashMap<String, Integer> donehash = new HashMap<String, Integer>();
            try {
                BufferedReader inread = new BufferedReader(new FileReader(saveblastname));
                String inline, donename;
                while ((inline = inread.readLine()) != null) {
                    if (inline.startsWith("# done for:")) {
                        donename = inline.substring(11).trim();
                        System.out.println(donename + " completed");
                        donehash.put(donename, null);
                    }
                }
                inread.close();
            } catch (IOException ioe) {
                System.err.println("IOERROR reading from '" + saveblastname + "'");
            }
            //now remove from the queryvec all those that are already done
            for (int i = queryvec.size(); --i >= 0;) {
                if (donehash.containsKey(queryvec.get(i).name)) {
                    queryvec.remove(i);
                }
            }//end for i
        }//end if readblast
        String dbfilename = "___placeholder___";
        if (queryvec.size() == 0) {
            System.out.println("All sequences accounted for, skipping blast search");
        } else {
            //first see whether I need to update the blastpath to include other things
            if (addblastvbparam == true) {//default: true!
                //i.e. if I want to check to see whether the -v or -b parameters were set in blastpath
                if (allseqnum > 250) {
                    //this is only of relevance if I have more than 250 sequences (default for blast returns)
                    if (isblastplus) {
                        //if I am using the new blastplus
                        if (blastpath.indexOf("-num_descriptions ") > -1 || blastpath.indexOf("-num_alignments ") > -1) {
                            //if this IS set, then don't do anything
                        } else {
                            //set these to allow at least one hit (on average) per sequence
                            this.blastpath += " -num_descriptions " + allseqnum + " -num_alignments " + allseqnum;
                        }
                    } else {
                        //if I am using legacy blast
                        if (blastpath.indexOf("-b ") > -1 || blastpath.indexOf("-v ") > -1) {
                            //if this IS set, then don't do anything
                        } else {
                            //set these to allow at least one hit (on average) per sequence
                            this.blastpath += " -b " + allseqnum + " -v " + allseqnum;
                        }
                    }
                }
            }
            //now write the sequences to a database for BLASTING
            dbfilename = new File(String.valueOf(System.currentTimeMillis() + ".dbfile")).getAbsolutePath();
            try {
                PrintWriter dbwrite = new PrintWriter(new BufferedWriter(new FileWriter(dbfilename)));
                for (int i = allseqnum; --i >= 0;) {
                    dbwrite.println(">" + queryseqs[i].name);
                    dbwrite.println(queryseqs[i].seq.replaceAll("-", ""));
                }//end for i
                dbwrite.close();
            } catch (IOException ioe) {
                System.err.println("IOERROR trying to write to blast database file '" + dbfilename + "'");
                System.err.println("Stopping search.");
                errbuff.append("-Error while creating blast database.\n");
                return null;
            }
            //now format the database for BLAST
            String formatcommand = "";
            if (isblastplus) {
                formatcommand = formatdbpath + " -in " + dbfilename;
            } else {
                formatcommand = formatdbpath + " -i " + dbfilename;//for the older blast version
            }
            if (cmd.length() > 0) {
                formatcommand = cmd + " " + formatcommand;
            }
            if (verbose > 1) {
                System.out.println("doing " + formatcommand);
            }
            StringBuffer errout = new StringBuffer();
            StringBuffer inout = new StringBuffer();
            threadstreamreader errread;
            threadstreamreader inread;
            Runtime rt = Runtime.getRuntime();
            try {
                Process p = rt.exec(formatcommand);
                errread = new threadstreamreader(new BufferedReader(new InputStreamReader(p.getErrorStream())), errout);
                inread = new threadstreamreader(new BufferedReader(new InputStreamReader(p.getInputStream())), inout);
                try {// wait for formatdb to finish
                    errread.start();
                    inread.start();
                    if (verbose > 1) {
                        System.out.println("waiting for " + formatcommand);
                    }
                    p.waitFor();
                } catch (InterruptedException e) {
                    System.err.println("ERROR Interrupted formatdb");
                }
                System.err.print(errout.toString());
            } catch (IOException ex2) {
                System.err.println("ERROR starting process " + formatcommand);
                errbuff.append("ERROR starting process " + formatcommand + "\n");
                ex2.printStackTrace();
            }
            if (verbose > 1) {
                System.out.println("done formatting blast database with :" + formatcommand);
            }
            //now define the temporary file to which to write the results
            try {
                outwrite = new PrintWriter(new BufferedWriter(new FileWriter(saveblastname, true)));//append to file
            } catch (IOException ioe) {
                System.err.println("IOEROR trying to open " + saveblastname + " for writing (appending)");
            }
            //now start the blast threads for blasting and get the results
            blastthread[] mythreads = new blastthread[cpu];
            for (int i = 0; i < cpu; i++) {
                mythreads[i] = new blastthread(i, blastblocksize, dbfilename);
            }//end for i
            //now that they are initialized, start them and check their status
            boolean alldone = false;
            while (alldone == false) {
                alldone = true;
                synchronized (syncdone) {
                    synchronized (queryvec) {
                        for (int i = cpu; --i >= 0;) {
                            if (mythreads[i].done == false) {
                                alldone = false;
                            } else if (queryvec.size() > 0) {
                                mythreads[i] = new blastthread(i, blastblocksize, dbfilename);
                                mythreads[i].start();
                            }
                        }//end for
                        if (queryvec.size() > 0) {
                            alldone = false;
                        }
                    }//end sync queryvec
                    if (alldone == false) {
                        try {
                            syncdone.wait();
                        } catch (InterruptedException ie) {
                            System.err.println("ERROR on wait in searchblastv2 main get method");
                        }
                    }
                }//end sync syncdone
            }//end while alldone==false
            outwrite.close();
        }//end else queryvec.size==0
        //now I have written all results to the temporary file.
        //next, read them back from the file and send them to CLANS
        System.out.println("reading all results from temporary file '" + saveblastname + "'");
        HashMap<String, MinimalHsp> allhash = new HashMap<String, MinimalHsp>();
        try {
            BufferedReader fileread = new BufferedReader(new FileReader(saveblastname));
            String inline = "", idline = "", valline = "", tmpstr = "";
            String[] tmparr;
            HashMap<String, MinimalHsp> roundhash = new HashMap<String, MinimalHsp>();
            //read until I hit the "#done for sequence" lines and then pass the current elements to the final hash containing the data
            MinimalHsp myhsp;
            boolean passrounds = false;
            HashMap<Integer, Integer> donehash = new HashMap<Integer, Integer>();
            int count = 0;
            try {
                while ((inline = fileread.readLine()) != null) {
                    if (inline.startsWith("hsp: ")) {
                        tmparr = idline.split(";");
                        if (tmparr.length >= 2) {
                            tmpstr = tmparr[0] + ";" + tmparr[1];
                            if (roundhash.containsKey(tmpstr)) {
                                myhsp = roundhash.get(tmpstr);
                                tmpstr = valline;
                                myhsp.addpval(Double.parseDouble(tmpstr));
                            } else {
                                myhsp = new MinimalHsp();
                                tmpstr = tmparr[0];
                                myhsp.query = Integer.parseInt(tmpstr);
                                tmpstr = tmparr[1];
                                myhsp.hit = Integer.parseInt(tmpstr);
                                tmpstr = valline;
                                myhsp.addpval(Double.parseDouble(tmpstr));
                                roundhash.put(myhsp.query + ";" + myhsp.hit, myhsp);
                            }
                        } else {
                            System.out.println("WARNING: hsp found with less than two elements on '" + idline + "'");
                        }
                        idline = inline.substring(5);
                        inline = fileread.readLine();
                        valline = inline.substring(inline.indexOf(":") + 1);
                    } else if (inline.startsWith("SECTION:")) {
                        if (passrounds == true) {
                            System.out.print(".");
                            if (count % 100 == 0) {
                                System.out.println(count);
                            }
                            count++;
                            //then move the data from the roundhash to the allhash for those sequences in donehash
                            MinimalHsp[] checkarr = roundhash.values().toArray(new MinimalHsp[0]);
                            for (int i = checkarr.length; --i >= 0;) {
                                myhsp = checkarr[i];
                                if (donehash.containsKey(new Integer(myhsp.query))) {
                                    allhash.put(myhsp.query + ";" + myhsp.hit, myhsp);
                                }
                            }//end for i
                            roundhash.clear();
                            donehash.clear();
                            passrounds = false;
                        }
                    } else if (inline.startsWith("#")) {
                        passrounds = true;
                        //add these elements from roundshash to allhash and overwrite previous ones if necessary
                        tmpstr = inline.substring(inline.indexOf(":") + 1);
                        donehash.put(nameshash.get(tmpstr), null);
                    } else {
                        if (inline.length() > 0) {
                            System.err.println("WARNING: unknown line '" + inline + "' found in temp file; skipping entry!");
                        }
                    }
                }//end while reading
                if (passrounds == true) {
                    System.out.print(".");
                    if (count % 100 == 0) {
                        System.out.println(count);
                    }
                    count++;
                    //then move the data from the roundhash to the allhash for those sequences in donehash
                    MinimalHsp[] checkarr = roundhash.values().toArray(new MinimalHsp[0]);
                    for (int i = checkarr.length; --i >= 0;) {
                        myhsp = checkarr[i];
                        if (donehash.containsKey(new Integer(myhsp.query))) {
                            allhash.put(myhsp.query + ";" + myhsp.hit, myhsp);
                        }
                    }//end for i
                    roundhash.clear();
                    donehash.clear();
                }
            } catch (NumberFormatException ne) {
                System.err.println("ERROR trying to parse number from '" + tmpstr + "'");
            }
            fileread.close();
        } catch (IOException ioe) {
            System.err.println("IOERROR reading from '" + saveblastname + "'");
            return null;
        }
        //now clean up the temporary files I used
        try {
            //delete the blast database and the files that formatdb created and the blast query files
            //files to check are: blastdbname; blastdbname.psd; blastdbname.phr; blastdbname.pin; blastdbname.pnd
            //blastdbname.pni; blastdbname.psi; blastdbname.psq
            File rmfile;
            for (int i = 0; i < cpu; i++) {
                rmfile = new File("tmp" + i + ".query");
                if (rmfile.exists()) {
                    rmfile.delete();
                }
                rmfile = new File("tmp" + i + ".chkpnt");
                if (rmfile.exists()) {
                    rmfile.delete();
                }
            }// end for i<cpu
            rmfile = new File(dbfilename);
            if (rmfile.exists()) {
                rmfile.delete();
            }
            rmfile = new File(dbfilename + ".psd");
            if (rmfile.exists()) {
                rmfile.delete();
            }
            rmfile = new File(dbfilename + ".phr");
            if (rmfile.exists()) {
                rmfile.delete();
            }
            rmfile = new File(dbfilename + ".pin");
            if (rmfile.exists()) {
                rmfile.delete();
            }
            rmfile = new File(dbfilename + ".pnd");
            if (rmfile.exists()) {
                rmfile.delete();
            }
            rmfile = new File(dbfilename + ".pni");
            if (rmfile.exists()) {
                rmfile.delete();
            }
            rmfile = new File(dbfilename + ".psi");
            if (rmfile.exists()) {
                rmfile.delete();
            }
            rmfile = new File(dbfilename + ".psq");
            if (rmfile.exists()) {
                rmfile.delete();
            }
            rmfile = new File("formatdb.log");
            if (rmfile.exists()) {
                rmfile.delete();
            }
            //don't delete the tmpblasthsp file by default, as the program micht also crash before the first saving of the data!!!
            //rmfile=new File("tmpblasthsp.txt");
            //if(rmfile.exists()){
            //    rmfile.delete();
            //}
        } catch (RuntimeException e) {
            System.err.println("unable to delete all temporary files for CLANS");
            e.printStackTrace();
        }
        //now see what values i read

        return allhash.values().toArray(new MinimalHsp[0]);
    }//end get

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //I need a thread that does the blast searches
    class blastthread extends java.lang.Thread {

        public blastthread(int threadnum, int blastblocksize, String dbfilename) {
            this.threadnum = threadnum;
            this.blastblock = blastblocksize;
            this.dbfilename = dbfilename;
            rt = Runtime.getRuntime();
        }// end blastthread
        Runtime rt;
        int threadnum;
        int blastblock = -1;
        String dbfilename;
        boolean done = true;
        boolean stopthread = false;
        StringBuffer errbuff = new StringBuffer();
        //--------------------------------------------------------------------------

        @Override
        public void run() {
            //check to see if there is sth assigned to the query
            //now get all of the queries you want to run the blast for in this block
            AminoAcidSequence[] myqueries;
            synchronized (syncdone) {
                done = false;
                synchronized (queryvec) {
                    if (queryvec.size() == 0) {
                        System.out.println("no more queries to do, skipping thread " + threadnum);
                        done = true;//I am already in a sync on syncdone
                        syncdone.notify();
                        return;
                    }
                    if (queryvec.size() > blastblock) {
                        myqueries = new AminoAcidSequence[blastblock];
                        for (int i = blastblock; --i >= 0;) {
                            myqueries[i] = queryvec.remove(0);
                        }//end for i
                    } else {
                        myqueries = new AminoAcidSequence[queryvec.size()];
                        for (int i = queryvec.size(); --i >= 0;) {
                            myqueries[i] = queryvec.remove(0);
                        }//end for i
                    }
                    System.out.println("start for thread " + threadnum + ": " + queryvec.size() + " sequences left to do");
                }//end synchronized queryvec
            }
            boolean allgood = true;
            stopthread = false;
            //String blastcommand="";
            String tmpfilestring = "tmp" + threadnum + ".query";
            //String tmpcheckfile = "tmp" + threadnum + ".chkpnt";
            String tmpoutfile = "tmp" + threadnum + ".outfile";
            //String dbstring;
            String[] cmdarr;
            //convert the String commands to a array of Strings (resolves some runtime.exec problems)
            Vector<String> tmpvec = new Vector<String>();
            if (isblastplus) {
                cmdarr = blastpath.split("\\s+");
                tmpvec.addElement("-outfmt");
                tmpvec.addElement("5");//xml output
                //tmpvec.addElement("-show_gis");//don't need that one
                tmpvec.addElement("-evalue");
                tmpvec.addElement(String.valueOf(eval));
                tmpvec.addElement("-query");
                tmpvec.addElement(new File(tmpfilestring).getAbsolutePath());
                tmpvec.addElement("-out");
                tmpvec.addElement(new File(tmpoutfile).getAbsolutePath());
                tmpvec.addElement("-db");
                tmpvec.addElement(dbfilename);
                for (int i = cmdarr.length; --i >= 0;) {
                    tmpvec.add(0, cmdarr[i]);
                }
                if (cmd.length() > 0) {
                    tmpvec.add(0, cmd);
                }
            } else {
                cmdarr = blastpath.split("\\s+");
                //tmpvec.addElement("-I");//don't need to show gi's as I have none!
                //tmpvec.addElement("T");
                tmpvec.addElement("-m");
                tmpvec.addElement("7");//xml output
                tmpvec.addElement("-e");
                tmpvec.addElement(String.valueOf(eval));
                tmpvec.addElement("-i");
                tmpvec.addElement(new File(tmpfilestring).getAbsolutePath());
                tmpvec.addElement("-o");
                tmpvec.addElement(new File(tmpoutfile).getAbsolutePath());
                tmpvec.addElement("-d");
                tmpvec.addElement(dbfilename);
                for (int i = cmdarr.length; i > 0; i--) {
                    tmpvec.add(0, cmdarr[i - 1]);
                }
                if (cmd.length() > 0) {
                    //blastcommand=cmd+" "+blastcommand;
                    tmpvec.add(0, cmd);
                }
            }
            cmdarr = new String[tmpvec.size()];
            tmpvec.copyInto(cmdarr);
            StringBuffer inout = new StringBuffer();
            StringBuffer errout = new StringBuffer();
            if (verbose > 2) {
                System.out.print(threadnum + " trying:");
                for (int i = 0; i < cmdarr.length; i++) {
                    System.out.print("'" + cmdarr[i] + "' ");
                }//end for i
                System.out.println();
            }
            //String myblast = "";
            try {
                //write the query to a file instead of using stdin; !stdin gives problems!
                PrintWriter outwrite = new PrintWriter(new BufferedWriter(new FileWriter(tmpfilestring)));
                synchronized (queryvec) {
                    AminoAcidSequence currquery;
                    for (int i = myqueries.length; --i >= 0;) {
                        currquery = myqueries[i];
                        outwrite.println(">" + currquery.name);
                        outwrite.println(currquery.seq.replaceAll("-", ""));
                    }//end for i
                }//end synchronized queryvec
                outwrite.close();
            } catch (IOException ioe) {
                System.out.println("IOERROR writing to query file " + tmpfilestring);
                errbuff.append("IOEROR writing to query file " + tmpfilestring);
                stopthread = true;
                allgood = false;
            }
            if (stopthread == false) {
                try {
                    BufferedReader perr;
                    BufferedReader pin;
                    //PrintWriter pout;
                    threadstreamreader perrread;
                    threadstreamreader pinread;
                    Process p = rt.exec(cmdarr);
                    perr = new BufferedReader(new InputStreamReader(p.getErrorStream()));//process error stream
                    pin = new BufferedReader(new InputStreamReader(p.getInputStream()));//output of process
                    //note: threadstreamreader closes the input bufferedreader
                    //threadstreamreaders are used because otherwise deadlock can occur if the output of the process is larger than the
                    //buffersize java allocates to the bufferedreaer.
                    perrread = new threadstreamreader(perr, errout);
                    pinread = new threadstreamreader(pin, inout);
                    try {
                        perrread.start();
                        pinread.start();
                        if (verbose > 2) {
                            System.out.print("Done setting sequence data for" + threadnum + " waiting for command:");
                            for (int i = 0; i < cmdarr.length; i++) {
                                System.out.print("'" + cmdarr[i] + "' ");
                            }//end for i
                            System.out.println(" to finish");
                        }
                        p.waitFor();
                        if (p.exitValue() != 0) {
                            System.err.println("non-perfect exit from blast for thread" + threadnum);
                            errbuff.append("non-perfect exit from blast for thread " + threadnum + "\n");
                            allgood = false;
                        }
                        if (verbose > 0) {
                            System.out.println(threadnum + ": blast done for " + blastblock + "queries:");
                            for (int i = myqueries.length; --i >= 0;) {
                                System.out.println("\t" + myqueries[i].name);
                            }
                        }
                        synchronized (inout) {
                            while (pinread.done == false) {
                                //System.out.println("waiting for process reader to finish");
                                try {
                                    //sleep(10l);
                                    inout.wait(10l);
                                } catch (InterruptedException e) {
                                    System.err.println("interrupted sleep in blastthread");
                                    e.printStackTrace();
                                }
                            }
                        }//end synchronized inout
                    } catch (InterruptedException e) {
                        System.out.print("Interrupted process for thread " + threadnum + " command:");
                        errbuff.append("Interrupted process for thread " + threadnum + " command:");
                        for (int i = 0; i < cmdarr.length; i++) {
                            System.out.print("'" + cmdarr[i] + "' ");
                            errbuff.append("'" + cmdarr[i] + "' ");
                        }//end for i
                        System.out.println();
                        stopthread = true;
                        allgood = false;
                    }
                    //myblast = inout.toString();
                    perrread = null;//.clear();
                    pinread = null;//.clear();
                } catch (IOException ioe) {
                    System.out.print("IOERROR for running thread " + threadnum + " command:");
                    errbuff.append("IOERROR running thread" + threadnum + " command\n");
                    for (int i = 0; i < cmdarr.length; i++) {
                        System.out.print("'" + cmdarr[i] + "' ");
                        errbuff.append("'" + cmdarr[i] + "'");
                    }//end for i
                    System.out.println();
                    stopthread = true;
                    allgood = false;
                }
            }//end if stopthread==false
            if (stopthread) {
                System.err.println("stopped thread " + threadnum);
                synchronized (syncdone) {
                    this.done = true;
                    syncdone.notify();
                }// end synchronized parent
                return;
            }// end if stopthread
            if (errout.indexOf("Error") > -1) {
                System.err.print(errout.toString());
                errbuff.append("ERRORS:\n" + errout.toString() + "\n");
            }
            if (verbose > 10) {
                System.out.println(inout.toString());
            }
            if (verbose > 2) {
                System.out.println("getting hsp's from thread " + threadnum);
            }
            inout = null;//free the memory from these stringbuffers (more or less)
            errout = null;
            ArrayList<hsp> hsplist = hspgetv2.get(tmpoutfile, eval, pval, coverage, scval, ident, verbose);
            if (hsplist == null) {
                System.err.println("ERROR trying to parse hsp's from blast output, parsing returned \"null\"");
                allgood = false;
            } else {
                //System.out.println("hsp's returned: "+hsplist.size());
                //now I have all the hsp's in an arrayList.
                //next, convert all of the hsp's to minhsp objects
                hsp myhsp;
                int qnum, hnum;
                double val;
                HashMap<String, MinimalHsp> hsphash = new HashMap<String, MinimalHsp>();
                for (int i = hsplist.size(); --i >= 0;) {
                    myhsp = hsplist.get(i);
                    if (nameshash.containsKey(myhsp.qname)) {
                        qnum = nameshash.get(myhsp.qname);
                    } else {
                        System.err.println("WARNING: unrecognized sequence name '" + myhsp.qname + "'");
                        allgood = false;
                        qnum = -1;
                    }
                    if (nameshash.containsKey(myhsp.hname)) {
                        hnum = nameshash.get(myhsp.hname);
                    } else {
                        System.err.println("WARNING: unrecognized sequence name '" + myhsp.hname + "'");
                        allgood = false;
                        hnum = -1;
                    }
                    val = myhsp.value;
                    if (qnum != hnum) {//I don't want the values of one sequence to itself
                        //System.out.println("doing hsp:"+qnum+";"+hnum);
                        if (hsphash.containsKey(qnum + ";" + hnum)) {
                            hsphash.get(qnum + ";" + hnum).addpval(val);
                        } else {
                            hsphash.put(qnum + ";" + hnum, new MinimalHsp(qnum, hnum, val));
                        }
                    }
                }//end for i
                //now I should have converted all of my data to minhsp's containing query, hit and the (possible multiple) hits they had
                MinimalHsp[] mymins = hsphash.values().toArray(new MinimalHsp[0]);
                //next synchronize on file output and write the results to file
                if (allgood) {
                    synchronized (filesync) {
                        outwrite.println("SECTION:");
                        MinimalHsp mymin;
                        for (int i = mymins.length; --i >= 0;) {
                            mymin = mymins[i];
                            for (int j = mymin.val.length; --j >= 0;) {
                                outwrite.println("hsp: " + mymin.query + ";" + mymin.hit + ";" + j);
                                outwrite.println("\tvalue: " + mymin.val[j]);
                                //System.out.println("hsp: "+mymin.query+";"+mymin.hit+";"+j);
                                //System.out.println("\tvalue: "+mymin.val[j]);

                            }//end for j
                        }//end for i
                        //now print out the sequence id's that were done
                        for (int i = myqueries.length; --i >= 0;) {
                            outwrite.println("# done for:" + myqueries[i].name);
                        }//end for i
                        outwrite.flush();
                    }//end synchronized filesync
                }else{
                	System.err.println("ERROR, NOT ALL GOOD in reading from BLAST");                	
                }
                //all done, notify parent to get the data
                if (verbose > 2) {
                    System.out.println("All done, notifying parent for thread" + threadnum);
                }
            }
            synchronized (syncdone) {
                this.done = true;
                syncdone.notify();
            }// end synchronized parent
        }// end run
    }
    //end class blastthread
}//end searchblastv2 class

