/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package clans;

import java.io.*;
import java.util.*;

/**
 *
 * @author  tancred
 */
public class hspget {

    /** Creates a new instance of hspget */
    public hspget() {
    }

    public static Vector get(boolean lowmem, String filename, Vector invec, double eval, double pval, float coverage, float scval, float ident, int verbose, HashMap nameshash, boolean useallrounds) {
        //simple fix, read blast file to ram
        System.out.println("lowmem, reading blast from " + filename);
        StringBuffer tmpstr = new StringBuffer();
        String inline;
        try {
            BufferedReader inread = new BufferedReader(new FileReader(filename));
            while ((inline = inread.readLine()) != null) {
                tmpstr.append(inline + "\n");
            }
            inread.close();
        } catch (IOException e) {
            System.err.println("IOError reading from " + filename);
        }
        return get(tmpstr.toString(), new Vector(), eval, pval, coverage, scval, ident, verbose, nameshash, useallrounds);
    }//end get for lowmem

    public static Vector get(String inreadstring, Vector invec, double eval, double pval, float coverage, float scval, float ident, int verbose, HashMap nameshash, boolean useallrounds) {
        //this should read the data from inread (a bufferedReader with blast output in html format)
        //to be able to filter by pval I need to search through the inread BufferedReader for the database size
        int seqnum = nameshash.size();
        double pvalfactor = getpvalfactor(inreadstring, seqnum);//get the factor i need to multiply the evalue by to get the pvalue
        //System.out.println("pvalfactor="+pvalfactor);
        if (pvalfactor == -1) {
            if (verbose > 0) {
                System.out.println("No hits found");
            }
            return invec;
        } else if (pvalfactor == 1) {
            String tmpfilename = System.currentTimeMillis() + "_blasterror.txt";
            System.err.println("Saving output to file " + tmpfilename);
            //System.err.println("tmpstr="+inreadstring);
            try {
                PrintWriter outwrite = new PrintWriter(new BufferedWriter(new FileWriter(tmpfilename)));
                outwrite.println(inreadstring);
                outwrite.close();
            } catch (IOException er) {
                System.err.println("IOError writing to " + tmpfilename);
            }
        }
        BufferedReader inread = new BufferedReader(new StringReader(inreadstring));
        boolean psiblastresult = false;
        String inline = "";
        String psiround = "";//how many psiblast rounds were done
        if (useallrounds == false) {
            try {
                //see if this file has a "Results from round" line before Sequences producing significant alignments:
                boolean readon = true;
                while (((inline = inread.readLine()) != null) && (readon == true)) {//while I can read
                    if (inline.startsWith("Results from round")) {
                        psiblastresult = true;
                        psiround = inline.trim();
                    } else if (inline.startsWith("Sequences producing significant") && (psiblastresult == false)) {
                        readon = false;
                    }
                }//end while reading
                inread.close();
            } catch (IOException e) {
                System.err.println("IOERROR reading from blast results: " + inreadstring);
            }
        }//end if useallrounds==false
        //now I should know wether I have a blast or psiblast output and if psiblast, then how many rounds were done.
        inread = new BufferedReader(new StringReader(inreadstring));
        //filter all hsp's by pvalue,coverage,scval and ident and return the valid hits
        boolean querydone = false;
        boolean readquery = false;
        boolean readname = false;
        boolean readseq = false;
        String queryname = "";
        String currname = "";
        String currquery = "";
        String currseq = "";
        String tmpstr = "";
        int querylength = 0;
        int subjectlength = 0;
        float myscval = -1;
        float myident = -1;
        double myeval = -1;
        double mypval = -1;
        double dblength = -1;
        double dblengtheff = -1;
        int qlength = -1;
        int qlengtheff = -1;
        hsp myhsp = new hsp();
        int index1;
        boolean foundlastround = false;
        boolean loopread = false;
        boolean readon = true;
        try {
            while (readon) {
                if (loopread) {//don't read next line
                    loopread = false;
                } else if ((inline = inread.readLine()) != null) {
                    //then I have update inline
                } else {
                    readon = false;
                    break;
                }
                inline = inline.trim();
                if (querydone == false) {
                    if (inline.indexOf("<b>Query=</b>") > -1) {
                        readquery = true;
                    }
                    int endql;
                    if ((endql = inline.indexOf("letters)")) > -1) {//old blast
                        int startql;
                        if ((startql = inline.indexOf("(")) > -1) {
                            try {
                                tmpstr = (inline.substring(startql + 1, endql)).trim();
                                tmpstr = tmpstr.replaceAll("\\D", "");//remove all non-numbers
                                querylength = Integer.parseInt(tmpstr);
                            } catch (NumberFormatException e) {
                                System.err.println("error in parsing query length:" + inline.substring(startql + 1, endql));
                            }
                        }// end if startql
                        querydone = true;
                        queryname = extractqueryname(queryname, nameshash);
                        continue;
                    }//end if endql
                    if ((endql = inline.indexOf("Length=")) > -1) {//blast+
                        int startql;
                        if ((startql = inline.indexOf("=", endql)) > -1) {
                            try {
                                tmpstr = (inline.substring(startql + 1)).trim();
                                tmpstr = tmpstr.replaceAll("\\D", "");//remove all non-numbers
                                querylength = Integer.parseInt(tmpstr);
                            } catch (NumberFormatException e) {
                                System.err.println("error in parsing query length:" + inline.substring(startql + 1, endql));
                            }
                        }// end if startql
                        querydone = true;
                        queryname = extractqueryname(queryname, nameshash);
                        continue;
                    }//end if endql
                    if (readquery) {
                        queryname = queryname + inline;
                        continue;
                    }
                } else {// end if querydone==false
                    if ((psiblastresult) && (foundlastround == false)) {//if I am reading a psiblast output
                        //skip to the last round
                        if (inline.startsWith(psiround)) {
                            foundlastround = true;
                        }
                        while (((inline = inread.readLine()) != null) && (foundlastround == false)) {
                            if (inline.startsWith(psiround)) {
                                foundlastround = true;
                            }
                        }
                        if (foundlastround == false) {
                            System.err.println("unable to find last " + psiround + " for:");
                            System.err.println(inreadstring.toString());
                            System.exit(0);
                        }
                        continue;
                    }//end if psiblast results
                    if (readseq == false) {
                        if ((inline.indexOf("<a name = ") > -1) || (inline.indexOf("<pre>&gt;<a name=") > -1)) {//if this substring exists; first for blast, second for psyblast
                            //try if((inline.indexOf("><a name = ")> -1)||(inline.indexOf("<pre>&gt;<a name=")>-1)){//if this substring exists; first for blast, second for psyblast
                            //start storing the name and prepare for the new block
                            currname = "";
                            readname = true;
                        } else if (inline.indexOf("<a name=") > -1) {//for the new blast+ if this substring exists; first for blast, second for psyblast
                            //try if((inline.indexOf("><a name = ")> -1)||(inline.indexOf("<pre>&gt;<a name=")>-1)){//if this substring exists; first for blast, second for psyblast
                            //start storing the name and prepare for the new block
                            currname = "";
                            readname = true;
                        }//end if ><a name =
                        if (readname == true) {
                            int strs;
                            if ((strs = inline.indexOf("Length =")) > -1 || (strs = inline.indexOf("Length=")) > -1) {//old blast & new blast+
                                try {
                                    //blast has the disturbing tendency to add a comma in numbers above 1000. remove that before parsing
                                    tmpstr = (inline.substring(strs + 8)).trim();
                                    tmpstr = tmpstr.replaceAll("\\D", "");
                                    subjectlength = Integer.parseInt(tmpstr);
                                } catch (NumberFormatException e) {
                                    System.err.println("unable to parse int from " + (inline.substring(strs + 8)).trim() + " for " + currname);
                                    subjectlength = 1;
                                }
                                readname = false;
                                currname = extractname(currname, nameshash);
                                continue;
                            }
                            //if this currstr has no "Length =" just add it to the current namestr.
                            currname = currname + inline;
                            continue;
                        }// end if readname
                        int strs;
                        if ((strs = inline.indexOf("Score =")) > -1) {
                            int endstrs;
                            if ((endstrs = inline.indexOf("bits", strs + 7)) > -1) {
                                String scorestring = (inline.substring(strs + 7, endstrs)).trim();
                                try {
                                    myscval = Float.parseFloat(scorestring);
                                } catch (NumberFormatException e) {
                                    System.err.println("Error parsing Score " + scorestring + " for " + currname);
                                    myscval = -1;
                                }
                            }
                        }// end if Score
                        //if((index1=inline.indexOf("Expect ="))>-1){//NOTE; this does not work for tblastx where Expect(\d) = may ocurr
                        if ((index1 = inline.indexOf("Expect")) > -1) {
                            //now I need to differentiate between blast and psiblast which differ in their outputs
                            if (inline.indexOf("Method:") > -1) {//in this case i'm using psiblast
                                String tmpeval = inline.substring(index1 + 8, inline.indexOf(",", index1)).trim();
                                if (tmpeval.charAt(0) == 'e') {
                                    tmpeval = "1" + tmpeval;
                                }
                                try {
                                    myeval = Double.parseDouble(tmpeval);
                                    mypval = myeval * pvalfactor;
                                } catch (NumberFormatException e) {
                                    System.err.println("Error parsing E-value " + tmpeval + " for " + currname);
                                    myeval = -1;
                                    mypval = -1;
                                }
                            } else {//i'm using blast
                                String[] tmparr = inline.split("\\s", 0);
                                if (tmparr[java.lang.reflect.Array.getLength(tmparr) - 1].charAt(0) == 'e') {
                                    tmparr[java.lang.reflect.Array.getLength(tmparr) - 1] = "1" + tmparr[java.lang.reflect.Array.getLength(tmparr) - 1];
                                }
                                try {
                                    myeval = Double.parseDouble(tmparr[java.lang.reflect.Array.getLength(tmparr) - 1]);
                                    mypval = myeval * pvalfactor;
                                //mypval=1-java.lang.Math.exp(-myeval);//NCBI formula doesn't work as well (may be correct but works less well)
                                } catch (NumberFormatException e) {
                                    System.err.println("Error parsing E-value " + tmparr[java.lang.reflect.Array.getLength(tmparr) - 1] + " for " + currname);
                                    myeval = -1;
                                    mypval = -1;
                                }
                            }
                        }// end if expect
                        if (inline.indexOf("Identities =") > -1) {
                            //the identities in percent are stored in the first () after IDENTITES=
                            int start;
                            int end;
                            if ((start = inline.indexOf("(")) > -1) {//start of first bracket
                                if ((end = inline.indexOf("%")) > -1) {//end of percent number
                                    try {
                                        myident = Integer.parseInt(inline.substring(start + 1, end));
                                    } catch (NumberFormatException e) {
                                        System.err.println("unable to parse int from " + inline.substring(start + 1, end) + " for ident");
                                        e.printStackTrace();
                                    }
                                }
                            }
                            readseq = true;
                            myident = myident / 100;
                            continue;
                        }// end if identities
                        if (inline.indexOf("Frame = ") > -1) {
                            continue;
                        }
                    } else { //if(readseq){
                        if (inline.indexOf("Query:") > -1 || inline.indexOf("Query  ") > -1) {
                            currquery = currquery + " " + inline;
                            continue;
                        }
                        if (inline.indexOf("Sbjct:") > -1 || inline.indexOf("Sbjct  ") > -1) {
                            currseq = currseq + " " + inline;
                            continue;
                        }
                        if ((inline.indexOf("</PRE>") > -1) || (inline.indexOf("</pre>") > -1)) {// if </PRE> is on the line /PRE for blast /pre for psyblast
                            //if I hit the end of a block put the current data into a vector
                            myhsp = new hsp();
                            myhsp.qname = queryname;
                            myhsp.hname = currname;
                            myhsp.qseq = currquery;
                            myhsp.hseq = currseq;
                            if (scval >= 0) {//see what the maximum length of a match could have been
                                //System.out.println("scval: qlength="+querylength+" slength="+subjectlength);
                                if (subjectlength > querylength) {
                                    //System.out.println("\tQlengthdiv");
                                    myhsp.value = -myscval / (float) querylength;
                                } else {
                                    //System.out.println("\t***Slengthdiv");
                                    myhsp.value = -myscval / (float) subjectlength;
                                }
                            } else {
                                myhsp.value = mypval;
                            }
                            readseq = false;
                            invec.addElement(myhsp);
                            currquery = "";
                            currseq = "";
                            myeval = -1;
                            mypval = -1;
                            myscval = -1;
                            myident = -1;
                            continue;
                        } else if ((inline.indexOf("><a name=") > -1)) {// if reading new sequence data in blast+ executable
                            //if I hit the end of a block put the current data into a vector
                            myhsp = new hsp();
                            myhsp.qname = queryname;
                            myhsp.hname = currname;
                            myhsp.qseq = currquery;
                            myhsp.hseq = currseq;
                            if (scval >= 0) {//see what the maximum length of a match could have been
                                //System.out.println("scval: qlength="+querylength+" slength="+subjectlength);
                                if (subjectlength > querylength) {
                                    //System.out.println("\tQlengthdiv");
                                    myhsp.value = -myscval / (float) querylength;
                                } else {
                                    //System.out.println("\t***Slengthdiv");
                                    myhsp.value = -myscval / (float) subjectlength;
                                }
                            } else {
                                myhsp.value = mypval;
                            }
                            readseq = false;
                            invec.addElement(myhsp);
                            currquery = "";
                            currseq = "";
                            myeval = -1;
                            mypval = -1;
                            myscval = -1;
                            myident = -1;
                            loopread = true;//do NOT read the next line, instead re-analyze this one
                            continue;
                        }else if(inline.indexOf("Score =")>-1){//blast+ multiple hit
                            //if I hit the end of a block put the current data into a vector
                            myhsp = new hsp();
                            myhsp.qname = queryname;
                            myhsp.hname = currname;
                            myhsp.qseq = currquery;
                            myhsp.hseq = currseq;
                            if (scval >= 0) {//see what the maximum length of a match could have been
                                //System.out.println("scval: qlength="+querylength+" slength="+subjectlength);
                                if (subjectlength > querylength) {
                                    //System.out.println("\tQlengthdiv");
                                    myhsp.value = -myscval / (float) querylength;
                                } else {
                                    //System.out.println("\t***Slengthdiv");
                                    myhsp.value = -myscval / (float) subjectlength;
                                }
                            } else {
                                myhsp.value = mypval;
                            }
                            readseq = false;
                            invec.addElement(myhsp);
                            currquery = "";
                            currseq = "";
                            myeval = -1;
                            mypval = -1;
                            myscval = -1;
                            myident = -1;
                            loopread = true;//do NOT read the next line, instead re-analyze this one
                            continue;
                        }// end if blast+ new sequence
                    }// end if readseq
                    //if i am at the end of the file reading data about the database size etc.
                    if (inline.indexOf("effective length of query:") > -1) {//if(inline.matches(".*effective length of query:.*")){
                        String[] tmparr = inline.split(":\\s*");
                        try {
                            tmparr[1] = tmparr[1].replaceAll(",", "");//get rid of the blast 1000 tics
                            qlengtheff = Integer.parseInt(tmparr[1]);
                        } catch (NumberFormatException e) {
                            System.err.println("unable to convert " + tmparr[1] + " to integer in qlengtheff");
                            qlengtheff = -1;
                        }
                    } else if (inline.indexOf("length of query:") > -1) {//inline.matches(".*length of query:.*")){
                        String[] tmparr = inline.split(":\\s*");
                        try {
                            tmparr[1] = tmparr[1].replaceAll(",", "");//get rid of the blast 1000 tics
                            qlength = Integer.parseInt(tmparr[1]);
                        } catch (NumberFormatException e) {
                            System.err.println("unable to convert " + tmparr[1] + " to integer in qlength");
                            qlength = -1;
                        }
                    } else if (inline.indexOf("effective length of database:") > -1) {//inline.matches(".*effective length of database:.*")){
                        String[] tmparr = inline.split(":\\s*");
                        try {
                            tmparr[1] = tmparr[1].replaceAll(",", "");//get rid of the blast 1000 tics
                            dblengtheff = Double.parseDouble(tmparr[1]);
                        } catch (NumberFormatException e) {
                            System.err.println("unable to convert " + tmparr[1] + " to double in dblengtheff");
                            dblengtheff = -1;
                        }
                    } else if (inline.indexOf("length of database:") > -1) {//inline.matches(".*length of database:.*")){
                        String[] tmparr = inline.split(":\\s*");
                        try {
                            tmparr[1] = tmparr[1].replaceAll(",", "");//get rid of the blast 1000 tics
                            dblength = Double.parseDouble(tmparr[1]);
                        } catch (NumberFormatException e) {
                            System.err.println("unable to convert " + tmparr[1] + " to double in dblength");
                            dblength = -1;
                        }
                    }
                }//end if querydone==true
            }// end while readline
        } catch (IOException e) {
            System.err.println("IOError in reading blast output");
            e.printStackTrace();
        }
        //now i have all the hsp's in invec
        //for all check if they comply to the coverage, scval,ident,pvalue cutoff
        //now loop through the Vector elements and splice query and subject out of the seqs
        //and get data for query and hit start and end
        for (int i = 0; i < invec.size(); i++) {
            invec.setElementAt(makefinal((hsp) (invec.elementAt(i))), i);
            if (verbose > 3) {
                System.out.println("hit to:'"+((hsp) invec.elementAt(i)).hname+"'");
            }
        }// end for i
        for (int i = 0; i < invec.size(); i++) {
            myhsp = (hsp) invec.elementAt(i);
            //myhsp.coverage=(myhsp.hend-myhsp.hstart)/querylength;
            //if this object passes the selection procedure, keep it.
            if (scval >= 0) {//then I want to use the score/column as attraction value (normalized to 0-1)
                if (myhsp.value <= -scval) {
                    invec.setElementAt(myhsp, i);
                } else {
                    invec.removeElementAt(i);
                    i--;
                }
            } else {
                if ((myhsp.value < pval)) {//&&(myhsp.coverage>=coverage)&&(myhsp.scval>=scval)&&(myhsp.ident>=ident)){
                    invec.setElementAt(myhsp, i);
                } else {
                    invec.removeElementAt(i);
                    i--;
                }
            }
        }// end for i
        if (((qlength == -1) || (qlengtheff == -1) || (dblength == -1) || (dblengtheff == -1)) && (inreadstring.toString().matches(".*\\*+\\s+No hits found\\s+\\*+.*"))) {
            System.err.println("not all data read");
            if (qlength == -1) {
                System.err.println("no qlength");
            }
            if (qlengtheff == -1) {
                System.err.println("no qlengtheff");
            }
            if (dblength == -1) {
                System.err.println("no dblength");
            }
            if (dblengtheff == -1) {
                System.err.println("no dblengtheff");
            }
            System.err.println(inreadstring);
        }
        return invec;
    }// end get

    //--------------------------------------------------------------------------
    static hsp makefinal(hsp inhsp) {
        //clean up qseq and hseq and get start and end for each seq.
        int qstart = -1;
        int hstart = -1;
        int qend = -1;
        int hend = -1;
        int currpos = 0;
        int endnum = 0;
        String currseq = "";
        //do the query
        String[] splitarr = inhsp.qseq.split(" ", 0);//split at spaces
        int elems = java.lang.reflect.Array.getLength(splitarr);
        int counter = 0;
        for (int i = 0; i < elems; i++) {
            if (splitarr[i].length() == 0) {//skip empty strings
                continue;
            }
            if (counter == 0) {//if I expect a Query:
                counter += 1;
                continue;
            } else if (counter == 1) {//if I expect a start number
                if (qstart == -1) {//if I have no qstart yet
                    try {
                        qstart = Integer.parseInt(splitarr[i]);
                    } catch (NumberFormatException e) {
                        System.err.println("unable to parse int from " + splitarr[i]);
                    }
                }
                counter += 1;
                continue;
            } else if (counter == 2) {//if this should be sequence data
                currseq += splitarr[i];
                counter += 1;
                continue;
            } else if (counter == 3) {//if this is end data
                try {
                    qend = Integer.parseInt(splitarr[i]);
                } catch (NumberFormatException e) {
                    System.err.println("unable to parse int from " + splitarr[i]);
                }
                counter = 0;
                continue;
            }
        }// end for i
        inhsp.qseq = currseq.toUpperCase();
        inhsp.qstart = qstart;
        inhsp.qend = qend;
        currseq = "";
        splitarr = inhsp.hseq.split(" ", 0);//split at spaces
        elems = java.lang.reflect.Array.getLength(splitarr);
        counter = 0;
        for (int i = 0; i < elems; i++) {
            if (splitarr[i].length() == 0) {//skip empty strings
                continue;
            }
            if (counter == 0) {//if I expect a Query:
                counter += 1;
                continue;
            } else if (counter == 1) {//if I expect a start number
                if (hstart == -1) {//if I have no qstart yet
                    try {
                        hstart = Integer.parseInt(splitarr[i]);
                    } catch (NumberFormatException e) {
                        System.err.println("unable to parse int from " + splitarr[i]);
                    }
                }
                counter += 1;
                continue;
            } else if (counter == 2) {//if this should be sequence data
                currseq += splitarr[i];
                counter += 1;
                continue;
            } else if (counter == 3) {//if this is end data
                try {
                    hend = Integer.parseInt(splitarr[i]);
                } catch (NumberFormatException e) {
                    System.err.println("unable to parse int from " + splitarr[i]);
                }
                counter = 0;
                continue;
            }
        }// end for i
        inhsp.hseq = currseq.toUpperCase();
        inhsp.hstart = hstart;
        inhsp.hend = hend;
        return inhsp;
    }// end makefinal

    //------------------------------------------------------------------------------
    static String extractqueryname(String instring, HashMap nameshash) {
        //get the name part of the query identifier block
        //see it it is present in nameshash.
        int begin;
        if ((begin = instring.indexOf("Query=</b>")) > -1) {//if this exists on this line
            instring = instring.substring(begin + 10).trim();
        }
        return instring;
    }//end extractqueryname

    //------------------------------------------------------------------------------
    static String extractname(String instring, HashMap nameshash) {
        //get the name part of the hitseq identifier block
        //see if it is present in nameshash
        int begin;
        int end;
        String tmpstr;
        if ((begin = instring.indexOf("</a>")) > -1) {
            tmpstr = instring.substring(begin + 4).trim();
            if ((end = tmpstr.indexOf(" ")) > -1) {
                instring = tmpstr.substring(0, end).trim();
            } else {
                instring = tmpstr.trim();
            }
        //if((end=instring.indexOf(" ",begin+5))>-1){
        //    instring=instring.substring(begin+4,end).trim();
        //}else{
        //    instring=instring.substring(begin+4).trim();
        //}
        }
        if (nameshash.containsKey(instring)) {
        } else {
            System.err.println("unable to parse name from " + instring);
        }
        return instring;
    }//end extractqueryname

    //------------------------------------------------------------------------------
    static double getpvalfactor(String inreadstring, int seqnum) {
        double retval = 1;
        try {
            BufferedReader inread = new BufferedReader(new StringReader(inreadstring));
            double dbsize = 0;
            double dbsizeeff = 0;
            boolean founddb = false;
            boolean founddbeff = false;
            double searchspace = -1;
            boolean nohits = false;
            // now try to read: length of database: xxx,xxx
            //and the line: effective length of database: xxx,xxx
            //use: (leffdb/ldb)*seqnum =seqnumeff
            //eval/seqnumeff=pval
            String currline;
            String[] tmp;
            while ((currline = inread.readLine()) != null) {
                if (currline.indexOf("No hits found") != -1) {
                    nohits = true;
                    break;
                }//
                if (currline.startsWith("length of database:") || currline.startsWith("Length of database:")) {
                    //get the number from this line
                    tmp = currline.split("database:\\s+", 2);
                    try {
                        dbsize = Double.parseDouble(tmp[1].replaceAll(",", ""));//get rid of the 1000-separators and convert to double
                        founddb = true;
                    } catch (NumberFormatException e) {
                        System.err.println("unable to parse double from " + tmp[1].replaceAll(",", ""));
                        founddb = false;
                    }
                }
                if (currline.startsWith("effective length of database:") || currline.startsWith("Effective length of database:")) {
                    //get the number from this line
                    tmp = currline.split("database:\\s+", 2);
                    try {
                        dbsizeeff = Double.parseDouble(tmp[1].replaceAll(",", ""));//get rid of the 1000-separators and convert to double
                        founddbeff = true;
                    } catch (NumberFormatException e) {
                        System.err.println("unable to parse double from " + tmp[1].replaceAll(",", ""));
                        founddbeff = false;
                    }
                }
                if (currline.startsWith("Effective search space used:")) {//for blast+
                    //get the number from this line
                    tmp = currline.split("used:\\s+", 2);
                    try {
                        searchspace = Double.parseDouble(tmp[1].replaceAll(",", ""));//get rid of the 1000-separators and convert to double
                    } catch (NumberFormatException e) {
                        System.err.println("unable to parse double from " + tmp[1].replaceAll(",", ""));
                    }
                }
            }//end while reading blast output
            if (founddb && founddbeff) {
                //compute retval so that evalue*retval=pval; assumes that each sequence has approx the same amount of info
                //retval=(dbsizeeff/dbsize)*(double)seqnum;
                //System.out.println("dbsize="+dbsize+" dbeff="+dbsizeeff+" seqnum="+seqnum);
                retval = dbsize / (dbsizeeff * (double) seqnum);
            } else if (nohits) {
                retval = -1;
            } else if (searchspace > -1) {//if I had blast+ output
                retval = 1/searchspace; //(E=P*searchspace, --> P=E/searchspace)
            } else {
                System.err.println("unable to read length of database and effective length of database from blast output");
                System.err.println("setting pvalue=evalue");
                retval = 1;
            }
        } catch (IOException e) {
            System.err.println("IOERROR in reading blast results; setting pvalue=evalue");
            retval = 1;
        }//end IOException
        return retval;
    }// end getpvalfactor
}
