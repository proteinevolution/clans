/*
 * enrichutils.java
 *
 * Created on October 29, 2003, 2:49 PM
 */
package clans;
import java.util.*;
import java.io.*;
/**
 *
 * @author  tancred
 */
public class enrichutils {
    
    /** Creates a new instance of enrichutils */
    public enrichutils() {
    }
    
    //--------------------------------------------------------------------------
    
    AminoAcidSequence[] enrich(AminoAcidSequence[] inseqs, String cmd, String blastpath, String formatdbpath, String[] referencedb, int cpu, double findeval, double rmsimeval, int maxseqs){
        //This "enriches" an input seq of sequences with similar but not too similar sequences from referencedb.
        //take each sequence in inseqs, do a blast/psiblast against the reference databases, and get all hits
        //better than findeval. NxN blast on this dataset should then give pairwise similarities. use these to
        //calculate the "enriched" set by removing all sequences more similar than rmsimeval or by keeping
        //the maxseqs most dissimilar sequences.
        Vector seqvec=new Vector();
        int inseqsnum=java.lang.reflect.Array.getLength(inseqs);
        for(int i=0;i<inseqsnum;i++){
            seqvec.add(inseqs[i]);
        }//end for i
        //now get all blast hits to the input sequences better than findeval but worse than rmsimeval
        System.out.println("getting blast hits for enrichment");
        String [] hitnames=getblasthits(inseqs,cmd,blastpath,referencedb,cpu,findeval,rmsimeval);
        System.out.println("done getting blast hits; hits="+java.lang.reflect.Array.getLength(hitnames));
        //now I have all sequences I want to reblast.
        //now get the full length sequences
        System.out.println("extracting sequences:");
        AminoAcidSequence[] tmpseqs=getsequences(hitnames,referencedb);
        System.out.println("done extracting; sequences="+java.lang.reflect.Array.getLength(tmpseqs));
        //add the inseqs to the tmpseqs
        for(int i=java.lang.reflect.Array.getLength(tmpseqs)-1;i>=0;i--){
            seqvec.add(tmpseqs[i]);
        }//end for i
        tmpseqs=new AminoAcidSequence[seqvec.size()];
        seqvec.copyInto(tmpseqs);
        seqvec.clear();
        //and then get the blast hit similarities matrix
        System.out.println("getting pairwise similarities");
        String tmpdbname=String.valueOf(System.currentTimeMillis())+".db";
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(tmpdbname)));
            for(int i=java.lang.reflect.Array.getLength(tmpseqs)-1;i>=0;i--){
                outwrite.println(">"+tmpseqs[i].name);
                outwrite.println(tmpseqs[i].seq);
            }//end for i
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR writing to "+tmpdbname+" in enrichutils");
            return tmpseqs;
        }//end writing blast db
        System.out.println("formatting new blast database "+tmpdbname);
        if(cmd.length()>0){
            seqvec.addElement(cmd);
        }
        String[] cmdarr=formatdbpath.split("\\s+");
        for(int i=java.lang.reflect.Array.getLength(cmdarr)-1;i>=0;i--){
            seqvec.add(0,cmdarr[i]);
        }//end for i
        seqvec.addElement("-i");
        seqvec.addElement(tmpdbname);
        cmdarr=new String[seqvec.size()];
        seqvec.copyInto(cmdarr);
        seqvec.clear();
        Runtime rt=Runtime.getRuntime();
        try{
            Process p=rt.exec(cmdarr);
            try{
                p.waitFor();
            }catch (InterruptedException ie){
                System.err.println("interrupted wait in formatdb for "+tmpdbname);
                return tmpseqs;
            }
        }catch (IOException e){
            System.err.println("IOError trying to execute formatdb for "+tmpdbname);
            for(int i=0;i<java.lang.reflect.Array.getLength(cmdarr);i++){
                System.out.println("'"+cmdarr[i]+"'");
            }
            System.err.println("skipping pairwise filtering");
            return tmpseqs;
        }
        System.out.println("Done formatdb");
        System.out.println("doing blast pairwise similarities");
        HashMap seqshash=new HashMap();
        int seqnum=java.lang.reflect.Array.getLength(tmpseqs);
        for(int i=0;i<seqnum;i++){
            seqshash.put(tmpseqs[i].name.substring(0,tmpseqs[i].name.indexOf(" ")),new Integer(i));
        }//end for i
        double[][] simmtx=getblastsim(cmd,blastpath,cpu,findeval,tmpdbname,seqshash,seqnum);//getblastsim(tmpseqs,blastpath,formatdbpath,cpu);(symmetric double[][])
        System.out.println("Done getting pairwise similarities");
        //now I have the pairwise similarities; next filter the sequences
        tmpseqs=filterseqs(inseqs,tmpseqs,rmsimeval,maxseqs,simmtx,seqshash);
        return tmpseqs;
    }//end enrich
    
    //--------------------------------------------------------------------------
    
    AminoAcidSequence[] filterseqs(AminoAcidSequence[] inseqs, AminoAcidSequence[] tmpseqs, double rmsimeval, int maxseqnum,double[][] simmtx,HashMap seqshash){
        //filter the sequences by maximum permissible similarity
        //however I want to keep the original input sequences!
        /*loop through the inseqs and remove all sequences from tmpseqs with similarities > rmsimeval.
         *start filtering with the inseqs.
         *take the best hit to that seq and filter all others for rmsimeval.
         *take the next best hit and re-filter remaining seqs
         *etc unit all hits are done. then take next input sequence and repeat.
         */
        int inseqsnum=java.lang.reflect.Array.getLength(inseqs);
        int currnum;
        int j;
        String tmpstr;
        int seqnum=java.lang.reflect.Array.getLength(tmpseqs);
        boolean[] keepseqs=new boolean[seqnum];
        for(int i=0;i<seqnum;i++){
            keepseqs[i]=true;
        }//end for i
        hitrank[] queryarr=new hitrank[inseqsnum];
        hitrank currhit;
        //filter for similarity to inseqs and remember the number and ranking of hits
        for(int i=0;i<inseqsnum;i++){
            currhit=new hitrank();
            tmpstr=inseqs[i].name.substring(0,inseqs[i].name.indexOf(" "));
            if(seqshash.get(tmpstr)!=null){
                currnum=((Integer)seqshash.get(tmpstr)).intValue();
                currhit.name=tmpstr;
                currhit.number=currnum;
                currhit.hits=0;
                currhit.hitarr=new hitvals[seqnum];
                for(j=0;j<seqnum;j++){
                    currhit.hitarr[j]=new hitvals();
                    currhit.hitarr[j].hitnum=j;
                    currhit.hitarr[j].hitval=simmtx[currnum][j];
                    if(simmtx[currnum][j]!=-1){
                        currhit.hits+=1;
                        if(keepseqs[j]==true){
                            if(simmtx[currnum][j]<=rmsimeval){
                                keepseqs[j]=false;
                            }
                        }
                    }
                }//end for j
            }else{//end if seqshash !=null
                System.err.println("unable to find input sequence "+inseqs[i].name+" in blasthits");
            }
            queryarr[i]=currhit;
        }//end for i
        //now filter for pairwise similarity in tmpseqs
        //start with inseq with the most hits (look through queryarr)
        //take the best hit to inseq 1 and filter all others
        //then take the next best hit, see if I still want to add it and if so, filter the rest
        //do this for all hits and all inseqs.
        //sort the queryarr by number of hits (most hits first)
        java.util.Arrays.sort(queryarr,new hitnumcomp());
        hitvalcomp myhitvalcomp=new hitvalcomp();
        int k;
        for(int i=0;i<inseqsnum;i++){
            java.util.Arrays.sort(queryarr[i].hitarr,myhitvalcomp);
            //now I have the hit values sorted by evalue (smallest first)
            for(j=0;j<seqnum;j++){
                if(keepseqs[queryarr[i].hitarr[j].hitnum]){//if this sequence is in the dataset
                    for(k=j+1;k<seqnum;k++){
                        if(keepseqs[k]){
                            if(simmtx[queryarr[i].hitarr[j].hitnum][k]>-1){
                                if(simmtx[queryarr[i].hitarr[j].hitnum][k]<=rmsimeval){
                                    keepseqs[k]=false;
                                }
                            }
                        }
                    }//end for k
                }//end if I wan to check this sequences hits
            }//end for j
        }//end for i
        //now make sure the input sequences are present
        for(int i=0;i<inseqsnum;i++){
            keepseqs[queryarr[i].number]=true;
        }//end for i
        Vector tmpvec=new Vector();
        for(int i=0;i<seqnum;i++){
            if(keepseqs[i]){
                tmpvec.addElement(tmpseqs[i]);
            }
        }//end for i
        AminoAcidSequence[] retarr=new AminoAcidSequence[tmpvec.size()];
        tmpvec.copyInto(retarr);
        return retarr;
    }//end fileterseqs
    
    class hitnumcomp implements Comparator{
        //inverse sort of hitranks
        public int compare(Object h1, Object h2){
            int i1=((hitrank)h1).hits;
            int i2=((hitrank)h2).hits;
            return(i1<i2 ? 1:(i1==i2 ? 0:-1));
        }//end compare
    }//en hitnumcomp
    
    class hitvalcomp implements Comparator{
        //sort of hitvals
        public int compare(Object h1, Object h2){
            double d1=((hitvals)h1).hitval;
            double d2=((hitvals)h1).hitval;
            return(d1<d2 ? -1:(d1==d2 ? 0:1));
        }//end compare
    }//end class hitvalcomp
    
    //--------------------------------------------------------------------------
    
    double[][] getblastsim(String cmd,String blastpath,int cpu,double findeval, String dbname, HashMap seqshash, int seqnum){
        //do an all against all blast search and parse the results for pairwise similarities
        Vector cmdvec=new Vector();
        if(cmd.length()>0){
            cmdvec.add(cmd);//pre-command (i.e. nice)
        }
        String[] cmdarr=blastpath.split("\\s+");
        for(int i=0;i<java.lang.reflect.Array.getLength(cmdarr);i++){
            cmdvec.addElement(cmdarr[i]);
        }
        cmdvec.addElement("-i");//input
        cmdvec.addElement(dbname);//input file
        cmdvec.addElement("-e");//maximum evalue
        cmdvec.addElement(String.valueOf(findeval));//maxeval
        cmdvec.addElement("-d");//databases to search
        cmdvec.addElement(dbname);//databases
        cmdvec.addElement("-b");
        cmdvec.addElement(String.valueOf(1));
        cmdvec.addElement("-v");//one line descriptions
        cmdvec.addElement(String.valueOf(10000));//get 10000 max
        cmdvec.addElement("-a");//multiprocessing?
        cmdvec.addElement(String.valueOf(cpu));//no. of cpu's to use
        cmdarr=new String[cmdvec.size()];
        cmdvec.copyInto(cmdarr);
        StringBuffer inout=new StringBuffer();
        StringBuffer errout=new StringBuffer();
        String myblast="";
        Runtime rt=Runtime.getRuntime();
        try{
            //write the query to a file instead of using stdin; !stdin gives problems!
            BufferedReader perr;
            BufferedReader pin;
            PrintWriter pout;
            threadstreamreader perrread;
            threadstreamreader pinread;
            Process p=rt.exec(cmdarr);
            perr=new BufferedReader(new InputStreamReader(p.getErrorStream()));//process error stream
            pin=new BufferedReader(new InputStreamReader(p.getInputStream()));//output of process
            //note: threadstreamreader closes the input bufferedreader
            //threadstreamreaders are used because otherwise deadlock can occur if the output of the process is larger than the
            //buffersize java allocates to the bufferedreaer.
            perrread=new threadstreamreader(perr,errout);
            pinread=new threadstreamreader(pin,inout);
            try{
                perrread.start();
                pinread.start();
                System.out.println("waiting for end of enrichment pairwise blast run");
                p.waitFor();
                if(p.exitValue()!=0){
                    System.err.println("non-perfect exit from blast for "+dbname);
                }
                while (pinread.done==false){
                    synchronized(inout){
                        try{
                            inout.wait(10l);//need to wait on inout else IllegalMonitorStateException (Thread does not own itself?)
                        }catch(InterruptedException e2){
                            System.err.println("interrupted sleep in blastthread");
                            e2.printStackTrace();
                        }
                    }
                }
            }catch (InterruptedException e){
                System.err.println("Interrupted pairwise enrichment process");
            }
            myblast=inout.toString();
            perrread=null;//.clear();
            pinread=null;//.clear();
        }catch (IOException ioe){
            System.err.println("IOError in enrich");
            for(int i=0;i<java.lang.reflect.Array.getLength(cmdarr);i++){
                System.err.println("'"+cmdarr[i]+"'");
            }
            System.out.println("inoutsize="+inout.length());
        }
        if(errout.indexOf("Error")>-1){
            System.err.print(errout.toString());
        }
        //now get the sequence names
        inout=null;//free the memory from these stringbuffers (more or less)
        errout=null;
        //now parse the myblast string for the evalues.
        System.out.println("parsing similarities");
        double[][] retarr=new double[seqnum][seqnum];
        int j;
        //init array with impossible values
        for(int i=0;i<seqnum;i++){
            retarr[i][i]=-1;
            for(j=i+1;j<seqnum;j++){
                retarr[i][j]=-1;
                retarr[j][i]=-1;
            }//end for j
        }//end for i
        //now parse the pairwise similarities from blast
        try{
            //System.out.println(myblast);
            BufferedReader inread=new BufferedReader(new StringReader(myblast));
            String inline;
            boolean doparse=false;
            String tmpname;
            String tmpnum;
            int hnum=-1;
            int qnum=-1;
            double thisnum;
            boolean startparse=false;
            boolean parsequery=false;
            while ((inline=inread.readLine())!=null){
                if(inline.startsWith("<b>Query=")){
                    tmpname=inline.substring(13,inline.indexOf(" ",14)).trim();
                    if(seqshash.get(tmpname)!=null){
                        qnum=((Integer)seqshash.get(tmpname)).intValue();
                    }else{
                        System.err.println("unable to resolve qnum from "+inline+" for '"+tmpname+"'");
                        qnum=-1;
                    }
                }else if(inline.startsWith("Query=")){
                    tmpname=inline.substring(7);
                    if(tmpname.indexOf(" ")>-1){
                        tmpname=tmpname.substring(0,inline.indexOf(" "));
                    }
                    if(seqshash.get(tmpname)!=null){
                        qnum=((Integer)seqshash.get(tmpname)).intValue();
                    }else{
                        System.err.println("unable to resolve qnum from "+inline+" for '"+tmpname+"'");
                        qnum=-1;
                    }
                }
                if(inline.startsWith("Sequences producing significant alignments:")){
                    startparse=true;
                    inline=inread.readLine();//skip the following line
                    continue;
                }
                if(startparse){
                    inline=inline.trim();
                    if((inline.length()<1)||(inline.equalsIgnoreCase("</PRE>"))){
                        //if I encounter an empty line or a <pre> this block is over!
                        startparse=false;
                        continue;
                    }
                    tmpname=inline.substring(0,inline.indexOf(" "));
                    tmpnum=inline.substring(inline.lastIndexOf(" ")+1);
                    if(tmpnum.startsWith("e")){
                        tmpnum="1"+tmpnum;
                    }
                    try{
                        thisnum=Double.parseDouble(tmpnum);
                    }catch (NumberFormatException ne){
                        System.err.println("ERROR parsing double from '"+tmpnum+"'");
                        continue;
                    }
                    if(qnum!=-1){
                        if(seqshash.get(tmpname)!=null){
                            hnum=((Integer)seqshash.get(tmpname)).intValue();
                        }else{
                            System.out.println("unable to get hnum from "+inline);
                            hnum=-1;
                            continue;
                        }
                    }else{
                        continue;
                    }
                    retarr[qnum][hnum]=thisnum;
                }//end if startparse
            }//end while
            inread.close();
        }catch (IOException e){
            System.err.println("IOERROR parsing similarities");
            return new double[0][0];
        }
        //now symmetrize the retarr
        for(int i=0;i<seqnum;i++){
            retarr[i][i]=-1;
            for(j=i+1;j<seqnum;j++){
                if((retarr[j][i]!=-1)&&(retarr[i][j]!=-1)){//if I don't have backvalidation
                    retarr[i][j]+=retarr[j][i];
                    retarr[i][j]=retarr[i][j]/2;
                    retarr[j][i]=retarr[i][j];
                }else{
                    retarr[j][i]=-1;
                    retarr[i][j]=-1;
                }
            }//end for j
        }//end for i
        return retarr;
    }//end getblastsim
    
    //--------------------------------------------------------------------------
    
    String[] getblasthits(AminoAcidSequence[] inseqs,String cmd, String blastpath, String[] referencedb, int cpu, double findeval, double rmsimeval){
        //search the databases in referencedb for all sequences with blast hits between rmsimeval and findeval.
        //extract those names and return them as an array.
        int seqnum=java.lang.reflect.Array.getLength(inseqs);
        //write the sequences to a start file
        String basename=String.valueOf(System.currentTimeMillis());
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(basename+".query")));
            for(int i=0;i<seqnum;i++){
                outwrite.println(">"+inseqs[i].name);
                outwrite.println(inseqs[i].seq);
            }//end for i
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOError writing to "+basename+".query");
            return new String[0];
        }
        //now run blast with the just created file as query
        //(necessary because blast dies on me in certain (rare) cases if I use standard input IO)
        Vector cmdvec=new Vector();
        String dbstring="";
        for(int i=0;i<java.lang.reflect.Array.getLength(referencedb);i++){
            dbstring+=referencedb[i]+" ";
        }//end for i
        if(cmd.length()>0){
            cmdvec.add(cmd);//pre-command (i.e. nice)
        }
        String[] cmdarr=blastpath.split("\\s+");
        for(int i=0;i<java.lang.reflect.Array.getLength(cmdarr);i++){
            cmdvec.addElement(cmdarr[i]);
        }
        cmdvec.addElement("-i");//input
        cmdvec.addElement(basename+".query");//input file
        cmdvec.addElement("-e");//maximum evalue
        cmdvec.addElement(String.valueOf(findeval));//maxeval
        cmdvec.addElement("-d");//databases to search
        cmdvec.addElement(dbstring.trim());//databases
        cmdvec.addElement("-b");
        cmdvec.addElement(String.valueOf(1));
        cmdvec.addElement("-v");//one line descriptions
        cmdvec.addElement(String.valueOf(10000));//get 10000 max
        cmdvec.addElement("-a");//multiprocessing?
        cmdvec.addElement(String.valueOf(cpu));//no. of cpu's to use
        cmdarr=new String[cmdvec.size()];
        cmdvec.copyInto(cmdarr);
        StringBuffer inout=new StringBuffer();
        StringBuffer errout=new StringBuffer();
        String myblast="";
        Runtime rt=Runtime.getRuntime();
        try{
            //write the query to a file instead of using stdin; !stdin gives problems!
            BufferedReader perr;
            BufferedReader pin;
            PrintWriter pout;
            threadstreamreader perrread;
            threadstreamreader pinread;
            Process p=rt.exec(cmdarr);
            perr=new BufferedReader(new InputStreamReader(p.getErrorStream()));//process error stream
            pin=new BufferedReader(new InputStreamReader(p.getInputStream()));//output of process
            //note: threadstreamreader closes the input bufferedreader
            //threadstreamreaders are used because otherwise deadlock can occur if the output of the process is larger than the
            //buffersize java allocates to the bufferedreaer.
            perrread=new threadstreamreader(perr,errout);
            pinread=new threadstreamreader(pin,inout);
            try{
                perrread.start();
                pinread.start();
                System.out.println("waiting for end of enrichment blast run");
                p.waitFor();
                if(p.exitValue()!=0){
                    System.err.println("non-perfect exit from blast for "+basename+".query");
                    for(int i=0;i<java.lang.reflect.Array.getLength(cmdarr);i++){
                        System.err.println("'"+cmdarr[i]+"'");
                    }
                }
                while (pinread.done==false){
                    try{
                        this.wait(10l);
                    }catch(InterruptedException e2){
                        System.err.println("interrupted sleep in blastthread");
                        e2.printStackTrace();
                    }
                }
            }catch (InterruptedException e){
                System.err.println("Interrupted process");
            }
            myblast=inout.toString();
            perrread=null;//.clear();
            pinread=null;//.clear();
        }catch (IOException ioe){
            System.err.println("IOError in enrich");
            for(int i=0;i<java.lang.reflect.Array.getLength(cmdarr);i++){
                System.err.println("'"+cmdarr[i]+"'");
            }
            System.out.println("inoutsize="+inout.length());
        }
        if(errout.indexOf("Error")>-1){
            System.err.print(errout.toString());
        }
        //now get the sequence names
        inout=null;//free the memory from these stringbuffers (more or less)
        errout=null;
        return getseqnames(myblast,findeval,rmsimeval);
    }//end getblasthits
    
    //--------------------------------------------------------------------------
    
    String[] getseqnames(String blaststring, double findeval, double rmsimeval){
        //read the sequence names from the blast results
        HashMap retmap=new HashMap();
        String[] tmpstr=new String[0];
        String noadd="noadd";
        String add="add";
        try{
            BufferedReader inread=new BufferedReader(new StringReader(blaststring));
            String inline;
            String tmpname;
            String tmpnum;
            double thisnum;
            boolean startparse=false;
            while ((inline=inread.readLine())!=null){
                if(inline.startsWith("Sequences producing significant alignments:")){
                    startparse=true;
                    inline=inread.readLine();//skip the following line
                    continue;
                }
                if(startparse){
                    inline=inline.trim();
                    if((inline.length()<1)||(inline.equalsIgnoreCase("</PRE>"))){
                        //if I encounter an empty line or a <pre> this block is over!
                        startparse=false;
                        continue;
                    }
                    tmpname=inline.substring(0,inline.indexOf(" "));
                    tmpnum=inline.substring(inline.lastIndexOf(" ")+1);
                    if(tmpnum.startsWith("e")){
                        tmpnum="1"+tmpnum;
                    }
                    try{
                        thisnum=Double.parseDouble(tmpnum);
                    }catch (NumberFormatException ne){
                        System.err.println("ERROR parsing double from '"+tmpnum+"'");
                        continue;
                    }
                    if(thisnum<=rmsimeval){
                        //I don't want to add it but want to remember it
                        retmap.put(tmpname,noadd);
                    }else if(thisnum<=findeval){
                        //see if I want to add it
                        if(retmap.get(tmpname)!=noadd){
                            retmap.put(tmpname,add);
                        }
                    }
                }//end if startparse
            }//end while reading
            inread.close();
        }catch (IOException e){
            System.err.println("IOERROR reading from blastresults in enrich");
            return new String[0];
        }
        //now go through the hashmap and see if the names should be added or not added
        Vector tmpvec=new Vector();
        tmpstr=(String[])(retmap.keySet().toArray(tmpstr));
        int namesnum=java.lang.reflect.Array.getLength(tmpstr);
        for(int i=0;i<namesnum;i++){
            if(retmap.get(tmpstr[i])==add){
                tmpvec.addElement(tmpstr[i]);
            }
        }//end for i
        String[] retarr=new String[tmpvec.size()];
        tmpvec.copyInto(retarr);
        return retarr;
    }//end getseqnames
    
    //--------------------------------------------------------------------------
    
    static AminoAcidSequence[] getsequences(String[] seqnames, String[] databases){
        //find the full length sequences from seqnames in databases and extract to aaseq array
        //databases are in fasta format, names are everything between > and first space char
        int seqnamesize=java.lang.reflect.Array.getLength(seqnames);
        //first make a hash for fast name lookup
        HashMap seqnameshash=new HashMap((int)(seqnamesize/0.75),0.8f);
        proghsp spacer=new proghsp(); //placeholder
        for(int i=0;i<seqnamesize;i++){
            seqnameshash.put(seqnames[i],spacer);
        }//end for i
        //now read from the databases and get the sequences matching the names
        int dbnum=java.lang.reflect.Array.getLength(databases);
        String instring;
        String tmpstr;
        AminoAcidSequence curraaseq;
        Vector tmpvec=new Vector();
        int counter=0;
        for(int i=0;i<dbnum;i++){
            try{
                BufferedReader inread=new BufferedReader(new FileReader(databases[i]));
                System.out.println("reading from "+databases[i]);
                while (((instring=inread.readLine())!=null)&&(counter<seqnamesize)){
                    if(instring.startsWith(">")){
                        tmpstr=instring.substring(1,instring.indexOf(" "));//get from > to first space
                        if(seqnameshash.get(tmpstr)!=null){
                            System.out.println("found "+tmpstr);
                            seqnameshash.remove(tmpstr);//change the mapping to null, don't re-extract again
                            counter++;
                            curraaseq=new AminoAcidSequence();
                            curraaseq.name=instring.substring(1);
                            curraaseq.seq="";
                            instring=inread.readLine();
                            while(((instring.startsWith(">")==false))&&(instring!=null)){
                                curraaseq.seq+=instring.trim();
                                instring=inread.readLine();
                            }//end while
                            tmpvec.addElement(curraaseq);
                        }//end if name exists in nameshash
                    }//end if starts with >
                }//end while
                inread.close();
            }catch (IOException e){
                System.err.println("IOERROR reading from "+databases[i]);
                return new AminoAcidSequence[0];
            }
        }//end for i
        AminoAcidSequence[] retarr=new AminoAcidSequence[tmpvec.size()];
        tmpvec.copyInto(retarr);
        return retarr;
    }//end getsequences
    
    //----------------------------------------------------------------------
    
    void testblastsim(double[][] inmtx){
        int num=java.lang.reflect.Array.getLength(inmtx);
        for(int i=0;i<num;i++){
            System.out.print(i+" "+inmtx[i][0]);
            for(int j=1;j<num;j++){
                System.out.print(";"+inmtx[i][j]);
            }//endf or j
            System.out.println();
        }//end for i
    }//end inmtx
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    
    class hitrank{
        
        public hitrank(){
            
        }
        
        String name="";
        int number=-1;
        int hits=-1;
        hitvals[] hitarr=new hitvals[0];
        
    }//end class hitrank
    
    class hitvals{
        public hitvals(){
        }
        
        int hitnum=-1;
        double hitval=-1;
        
    }//end class hitvals
    
}//end class enrichutils
