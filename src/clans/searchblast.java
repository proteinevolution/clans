package clans;

import java.util.*;
import java.io.*;

public class searchblast {

    public searchblast(StringBuffer errbuff,boolean addblastvbparam) {
        this.errbuff=errbuff;
        this.addblastvbparam=addblastvbparam;
    }
    
    StringBuffer errbuff;
    boolean addblastvbparam;

    /**
     * get the blast hits for the new sequences (newaln) agains a concatenated database of old and new seqs
     * 
     * @param oldaln
     * @param oldblasthits
     * @param newaln
     * @param cmd
     * @param formatdbpath
     * @param blastpath
     * @param cpu
     * @param eval
     * @param pval
     * @param coverage
     * @param scval
     * @param ident
     * @param verbose
     * @param nameshash
     * @param useallrounds
     * @param lowmem
     * @param referencedb
     * @param exhaustive
     * @param readblast
     * @param newblast
     * @return
     */
    public MinimalHsp[] gethits(AminoAcidSequence[] oldaln, MinimalHsp[] oldblasthits, AminoAcidSequence[] newaln,
            String cmd, String formatdbpath, String blastpath, int cpu, double eval, double pval, float coverage,
            float scval, float ident, int verbose, HashMap<String, Integer> nameshash, boolean useallrounds, boolean lowmem,
            String[] referencedb, int exhaustive, boolean readblast, boolean newblast) {
        System.out.println("doing searchblast for new sequences");
        int oldelements=oldaln.length;
        int newelements=newaln.length;
        int allelements=oldelements+newelements;
        double cutoff=pval;
        if(scval>=0){
            cutoff=scval;
        }
        AminoAcidSequence[] allaln=new AminoAcidSequence[allelements];
        for(int i=0;i<oldelements;i++){
            allaln[i]=oldaln[i];
        }
        for(int i=0;i<newelements;i++){
            allaln[oldelements+i]=newaln[i];
        }
        
        new File("tmpblasthsp.txt"); // this might not be necessary
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter("tmpblasthsp.txt")));
            String blastdbname=new File(String.valueOf(System.currentTimeMillis())).getAbsolutePath();//set the name for the blast database to search against
            int seqnum=allelements;
            MinimalHsp[] retarr=new MinimalHsp[0];

            if(newblast==true){
                outwrite.println(seqnum+" sequences");
            }

            //write the database to blast against to file
            try{
                PrintWriter blastdb=new PrintWriter(new BufferedWriter(new FileWriter(blastdbname)));
                for(int i=0;i<seqnum;i++){
                    allaln[i].seq=(allaln[i].seq.replaceAll("-","")).toUpperCase();
                    blastdb.println(">"+allaln[i].name);
                    blastdb.println(allaln[i].seq);//remove all gaps before printing
                }
                blastdb.close();
            }catch (IOException e){
                System.err.println("Error while creating blast database.");
                errbuff.append("-Error while creating blast database.\n");
                e.printStackTrace();
            }
            //now format the database
            if(verbose>0){
                System.out.println("formating blast database");
            }
            String formatcommand=formatdbpath+" -i "+blastdbname;//default: for the older blast version
            //System.out.println("formatdbpath='"+formatdbpath+"'");
            if(formatdbpath.indexOf("makeblastdb")>-1){//if I am using the new blast+ version of formatdb
                formatcommand=formatdbpath+" -in "+blastdbname;
            }//else{
             //   System.out.println("No makeblastdb; using old version style");
            //}
            if(cmd.length()>0){
                formatcommand=cmd+" "+formatcommand;
            }
            if(verbose>1){
                System.out.println("doing "+formatcommand);
            }
            StringBuffer errout=new StringBuffer();
            StringBuffer inout=new StringBuffer();
            threadstreamreader errread;
            threadstreamreader inread;
            Runtime rt=Runtime.getRuntime();
            try{
                Process p=rt.exec(formatcommand);
                errread=new threadstreamreader(new BufferedReader(new InputStreamReader(p.getErrorStream())),errout);
                inread=new threadstreamreader(new BufferedReader(new InputStreamReader(p.getInputStream())),inout);
                try{// wait for formatdb to finish
                    errread.start();
                    inread.start();
                    if(verbose>1){
                        System.out.println("waiting for "+formatcommand);
                    }
                    p.waitFor();
                }catch (InterruptedException e){
                    System.err.println("ERROR Interrupted formatdb");
                }
                System.err.print(errout.toString());
            }catch (IOException ex2){
                System.err.println("ERROR starting process "+formatcommand);
                errbuff.append("-ERROR starting process "+formatcommand+"\n");
                ex2.printStackTrace();
            }
            if(verbose>1){
                System.out.println("done formatting blast database with :"+formatcommand);
            }
            if(verbose>0){
                System.out.println("starting blast runs");
            }
            //I have a formatted database, so now start -cpu threads that do the blast searches.
            blastthread[] mythreads=new blastthread[cpu];
            int alldone=0;
            int allstart=0;
            if(exhaustive<2){
                alldone=oldelements;
                allstart=oldelements;
            }
            for(int i=0;(i<cpu)&&(i<newelements);i++){
                mythreads[i] = new blastthread(this, seqnum, cmd, blastpath, blastdbname, allaln[allstart], eval, pval,
                        coverage, scval, ident, new Vector[seqnum], nameshash, verbose, rt, i, useallrounds,
                        referencedb, lowmem);
                mythreads[i].start();
                allstart++;
            }
            if(verbose>1){
                System.out.println("init-blastthreads started");
            }
            int freethread;
            //now for the lowmem case:
            while (allstart<seqnum){
                freethread=-1;
                try{
                    synchronized (this){
                        while (freethread==-1){
                            for(int i=0;(i<cpu)&&(i<seqnum);i++){//see if any threads have finished
                                if(mythreads[i].done){
                                    freethread=i;
                                    break;
                                }
                            }// end for
                            if(freethread==-1){//if I don't have a free thread
                                if(verbose>1){
                                    System.out.println("Waiting(loop).........."+allstart);
                                }
                                this.wait();//release lock on this and wait
                            }
                        }// end while freethread==-1
                    }//end synchronized this
                }catch (InterruptedException e){
                    System.err.println("Interrupted wait in searchblast");
                    errbuff.append("-Interrupted wait in searchblast\n");
                    e.printStackTrace();
                }

                //write the data to disk
                printout.saveblastappend(outwrite,mythreads[freethread].retarr,nameshash.get(mythreads[freethread].query.name).intValue());
                if(mythreads[freethread].errbuff.length()>0){
                    errbuff.append(mythreads[freethread].errbuff);
                }
                mythreads[freethread] = new blastthread(this, seqnum, cmd, blastpath, blastdbname, allaln[allstart],
                        eval, pval, coverage, scval, ident, new Vector[seqnum], nameshash, verbose, rt, freethread,
                        useallrounds, referencedb, lowmem);
                alldone++;
                mythreads[freethread].start();
                allstart++;
            }// end while allstart
            //now everything has been started, now wait for the last threads to finish.
            if(verbose>1){
                System.out.println("waiting for all to finish");
            }
            while (alldone<seqnum){
                freethread=-1;
                while (freethread==-1){
                    try{
                        synchronized(this){
                            for(int i=0;(i<cpu)&&(i<seqnum);i++){
                                if(mythreads[i].done){
                                    freethread=i;
                                    break;
                                }
                            }
                            if(freethread==-1){
                                if(verbose>1){
                                    System.out.println("Waiting(end)..."+alldone);
                                }
                                this.wait();
                            }
                        }// end synchronized this
                    }catch(InterruptedException e){
                        System.err.println("interrupted wait in searchblast");
                        errbuff.append("interrupted wait in searchblast\n");
                        e.printStackTrace();
                    }
                }// end while freethread==-1

                //write the data to disk and clear that array element
                printout.saveblastappend(outwrite,mythreads[freethread].retarr,nameshash.get(mythreads[freethread].query.name).intValue());
                if(mythreads[freethread].errbuff.length()>0){
                    errbuff.append(mythreads[freethread].errbuff);
                }
                alldone++;
                mythreads[freethread].done=false;
                if((exhaustive!=1)){
                    outwrite.println("DONE");
                }
                if(exhaustive!=1){
                    //if exhaustive==1 I still want to write to the outputfile
                    outwrite.close();
                }else{
                    outwrite.flush();
                }
            }// end while alldone
            
            if(exhaustive==1){
                //if non-exhaustive searching was done and I want to back-validate the hits I collected for the new sequences
                //repeat the search from above for the sequences with hits to the new seqs
                //first get the sequences with hits
                retarr=readsave.blast("tmpblasthsp.txt",cutoff);
                Vector<Integer> tmpvec = new Vector<Integer>();
                boolean[] redoseqs=new boolean[oldelements];
                for(int i=0;i<oldelements;i++){
                    redoseqs[i]=false;
                }
                for(int i=retarr.length;--i>=0;){
                    if((retarr[i].query>=oldelements)&&(retarr[i].hit<oldelements)){
                        if(redoseqs[retarr[i].hit]==false){
                            tmpvec.addElement(new Integer(retarr[i].hit));
                        }
                        redoseqs[retarr[i].hit]=true;
                    }
                }
                int backvalnum=tmpvec.size();
                int[] backvalseqs=new int[backvalnum];
                for(int i=0;i<backvalnum;i++){
                    backvalseqs[i]=tmpvec.elementAt(i).intValue();
                }
                tmpvec=null;
                if(verbose>0){
                    System.out.println("Backvalidation required for "+backvalnum+" sequences");
                }
                //retarr=null;//clear a large chunk of memory
                //now I should have an array of sequence numbers for which I need to repeat the blast runs
                mythreads=new blastthread[cpu];
                alldone=0;
                allstart=0;
                for(int i=0;(i<cpu)&&(i<backvalnum);i++){
                    mythreads[i] = new blastthread(this, seqnum, cmd, blastpath, blastdbname,
                            allaln[backvalseqs[allstart]], eval, pval, coverage, scval, ident, new Vector[seqnum],
                            nameshash, verbose, rt, i, useallrounds, referencedb, lowmem);
                    mythreads[i].start();
                    allstart++;
                }
                if(verbose>1){
                    System.out.println("init-backvalidation started");
                }

                //now for the lowmem case:
                while (allstart<backvalnum){
                    freethread=-1;
                    try{
                        synchronized (this){
                            while (freethread==-1){
                                for(int i=0;(i<cpu)&&(i<backvalnum);i++){//see if any threads have finished
                                    if(mythreads[i].done){
                                        freethread=i;
                                        break;
                                    }
                                }// end for
                                if(freethread==-1){//if I don't have a free thread
                                    if(verbose>1){
                                        System.out.println("Waiting(backvalidation loop).........."+allstart);
                                    }
                                    this.wait();//release lock on this and wait
                                }
                            }// end while freethread==-1
                        }//end synchronized this
                    }catch (InterruptedException e){
                        System.err.println("Interrupted wait in searchblast");
                        errbuff.append("interrupted wait in searchblast\n");
                        e.printStackTrace();
                    }

                    //write the data to disk and clear that array element
                    printout.saveblastappend(outwrite,mythreads[freethread].retarr,nameshash.get(mythreads[freethread].query.name).intValue());
                    if(mythreads[freethread].errbuff.length()>0){
                        errbuff.append(mythreads[freethread].errbuff);
                    }
                    mythreads[freethread] = new blastthread(this, seqnum, cmd, blastpath, blastdbname,
                            allaln[backvalseqs[allstart]], eval, pval, coverage, scval, ident, new Vector[seqnum],
                            nameshash, verbose, rt, freethread, useallrounds, referencedb, lowmem);
                    alldone++;
                    mythreads[freethread].start();
                    allstart++;
                }// end while allstart
                //now everything has been started, now wait for the last threads to finish.
                if(verbose>1){
                    System.out.println("waiting for backvalidation to finish");
                }
                while (alldone<backvalnum){
                    freethread=-1;
                    while (freethread==-1){
                        try{
                            synchronized(this){
                                for(int i=0;(i<cpu)&&(i<seqnum);i++){
                                    if(mythreads[i].done){
                                        freethread=i;
                                        break;
                                    }
                                }
                                if(freethread==-1){
                                    if(verbose>1){
                                        System.out.println("Waiting(backvalidation end)..."+alldone);
                                    }
                                    this.wait();
                                }
                            }// end synchronized this
                        }catch(InterruptedException e){
                            System.err.println("interrupted wait in searchblast");
                            errbuff.append("interrupted wait in searchblast\n");
                            e.printStackTrace();
                        }
                    }// end while freethread==-1

                    //write the data to disk and clear that array element
                    printout.saveblastappend(outwrite,mythreads[freethread].retarr,nameshash.get(mythreads[freethread].query.name).intValue());
                    if(mythreads[freethread].errbuff.length()>0){
                        errbuff.append(mythreads[freethread].errbuff);
                    }
                    alldone++;
                    mythreads[freethread].done=false;

                    outwrite.println("#DONE all backvalidation blast runs");
                }

                retarr=null;
                outwrite.println("DONE");
                outwrite.close();
            }//end if exhaustive==1
            System.out.println("re-reading hsp's");
            retarr=readsave.blast("tmpblasthsp.txt",cutoff);
            if(exhaustive==2){
                //then I redid all blast runs and have all data in retarr
            }else{
                //I need to look at what parts to keep and which to add
                Vector<MinimalHsp> tmpvec = new Vector<MinimalHsp>();
                if(exhaustive<2){
                    //if I didn't redo the old blast runs; transfer the old info
                    for(int i=retarr.length-1;i>=0;i--){
                        tmpvec.addElement(retarr[i]);
                    }
                    for(int i=0;i<oldelements;i++){
                        tmpvec.addElement(oldblasthits[i]);
                    }
                }//end if exhaustive<2
                //now for a fast version just mirror the blast values for the new sequences
                if(exhaustive==0){
                    //if I want to make a->b == b->a (one way search)
                    for(int i=retarr.length-1;i>=0;i--){
                        tmpvec.addElement(new MinimalHsp(retarr[i].hit,retarr[i].query,retarr[i].val));//inverted query and hit info
                    }
                }//end if exhaustive==0
            }
            //all threads started and done. the values are stored in retarr.
            try{
                //delete the blast database and the files that formatdb created and the blast query files
                //files to check are: blastdbname; blastdbname.psd; blastdbname.phr; blastdbname.pin; blastdbname.pnd
                //blastdbname.pni; blastdbname.psi; blastdbname.psq
                File rmfile;
                for(int i=0;i<cpu;i++){
                    rmfile=new File("tmp"+i+".query");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    rmfile=new File("tmp"+i+".chkpnt");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                }
                rmfile=new File(blastdbname);
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".psd");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".phr");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".pin");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".pnd");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".pni");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".psi");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".psq");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File("formatdb.log");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File("tmpblasthsp.txt");
                if(rmfile.exists()){
                    rmfile.delete();
                }
            }catch (RuntimeException e){
                System.err.println("unable to delete formatdb files of "+blastdbname);
                e.printStackTrace();
            }
            return retarr;
        }catch (IOException ioe){
            System.err.println("IOError writing to tmpblasthsp.txt");
            errbuff.append("-IOERROR writing to tmpblasthsp.txt\n");
            ioe.printStackTrace();
            return new MinimalHsp[0];
        }
    }

    /**
     * run an all against all blast search and save the hsp's in vector[][]
     * 
     * @param inarr
     * @param cmd
     * @param formatdbpath
     * @param blastpath
     * @param cpu
     * @param eval
     * @param pval
     * @param coverage
     * @param scval
     * @param ident
     * @param verbose
     * @param nameshash
     * @param useallrounds
     * @param lowmem
     * @param referencedb
     * @param readblast
     * @param newblast
     * @return
     */
    public MinimalHsp[] gethits(AminoAcidSequence[] inarr, String cmd, String formatdbpath, String blastpath, int cpu,
            double eval, double pval, float coverage, float scval, float ident, int verbose,
            HashMap<String, Integer> nameshash, boolean useallrounds, boolean lowmem, String[] referencedb,
            boolean readblast, boolean newblast) {
        System.out.println("doing searchblast");
        File tmpfile=new File("tmpblasthsp.txt");
        if(tmpfile.canRead()&&readblast){
            System.out.println("reading former blastruns from tmpblasthsp.txt");
            MinimalHsp[] retarr=readsave.blast("tmpblasthsp.txt",pval);
            return retarr;
        }else{
            try{
                PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter("tmpblasthsp.txt")));
                String blastdbname=new File(String.valueOf(System.currentTimeMillis())).getAbsolutePath();//set the name for the blast database to search against
                int seqnum=inarr.length;
                //if(lowmem){
                if(newblast==true){
                    outwrite.println(seqnum+" sequences");
                }
                //}
                //write the database to blast against to file
                try{
                    PrintWriter blastdb=new PrintWriter(new BufferedWriter(new FileWriter(blastdbname)));
                    //if(lowmem){
                    for(int i=0;i<seqnum;i++){
                        inarr[i].seq=(inarr[i].seq.replaceAll("-","")).toUpperCase();
                        blastdb.println(">"+inarr[i].name);
                        blastdb.println(inarr[i].seq);//remove all gaps before printing
                    }
                    //}else{
                    //    for(int i=0;i<seqnum;i++){
                    //        blastdb.println(">"+inarr[i].name);
                    //        blastdb.println((inarr[i].seq.replaceAll("-","")).toUpperCase());//remove all gaps before printing
                    //    }
                    //}
                    blastdb.close();
                }catch (IOException e){
                    System.err.println("Error while creating blast database.");
                    errbuff.append("Error while creating blast database.\n");
                    e.printStackTrace();
                }
                //now format the database
                if(verbose>0){
                    System.out.println("formating blast database");
                }
                String formatcommand=formatdbpath+" -i "+blastdbname;//default: for the older blast version
                //System.out.println("formatdbpath='"+formatdbpath+"'");
                if(formatdbpath.indexOf("makeblastdb")>-1){//if I am using the new blast+ version of formatdb
                    formatcommand=formatdbpath+" -in "+blastdbname;
                }//else{
                 //   System.out.println("No makeblastdb; using old version style");
                //}
                if(cmd.length()>0){
                    formatcommand=cmd+" "+formatcommand;
                }
                if(verbose>1){
                    System.out.println("doing "+formatcommand);
                }
                StringBuffer errout=new StringBuffer();
                StringBuffer inout=new StringBuffer();
                threadstreamreader errread;
                threadstreamreader inread;
                Runtime rt=Runtime.getRuntime();
                try{
                    Process p=rt.exec(formatcommand);
                    errread=new threadstreamreader(new BufferedReader(new InputStreamReader(p.getErrorStream())),errout);
                    inread=new threadstreamreader(new BufferedReader(new InputStreamReader(p.getInputStream())),inout);
                    try{// wait for formatdb to finish
                        errread.start();
                        inread.start();
                        if(verbose>1){
                            System.out.println("waiting for "+formatcommand);
                        }
                        p.waitFor();
                    }catch (InterruptedException e){
                        System.err.println("ERROR Interrupted formatdb");
                    }
                    System.err.print(errout.toString());
                }catch (IOException ex2){
                    System.err.println("ERROR starting process "+formatcommand);
                    errbuff.append("-ERROR starting process "+formatcommand+"\n");
                    ex2.printStackTrace();
                }
                if(verbose>1){
                    System.out.println("done formatting blast database with :"+formatcommand);
                }
                if(verbose>0){
                    System.out.println("starting blast runs");
                }

                //I have a formatted database, so now start -cpu threads that do the blast searches.
                if(verbose>1){
                    System.out.println("done Vector initialization");
                }
                blastthread[] mythreads=new blastthread[cpu];
                int alldone=0;
                int allstart=0;
                for(int i=0;(i<cpu)&&(i<seqnum);i++){
                    mythreads[i] = new blastthread(this, seqnum, cmd, blastpath, blastdbname, inarr[allstart], eval,
                            pval, coverage, scval, ident, new Vector[seqnum], nameshash, verbose, rt, i, useallrounds,
                            referencedb, lowmem);
                    mythreads[i].start();
                    allstart++;
                }
                
                if(verbose>1){
                    System.out.println("init-blastthreads started");
                }
                int freethread;
                //now for the lowmem case:
                while (allstart<seqnum){
                    freethread=-1;
                    try{
                        synchronized (this){
                            while (freethread==-1){
                                for(int i=0;(i<cpu)&&(i<seqnum);i++){//see if any threads have finished
                                    if(mythreads[i].done){
                                        freethread=i;
                                        break;
                                    }
                                }// end for
                                if(freethread==-1){//if I don't have a free thread
                                    if(verbose>1){
                                        System.out.println("Waiting(loop).........."+allstart);
                                    }
                                    this.wait();//release lock on this and wait
                                }
                            }// end while freethread==-1
                        }//end synchronized this
                    }catch (InterruptedException e){
                        System.err.println("Interrupted wait in searchblast");
                        errbuff.append("Interrupted wait in searchblast\n");
                        e.printStackTrace();
                    }
                    //write the data to disk and clear that array element
                    printout.saveblastappend(outwrite,mythreads[freethread].retarr,nameshash.get(mythreads[freethread].query.name).intValue());

                    if(mythreads[freethread].errbuff.length()>0){
                        errbuff.append(mythreads[freethread].errbuff);
                    }
                    mythreads[freethread] = new blastthread(this, seqnum, cmd, blastpath, blastdbname, inarr[allstart],
                            eval, pval, coverage, scval, ident, new Vector[seqnum], nameshash, verbose, rt, freethread,
                            useallrounds, referencedb, lowmem);
                    alldone++;
                    mythreads[freethread].start();
                    allstart++;
                }
                
                //now everything has been started, now wait for the last threads to finish.
                if(verbose>1){
                    System.out.println("waiting for all to finish");
                }
                while (alldone<seqnum){
                    freethread=-1;
                    while (freethread==-1){
                        try{
                            synchronized(this){
                                for(int i=0;(i<cpu)&&(i<seqnum);i++){
                                    if(mythreads[i].done){
                                        freethread=i;
                                        break;
                                    }
                                }
                                if(freethread==-1){
                                    if(verbose>1){
                                        System.out.println("Waiting(end)..."+alldone);
                                    }
                                    this.wait();
                                }
                            }// end synchronized this
                        }catch(InterruptedException e){
                            System.err.println("interrupted wait in searchblast");
                            errbuff.append("Interrupted wait in searchblast\n");
                            e.printStackTrace();
                        }
                    }// end while freethread==-1
                    //if(lowmem){
                    //write the data to disk and clear that array element
                    printout.saveblastappend(outwrite,mythreads[freethread].retarr,nameshash.get(mythreads[freethread].query.name).intValue());
                    //retarr[((Integer)nameshash.get(mythreads[freethread].query.name)).intValue()]=new Vector[seqnum];
                    //}else{
                    //    retarr[((Integer)nameshash.get(mythreads[freethread].query.name)).intValue()]=mythreads[freethread].retarr;//get the array with the computed values from that thread
                    //}
                    if(mythreads[freethread].errbuff.length()>0){
                        errbuff.append(mythreads[freethread].errbuff);
                    }
                    alldone++;
                    mythreads[freethread].done=false;
                    outwrite.println("DONE");
                    outwrite.close();
                }// end while alldone
                //if(lowmem){
                System.out.println("re-reading hsp's");
                MinimalHsp[] retarr=readsave.blast("tmpblasthsp.txt",pval);
                //}
                //all threads started and done. the values are stored in retarr.
                try{
                    //delete the blast database and the files that formatdb created and the blast query files
                    //files to check are: blastdbname; blastdbname.psd; blastdbname.phr; blastdbname.pin; blastdbname.pnd
                    //blastdbname.pni; blastdbname.psi; blastdbname.psq
                    File rmfile;
                    for(int i=0;i<cpu;i++){
                        rmfile=new File("tmp"+i+".query");
                        if(rmfile.exists()){
                            rmfile.delete();
                        }
                        rmfile=new File("tmp"+i+".chkpnt");
                        if(rmfile.exists()){
                            rmfile.delete();
                        }
                    }
                    rmfile=new File(blastdbname);
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    rmfile=new File(blastdbname+".psd");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    rmfile=new File(blastdbname+".phr");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    rmfile=new File(blastdbname+".pin");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    rmfile=new File(blastdbname+".pnd");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    rmfile=new File(blastdbname+".pni");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    rmfile=new File(blastdbname+".psi");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    rmfile=new File(blastdbname+".psq");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    rmfile=new File("formatdb.log");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    //rmfile=new File("tmpblasthsp.txt");
                    //if(rmfile.exists()){
                    //    rmfile.delete();
                    //}
                }catch (RuntimeException e){
                    System.err.println("unable to delete formatdb files of "+blastdbname);
                    e.printStackTrace();
                }
                return retarr;
            }catch (IOException ioe){
                System.err.println("IOError writing to tmpblasthsp.txt");
                errbuff.append("-IOERROR writing to tmpblasthsp.txt");
                ioe.printStackTrace();
                return new MinimalHsp[0];
            }
        }
    }

    /**
     * do all blast runs for sequences present in seqstodo NOT FOR THE OTHERS (should only happen for interrupted runs)
     * 
     * @param blasthits
     * @param inarr
     * @param cmd
     * @param formatdbpath
     * @param blastpath
     * @param cpu
     * @param eval
     * @param pval
     * @param coverage
     * @param scval
     * @param ident
     * @param verbose
     * @param nameshash
     * @param useallrounds
     * @param lowmem
     * @param referencedb
     * @param seqstodo
     * @param readblast
     * @param newblast
     * @return
     */
    public MinimalHsp[] gethits(MinimalHsp[] blasthits, AminoAcidSequence[] inarr, String cmd, String formatdbpath,
            String blastpath, int cpu, double eval, double pval, float coverage, float scval, float ident, int verbose,
            HashMap<String, Integer> nameshash, boolean useallrounds, boolean lowmem, String[] referencedb,
            int[] seqstodo, boolean readblast, boolean newblast) {

        System.out.println("continuing searchblast");
        File tmpfile=new File("tmpblasthsp.txt");
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(tmpfile,true)));//append to the end of the file
            String blastdbname=new File(String.valueOf(System.currentTimeMillis())).getAbsolutePath();//set the name for the blast database to search against
            int seqnum=inarr.length;
            //if(lowmem){
            if(newblast==true){
                outwrite.println(seqnum+" sequences");
            }
            //}
            //write the database to blast against to file
            try{
                PrintWriter blastdb=new PrintWriter(new BufferedWriter(new FileWriter(blastdbname)));
                //if(lowmem){
                for(int i=0;i<seqnum;i++){
                    inarr[i].seq=(inarr[i].seq.replaceAll("-","")).toUpperCase();
                    blastdb.println(">"+inarr[i].name);
                    blastdb.println(inarr[i].seq);//remove all gaps before printing
                }
                //}else{
                //    for(int i=0;i<seqnum;i++){
                //        blastdb.println(">"+inarr[i].name);
                //        blastdb.println((inarr[i].seq.replaceAll("-","")).toUpperCase());//remove all gaps before printing
                //    }
                //}
                blastdb.close();
            }catch (IOException e){
                System.err.println("Error while creating blast database.");
                errbuff.append("-Error while creating blast database.\n");
                e.printStackTrace();
            }
            //now format the database
            if(verbose>0){
                System.out.println("formatting blast database");
            }
            String formatcommand=formatdbpath+" -i "+blastdbname;//default: for the older blast version
            //System.out.println("formatdbpath='"+formatdbpath+"'");
            if(formatdbpath.indexOf("makeblastdb")>-1){//if I am using the new blast+ version of formatdb
                formatcommand=formatdbpath+" -in "+blastdbname;
            }//else{
             //   System.out.println("No makeblastdb; using old version style");
            //}
            if(cmd.length()>0){
                formatcommand=cmd+" "+formatcommand;
            }
            if(verbose>1){
                System.out.println("doing "+formatcommand);
            }
            StringBuffer errout=new StringBuffer();
            StringBuffer inout=new StringBuffer();
            threadstreamreader errread;
            threadstreamreader inread;
            Runtime rt=Runtime.getRuntime();
            try{
                Process p=rt.exec(formatcommand);
                errread=new threadstreamreader(new BufferedReader(new InputStreamReader(p.getErrorStream())),errout);
                inread=new threadstreamreader(new BufferedReader(new InputStreamReader(p.getInputStream())),inout);
                try{// wait for formatdb to finish
                    errread.start();
                    inread.start();
                    if(verbose>1){
                        System.out.println("waiting for "+formatcommand);
                    }
                    p.waitFor();
                }catch (InterruptedException e){
                    System.err.println("ERROR Interrupted formatdb");
                }
                System.err.print(errout.toString());
            }catch (IOException ex2){
                System.err.println("ERROR starting process "+formatcommand);
                errbuff.append("ERROR starting process "+formatcommand+"\n");
                ex2.printStackTrace();
            }
            if(verbose>1){
                System.out.println("done formatting blast database with :"+formatcommand);
            }
            //I have a formatted database, so now start -cpu threads that do the blast searches.
            blastthread[] mythreads=new blastthread[cpu];
            int seqstodonum=seqstodo.length;
            int alldone=0;
            int allstart=0;
            for(int i=0;(i<cpu)&&(i<seqstodonum);i++){
                mythreads[i] = new blastthread(this, seqnum, cmd, blastpath, blastdbname, inarr[seqstodo[allstart]],
                        eval, pval, coverage, scval, ident, new Vector[seqnum], nameshash, verbose, rt, i,
                        useallrounds, referencedb, lowmem);
                mythreads[i].start();
                allstart++;
            }
            if(verbose>1){
                System.out.println("init-blastthreads started");
            }
            int freethread;
            //now for the lowmem case:
            while ((allstart<seqnum)&&(allstart<seqstodonum)){
                freethread=-1;
                try{
                    synchronized (this){
                        while (freethread==-1){
                            for(int i=0;(i<cpu)&&(i<seqnum);i++){//see if any threads have finished
                                if(mythreads[i].done){
                                    freethread=i;
                                    break;
                                }
                            }// end for
                            if(freethread==-1){//if I don't have a free thread
                                if(verbose>1){
                                    System.out.println("Waiting(loop).........."+allstart);
                                }
                                this.wait();//release lock on this and wait
                            }
                        }// end while freethread==-1
                    }//end synchronized this
                }catch (InterruptedException e){
                    System.err.println("Interrupted wait in searchblast");
                    errbuff.append("-Interrupted wait in searchblast\n");
                    e.printStackTrace();
                }

                //write the data to disk and clear that array element
                printout.saveblastappend(outwrite,mythreads[freethread].retarr,nameshash.get(mythreads[freethread].query.name).intValue());
                if(mythreads[freethread].errbuff.length()>0){
                    errbuff.append(mythreads[freethread].errbuff);
                }
                mythreads[freethread] = new blastthread(this, seqnum, cmd, blastpath, blastdbname,
                        inarr[seqstodo[allstart]], eval, pval, coverage, scval, ident, new Vector[seqnum], nameshash,
                        verbose, rt, freethread, useallrounds, referencedb, lowmem);
                alldone++;
                mythreads[freethread].start();
                allstart++;
            }// end while allstart
            
            //now everything has been started, now wait for the last threads to finish.
            if(verbose>1){
                System.out.println("waiting for all to finish");
            }
            while ((alldone<seqnum)&&(alldone<seqstodonum)){
                freethread=-1;
                while (freethread==-1){
                    try{
                        synchronized(this){
                            for(int i=0;(i<cpu)&&(i<seqnum);i++){
                                if(mythreads[i].done){
                                    freethread=i;
                                    break;
                                }
                            }
                            if(freethread==-1){
                                if(verbose>1){
                                    System.out.println("Waiting(end)..."+alldone);
                                }
                                this.wait();
                            }
                        }// end synchronized this
                    }catch(InterruptedException e){
                        System.err.println("interrupted wait in searchblast");
                        errbuff.append("-Interrupted wait in searchblast\n");
                        e.printStackTrace();
                    }
                }

                //write the data to disk and clear that array element
                printout.saveblastappend(outwrite,mythreads[freethread].retarr,nameshash.get(mythreads[freethread].query.name).intValue());
                if(mythreads[freethread].errbuff.length()>0){
                    errbuff.append(mythreads[freethread].errbuff);
                }
                alldone++;
                mythreads[freethread].done=false;
                outwrite.println("DONE");
                outwrite.close();
            }

            System.out.println("re-reading hsp's");
            blasthits=readsave.blast("tmpblasthsp.txt",pval);
            //}
            //all threads started and done. the values are stored in retarr.
            try{
                //delete the blast database and the files that formatdb created and the blast query files
                //files to check are: blastdbname; blastdbname.psd; blastdbname.phr; blastdbname.pin; blastdbname.pnd
                //blastdbname.pni; blastdbname.psi; blastdbname.psq
                File rmfile;
                for(int i=0;i<cpu;i++){
                    rmfile=new File("tmp"+i+".query");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                    rmfile=new File("tmp"+i+".chkpnt");
                    if(rmfile.exists()){
                        rmfile.delete();
                    }
                }
                rmfile=new File(blastdbname);
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".psd");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".phr");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".pin");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".pnd");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".pni");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".psi");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File(blastdbname+".psq");
                if(rmfile.exists()){
                    rmfile.delete();
                }
                rmfile=new File("formatdb.log");
                if(rmfile.exists()){
                    rmfile.delete();
                }
            }catch (RuntimeException e){
                System.err.println("unable to delete formatdb files of "+blastdbname);
                e.printStackTrace();
            }
            return blasthits;
        }catch (IOException ioe){
            System.err.println("IOError writing to tmpblasthsp.txt");
            errbuff.append("IOERROR writing to tmpblasthsp.txt\n");
            ioe.printStackTrace();
            return new MinimalHsp[0];
        }
    }// end gethits redoing for missing sequences
    
    //--------------------------------------------------------------------------
    
}

class blastthread extends java.lang.Thread{
    //this class should perform the blast run and filter the output according to coverage,scval,ident,etc.
    
    /**
     * 
     * @param parent
     * @param allseqnum
     * @param cmd
     * @param blastpath
     * @param blastdbname
     * @param query
     * @param eval
     * @param pval
     * @param coverage
     * @param scval
     * @param ident
     * @param retvecarr
     * @param nameshash
     * @param verbose
     * @param rt
     * @param threadnum
     * @param useallrounds
     * @param referencedb
     * @param lowmem
     */
    public blastthread(searchblast parent, int allseqnum, String cmd, String blastpath, String blastdbname,
            AminoAcidSequence query, double eval, double pval, float coverage, float scval, float ident,
            Vector<hsp>[] retvecarr, HashMap<String, Integer> nameshash, int verbose, Runtime rt, int threadnum,
            boolean useallrounds, String[] referencedb, boolean lowmem) {
        this.parent=parent;
        this.allseqnum=allseqnum;
        this.cmd=cmd;
        this.blastpath=blastpath;
        this.blastdbname=blastdbname;
        this.query=query;
        this.eval=eval;
        this.pval=pval;
        this.coverage=coverage;
        this.scval=scval;
        this.ident=ident;
        this.verbose=verbose;
        this.retarr=retvecarr;
        this.nameshash=nameshash;
        this.done=false;
        this.rt=rt;
        this.threadnum=threadnum;
        this.useallrounds=useallrounds;
        this.referencedb=referencedb;
        this.lowmem=lowmem;
        if(parent.addblastvbparam==true){//default: true!
            //System.out.println("checking addblastvb");
            //i.e. if I want to check to see whether the -v or -b parameters were set in blastpath
            if(allseqnum>250){
                if(blastpath.indexOf("blastall")>-1 || blastpath.indexOf("blastpgp")>-1){//if using the old blast version
                    //System.out.println("allseqnum>250");
                    //if I have more sequences than BLAST would normally return hits for
                    //see whether the -v and -b settings are set
                    if(blastpath.indexOf("-b ")>-1 || blastpath.indexOf("-v ")>-1){
                        //System.out.println("blastpath '"+blastpath+"' contains -v or -b entry; skipping");
                        //if this IS set, then don't do anything
                    }else{
                        //set these to allow at least one hit (on average) per sequence
                        this.blastpath+=" -b "+allseqnum+" -v "+allseqnum;
                        //System.out.println("setting blastpath to '"+this.blastpath+"'");
                    }
                }else{//if using the newer blast+ version
                    if(blastpath.indexOf("-num_descriptions ")>-1 || blastpath.indexOf("-num_alignments ")>-1){
                        //System.out.println("blastpath '"+blastpath+"' contains -v or -b entry; skipping");
                        //if this IS set, then don't do anything
                    }else{
                        //set these to allow at least one hit (on average) per sequence
                        this.blastpath+=" -num_descriptions "+allseqnum+" -num_alignments "+allseqnum;
                        //System.out.println("setting blastpath to '"+this.blastpath+"'");
                    }
                }
            }
        }
    }// end blastthread
    
    String cmd;
    String blastpath;
    String blastdbname;
    String[] referencedb;
    AminoAcidSequence query;
    double eval;
    double pval;
    float coverage;
    float scval;
    float ident;
    int verbose;
    int allseqnum=0;
    Vector<hsp>[] retarr;
    HashMap<String, Integer> nameshash;
    searchblast parent;
    boolean done=false;
    boolean stopthread=false;
    boolean useallrounds=false;//use all psiblast rounds to get hsp's?
    boolean lowmem=false;
    Runtime rt;
    int threadnum;
    StringBuffer errbuff=new StringBuffer();
    //--------------------------------------------------------------------------
    
    public void run(){
        //check to see if there is sth assigned to the query
        done=false;
        stopthread=false;
        //String blastcommand="";
        String tmpfilestring="tmp"+threadnum+".query";
        String tmpcheckfile="tmp"+threadnum+".chkpnt";
        String tmpoutfile="tmp"+threadnum+".outfile";
        String dbstring;
        String[] cmdarr;
        //convert the String commands to a array of Strings (resolves some runtime.exec problems)
        Vector<String> tmpvec = new Vector<String>();
        if(blastpath.indexOf("blastall")>-1 || blastpath.indexOf("blastpgp")>-1){//for the old blast version
            cmdarr=blastpath.split("\\s+");
            tmpvec.addElement("-T");
            tmpvec.addElement("T");
            tmpvec.addElement("-e");
            tmpvec.addElement(String.valueOf(eval));
            tmpvec.addElement("-i");
            tmpvec.addElement(new File(tmpfilestring).getAbsolutePath());
            tmpvec.addElement("-d");
            tmpvec.addElement(blastdbname);
            for(int i=cmdarr.length;i>0;i--){
                tmpvec.add(0,cmdarr[i-1]);
            }
            if(cmd.length()>0){
                //blastcommand=cmd+" "+blastcommand;
                tmpvec.add(0, cmd);
            }
        }else{//for the new blast+ executables
            cmdarr=blastpath.split("\\s+");
            tmpvec.addElement("-html");
            tmpvec.addElement("-show_gis");
            tmpvec.addElement("-evalue");
            tmpvec.addElement(String.valueOf(eval));
            tmpvec.addElement("-query");
            tmpvec.addElement(new File(tmpfilestring).getAbsolutePath());
            tmpvec.addElement("-db");
            tmpvec.addElement(blastdbname);
            for(int i=cmdarr.length;i>0;i--){
                tmpvec.add(0,cmdarr[i-1]);
            }
            if(cmd.length()>0){
                //blastcommand=cmd+" "+blastcommand;
                tmpvec.add(0, cmd);
            }
        }
        
        if((referencedb.length>0)&&(blastpath.indexOf("blastpgp")>-1)){
            if(verbose>1){
                System.out.println("profile psiblast run on "+threadnum);
            }
            
            @SuppressWarnings("unchecked")
            Vector<String> tmpvec2 = (Vector<String>) tmpvec.clone();// supposed to hold the command parameters
            
            //first do the psiblast runs agains all specified databases
            //first put all database names into one string
            dbstring=blastdbname;
            for(int i=0;i<referencedb.length;i++){
                dbstring+=" "+referencedb[i];
            }
            //blastcommand=blastpath+" -T T -e "+eval+" -i "+tmpfilestring+" -C "+tmpcheckfile+" -d "+dbstring;
            //cmdarr=blastpath.split("\\s+");
            //tmpvec.addElement("-T");
            //tmpvec.addElement("T");
            //tmpvec.addElement("-e");
            //tmpvec.addElement(String.valueOf(eval));
            //tmpvec.addElement("-i");
            //tmpvec.addElement(tmpfilestring);
            tmpvec2.addElement("-C");
            tmpvec2.addElement(tmpcheckfile);
            tmpvec2.addElement("-d");
            tmpvec2.addElement(dbstring);
            if(lowmem){
                tmpvec2.addElement("-o");
                tmpvec2.addElement(tmpoutfile);
            }
            //for(int i=cmdarr.length;i>0;i--){
            //    tmpvec.add(0,cmdarr[i-1]);
            //}
            //if(cmd.length()>0){
            //    blastcommand=cmd+" "+blastcommand;
            //    tmpvec.add(0, cmd);
            //}
            cmdarr=new String[tmpvec2.size()];
            tmpvec2.copyInto(cmdarr);
            StringBuffer inout=new StringBuffer();
            StringBuffer errout=new StringBuffer();
            if(verbose>2){
                //System.out.println("for "+query.name+" trying command="+blastcommand);
                System.out.print(query.name+" trying:");
                for(int i=0;i<cmdarr.length;i++){
                    System.out.print("'"+cmdarr[i]+"' ");
                }
                System.out.println();
            }
            try{
                //write the query to a file instead of using stdin; !stdin gives problems!
                PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(tmpfilestring)));
                outwrite.println(">"+query.name);
                outwrite.println(query.seq.replaceAll("-",""));
                outwrite.close();
                //System.out.println("done writing for "+tmpfilestring);
                BufferedReader perr;
                BufferedReader pin;
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
                    if(verbose>2){
                        System.out.println("feeding seqdata >"+query.name);
                        System.out.println("seq: "+query.seq.replaceAll("-",""));
                    }
                    if(verbose>2){
                        //System.out.println("Done feeding, waiting for "+blastcommand+" query="+query.name);
                        System.out.print("Done feeding, waiting for:");
                        for(int i=0;i<cmdarr.length;i++){
                            System.out.print("'"+cmdarr[i]+"' ");
                        }
                        System.out.println(" query="+query.name);
                    }
                    p.waitFor();
                    if(p.exitValue()!=0){
                        System.err.println("non-perfect exit from blast for "+query.name);
                        errbuff.append("non-perfect exit from blast for "+query.name+"\n");
                    }
                    synchronized(inout){//that's also what the threadstreamreader syncs by
                        while (pinread.done==false){
                            //System.out.println("waiting for process reader to finish");
                            try{
                                //sleep(10l);
                                inout.wait(10l);
                            }catch(InterruptedException e){
                                System.err.println("interrupted sleep in blastthread");
                                e.printStackTrace();
                            }
                        }
                    }//end synchronized inout
                }catch (InterruptedException e){
                    //System.err.println("Interrupted process "+blastcommand+" for "+query.name);
                    System.out.print("Interrupted process for:"+query.name+"; command=");
                    for(int i=0;i<cmdarr.length;i++){
                        System.out.print("'"+cmdarr[i]+"' ");
                    }
                    System.out.println();
                    stopthread=true;
                }
                perrread=null;//.clear();
                pinread=null;//.clear();
            }catch (IOException ioe){
                //System.err.println("IOError in "+blastcommand+" for "+query.name);
                System.err.print("IOERROR for "+query.name+" command:");
                errbuff.append("IOERROR for "+query.name+" command\n");
                for(int i=0;i<cmdarr.length;i++){
                    System.err.print("'"+cmdarr[i]+"';");
                }
                System.err.println();
                //for(int i=0;i<cmdarr.length;i++){
                //    System.err.println(cmdarr[i]);
                //}
                System.out.println("inoutsize="+inout.length());
                stopthread=true;
            }
            if(stopthread){
                System.err.println("stopped thread "+query.name);
                synchronized (parent){
                    this.done=true;
                    parent.notify();
                }// end synchronized parent
                return;
            }// end if stopthread
            System.err.print(errout.toString());
            //now the psiblast run should have generated a checkpoint file
            //then do the last round using the psiblast checkpoint file against blastdb
            if(verbose>1){
                System.out.println("actual psiblast run on "+threadnum);
            }

            tmpvec.addElement("-d");
            tmpvec.addElement(blastdbname);
            tmpvec.addElement("-R");
            tmpvec.addElement(tmpcheckfile);
        }
        cmdarr=new String[tmpvec.size()];
        tmpvec.copyInto(cmdarr);
        //if(cmd.length()>0){
        //    blastcommand=cmd+" "+blastcommand;
        //}
        StringBuffer inout=new StringBuffer();
        StringBuffer errout=new StringBuffer();
        if(verbose>2){
            //System.out.println("for "+query.name+" trying command="+blastcommand);
            System.out.print(query.name+" trying:");
            for(int i=0;i<cmdarr.length;i++){
                System.out.print("'"+cmdarr[i]+"' ");
            }
            System.out.println();
        }
        String myblast="";
        try{
            //write the query to a file instead of using stdin; !stdin gives problems!
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(tmpfilestring)));
            outwrite.println(">"+query.name);
            outwrite.println(query.seq.replaceAll("-",""));
            outwrite.close();
        }catch (IOException ioe){
            System.out.println("IOERROR writing to query file "+tmpfilestring);
        }
        //System.out.println("donw writing 2 for "+tmpfilestring);
        try{
            BufferedReader perr;
            BufferedReader pin;
            threadstreamreader perrread;
            threadstreamreader pinread;
            //Process p=rt.exec(blastcommand);
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
                if(verbose>2){
                    System.out.println("feeding seqdata >"+query.name);
                    System.out.println("seq: "+query.seq.replaceAll("-",""));
                }
                if(verbose>2){
                    //System.out.println("Done feeding, waiting for "+blastcommand+" query="+query.name);
                    System.out.print("Done feeding for "+query.name+" waiting for command:");
                    for(int i=0;i<cmdarr.length;i++){
                        System.out.print("'"+cmdarr[i]+"' ");
                    }
                    System.out.println();
                }
                p.waitFor();
                if(p.exitValue()!=0){
                    System.err.println("non-perfect exit from blast for "+query.name);
                    errbuff.append("non-perfect exit from blast for "+query.name+"\n");
                }
                if(verbose>0){
                    System.out.println(threadnum+": blast done for "+query.name);
                }
                synchronized(inout){
                    while (pinread.done==false){
                        //System.out.println("waiting for process reader to finish");
                        try{
                            //sleep(10l);
                            inout.wait(10l);
                        }catch(InterruptedException e){
                            System.err.println("interrupted sleep in blastthread");
                            e.printStackTrace();
                        }
                    }
                }//end synchronized inout
            }catch (InterruptedException e){
                //System.err.println("Interrupted process "+blastcommand+" for "+query.name);
                System.out.print("Interrupted process for "+query.name+" command:");
                errbuff.append("Interrupted process for "+query.name+" command:");
                for(int i=0;i<cmdarr.length;i++){
                    System.out.print("'"+cmdarr[i]+"' ");
                }
                System.out.println();
                stopthread=true;
            }
            myblast=inout.toString();
            perrread=null;//.clear();
            pinread=null;//.clear();
        }catch (IOException ioe){
            //System.err.println("IOError in "+blastcommand+" for "+query.name);
            //for(int i=0;i<cmdarr.length;i++){
            //    System.err.println(cmdarr[i]);
            //}
            System.out.print("IOERROR for running "+query.name+" command:");
            errbuff.append("IOERROR running "+query.name+" command\n");
            for(int i=0;i<cmdarr.length;i++){
                System.out.print("'"+cmdarr[i]+"' ");
            }
            System.out.println();
            stopthread=true;
        }
        if(stopthread){
            System.err.println("stopped thread "+query.name);
            synchronized (parent){
                this.done=true;
                parent.notify();
            }// end synchronized parent
            return;
        }// end if stopthread
        if(errout.indexOf("Error")>-1){
            System.err.print(errout.toString());
            errbuff.append("ERRORS:\n"+errout.toString()+"\n");
        }
        if(verbose>10){
            System.out.println(inout.toString());
        }
        if(verbose>2){
            System.out.println("getting hsp's from "+query.name);
        }
        Vector<hsp> hspvec;
        if (lowmem) {
            hspvec = hspget.get(lowmem, tmpoutfile, new Vector<hsp>(), eval, pval, coverage, scval, ident, verbose,
                    nameshash, useallrounds);// here the blast output is filtered for valid hsp's
        } else {
            hspvec = hspget.get(myblast, new Vector<hsp>(), eval, pval, coverage, scval, ident, verbose, nameshash,
                    useallrounds);// here the blast output is filtered for valid hsp's
        }
        
        if(verbose>1){
            System.out.println("hsp's="+hspvec.size());
        }
        inout=null;//free the memory from these stringbuffers (more or less)
        errout=null;
        //now I have all the hsp's in a vector.
        //next assign the hsp's by name of hit sequence to a position in the retarr
        int hspnum=hspvec.size();
        int arrsize=retarr.length;
        for(int i=0;i<arrsize;i++){
            retarr[i]=new Vector<hsp>();
        }
        int namenum=-1;
        for(int i=0;i<hspnum;i++){
            if(nameshash.containsKey(hspvec.elementAt(i).hname)){
                namenum=nameshash.get(hspvec.elementAt(i).hname).intValue();
                retarr[namenum].addElement(hspvec.elementAt(i));
            }else{
                System.err.println("unknown name for "+hspvec.elementAt(i).hname);
            }
        }
        //all done, notify parent to get the data
        if(verbose>2){
            System.out.println("notifying parent for "+query.name);
        }
        synchronized(parent){
            this.done=true;
            parent.notify();
        }// end synchronized parent
    }// end run
    
    //--------------------------------------------------------------------------
    
}
