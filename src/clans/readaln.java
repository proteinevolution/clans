/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clans;
import java.io.*;
import java.util.*;
/**
 *
 * @author tancred
 */
public class readaln {
    //this class should contain everything necessary to read alignment files and
    //output them as an array of aaseq objects.
    
    /** Creates a new instance of readaln */
    public readaln() {
    }
    
    //--------------------------------------------------------------------------
    
    public static aaseq[] read(String instring){
        //this converts a string to a filename and tries to read the seq alignment
        File infile=new File(instring);
        return read(infile);
    }//end read string
    
    public static aaseq[] read(File infile){
        try{
            BufferedReader inread;
            if(infile.getName().equalsIgnoreCase("STDIN")){
                inread=new BufferedReader(new InputStreamReader(System.in));
            }else{
                inread=new BufferedReader(new FileReader(infile));
            }
            //now get the first non space line. this should tell if fasta or clustal or stockholm
            String inline;
            while((inline=inread.readLine())!=null){
                inline=inline.trim();
                if(inline.length()>0){//if there is sth. on this line
                    if(inline.startsWith(">")){//if this is fasta format
                        inread.close();
                        return fastaread(infile);
                    }else if((inline.length()>6)&&(inline.substring(0,7).equalsIgnoreCase("clustal"))){
                        inread.close();
                        return clustalread(infile);
                    }else if((inline.length()>8)&&(inline.toUpperCase().matches(".*STOCKHOLM.*"))){//if this is stockholm format
                        inread.close();
                        return stockholmread(infile);
                    }else{//see if the line read has one or two numbers (if neither :unknown format)
                        String[] tmparr=inline.split("\\s+",0);
                        int arrsize=java.lang.reflect.Array.getLength(tmparr);
                        boolean didconv=true;
                        try{//try to convert the string to numbers
                            for(int i=0;i<arrsize;i++){
                                int tmpint=Integer.parseInt(tmparr[i]);
                            }
                        }catch (NumberFormatException e){
                            didconv=false;
                        }
                        if(didconv){
                            if(arrsize==1){//treecon
                                inread.close();
                                return treeconread(infile);
                            }else if(arrsize==2){//phylip
                                inread.close();
                                return phylipread(infile);
                            }
                        }
                        System.err.println("unknown file format for "+infile.getAbsolutePath());
                        return new aaseq[0];
                    }
                }
            }// end while
            //if I get here the file was empty
        }catch (IOException e){
            System.err.println("IOError reading "+infile.getAbsolutePath());
            e.printStackTrace();
        }
        System.err.println("empty file for "+infile.getAbsolutePath());
        return new aaseq[0];
    }//end read file
    
    //--------------------------------------------------------------------------
    
    public static aaseq[] fastaread(String infilename){
        File tmpfile=new File(infilename);
        return fastaread(tmpfile);
    }
    
    public static aaseq[] fastaread(File infile){
        //this should read fasta format data from infile and return it as an array
        //of aasseq objects (name,seq)
        Vector tmpvec=new Vector();
        try{
            BufferedReader inread;
            if(infile.getName().equalsIgnoreCase("STDIN")){
                System.out.println("reading sequence data from STDIN");
                inread=new BufferedReader(new InputStreamReader(System.in));
            }else{
                inread=new BufferedReader(new FileReader(infile));
            }
            String inline;
            StringBuffer seqbuff=new StringBuffer();
            aaseq myaaseq=new aaseq();
            while((inline=inread.readLine())!=null){
                inline=inline.trim();
                if(inline.startsWith(">")){//if this is a fasta name
                    if(seqbuff.length()>0){//if I have read a seq
                        myaaseq.seq=seqbuff.toString().toUpperCase();
                        seqbuff.setLength(0);
                        tmpvec.addElement(myaaseq);
                    }
                    myaaseq=new aaseq();
                    myaaseq.name=inline.substring(1);//skip the >
                }else{
                    seqbuff.append(inline);
                }
            }// end while
            if(seqbuff.length()>0){//if I have read a seq
                myaaseq.seq=seqbuff.toString();
                seqbuff.setLength(0);
                tmpvec.addElement(myaaseq);
            }
            inread.close();
        }catch (IOException e){
            System.err.println("IOError reading from "+infile.getAbsolutePath());
            e.printStackTrace();
            return new aaseq[0];
        }
        aaseq[] retarr=new aaseq[tmpvec.size()];
        tmpvec.copyInto(retarr);
        return retarr;
    }// end fastaread
    
    //--------------------------------------------------------------------------
    
    public static aaseq[] clustalread(File infile){
        //read a clustal alignment and convert to array of aaseq objects
        //this format is: CLUSTAL in the first line
        //name 'space(s)' seq on a line
        //name2 'space(s)' seq2 on one line
        //maybe a line with *. etc but no name
        //maybe an empty line
        //name 'space(s)' seq part 2 on a line
        //etc...
        HashMap myhash=new HashMap();
        Vector alnnames=new Vector();
        try{
            BufferedReader inread;
            if(infile.getName().equalsIgnoreCase("STDIN")){
                inread=new BufferedReader(new InputStreamReader(System.in));
            }else{
                inread=new BufferedReader(new FileReader(infile));
            }
            String inline;
            String inname;
            StringBuffer seqbuff=new StringBuffer();
            int namespace;
            while((inline=inread.readLine())!=null){
                inline=inline.trim();
                //skip the clustal file
                if(inline.length()==0){//if the line is empty
                    continue;
                }else if((inline.length()>6)&&(inline.substring(0,7).equalsIgnoreCase("clustal"))){//if it is the clustal line
                    continue;
                }else{
                    //here I read the seqnames and seqs (and opt. lines without AA or -), skip those
                    //name is beginning of line to first space char
                    if(inline.matches(".*[a-zA-Z-]+.*")==false){//if it has no a-z or gap "-"
                        continue;
                    }
                    if((namespace=inline.indexOf(" "))>-1){//if I can split this line in two at the first space char
                        inname=inline.substring(0,namespace);
                        //System.out.println("dataline! name="+inname);
                        if(myhash.containsKey(inname)){
                            ((StringBuffer)myhash.get(inname)).append((inline.substring(namespace)).trim());
                        }else{
                            myhash.put(inname,new StringBuffer((inline.substring(namespace)).trim()));
                            alnnames.addElement(inname);
                        }
                    }// end if I can split this seq at a space
                }// end else
            }// end while
            inread.close();
        }catch (IOException e){
            System.err.println("IOError reading from "+infile.getAbsolutePath());
            e.printStackTrace();
            return new aaseq[0];
        }
        //now take the elements of myhash and convert them to aaseq[]
        int keysnum=myhash.size();
        if(keysnum!=alnnames.size()){
            System.err.println("Unequal sequence number in readclustal; exiting");
            System.exit(0);
        }
        aaseq[] retarr=new aaseq[keysnum];
        String[] keys=(String[])(myhash.keySet().toArray(new String[0]));
        for(int i=0;i<keysnum;i++){
            retarr[i]=new aaseq();
            retarr[i].name=(String)alnnames.elementAt(i);
            retarr[i].seq=((StringBuffer)myhash.get(retarr[i].name)).toString().toUpperCase();
        }
        return retarr;
    }// end clustalread
    
    //--------------------------------------------------------------------------
    
    public static aaseq[] treeconread(File filenamein)
    //read in a treecon-format alignment.
    //format: one number on the first line
    //next a line with the name
    //following lines: the sequence
    //empty lines (opt.)
    //next name...
    throws IOException{
        String aaseqname;
        Vector seqvector = new Vector();
        aaseq seqobj;
        int seqlength, iaa;
        int count=0;
        char seqaa;
        int seqnumber=0;
        BufferedReader alignread;
        if(filenamein.getName().equalsIgnoreCase("STDIN")){
            alignread=new BufferedReader(new InputStreamReader(System.in));
        }else{
            alignread=new BufferedReader(new FileReader(filenamein));
        }
        try{
            seqlength = Integer.parseInt(alignread.readLine());
        }catch (NumberFormatException e){
            System.err.println("Invalid Alignment");
            aaseq[] errorseq = new aaseq[0];
            return (errorseq);
        }//return a dummy aaseq array
        while (((aaseqname=(alignread.readLine()))!=null)){
            if (((aaseqname.length())!=0)){
                count=0;
                seqnumber+=1;
                char [] seqseq = new char[seqlength];
                while (count<(seqlength)){
                    iaa=(alignread.read());
                    seqaa=(char)iaa;//was :convert.tochar(iaa);
                    if (Character.isJavaIdentifierPart(seqaa)||(seqaa=='-')||(seqaa=='?')||(seqaa=='_')){
                        if((Character.isLetter(seqaa))==false){
                            seqaa='-';
                        }
                        seqseq [count]= seqaa;
                        count += 1;
                    }//end of ifcharacterisjavaidentifierpart
                    if ((Character.isJavaIdentifierPart(seqaa)||(seqaa=='-')||(seqaa=='?')||(seqaa=='_')||(Character.isWhitespace(seqaa)))==false){
                        System.err.println("false sequence length or illegal character:"+seqaa+" at sequence "+aaseqname+" position="+count);
                    }
                }//end of while count<seqlenght
                alignread.mark(2);
                seqaa=(char)alignread.read();
                while((Character.isWhitespace(seqaa))){//if there was a newline or a space
                    alignread.mark(2);//save the current position
                    seqaa=(char)alignread.read();
                }//if another newline repeat loop
                alignread.reset();//if not, start reading after the last mark;
                aaseq currentread = new aaseq();
                currentread.name=aaseqname;
                currentread.seq=new String(seqseq);
                seqvector.addElement(currentread);
            }//end of if Stringlength!=0
        }//end of while aaseq!=null
        aaseq[] seqarrayout = new aaseq[seqvector.size()];
        seqvector.copyInto(seqarrayout);
        alignread.close();
        return seqarrayout;
    }//end of treeconread()
    
    //--------------------------------------------------------------------------
    
    public static aaseq[] phylipread(File infile){
        //this needs to decide wether I am reading phylip sequential or interleaved
        //physeq: first line : number of species + sequence length
        //opt. empty lines
        //name[10],seq[11-end].
        //opt. empty lines
        //next seq
        
        //phy interleaved
        //first line : number of spec + seqlength
        //opt. empty lines
        //name[10](opt. 1 space),seq[11-end];
        //next name[10],seq
        //... for rest of species
        //opt. one empty line
        //(opt. spaces[10 or 11]) seq1
        //seq2 etc... to end
        
        //the way I differentiate between phylip sequential and interleaved is to read in both in parallel
        //and see which one exits with an error first. if both read to end, return physeq; if neither reads
        //to end return aaseq[0] and print error
        
        //NOTE!!!: remove any space chars from within the sequence
        
        boolean physeq=true;
        boolean phyint=true;
        int specnum=0;
        int seqlength=0;
        Vector phyintvec=new Vector();
        Vector physeqvec=new Vector();
        aaseq[] retarr=new aaseq[0];
        try{
            BufferedReader inread;
            if(infile.getName().equalsIgnoreCase("STDIN")){
                inread=new BufferedReader(new InputStreamReader(System.in));
            }else{
                inread=new BufferedReader(new FileReader(infile));
            }
            String inline="";
            //read the first line and parse the numbers
            boolean readon=true;
            String[] tmparr=new String[0];
            while((readon)&&((inline=inread.readLine())!=null)){//while I can read and have not read the numbers
                inline=inline.trim();
                if(inline.length()>0){
                    tmparr=inline.split("\\s+",0);
                    readon=false;
                }
            }
            try{
                if(java.lang.reflect.Array.getLength(tmparr)==2){
                    specnum=Integer.parseInt(tmparr[0]);
                    seqlength=Integer.parseInt(tmparr[1]);
                }else{
                    System.err.println("Could not generate two numbers from "+inline);
                    return retarr;
                }
            }catch (NumberFormatException e){
                System.err.println("unable to parse species and seqlength numbers from phylip alignment");
                return retarr;
            }
            //if I get here I have successfully parsed specnum and seqlength from the first line of the file
            //now skip all empty lines and read the seqdata
            for(int i=0;i<specnum;i++){
                phyintvec.addElement(new aaseq());
                physeqvec.addElement(new aaseq());
            }
            int seqspecs=-1;
            int intspecs=-1;
            while(((inline=inread.readLine())!=null)&&(physeq||phyint)){//while I am reading from the file and either seq or int are still ok
                inline=inline.trim();
                if(inline.length()==0){
                    continue;
                }
                //I have a non-empty line and I know how many specs I have read
                //try sequential read
                if(physeq){
                    //I know which species I am at and I know what sequence position I am at
                    if(seqspecs>=specnum){
                        System.out.println("PHYSEQ=FALSE seqspecs>specnum");
                        physeq=false;
                    }else{
                        if(seqspecs>-1){
                            int currseqpos=((aaseq)physeqvec.elementAt(seqspecs)).seq.length();
                            if(currseqpos<seqlength){//if I have not read the full sequence yet
                                ((aaseq)physeqvec.elementAt(seqspecs)).seq+=inline;
                            }
                            if(currseqpos>seqlength){
                                System.err.println("sequence "+((aaseq)physeqvec.elementAt(seqspecs)).seq+" is longer than "+seqlength);
                                System.out.println("PHYSEQ=FALSE");
                                physeq=false;
                            }
                            if(currseqpos==seqlength){
                                //in this case the last sequence is fully accounted for and this line contains a species name in the first 10 characters
                                seqspecs++;
                                ((aaseq)physeqvec.elementAt(seqspecs)).name=inline.substring(0,10).trim();
                                ((aaseq)physeqvec.elementAt(seqspecs)).seq=(inline.substring(10).trim()).replaceAll(" ","");
                            }
                        }else{//if this is the first seqpart I read:
                            seqspecs++;
                            ((aaseq)physeqvec.elementAt(seqspecs)).name=inline.substring(0,10).trim();
                            ((aaseq)physeqvec.elementAt(seqspecs)).seq=(inline.substring(10).trim()).replaceAll(" ","");
                        }
                    }
                }
                //try interleaved read
                if(phyint){
                    // I know what species I am at and where in the sequence I am.
                    if(intspecs<(specnum-1)){
                        //if I am reading the first block
                        intspecs++;
                        ((aaseq) phyintvec.elementAt(intspecs)).name=inline.substring(0,10).trim();
                        ((aaseq) phyintvec.elementAt(intspecs)).seq=(inline.substring(10).trim()).replaceAll(" ","");
                    }else{
                        //if i am reading a subsequent block
                        //do not read names (no need to)
                        //get the real species number
                        intspecs++;
                        int currspec=intspecs;
                        while (currspec>=specnum){
                            currspec-=specnum;
                        }
                        ((aaseq) phyintvec.elementAt(currspec)).seq+=(inline.trim()).replaceAll(" ","");
                        if(((aaseq) phyintvec.elementAt(currspec)).seq.length()>seqlength){
                            System.err.println("sequence "+((aaseq) phyintvec.elementAt(currspec)).name+" is longer than "+seqlength);
                            System.err.println(((aaseq) phyintvec.elementAt(currspec)).seq);
                            System.out.println("PHYINT=FALSE");
                            phyint=false;
                        }
                    }
                }//end phyint
                //done both, read the next line
            }// end while readline
            
        }catch (IOException e){
            System.err.println("IOError in reading from "+infile.getName());
            return retarr;
        }
        if(physeq){//if phylip sequential is still true at the end of the file
            retarr=new aaseq[physeqvec.size()];
            physeqvec.copyInto(retarr);
            return retarr;
        }else if(phyint){
            retarr=new aaseq[phyintvec.size()];
            phyintvec.copyInto(retarr);
            return retarr;
        }else{
            System.err.println("unable to read Phylip alignment from "+infile.getName());
        }
        return retarr;
    }// end phylipread
    
    //--------------------------------------------------------------------------
    
    public static aaseq[] stockholmread(File infile){
        aaseq[] retarr=new aaseq[0];
        Vector seqvec=new Vector();
        boolean doneall=false;
        try{
            BufferedReader inread;
            if(infile.getName().equalsIgnoreCase("STDIN")){
                inread=new BufferedReader(new InputStreamReader(System.in));
            }else{
                inread=new BufferedReader(new FileReader(infile));
            }
            String inline;
            aaseq curr;
            while (((inline=inread.readLine())!=null)&&(doneall==false)){
                //any line with # is not of interest
                //the non-# lines (except for the terminating //) contain the sequences I want to read
                //name[29] seq[to end of line]
                if(inline.startsWith("#")){
                    continue;
                }else{
                    if(inline.startsWith("//")){
                        //if this is the termination line
                        doneall=true;
                        continue;
                    }
                    curr=new aaseq();
                    curr.name=inline.substring(0,29).trim();
                    curr.seq=inline.substring(29).trim();
                    seqvec.addElement(curr);
                }
            }// end while readline
            inread.close();
        }catch(IOException e){
            System.err.println("IOError in stockholmread for "+infile.getName());
            return retarr;
        }
        retarr=new aaseq[seqvec.size()];
        seqvec.copyInto(retarr);
        return retarr;
    }//end stockholmread
    
    //--------------------------------------------------------------------------
    
}
