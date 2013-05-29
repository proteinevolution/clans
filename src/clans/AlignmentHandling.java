package clans;

import java.io.*;
import java.util.*;

/**
 * this class should contain everything necessary to read alignment files and output them as an array of aaseq objects.
 */
public class AlignmentHandling {

	public static AminoAcidSequence[] read(String instring) {
		File infile = new File(instring);
		return read(infile);
	}

	/**
	 * 
	 * @param infile
	 * @return
	 */
    public static AminoAcidSequence[] read(File infile){
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
                        return parse_fasta_format(infile);
                    }else if((inline.length()>6)&&(inline.substring(0,7).equalsIgnoreCase("clustal"))){
                        inread.close();
                        return parse_clustal_format(infile);
                    }else if((inline.length()>8)&&(inline.toUpperCase().matches(".*STOCKHOLM.*"))){//if this is stockholm format
                        inread.close();
                        return parse_stockholm_format(infile);
                    }else{//see if the line read has one or two numbers (if neither :unknown format)
                        String[] tmparr=inline.split("\\s+",0);
                        int arrsize=tmparr.length;
                        boolean didconv=true;
                        try{//try to convert the string to numbers
                            for(int i=0;i<arrsize;i++){
                            }
                        }catch (NumberFormatException e){
                            didconv=false;
                        }
                        if(didconv){
                            if(arrsize==1){//treecon
                                inread.close();
                                return parse_treecon_format(infile);
                            }else if(arrsize==2){//phylip
                                inread.close();
                                return parse_phylip_format(infile);
                            }
                        }
                        System.err.println("unknown file format for "+infile.getAbsolutePath());
                        return new AminoAcidSequence[0];
                    }
                }
            }// end while
            //if I get here the file was empty
        }catch (IOException e){
            System.err.println("IOError reading "+infile.getAbsolutePath());
            e.printStackTrace();
        }
        System.err.println("empty file for "+infile.getAbsolutePath());
        return new AminoAcidSequence[0];
    }
    
    /**
     * 
     * @param infilename
     * @return
     */
    public static AminoAcidSequence[] parse_fasta_format(String infilename) {
        File tmpfile = new File(infilename);
        return parse_fasta_format(tmpfile);
    }
    
    /**
     * this should read fasta format data from infile and return it as an array of aasseq objects (name,seq)
     * 
     * @param infile
     * @return
     */
    public static AminoAcidSequence[] parse_fasta_format(File infile) {
		Vector<AminoAcidSequence> tmpvec = new Vector<AminoAcidSequence>();
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
            AminoAcidSequence myaaseq=new AminoAcidSequence();
            while((inline=inread.readLine())!=null){
                inline=inline.trim();
                if(inline.startsWith(">")){//if this is a fasta name
                    if(seqbuff.length()>0){//if I have read a seq
                        myaaseq.seq=seqbuff.toString().toUpperCase();
                        seqbuff.setLength(0);
                        tmpvec.addElement(myaaseq);
                    }
                    myaaseq=new AminoAcidSequence();
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
            return new AminoAcidSequence[0];
        }
        AminoAcidSequence[] retarr=new AminoAcidSequence[tmpvec.size()];
        tmpvec.copyInto(retarr);
        return retarr;
    }

    /**
     * read a clustal alignment and convert to array of aaseq objects
     * 
     * @param infile
     * @return
     */
    public static AminoAcidSequence[] parse_clustal_format(File infile) {
        //this format is: CLUSTAL in the first line
        //name 'space(s)' seq on a line
        //name2 'space(s)' seq2 on one line
        //maybe a line with *. etc but no name
        //maybe an empty line
        //name 'space(s)' seq part 2 on a line
        //etc...
        
		HashMap<String, StringBuffer> myhash = new HashMap<String, StringBuffer>();
		Vector<String> alnnames = new Vector<String>();
        try{
            BufferedReader inread;
            if(infile.getName().equalsIgnoreCase("STDIN")){
                inread=new BufferedReader(new InputStreamReader(System.in));
            }else{
                inread=new BufferedReader(new FileReader(infile));
            }
            String inline;
            String inname;
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
            return new AminoAcidSequence[0];
        }
        //now take the elements of myhash and convert them to aaseq[]
        int keysnum=myhash.size();
        if(keysnum!=alnnames.size()){
            System.err.println("Unequal sequence number in readclustal; exiting");
            System.exit(0);
        }
        AminoAcidSequence[] retarr=new AminoAcidSequence[keysnum];
        for(int i=0;i<keysnum;i++){
            retarr[i]=new AminoAcidSequence();
            retarr[i].name=(String)alnnames.elementAt(i);
            retarr[i].seq=((StringBuffer)myhash.get(retarr[i].name)).toString().toUpperCase();
        }
        return retarr;
    }

    /**
     * read in a treecon-format alignment.
     * 
     * @param filenamein
     * @return
     * @throws IOException
     */
    public static AminoAcidSequence[] parse_treecon_format(File filenamein) throws IOException {
        //format: one number on the first line
        //next a line with the name
        //following lines: the sequence
        //empty lines (opt.)
        //next name...
        String aaseqname;
		Vector<AminoAcidSequence> seqvector = new Vector<AminoAcidSequence>();
        int seqlength, iaa;
        int count=0;
        char seqaa;
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
            AminoAcidSequence[] errorseq = new AminoAcidSequence[0];
            return (errorseq);
        }//return a dummy aaseq array
        while (((aaseqname=(alignread.readLine()))!=null)){
            if (((aaseqname.length())!=0)){
                count=0;
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
                AminoAcidSequence currentread = new AminoAcidSequence();
                currentread.name=aaseqname;
                currentread.seq=new String(seqseq);
                seqvector.addElement(currentread);
            }//end of if Stringlength!=0
        }//end of while aaseq!=null
        AminoAcidSequence[] seqarrayout = new AminoAcidSequence[seqvector.size()];
        seqvector.copyInto(seqarrayout);
        alignread.close();
        return seqarrayout;
    }

    /**
     * Parse a Phylip sequential or interleaved file
     * 
     * @param infile
     * @return
     */
    public static AminoAcidSequence[] parse_phylip_format(File infile) {
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
		Vector<AminoAcidSequence> phyintvec = new Vector<AminoAcidSequence>();
		Vector<AminoAcidSequence> physeqvec = new Vector<AminoAcidSequence>();
        AminoAcidSequence[] retarr=new AminoAcidSequence[0];
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
                if(tmparr.length==2){
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
                phyintvec.addElement(new AminoAcidSequence());
                physeqvec.addElement(new AminoAcidSequence());
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
                            int currseqpos=((AminoAcidSequence)physeqvec.elementAt(seqspecs)).seq.length();
                            if(currseqpos<seqlength){//if I have not read the full sequence yet
                                ((AminoAcidSequence)physeqvec.elementAt(seqspecs)).seq+=inline;
                            }
                            if(currseqpos>seqlength){
                                System.err.println("sequence "+((AminoAcidSequence)physeqvec.elementAt(seqspecs)).seq+" is longer than "+seqlength);
                                System.out.println("PHYSEQ=FALSE");
                                physeq=false;
                            }
                            if(currseqpos==seqlength){
                                //in this case the last sequence is fully accounted for and this line contains a species name in the first 10 characters
                                seqspecs++;
                                ((AminoAcidSequence)physeqvec.elementAt(seqspecs)).name=inline.substring(0,10).trim();
                                ((AminoAcidSequence)physeqvec.elementAt(seqspecs)).seq=(inline.substring(10).trim()).replaceAll(" ","");
                            }
                        }else{//if this is the first seqpart I read:
                            seqspecs++;
                            ((AminoAcidSequence)physeqvec.elementAt(seqspecs)).name=inline.substring(0,10).trim();
                            ((AminoAcidSequence)physeqvec.elementAt(seqspecs)).seq=(inline.substring(10).trim()).replaceAll(" ","");
                        }
                    }
                }
                //try interleaved read
                if(phyint){
                    // I know what species I am at and where in the sequence I am.
                    if(intspecs<(specnum-1)){
                        //if I am reading the first block
                        intspecs++;
                        ((AminoAcidSequence) phyintvec.elementAt(intspecs)).name=inline.substring(0,10).trim();
                        ((AminoAcidSequence) phyintvec.elementAt(intspecs)).seq=(inline.substring(10).trim()).replaceAll(" ","");
                    }else{
                        //if i am reading a subsequent block
                        //do not read names (no need to)
                        //get the real species number
                        intspecs++;
                        int currspec=intspecs;
                        while (currspec>=specnum){
                            currspec-=specnum;
                        }
                        ((AminoAcidSequence) phyintvec.elementAt(currspec)).seq+=(inline.trim()).replaceAll(" ","");
                        if(((AminoAcidSequence) phyintvec.elementAt(currspec)).seq.length()>seqlength){
                            System.err.println("sequence "+((AminoAcidSequence) phyintvec.elementAt(currspec)).name+" is longer than "+seqlength);
                            System.err.println(((AminoAcidSequence) phyintvec.elementAt(currspec)).seq);
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
            retarr=new AminoAcidSequence[physeqvec.size()];
            physeqvec.copyInto(retarr);
            return retarr;
        }else if(phyint){
            retarr=new AminoAcidSequence[phyintvec.size()];
            phyintvec.copyInto(retarr);
            return retarr;
        }else{
            System.err.println("unable to read Phylip alignment from "+infile.getName());
        }
        return retarr;
    }

    /**
     * Parse a file in Stockholm format.
     * @param infile
     * @return
     */
    public static AminoAcidSequence[] parse_stockholm_format(File infile){
        AminoAcidSequence[] retarr=new AminoAcidSequence[0];
		Vector<AminoAcidSequence> seqvec = new Vector<AminoAcidSequence>();
        boolean doneall=false;
        try{
            BufferedReader inread;
            if(infile.getName().equalsIgnoreCase("STDIN")){
                inread=new BufferedReader(new InputStreamReader(System.in));
            }else{
                inread=new BufferedReader(new FileReader(infile));
            }
            String inline;
            AminoAcidSequence curr;
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
                    curr=new AminoAcidSequence();
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
        retarr=new AminoAcidSequence[seqvec.size()];
        seqvec.copyInto(retarr);
        return retarr;
    }
}