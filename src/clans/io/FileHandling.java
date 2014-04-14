package clans.io;

import java.io.*;
import java.util.*;

import clans.model.proteins.AminoAcidSequence;
import clans.model.proteins.HighScoringSegmentPair;

public class FileHandling {

    /**
     * print the working alignment to file outfilename in the specified format
     * 
     * @param inaln
     * @param outfilename
     * @param oformat
     * @param blocksize
     */
    static void printaln(AminoAcidSequence[] inaln,String outfilename, String oformat, int blocksize){
        if(oformat.equalsIgnoreCase("clustal")){
            printclustal(inaln,outfilename,blocksize);
        }else{
            printfasta(inaln,outfilename);
        }
    }
    
    /**
     * print the position specific confidences for each sequence
     * 
     * @param inaln
     * @param confarr
     * @param outfilename
     * @param oformat
     * @param blocksize
     */
    static void printconf(AminoAcidSequence[] inaln, double[][] confarr, String outfilename, String oformat,
            int blocksize) {
       if(oformat.equalsIgnoreCase("clustal")){
            printclustal(inaln,confarr,outfilename,blocksize);
        }else{
            printfasta(inaln,confarr,outfilename);
        }
    }

    /**
     * print the input alignment filtered for confidence. replace low confidence residues by dashes
     * 
     * @param inaln
     * @param confarr
     * @param outfilename
     * @param oformat
     * @param blocksize
     * @param filtercutoff
     */
    static void printfiltered(AminoAcidSequence[] inaln, double[][] confarr, String outfilename, String oformat,
            int blocksize, float filtercutoff) {
        int elements=inaln.length;
        int seqlength=inaln[0].seq.length();
        AminoAcidSequence[] tmpaln=new AminoAcidSequence[elements];
        for(int i=0;i<elements;i++){
            tmpaln[i]=new AminoAcidSequence();
            tmpaln[i].name=inaln[i].name;
            tmpaln[i].seq=inaln[i].seq;
        }//end for i
        int j;
        StringBuffer tmpstr=new StringBuffer(blocksize+3);
        for(int i=0;i<elements;i++){
            tmpstr.setLength(0);
            for(j=0;j<seqlength;j++){
                if(confarr[i][j]>=filtercutoff){
                    tmpstr.append(tmpaln[i].seq.charAt(j));
                }else{
                    tmpstr.append("-");
                }
            }//end for j
            tmpaln[i].seq=tmpstr.toString();
        }//end for i
        if(oformat.equalsIgnoreCase("clustal")){
            printclustal(tmpaln,outfilename,blocksize);
        }else{
            printfasta(tmpaln,outfilename);
        }
    }
    
    /**
     * print the alignment in fasta format
     * 
     * @param inaln
     * @param outfilename
     */
    static void printfasta(AminoAcidSequence[] inaln, String outfilename) {
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfilename)));
            int elements=inaln.length;
            for(int i=0;i<elements;i++){
                outwrite.println(">"+inaln[i].name);
                outwrite.println(inaln[i].seq);
            }
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR in printfiltered");
            e.printStackTrace();
        }
    }
    
    /**
     * print the alignment in fasta format
     * 
     * @param inaln
     * @param outfile
     */
    public static void printfasta(AminoAcidSequence[] inaln, File outfile) {
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
            int elements=inaln.length;
            for(int i=0;i<elements;i++){
                outwrite.println(">"+inaln[i].name);
                outwrite.println(inaln[i].seq);
            }
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR in printfiltered");
            e.printStackTrace();
        }
    }

    /**
     * print the confarr in fasta format
     * 
     * @param inaln
     * @param confarr
     * @param outfilename
     */
    static void printfasta(AminoAcidSequence[] inaln, double[][] confarr, String outfilename) {
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfilename)));
            int elements=inaln.length;
            if(elements<1){
                System.err.println("ERROR: empty alignment");
                outwrite.close();
                return;
            }
            int seqlength=inaln[0].seq.length();
            int j;

            for(int i=0;i<elements;i++){
                outwrite.println(">"+inaln[i].name);
                for(j=0;j<seqlength;j++){
                    outwrite.print(((float)(int)(confarr[i][j]*100))/100+";");
                }//end for j
                outwrite.println();
            }//end for i
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR in printfiltered");
            e.printStackTrace();
        }
    }
    
    /**
     * print the alignment in clustal format
     * 
     * @param inaln
     * @param outfilename
     * @param blocksize
     */
    static void printclustal(AminoAcidSequence[] inaln, String outfilename, int blocksize){
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfilename)));
            int elements=inaln.length;
            if(elements<1){
                System.err.println("ERROR: empty alignment");
                outwrite.close();
                return;
            }
            int position=0;
            int seqlength=inaln[0].seq.length();
            int maxnamelength=10;
            int i;
            for(i=0;i<elements;i++){
                if(inaln[i].name.length()>maxnamelength){
                    maxnamelength=inaln[i].name.length();
                }
            }//end for i
            StringBuffer tmpstr=new StringBuffer(maxnamelength+3);
            String[] namearr=new String[elements];
            for(i=0;i<elements;i++){
                tmpstr.setLength(0);
                tmpstr.append(inaln[i].name);
                while (tmpstr.length()<maxnamelength){
                    tmpstr.append(" ");
                }
                namearr[i]=tmpstr.toString();
            }//end for i
            outwrite.println("CLUSTAL");
            outwrite.println();
            outwrite.println();
            while ((position+blocksize)<=seqlength){
                for(i=0;i<elements;i++){
                    outwrite.println(namearr[i]+"  "+inaln[i].seq.substring(position,position+blocksize));
                }//end for i
                outwrite.println();
                position+=blocksize;
            }//end while
            for(i=0;i<elements;i++){
                tmpstr.setLength(0);
                tmpstr.append(inaln[i].name);
                while (tmpstr.length()<maxnamelength){
                    tmpstr.append(" ");
                }
                outwrite.println(tmpstr.toString()+"  "+inaln[i].seq.substring(position));
            }//end for i
            outwrite.println();
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR in printfiltered");
            e.printStackTrace();
        }
    }
    
    /**
     * print the confarr in clustal format
     * 
     * @param inaln
     * @param confarr
     * @param outfilename
     * @param blocksize
     */
    static void printclustal(AminoAcidSequence[] inaln, double[][] confarr, String outfilename, int blocksize) {

        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfilename)));
            int elements=inaln.length;
            if(elements<1){
                System.err.println("ERROR: empty alignment");
                outwrite.close();
                return;
            }
            int position=0;
            int seqlength=inaln[0].seq.length();
            int maxnamelength=10;
            int i,j;
            for(i=0;i<elements;i++){
                if(inaln[i].name.length()>maxnamelength){
                    maxnamelength=inaln[i].name.length();
                }
            }//end for i
            StringBuffer tmpstr=new StringBuffer(maxnamelength+3);
            String[] namearr=new String[elements];
            for(i=0;i<elements;i++){
                tmpstr.setLength(0);
                tmpstr.append(inaln[i].name);
                while (tmpstr.length()<maxnamelength){
                    tmpstr.append(" ");
                }
                namearr[i]=tmpstr.toString();
            }//end for i
            outwrite.println("CLUSTAL");
            outwrite.println();
            outwrite.println();
            while ((position+blocksize)<=seqlength){
                for(i=0;i<elements;i++){
                    outwrite.print(namearr[i]+"  ");
                    for(j=position;j<position+blocksize;j++){
                        outwrite.print(((float)(int)(confarr[i][j]*100))/100+";");
                    }//end for j
                    outwrite.println();
                }//end for i
                outwrite.println();
                position+=blocksize;
            }//end while
            for(i=0;i<elements;i++){
                outwrite.print(namearr[i]+"  ");
                for(j=position;j<seqlength;j++){
                    outwrite.print(confarr[i][j]+";");
                }//end for j
                outwrite.println();
            }//end for i
            outwrite.println();
            outwrite.close();
        }catch (IOException e){
            System.err.println("IOERROR in printfiltered");
            e.printStackTrace();
        }
    }
    
    /**
     * 
     * @param outfilename
     * @param blastvals
     */
    static void saveblast(String outfilename, Vector<HighScoringSegmentPair>[][] blastvals){
        try{
            PrintWriter outwrite=new PrintWriter(new BufferedWriter(new FileWriter(outfilename)));
            int vecsize=blastvals.length;
            int i,j,k;
            int tmpvecsize;
            Vector<HighScoringSegmentPair> tmpvec;
            HighScoringSegmentPair tmphsp;
            outwrite.println(vecsize+" sequences");
            for(i=0;i<vecsize;i++){
                for(j=0;j<vecsize;j++){
                    //print the hsp data for this sequence comparison
                    tmpvec=blastvals[i][j];
                    tmpvecsize=tmpvec.size();
                    for(k=0;k<tmpvecsize;k++){
                        tmphsp = tmpvec.elementAt(k);
                        outwrite.println("hsp: "+i+";"+j+";"+k);
                        outwrite.println("\tvalue: "+tmphsp.value);
                    }
                }
            }
            outwrite.close();
            
        }catch (IOException e){
            System.out.println("Error printing blast results to "+outfilename);
        }
    }
    
    /**
     * 
     * @param outwrite
     * @param blastvals
     * @param i
     */
    static void saveblastappend(PrintWriter outwrite, Vector<HighScoringSegmentPair>[] blastvals, int i) {
        int vecsize=blastvals.length;
        int j,k;
        int tmpvecsize;
        Vector<HighScoringSegmentPair> tmpvec;
        HighScoringSegmentPair tmphsp;
        for(j=0;j<vecsize;j++){
            //print the hsp data for this sequence comparison
            tmpvec=blastvals[j];
            if(tmpvec!=null){
                tmpvecsize=tmpvec.size();
                for(k=0;k<tmpvecsize;k++){
                    tmphsp = tmpvec.elementAt(k);
                    outwrite.println("hsp: "+i+";"+j+";"+k);
                    outwrite.println("\tvalue: "+tmphsp.value);
                }
            }
        }
    }

    /**
     * 
     * @param outwrite
     * @param positions
     */
    public static void save_semicolon_delimited_positions(PrintWriter outwrite, float[][] positions) {
        outwrite.println("positions:");
        for (int i = 0; i < positions.length; i++) {
            outwrite.println(i + ";" + positions[i][0] + ";" + positions[i][1] + ";" + positions[i][2]);
        }
        outwrite.flush();
    }
}