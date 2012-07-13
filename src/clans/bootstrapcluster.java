/*
 * bootstrapcluster.java
 *
 * Created on May 17, 2004, 5:00 PM
 */
package clans;
import java.util.*;
/**
 *
 * @author  tancred
 */
public class bootstrapcluster {
    
    /** Creates a new instance of bootstrapcluster */
    public bootstrapcluster() {
    }
    
    static boolean bootstrapconvex(minattvals[] dataarr, Vector clustervec, String clustermethod, int replicates, float remove, float sigmafac, int minseqnum, int elements){
        //remove "remove" of data and then recluster for "replicates" replicates.
        //then get the cluster confidences AND sequence to cluster confidences
        Vector[] replicate=new Vector[replicates];
        int attnum=java.lang.reflect.Array.getLength(dataarr);
        Vector tmpdata=new Vector();
        Random rand=new Random(System.currentTimeMillis());
        minattvals[] tmparr;
        for(int i=0;i<replicates;i++){
            System.out.println("replicate "+i);
            //first remove "remove" fraction of data
            for(int j=0;j<attnum;j++){
                if(rand.nextFloat()>=remove){
                tmpdata.addElement(dataarr[j]);
                }
            }//end for j
            tmparr=new minattvals[tmpdata.size()];
            tmpdata.copyInto(tmparr);
            tmpdata.clear();
            //now get the cluster for this data
            replicate[i]=findclusters.getconvex(tmparr,sigmafac,minseqnum,elements);
        }//end for i
        //now compare the replicate clusters to the original
        checkreplicates(clustervec,replicate);
        return true;
    }//end bootstrap
    
    //--------------------------------------------------------------------------
    
    static boolean bootstraplinkage(minattvals[] dataarr, Vector clustervec, String clustermethod, int replicates, float remove, int minlinks, int minseqnum,int elements){
        //remove "remove" of data and then recluster for "replicates" replicates.
        //then get the cluster confidences AND sequence to cluster confidences
        Vector[] replicate=new Vector[replicates];
        int attnum=java.lang.reflect.Array.getLength(dataarr);
        Vector tmpdata=new Vector();
        Random rand=new Random(System.currentTimeMillis());
        minattvals[] tmparr;
        for(int i=0;i<replicates;i++){
            System.out.println("replicate "+i);
            //first remove "remove" fraction of data
            for(int j=0;j<attnum;j++){
                if(rand.nextFloat()>=remove){
                tmpdata.addElement(dataarr[j]);
                }
            }//end for j
            tmparr=new minattvals[tmpdata.size()];
            tmpdata.copyInto(tmparr);
            tmpdata.clear();
            //now get the cluster for this data
            replicate[i]=findclusters.getlinkage(tmparr,minlinks,minseqnum,elements);
        }//end for i
        //now compare the replicate clusters to the original
        checkreplicates(clustervec,replicate);
        return true;
    }//end bootstrap
    
    //--------------------------------------------------------------------------
    
    static boolean bootstrapnetwork(minattvals[] dataarr, Vector <cluster> clustervec, String clustermethod, int replicates, float remove, int minseqnum, boolean dooffset,boolean globalaverage, int elements,int maxrounds){
        //remove "remove" of data and then recluster for "replicates" replicates.
        //then get the cluster confidences AND sequence to cluster confidences
        Vector[] replicate=new Vector[replicates];
        int attnum=java.lang.reflect.Array.getLength(dataarr);
        Vector tmpdata=new Vector();
        Random rand=new Random(System.currentTimeMillis());
        minattvals[] tmparr;
        for(int i=0;i<replicates;i++){
            System.out.println("replicate "+i);
            //first remove "remove" fraction of data
            for(int j=0;j<attnum;j++){
                if(rand.nextFloat()>=remove){
                tmpdata.addElement(dataarr[j]);
                }
            }//end for j
            tmparr=new minattvals[tmpdata.size()];
            tmpdata.copyInto(tmparr);
            tmpdata.clear();
            //now get the cluster for this data
            replicate[i]=findclusters.getnetwork(tmparr,minseqnum,dooffset,globalaverage,elements,maxrounds);
            System.out.println("clusternum="+replicate[i].size());
            for(int j=replicate[i].size()-1;j>=0;j--){
                System.out.println("\t seqs="+java.lang.reflect.Array.getLength(((cluster)replicate[i].elementAt(j)).member));
            }//end for j
        }//end for i
        //now compare the replicate clusters to the original
        checkreplicates(clustervec,replicate);
        return true;
    }//end bootstrap
    
    //--------------------------------------------------------------------------
    
    static void checkreplicates(Vector clustervec, Vector[] replicate){
        System.out.println("comparing the replicates");
        //for this just check the clusters for how often they appear EXACTLY in the replicates
        int replicates=java.lang.reflect.Array.getLength(replicate);
        int clusternum=clustervec.size();
        cluster currcluster,checkcluster;
        int replicatessize,i,j,k,l,m;
        boolean foundcluster;
        String checkstring;
        //first sort all the members in all clusters in the same fashion
        for(i=0;i<replicates;i++){
            replicatessize=replicate[i].size();
            for(j=0;j<replicatessize;j++){
                java.util.Arrays.sort(((cluster)replicate[i].elementAt(j)).member);
            }//end for j
        }//end for i
        for(i=0;i<clusternum;i++){
            currcluster=(cluster)clustervec.elementAt(i);
            java.util.Arrays.sort(currcluster.member);
            currcluster.clusterconfidence=0;
            checkstring=makecheckstring(currcluster);
            for(j=0;j<replicates;j++){
                foundcluster=false;
                replicatessize=replicate[j].size();
                for(k=0;((k<replicatessize)&&(foundcluster==false));k++){
                    //now see if any of these clusters matches the currcluster elements exactly
                    checkcluster=(cluster)replicate[j].elementAt(k);
                    if(checkstring.equals(makecheckstring(checkcluster))){
                        foundcluster=true;
                    }
                }//end for k
                if(foundcluster){
                    currcluster.clusterconfidence++;
                }
            }//end for j
            currcluster.clusterconfidence/=replicates;
        }//end for i
        getseqconfidences(clustervec,replicate);
    }//end checkreplicates
    
    //--------------------------------------------------------------------------
    
    static String makecheckstring(cluster incluster){
        //convert the cluster elements into a sorted string
        StringBuffer tmp=new StringBuffer();
        int elements=java.lang.reflect.Array.getLength(incluster.member);
        for(int i=0;i<elements;i++){
            tmp.append(incluster.member[i]+";");
        }//end for i
        return tmp.toString();
    }//end makecheckstring
    
    //--------------------------------------------------------------------------
    
    static void getseqconfidences(Vector clustervec, Vector[] replicate){
        System.out.println("calculating sequence confidences");
        //get the confidence with which each sequence is assigned to each cluster
        //sum how often each cluster sequence pair is present.
        int clusternum=clustervec.size();
        int replicates=java.lang.reflect.Array.getLength(replicate);
        int elements,seqnum,checkseqnum,i,j,k,l,m,n;
        float sharednum,nonsharednum;
        cluster currcluster,checkcluster;
        boolean hasseq;
        for(i=0;i<clusternum;i++){//for each cluster
            currcluster=(cluster)clustervec.elementAt(i);
            seqnum=java.lang.reflect.Array.getLength(currcluster.member);
            currcluster.seqconfidence=new float[seqnum];
            for(j=0;j<seqnum;j++){
                currcluster.seqconfidence[j]=0;
            }//end for j
            //now see how often each sequence in this cluster is recovered as belonging to the same group in the replicates
            for(j=0;j<replicates;j++){
                for(k=replicate[j].size()-1;k>=0;k--){
                    sharednum=0;
                    nonsharednum=0;
                    checkcluster=(cluster)replicate[j].elementAt(k);
                    checkseqnum=java.lang.reflect.Array.getLength(checkcluster.member);
                    //now see which and how many sequences these clusters share
                    for(l=0;l<seqnum;l++){
                        hasseq=false;
                        for(m=0;m<checkseqnum;m++){
                            if(checkcluster.member[m]==currcluster.member[l]){//if they share this sequence
                                hasseq=true;
                                break;
                            }
                        }//end for m
                        if(hasseq){
                            sharednum++;
                        }else{
                            nonsharednum++;
                        }
                    }//end for l
                    //now I know how many sequences these two clusters share, now add the number of shareds to the seqs present in both.
                    if(sharednum>0){
                        for(l=0;l<seqnum;l++){
                            for(m=0;m<checkseqnum;m++){
                                if(checkcluster.member[m]==currcluster.member[l]){//if they share this sequence
                                    currcluster.seqconfidence[l]+=((sharednum-1)/seqnum);//fraction of sequences they shared
                                    currcluster.seqconfidence[l]-=(nonsharednum/checkseqnum);//remove fraction of sequences they did not share
                                }
                            }//end for m
                        }//end for l
                    }
                }//end for k
            }//end for j
            System.out.println("cluster "+i);
            for(j=0;j<seqnum;j++){
                currcluster.seqconfidence[j]/=replicates;
                System.out.println("\t"+j+" "+currcluster.member[j]+" "+currcluster.seqconfidence[j]);
            }//end for j
        }//end for i
    }//end getseqconfidences
    
    //--------------------------------------------------------------------------
    
    static void checkreplicatesold(Vector clustervec, Vector[] replicate){
        //check the bootstrap replicates
        int replicates=java.lang.reflect.Array.getLength(replicate);
        int clusternum=clustervec.size();
        int replicateclusternum;
        HashMap clusterseqshash=new HashMap();
        int sequences,tmpint,j,k,l,m,n;
        int[] tmpintarr;
        float[] tmpfloatarr;
        String tmpname;
        checkobject mycheck;
        //now make a quick lookup for which sequences were defined in which cluster
        for(int i=0;i<clusternum;i++){
            sequences=java.lang.reflect.Array.getLength(((cluster)clustervec.elementAt(i)).member);
            for(j=0;j<sequences;j++){
                tmpname=String.valueOf(((cluster)clustervec.elementAt(i)).member[j]);
                if(clusterseqshash.containsKey(tmpname)){
                    mycheck=(checkobject)clusterseqshash.get(tmpname);
                    tmpintarr=mycheck.clusters;
                    tmpfloatarr=mycheck.sumvals;
                    tmpint=java.lang.reflect.Array.getLength(tmpintarr);
                    mycheck.clusters=new int[tmpint+1];
                    mycheck.sumvals=new float[tmpint+1];
                    mycheck.entries=new int[tmpint+1];
                    for(k=0;k<tmpint;k++){
                        mycheck.clusters[k]=tmpintarr[k];
                        mycheck.sumvals[k]=tmpfloatarr[k];
                        mycheck.entries[k]=0;
                    }//end for k
                    mycheck.clusters[tmpint]=i;
                    mycheck.sumvals[tmpint]=0;
                    mycheck.entries[tmpint]=0;
                    clusterseqshash.put(tmpname,mycheck);
                }else{
                    mycheck=new checkobject();
                    mycheck.clusters=new int[1];
                    mycheck.clusters[0]=i;
                    mycheck.sumvals=new float[1];
                    mycheck.sumvals[0]=0;
                    mycheck.entries=new int[1];
                    mycheck.entries[0]=0;
                    clusterseqshash.put(tmpname,mycheck);
                }
            }//end for i
        }//end for i
        //next, see how well these definitions fit with what the replicate clusters say
        cluster checkcluster;
        int[] myclusternums;
        int[] checkclusternums;
        boolean sharednum;
        float[] seqsums;
        int[] entries;
        for(int i=0;i<replicates;i++){
            System.out.println("r"+i);
            for(j=replicate[i].size()-1;j>=0;j--){
                System.out.print(".");
                checkcluster=(cluster)replicate[i].elementAt(j);
                //now see how many sequences are shared between this cluster and the reference
                sequences=java.lang.reflect.Array.getLength(checkcluster.member);
                //now find how many sequences it shares with the reference cluster
                for(k=0;k<sequences;k++){
                    if(clusterseqshash.containsKey(String.valueOf(checkcluster.member[k]))){//if this sequence is known as part of a cluster
                        myclusternums=((checkobject)clusterseqshash.get(String.valueOf(checkcluster.member[k]))).clusters;
                        seqsums=((checkobject)clusterseqshash.get(String.valueOf(checkcluster.member[k]))).sumvals;
                        entries=((checkobject)clusterseqshash.get(String.valueOf(checkcluster.member[k]))).entries;
                        //need the for loops as multiple cluster assignments are possible
                        for(l=0;l<sequences;l++){//see if these sequences shared a cluster in the references
                            if(l==k){
                                continue;
                            }
                            if(clusterseqshash.containsKey(String.valueOf(checkcluster.member[l]))){
                                checkclusternums=((checkobject)clusterseqshash.get(String.valueOf(checkcluster.member[l]))).clusters;
                                //now see if myclusternums and checkclusternums share a number
                                for(m=java.lang.reflect.Array.getLength(myclusternums)-1;m>=0;m--){
                                    sharednum=false;
                                    for(n=java.lang.reflect.Array.getLength(checkclusternums)-1;n>=0;n--){
                                        if(myclusternums[m]==checkclusternums[n]){
                                            sharednum=true;
                                            break;//stop this for loop
                                        }
                                    }//end for n
                                    if(sharednum){//if they do share a number, then this increases confidence for this sequence cluster
                                        seqsums[m]++;
                                    }else{//else it reduces it (seqs were in different clusters)
                                        seqsums[m]--;
                                    }
                                    entries[m]++;
                                }//end for m
                            }else{
                                //it the second sequence was not assigned a cluster
                                //add nothing, subtract nothing
                                for(m=java.lang.reflect.Array.getLength(myclusternums)-1;m>=0;m--){
                                    entries[m]++;
                                }//end for m
                            }
                        }//end for l
                        for(l=java.lang.reflect.Array.getLength(myclusternums)-1;l>=0;l--){//now get the average value for each
                            //    ((checkobject)clusterseqshash.get(String.valueOf(checkcluster.member[k]))).sumvals[l]+=seqsums[l]/((float)((sequences-1)*sequences));
                            ((checkobject)clusterseqshash.get(String.valueOf(checkcluster.member[k]))).sumvals[l]+=seqsums[l];
                            ((checkobject)clusterseqshash.get(String.valueOf(checkcluster.member[k]))).entries[l]+=entries[l];
                        }
                    }else{
                        //this sequence was not assigned to a cluster in the references
                        //???then reduce the sumvals for all sequences in that cluster???
                        //???or just let it be???
                        //actually, this should be treated as a null entry, but should add to the "entries" number for each cluster
                    }
                }//end for k
            }//end for j
        }//end for i
        //now average the confidence values for the different sequences in the different clusters
        String[] keysarr=(String[])clusterseqshash.keySet().toArray(new String[0]);
        int seqnamesnum=java.lang.reflect.Array.getLength(keysarr);
        int[] clusterconf=new int[clusternum];
        int[] clusterentries=new int[clusternum];
        int[] seqarr;
        for(int i=0;i<clusternum;i++){
            clusterconf[i]=0;
            clusterentries[i]=0;
            //I haven't created the seqconfidences array yet
            ((cluster)clustervec.elementAt(i)).seqconfidence=new float[java.lang.reflect.Array.getLength(((cluster)clustervec.elementAt(i)).member)];
        }//end for i
        for(int i=0;i<seqnamesnum;i++){
            mycheck=(checkobject)clusterseqshash.get(keysarr[i]);
            for(j=java.lang.reflect.Array.getLength(mycheck.clusters)-1;j>=0;j--){
                clusterconf[mycheck.clusters[j]]+=mycheck.sumvals[j];
                clusterentries[mycheck.clusters[j]]+=mycheck.entries[j];
                seqarr=((cluster)clustervec.elementAt(mycheck.clusters[j])).member;
                tmpint=-1;
                for(k=java.lang.reflect.Array.getLength(seqarr)-1;k>=0;k--){
                    if(String.valueOf(seqarr[k]).equals(keysarr[i])){
                        tmpint=k;
                        break;
                    }
                }//end for k
                if(tmpint>-1){
                    ((cluster)clustervec.elementAt(mycheck.clusters[j])).seqconfidence[tmpint]=mycheck.sumvals[j]/mycheck.entries[j];
                }
            }//end for j
        }//end for i
        for(int i=0;i<clusternum;i++){
            ((cluster)clustervec.elementAt(i)).clusterconfidence=(float)clusterconf[i]/(float)clusterentries[i];
        }//end for i
        testconfidences(clustervec);
    }//end checkreplicates
    
    //--------------------------------------------------------------------------
    
    static void testconfidences(Vector invec){
        int clusternum=invec.size();
        cluster currcluster;
        for(int i=0;i<clusternum;i++){
            currcluster=(cluster)invec.elementAt(i);
            System.out.println("Cluster:"+i);
            System.out.println("\tbootstrap:"+currcluster.clusterconfidence);
            for(int j=java.lang.reflect.Array.getLength(currcluster.member)-1;j>=0;j--){
                System.out.print(" "+currcluster.member[j]+":"+currcluster.seqconfidence[j]+";");
            }//end for j
            System.out.println();
        }//end for i
    }//end testconfidences
    
    //--------------------------------------------------------------------------
    
    static void checkreplicates_old(Vector clustervec, Vector[] replicate){
        //assign confidence values to each cluster of sequences (and to each sequence within the cluster)
        //for each sequence in each cluster see in how many of the other clusters it appears with each of its
        //neighbors in a cluster. if so, ++ else -- if an additional non-known member, do nothing.
        //finally the sum over all clusters should give a "confidence" value
        int vecnum=clustervec.size();
        int replicates=java.lang.reflect.Array.getLength(replicate);
        cluster currcluster;
        HashMap nameclusterhash=new HashMap();//quick lookup of which cluster each sequence is assigned to
        int seqnum,j,k,l,m;
        boolean foundseq;
        for(int i=0;i<vecnum;i++){
            currcluster=(cluster) clustervec.elementAt(i);
            seqnum=java.lang.reflect.Array.getLength(currcluster.member);
            for(j=0;j<seqnum;j++){
                nameclusterhash.put(String.valueOf(currcluster.member[j]),new Integer(i));
            }//end for j
        }//end for i
        //now I know which sequence was assigned to which cluster
        double clustersum;
        double[] seqsums;
        int mycluster;
        cluster currcluster2;
        int seqnum2;
        for(int i=0;i<vecnum;i++){//for each cluster I need to bootstrap
            currcluster=(cluster)clustervec.elementAt(i);
            clustersum=0;
            seqnum=java.lang.reflect.Array.getLength(currcluster.member);
            seqsums=new double[seqnum];
            mycluster=i;
            for(j=0;j<seqnum;j++){//for each sequence in that cluster
                //check how this sequence is assigned in the other clusters
                //see how this sequence matches it's relatives in this cluster
                for(k=0;k<replicates;k++){//for each of the bootstrap replicates
                    for(l=replicate[k].size()-1;l>=0;l--){//for each of the clusters in that replicate
                        //see if this cluster has the sequence I look for
                        foundseq=false;
                        //for each sequence in each of the clusters in each replicate
                        for(m=java.lang.reflect.Array.getLength(((cluster)replicate[k].elementAt(l)).member)-1;m>=0;m--){
                            if(((cluster)replicate[k].elementAt(l)).member[m]==currcluster.member[j]){
                                foundseq=true;
                                break;
                            }
                        }
                        if(foundseq){//if this cluster does contain the sequence,
                            //check how many other sequences correspond in that cluster
                            currcluster2=(cluster)replicate[k].elementAt(l);
                            seqnum2=java.lang.reflect.Array.getLength(currcluster.member);
                            for(m=0;m<seqnum2;m++){
                                //if this sequence is also part of the same cluster in the main vector
                                if(nameclusterhash.containsKey(String.valueOf(currcluster.member[m]))){
                                    if(((Integer)nameclusterhash.get(String.valueOf(currcluster.member[m]))).intValue()==i){
                                        seqsums[j]++;
                                    }else{
                                        seqsums[j]--;
                                    }
                                }
                            }
                        }
                    }//end for l
                }//end for k
            }//end for j
            currcluster.seqconfidence=new float[seqnum];
            for(j=0;j<seqnum;j++){
                currcluster.seqconfidence[j]=(float)(seqsums[j]/(replicates*seqnum));
                clustersum+=currcluster.seqconfidence[j];
            }//end for j
            currcluster.clusterconfidence=(float)(clustersum/(float)seqnum);
            clustervec.setElementAt(currcluster,i);//should be unnecessary but doesn't hurt and is more readable
        }//end for i
    }//end checkreplicate
    
    
}//end class

class checkobject{
    
    public checkobject(){}
    
    int[] clusters=null;
    float[] sumvals=null;
    int[] entries=null;
    
}//end class checkobject
