package clans.algorithms;

import java.util.*;

import clans.model.SequenceCluster;
import clans.model.proteins.MinimalAttractionValue;

public class ClusterDetectionBootstrapping {
    
    public static boolean bootstrapconvex(MinimalAttractionValue[] dataarr, Vector<SequenceCluster> clustervec, String clustermethod, int replicates,
            float remove, float sigmafac, int minseqnum, int elements) {
        //remove "remove" of data and then recluster for "replicates" replicates.
        //then get the cluster confidences AND sequence to cluster confidences
        Vector<SequenceCluster>[] replicate = new Vector[replicates];
        int attnum=dataarr.length;
        Vector<MinimalAttractionValue> tmpdata = new Vector<MinimalAttractionValue>();
        Random rand=new Random(System.currentTimeMillis());
        MinimalAttractionValue[] tmparr;
        for(int i=0;i<replicates;i++){
            System.out.println("replicate "+i);
            //first remove "remove" fraction of data
            for(int j=0;j<attnum;j++){
                if(rand.nextFloat()>=remove){
                tmpdata.addElement(dataarr[j]);
                }
            }//end for j
            tmparr=new MinimalAttractionValue[tmpdata.size()];
            tmpdata.copyInto(tmparr);
            tmpdata.clear();
            //now get the cluster for this data
            replicate[i]=ClusterDetection.getconvex(tmparr,sigmafac,minseqnum,elements);
        }//end for i
        //now compare the replicate clusters to the original
        checkreplicates(clustervec, replicate);
        return true;
    }//end bootstrap
    
    //--------------------------------------------------------------------------
    
    public static boolean bootstraplinkage(MinimalAttractionValue[] dataarr, Vector<SequenceCluster> clustervec, String clustermethod,
            int replicates, float remove, int minlinks, int minseqnum, int elements) {
        //remove "remove" of data and then recluster for "replicates" replicates.
        //then get the cluster confidences AND sequence to cluster confidences
        Vector<SequenceCluster>[] replicate = new Vector[replicates];
        int attnum=dataarr.length;
        Vector<MinimalAttractionValue> tmpdata=new Vector<MinimalAttractionValue>();
        Random rand=new Random(System.currentTimeMillis());
        MinimalAttractionValue[] tmparr;
        for(int i=0;i<replicates;i++){
            System.out.println("replicate "+i);
            //first remove "remove" fraction of data
            for(int j=0;j<attnum;j++){
                if(rand.nextFloat()>=remove){
                tmpdata.addElement(dataarr[j]);
                }
            }//end for j
            tmparr=new MinimalAttractionValue[tmpdata.size()];
            tmpdata.copyInto(tmparr);
            tmpdata.clear();
            //now get the cluster for this data
            replicate[i]=ClusterDetection.getlinkage(tmparr,minlinks,minseqnum,elements);
        }//end for i
        //now compare the replicate clusters to the original
        checkreplicates(clustervec,replicate);
        return true;
    }//end bootstrap
    
    //--------------------------------------------------------------------------
    
    public static boolean bootstrapnetwork(MinimalAttractionValue[] dataarr, Vector <SequenceCluster> clustervec, String clustermethod, int replicates, float remove, int minseqnum, boolean dooffset,boolean globalaverage, int elements,int maxrounds){
        //remove "remove" of data and then recluster for "replicates" replicates.
        //then get the cluster confidences AND sequence to cluster confidences
        Vector<SequenceCluster>[] replicate = new Vector[replicates];
        int attnum=dataarr.length;
        Vector<MinimalAttractionValue> tmpdata=new Vector<MinimalAttractionValue>();
        Random rand=new Random(System.currentTimeMillis());
        MinimalAttractionValue[] tmparr;
        for(int i=0;i<replicates;i++){
            System.out.println("replicate "+i);
            //first remove "remove" fraction of data
            for(int j=0;j<attnum;j++){
                if(rand.nextFloat()>=remove){
                tmpdata.addElement(dataarr[j]);
                }
            }//end for j
            tmparr=new MinimalAttractionValue[tmpdata.size()];
            tmpdata.copyInto(tmparr);
            tmpdata.clear();
            //now get the cluster for this data
            replicate[i]=ClusterDetection.getnetwork(tmparr,minseqnum,dooffset,globalaverage,elements,maxrounds);
            System.out.println("clusternum="+replicate[i].size());
            for(int j=replicate[i].size()-1;j>=0;j--){
                System.out.println("\t seqs=" + replicate[i].elementAt(j).member.length);
            }//end for j
        }//end for i
        //now compare the replicate clusters to the original
        checkreplicates(clustervec,replicate);
        return true;
    }

    /**
     * 
     * @param clustervec
     * @param replicate
     */
    static void checkreplicates(Vector<SequenceCluster> clustervec, Vector<SequenceCluster>[] replicate) {
        System.out.println("comparing the replicates");
        //for this just check the clusters for how often they appear EXACTLY in the replicates
        int replicates=replicate.length;
        int clusternum=clustervec.size();
        SequenceCluster currcluster,checkcluster;
        int replicatessize, i, j, k;
        boolean foundcluster;
        String checkstring;
        //first sort all the members in all clusters in the same fashion
        for(i=0;i<replicates;i++){
            replicatessize=replicate[i].size();
            for(j=0;j<replicatessize;j++){
                java.util.Arrays.sort(replicate[i].elementAt(j).member);
            }//end for j
        }//end for i
        for(i=0;i<clusternum;i++){
            currcluster=clustervec.elementAt(i);
            java.util.Arrays.sort(currcluster.member);
            currcluster.clusterconfidence=0;
            checkstring=makecheckstring(currcluster);
            for(j=0;j<replicates;j++){
                foundcluster=false;
                replicatessize=replicate[j].size();
                for(k=0;((k<replicatessize)&&(foundcluster==false));k++){
                    //now see if any of these clusters matches the currcluster elements exactly
                    checkcluster=(SequenceCluster)replicate[j].elementAt(k);
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
    
    static String makecheckstring(SequenceCluster incluster){
        //convert the cluster elements into a sorted string
        StringBuffer tmp=new StringBuffer();
        int elements=incluster.member.length;
        for(int i=0;i<elements;i++){
            tmp.append(incluster.member[i]+";");
        }//end for i
        return tmp.toString();
    }

    /**
     * 
     * @param clustervec
     * @param replicate
     */
    static void getseqconfidences(Vector<SequenceCluster> clustervec, Vector<SequenceCluster>[] replicate){
        System.out.println("calculating sequence confidences");
        //get the confidence with which each sequence is assigned to each cluster
        //sum how often each cluster sequence pair is present.
        int clusternum=clustervec.size();
        int replicates=replicate.length;
        int seqnum, checkseqnum, i, j, k, l, m;
        float sharednum,nonsharednum;
        SequenceCluster currcluster,checkcluster;
        boolean hasseq;
        for(i=0;i<clusternum;i++){//for each cluster
            currcluster=clustervec.elementAt(i);
            seqnum=currcluster.member.length;
            currcluster.seqconfidence=new float[seqnum];
            for(j=0;j<seqnum;j++){
                currcluster.seqconfidence[j]=0;
            }//end for j
            //now see how often each sequence in this cluster is recovered as belonging to the same group in the replicates
            for(j=0;j<replicates;j++){
                for(k=replicate[j].size()-1;k>=0;k--){
                    sharednum=0;
                    nonsharednum=0;
                    checkcluster = replicate[j].elementAt(k);
                    checkseqnum=checkcluster.member.length;
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
    
    static void checkreplicatesold(Vector<SequenceCluster> clustervec, Vector<SequenceCluster>[] replicate){
        //check the bootstrap replicates
        int replicates=replicate.length;
        int clusternum=clustervec.size();
        HashMap<String, checkobject> clusterseqshash=new HashMap<String, checkobject>();
        int sequences,tmpint,j,k,l,m,n;
        int[] tmpintarr;
        float[] tmpfloatarr;
        String tmpname;
        checkobject mycheck;
        //now make a quick lookup for which sequences were defined in which cluster
        for(int i=0;i<clusternum;i++){
            sequences = clustervec.elementAt(i).member.length;
            for(j=0;j<sequences;j++){
                tmpname = String.valueOf(clustervec.elementAt(i).member[j]);
                if(clusterseqshash.containsKey(tmpname)){
                    mycheck=clusterseqshash.get(tmpname);
                    tmpintarr=mycheck.clusters;
                    tmpfloatarr=mycheck.sumvals;
                    tmpint=tmpintarr.length;
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
        SequenceCluster checkcluster;
        int[] myclusternums;
        int[] checkclusternums;
        boolean sharednum;
        float[] seqsums;
        int[] entries;
        for(int i=0;i<replicates;i++){
            System.out.println("r"+i);
            for(j=replicate[i].size()-1;j>=0;j--){
                System.out.print(".");
                checkcluster = replicate[i].elementAt(j);
                
                //now see how many sequences are shared between this cluster and the reference
                sequences=checkcluster.member.length;
                //now find how many sequences it shares with the reference cluster
                for(k=0;k<sequences;k++){
                    if(clusterseqshash.containsKey(String.valueOf(checkcluster.member[k]))){//if this sequence is known as part of a cluster
                        myclusternums=clusterseqshash.get(String.valueOf(checkcluster.member[k])).clusters;
                        seqsums=clusterseqshash.get(String.valueOf(checkcluster.member[k])).sumvals;
                        entries=clusterseqshash.get(String.valueOf(checkcluster.member[k])).entries;
                        //need the for loops as multiple cluster assignments are possible
                        for(l=0;l<sequences;l++){//see if these sequences shared a cluster in the references
                            if(l==k){
                                continue;
                            }
                            if(clusterseqshash.containsKey(String.valueOf(checkcluster.member[l]))){
                                checkclusternums=clusterseqshash.get(String.valueOf(checkcluster.member[l])).clusters;
                                //now see if myclusternums and checkclusternums share a number
                                for(m=myclusternums.length-1;m>=0;m--){
                                    sharednum=false;
                                    for(n=checkclusternums.length-1;n>=0;n--){
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
                                for(m=myclusternums.length-1;m>=0;m--){
                                    entries[m]++;
                                }//end for m
                            }
                        }//end for l
                        for(l=myclusternums.length-1;l>=0;l--){//now get the average value for each
                            //    ((checkobject)clusterseqshash.get(String.valueOf(checkcluster.member[k]))).sumvals[l]+=seqsums[l]/((float)((sequences-1)*sequences));
                            clusterseqshash.get(String.valueOf(checkcluster.member[k])).sumvals[l]+=seqsums[l];
                            clusterseqshash.get(String.valueOf(checkcluster.member[k])).entries[l]+=entries[l];
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
        String[] keysarr=clusterseqshash.keySet().toArray(new String[0]);
        int seqnamesnum=keysarr.length;
        int[] clusterconf=new int[clusternum];
        int[] clusterentries=new int[clusternum];
        int[] seqarr;
        for(int i=0;i<clusternum;i++){
            clusterconf[i]=0;
            clusterentries[i]=0;
            //I haven't created the seqconfidences array yet
            clustervec.elementAt(i).seqconfidence = new float[clustervec.elementAt(i).member.length];
        }//end for i
        for(int i=0;i<seqnamesnum;i++){
            mycheck=clusterseqshash.get(keysarr[i]);
            for(j=mycheck.clusters.length-1;j>=0;j--){
                clusterconf[mycheck.clusters[j]]+=mycheck.sumvals[j];
                clusterentries[mycheck.clusters[j]]+=mycheck.entries[j];
                seqarr=((SequenceCluster)clustervec.elementAt(mycheck.clusters[j])).member;
                tmpint=-1;
                for(k=seqarr.length-1;k>=0;k--){
                    if(String.valueOf(seqarr[k]).equals(keysarr[i])){
                        tmpint=k;
                        break;
                    }
                }//end for k
                if(tmpint>-1){
                    ((SequenceCluster)clustervec.elementAt(mycheck.clusters[j])).seqconfidence[tmpint]=mycheck.sumvals[j]/mycheck.entries[j];
                }
            }//end for j
        }//end for i
        for(int i=0;i<clusternum;i++){
            ((SequenceCluster)clustervec.elementAt(i)).clusterconfidence=(float)clusterconf[i]/(float)clusterentries[i];
        }//end for i
        testconfidences(clustervec);
    }

    /**
     * 
     * @param invec
     */
    static void testconfidences(Vector<SequenceCluster> invec){
        int clusternum=invec.size();
        SequenceCluster currcluster;
        for(int i=0;i<clusternum;i++){
            currcluster=(SequenceCluster)invec.elementAt(i);
            System.out.println("Cluster:"+i);
            System.out.println("\tbootstrap:"+currcluster.clusterconfidence);
            for(int j=currcluster.member.length-1;j>=0;j--){
                System.out.print(" "+currcluster.member[j]+":"+currcluster.seqconfidence[j]+";");
            }
            System.out.println();
        }
    }
    

    /**
     * assign confidence values to each cluster of sequences (and to each sequence within the cluster) for each sequence
     * in each cluster see in how many of the other clusters it appears with each of its neighbors in a cluster. if so,
     * ++ else -- if an additional non-known member, do nothing. finally the sum over all clusters should give a
     * "confidence" value
     * 
     * @param clustervec
     * @param replicate
     */
    static void checkreplicates_old(Vector<SequenceCluster> clustervec, Vector<SequenceCluster>[] replicate) {
        int vecnum = clustervec.size();   
        int replicates=replicate.length;
        SequenceCluster currcluster;
        HashMap<String, Integer> nameclusterhash=new HashMap<String, Integer>();//quick lookup of which cluster each sequence is assigned to
        int seqnum,j,k,l,m;
        boolean foundseq;
        for(int i=0;i<vecnum;i++){
            currcluster=clustervec.elementAt(i);
            seqnum=currcluster.member.length;
            for(j=0;j<seqnum;j++){
                nameclusterhash.put(String.valueOf(currcluster.member[j]),new Integer(i));
            }//end for j
        }//end for i
        //now I know which sequence was assigned to which cluster
        double clustersum;
        double[] seqsums;
        int seqnum2;
        for(int i=0;i<vecnum;i++){//for each cluster I need to bootstrap
            currcluster=clustervec.elementAt(i);
            clustersum=0;
            seqnum=currcluster.member.length;
            seqsums=new double[seqnum];
            for(j=0;j<seqnum;j++){//for each sequence in that cluster
                //check how this sequence is assigned in the other clusters
                //see how this sequence matches it's relatives in this cluster
                for(k=0;k<replicates;k++){//for each of the bootstrap replicates
                    for(l=replicate[k].size()-1;l>=0;l--){//for each of the clusters in that replicate
                        //see if this cluster has the sequence I look for
                        foundseq=false;
                        //for each sequence in each of the clusters in each replicate
                        for (m = replicate[k].elementAt(l).member.length - 1; m >= 0; m--) {
                            if(((SequenceCluster)replicate[k].elementAt(l)).member[m]==currcluster.member[j]){
                                foundseq=true;
                                break;
                            }
                        }
                        if(foundseq){//if this cluster does contain the sequence,
                            //check how many other sequences correspond in that cluster
                            seqnum2=currcluster.member.length;
                            for(m=0;m<seqnum2;m++){
                                //if this sequence is also part of the same cluster in the main vector
                                if(nameclusterhash.containsKey(String.valueOf(currcluster.member[m]))){
                                    if(nameclusterhash.get(String.valueOf(currcluster.member[m])).intValue()==i){
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
