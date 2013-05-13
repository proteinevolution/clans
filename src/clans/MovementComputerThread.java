package clans;

import java.util.HashMap;

/**
 *
 * @author tancred
 */
public class MovementComputerThread extends java.lang.Thread {

    public MovementComputerThread(float[][] myposarr, minattvals[] myattvals, float[][] mymovearr, int thread_id, int cpu,
            HashMap<String, Integer> selectnamehash, int[] selectnames, String syncon, ClusterData parent) {
        this.done = false;
        this.posarr = myposarr;
        this.attvals = myattvals;
        this.movearr = mymovearr;
        this.thread_id = thread_id;
        this.cpu = cpu;
        this.doselected = true;
        this.tmphash = selectnamehash;
        this.parent = parent;
        this.syncon = syncon;
        this.selectnames = selectnames;
    }

    public MovementComputerThread(float[][] myposarr, minattvals[] myattvals, float[][] mymovearr, int thread_id, int cpu,
            String syncon, ClusterData parent) {
        this.done = false;
        this.posarr = myposarr;
        this.attvals = myattvals;
        this.movearr = mymovearr;
        this.thread_id = thread_id;
        this.cpu = cpu;
        this.doselected = false;
        this.tmphash = null;
        this.parent = parent;
        this.syncon = syncon;
    }

    boolean done;
    boolean doselected;
    float[][] posarr;
    minattvals[] attvals;
    float[][] movearr;
    int thread_id;
    int cpu;
    HashMap<String, Integer> tmphash = null;
    ClusterData parent;
    String syncon;
    int[] selectnames;// a local copy of parent.selectednames, as that may change during a calculation

    /**
     * use the positions of all elements and their attraction/repulsion values to calculate a movement vector for each
     * (take into account the last movement). repulsion doesn't have a specific value as all evalues below a certain
     * point are simply regarded as insignificant. therefore use a one formula for all to compute the repulsive forces
     * (see getrepulse)
     */
    @Override
    public void run() {

        int elements = posarr.length;
        double[] currmoverep = new double[3];
        double[] currmoveatt = new double[3];
        double totaldist = 0;

        int attnum = attvals.length;
        double minattract = parent.minattract;
        int repvalpow = parent.repvalpow;
        int attvalpow = parent.attvalpow;
        float repfactor = parent.repfactor;
        float attfactor = parent.attfactor;
        float maxmove = parent.maxmove;

        if (doselected) {
            int selectnames_count = selectnames.length;
            //now get from where to where I should do my calculations
            int start, end;
            if (thread_id == (cpu - 1)) {
                start = thread_id * (int) (selectnames_count / cpu);
                end = selectnames_count;
            } else {
                start = thread_id * (int) (selectnames_count / cpu);
                end = (thread_id + 1) * (int) (selectnames_count / cpu);
            }
            
            if (parent.cluster2d == true) {
                //System.out.println("clusterselected2D");
                //cluster only the selected sequences in 2D
                for (int i = start; i < end; i++) {

                    ClusterMethods.add_attraction_towards_origin(posarr[selectnames[i]], currmoveatt, minattract);
                    movearr[selectnames[i]][0] += currmoveatt[0];
                    movearr[selectnames[i]][1] += currmoveatt[1];
                    for (int j = elements; --j >= 0;) {
                        //currmoverep=getrepulse2d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow,repfactor,rand);
                        ClusterMethods.getrepulse2d(posarr[selectnames[i]], posarr[j], currmoverep, repvalpow, repfactor, ClusterMethods.rand);
                        movearr[selectnames[i]][0] += currmoverep[0];
                        movearr[selectnames[i]][1] += currmoverep[1];
                    }//end for j
                }//end for i
                //now add the attraction values, but only for the query or hit sequences in my part of the selectnames array (assigned in recluster3d)
                for (int i = attnum; --i >= 0;) {
                    if ((tmphash.containsKey(String.valueOf(attvals[i].query))) && (tmphash.get(String.valueOf(attvals[i].query)).intValue() == thread_id)) {
                        //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                        ClusterMethods.getattract2d(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt, attvalpow, attfactor);
                        movearr[attvals[i].query][0] += currmoveatt[0];
                        movearr[attvals[i].query][1] += currmoveatt[1];
                        //movement[attvals[i].query][2]+=currmoveatt[2];
                        if ((tmphash.containsKey(String.valueOf(attvals[i].hit))) && (tmphash.get(String.valueOf(attvals[i].hit)).intValue() == thread_id)) {
                            movearr[attvals[i].hit][0] -= currmoveatt[0];
                            movearr[attvals[i].hit][1] -= currmoveatt[1];
                            //movement[attvals[i].hit][2]-=currmoveatt[2];
                        }
                    } else if ((tmphash.containsKey(String.valueOf(attvals[i].hit))) && (tmphash.get(String.valueOf(attvals[i].hit)).intValue() == thread_id)) {
                        //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                        ClusterMethods.getattract2d(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt, attvalpow, attfactor);
                        movearr[attvals[i].hit][0] -= currmoveatt[0];
                        movearr[attvals[i].hit][1] -= currmoveatt[1];
                        //movement[attvals[i].hit][2]-=currmoveatt[2];
                    }
                }//end for i
                //double totaldistsq;
                for (int i = start; i < end; i++) {
                    movearr[selectnames[i]][0] /= elements;
                    movearr[selectnames[i]][1] /= elements;
                    movearr[selectnames[i]][2] = 0;
                    totaldist = java.lang.Math.sqrt((movearr[selectnames[i]][0] * movearr[selectnames[i]][0]) + (movearr[selectnames[i]][1] * movearr[selectnames[i]][1]));
                    if (totaldist > maxmove) {
                        movearr[selectnames[i]][0] *= maxmove / totaldist;
                        movearr[selectnames[i]][1] *= maxmove / totaldist;
                    }
                }//end for i
            } else {
                //cluster only the selected sequences in 3D
                //System.out.println("clusterselected3D");
                for (int i = start; i < end; i++) {

                    ClusterMethods.add_attraction_towards_origin(posarr[selectnames[i]], currmoveatt, minattract);
                    movearr[selectnames[i]][0] += currmoveatt[0];
                    movearr[selectnames[i]][1] += currmoveatt[1];
                    movearr[selectnames[i]][2] += currmoveatt[2];
                    for (int j = elements; --j >= 0;) {
                        //currmoverep=getrepulse3d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow,repfactor,rand);
                        ClusterMethods.getrepulse3d(posarr[selectnames[i]], posarr[j], currmoverep, repvalpow, repfactor, ClusterMethods.rand);
                        movearr[selectnames[i]][0] += currmoverep[0];
                        movearr[selectnames[i]][1] += currmoverep[1];
                        movearr[selectnames[i]][2] += currmoverep[2];
                    }//end for j
                }//end for i
                //now add the attraction values, but only for the query or hit sequences in my part of the selectnames array (assigned in recluster3d)
                for (int i = attnum; --i >= 0;) {
                    if ((tmphash.containsKey(String.valueOf(attvals[i].query))) && (tmphash.get(String.valueOf(attvals[i].query)).intValue() == thread_id)) {
                        //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                        ClusterMethods.getattract3d(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt, attvalpow, attfactor);
                        movearr[attvals[i].query][0] += currmoveatt[0];
                        movearr[attvals[i].query][1] += currmoveatt[1];
                        movearr[attvals[i].query][2] += currmoveatt[2];
                        if ((tmphash.containsKey(String.valueOf(attvals[i].hit))) && (tmphash.get(String.valueOf(attvals[i].hit)).intValue() == thread_id)) {
                            movearr[attvals[i].hit][0] -= currmoveatt[0];
                            movearr[attvals[i].hit][1] -= currmoveatt[1];
                            movearr[attvals[i].hit][2] -= currmoveatt[2];
                        }
                    } else if ((tmphash.containsKey(String.valueOf(attvals[i].hit))) && (tmphash.get(String.valueOf(attvals[i].hit)).intValue() == thread_id)) {
                        //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                        ClusterMethods.getattract3d(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt, attvalpow, attfactor);
                        movearr[attvals[i].hit][0] -= currmoveatt[0];
                        movearr[attvals[i].hit][1] -= currmoveatt[1];
                        movearr[attvals[i].hit][2] -= currmoveatt[2];
                    }
                }//end for i
                //double totaldistsq;
                for (int i = start; i < end; i++) {
                    movearr[selectnames[i]][0] /= elements;
                    movearr[selectnames[i]][1] /= elements;
                    movearr[selectnames[i]][2] /= elements;
                    totaldist = java.lang.Math.sqrt((movearr[selectnames[i]][0] * movearr[selectnames[i]][0]) + (movearr[selectnames[i]][1] * movearr[selectnames[i]][1]) + (movearr[selectnames[i]][2] * movearr[selectnames[i]][2]));
                    if (totaldist > maxmove) {
                        movearr[selectnames[i]][0] *= maxmove / totaldist;
                        movearr[selectnames[i]][1] *= maxmove / totaldist;
                        movearr[selectnames[i]][2] *= maxmove / totaldist;
                    }
                }//end for i
            }
        } else { // if no sequences were selected or all should be used

            // get start and end index of data concerning this thread
            int start, end;
            if (thread_id == (cpu - 1)) {
                // do everything from here to the last element to avoid rounding errors
                start = thread_id * (int) (elements / cpu);
                end = elements;
            } else {
                start = thread_id * (int) (elements / cpu);
                end = (thread_id + 1) * (int) (elements / cpu);
            }
            
            // add repulsive movements
            for (int i = start; i < end; i++) {
                ClusterMethods.add_attraction_towards_origin(posarr[i], currmoveatt, minattract);
                
                for (int j = 0; j < movearr[i].length; j ++) {
                    movearr[i][j] += currmoveatt[j];
                }

                for (int j = 0; j < elements; j++) {
                    if (i == j) {
                        continue;
                    } 

                    ClusterMethods.getrepulse3d(posarr[i], posarr[j], currmoverep, repvalpow, repfactor, ClusterMethods.rand);
                    
                    for (int k = 0; k < movearr[i].length; k ++) {
                        movearr[i][k] += currmoverep[k];
                    }
                }
            }
            
            // add attractive movements
            for (int i = attnum; --i >= 0;) {
                if (attvals[i].query >= start && attvals[i].query < end) {
                    ClusterMethods.getattract3d(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt, attvalpow, attfactor);
                    
                    for (int j = 0; j < movearr[attvals[i].query].length; j ++) { 
                        movearr[attvals[i].query][j] += currmoveatt[j];
                    }
                }
                
                if (attvals[i].hit >= start && attvals[i].hit < end) {
                    ClusterMethods.getattract3d(posarr[attvals[i].hit], posarr[attvals[i].query], attvals[i].att, currmoveatt, attvalpow, attfactor);

                    for (int j = 0; j < movearr[attvals[i].hit].length; j ++) {
                        movearr[attvals[i].hit][j] += currmoveatt[j];
                    }
                }
            }

            // restrict to the maximal allowed move distance
            for (int i = start; i < end; i++) {
                for (int j = 0; j < movearr[i].length; j ++) {
                    movearr[i][j] /= elements;
                }
                
                totaldist = java.lang.Math.sqrt((movearr[i][0] * movearr[i][0]) + (movearr[i][1] * movearr[i][1])
                        + (movearr[i][2] * movearr[i][2]));
                
                if (totaldist > maxmove) {
                    for (int j = 0; j < movearr[i].length; j ++) {
                        movearr[i][j] *= maxmove / totaldist;
                    }
                }
            }
        }

        this.done = true;
        synchronized (syncon) {
            syncon.notify();
        }
    }
}

