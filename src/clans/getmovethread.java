/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package clans;

/**
 *
 * @author tancred
 */
public class getmovethread extends java.lang.Thread {

    public getmovethread(float[][] myposarr, minattvals[] myattvals, float[][] mymovearr, int myi, int cpu, java.util.HashMap selectnamehash, int[] selectnames, String syncon, clusterdata parent) {
        this.done = false;
        this.posarr = myposarr;
        this.attvals = myattvals;
        this.movearr = mymovearr;
        this.myi = myi;
        this.cpu = cpu;
        this.doselected = true;
        this.tmphash = selectnamehash;
        this.parent = parent;
        this.syncon = syncon;
        this.selectnames = selectnames;
    }

    public getmovethread(float[][] myposarr, minattvals[] myattvals, float[][] mymovearr, int myi, int cpu, String syncon, clusterdata parent) {
        this.done = false;
        this.posarr = myposarr;
        this.attvals = myattvals;
        this.movearr = mymovearr;
        this.myi = myi;
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
    int myi;
    int cpu;
    java.util.HashMap tmphash = null;
    clusterdata parent;
    String syncon;
    int[] selectnames;//a local copy of parent.selectednames, as that may change during a calculation

    @Override
    public void run() {
        //System.out.println("start="+start);
        //use the positions of all elements and their attraction/repulsion values to
        //calculate a movement vector for each (take into account the last movement).
        //repulsion doesn't have a specific value as all evalues below a certain point
        //are simply regarded as insignificant. therefore use a one formula for all
        //to compute the repulsive forces (see getrepulse)
        int j;
        int elements = java.lang.reflect.Array.getLength(posarr);
        double[] currmoverep = new double[3];
        double[] currmoveatt = new double[3];
        double totaldist = 0;
        //double totalmove=0;
        //double avgmove=0;
        //double minattelems=minattract*elements;
        int attnum = java.lang.reflect.Array.getLength(attvals);
        double minattract = parent.minattract;
        int repvalpow = parent.repvalpow;
        int attvalpow = parent.attvalpow;
        float repfactor = parent.repfactor;
        float attfactor = parent.attfactor;
        float maxmove = parent.maxmove;

        if (doselected) {
            //System.out.println("getmoveselected");
            int selectnamesnum = java.lang.reflect.Array.getLength(selectnames);
            //now get from where to where I should do my calculations
            int start, end;//,attstart,attend;
            if (myi == (cpu - 1)) {
                start = myi * (int) (selectnamesnum / cpu);
                end = selectnamesnum;
                //attstart=myi*(int)(attnum/cpu);
                //attend=attnum;
            } else {
                start = myi * (int) (selectnamesnum / cpu);
                end = (myi + 1) * (int) (selectnamesnum / cpu);
                //attstart=myi*(int)(attnum/cpu);
                //attend=(myi+1)*(int)(attnum/cpu);
            }
            if (parent.cluster2d == true) {
                //System.out.println("clusterselected2D");
                //cluster only the selected sequences in 2D
                for (int i = start; i < end; i++) {
                    //currmoveatt=getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                    ClusterMethods.getminattract(posarr[selectnames[i]], currmoveatt, minattract);
                    movearr[selectnames[i]][0] += currmoveatt[0];
                    movearr[selectnames[i]][1] += currmoveatt[1];
                    for (j = elements; --j >= 0;) {
                        //currmoverep=getrepulse2d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow,repfactor,rand);
                        ClusterMethods.getrepulse2d(posarr[selectnames[i]], posarr[j], currmoverep, repvalpow, repfactor, ClusterMethods.rand);
                        movearr[selectnames[i]][0] += currmoverep[0];
                        movearr[selectnames[i]][1] += currmoverep[1];
                    }//end for j
                }//end for i
                //now add the attraction values, but only for the query or hit sequences in my part of the selectnames array (assigned in recluster3d)
                for (int i = attnum; --i >= 0;) {
                    if ((tmphash.containsKey(String.valueOf(attvals[i].query))) && (((Integer) tmphash.get(String.valueOf(attvals[i].query))).intValue() == myi)) {
                        //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                        ClusterMethods.getattract2d(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt, attvalpow, attfactor);
                        movearr[attvals[i].query][0] += currmoveatt[0];
                        movearr[attvals[i].query][1] += currmoveatt[1];
                        //movement[attvals[i].query][2]+=currmoveatt[2];
                        if ((tmphash.containsKey(String.valueOf(attvals[i].hit))) && (((Integer) tmphash.get(String.valueOf(attvals[i].hit))).intValue() == myi)) {
                            movearr[attvals[i].hit][0] -= currmoveatt[0];
                            movearr[attvals[i].hit][1] -= currmoveatt[1];
                            //movement[attvals[i].hit][2]-=currmoveatt[2];
                        }
                    } else if ((tmphash.containsKey(String.valueOf(attvals[i].hit))) && (((Integer) tmphash.get(String.valueOf(attvals[i].hit))).intValue() == myi)) {
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
                    //currmoveatt=getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                    ClusterMethods.getminattract(posarr[selectnames[i]], currmoveatt, minattract);
                    movearr[selectnames[i]][0] += currmoveatt[0];
                    movearr[selectnames[i]][1] += currmoveatt[1];
                    movearr[selectnames[i]][2] += currmoveatt[2];
                    for (j = elements; --j >= 0;) {
                        //currmoverep=getrepulse3d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow,repfactor,rand);
                        ClusterMethods.getrepulse3d(posarr[selectnames[i]], posarr[j], currmoverep, repvalpow, repfactor, ClusterMethods.rand);
                        movearr[selectnames[i]][0] += currmoverep[0];
                        movearr[selectnames[i]][1] += currmoverep[1];
                        movearr[selectnames[i]][2] += currmoverep[2];
                    }//end for j
                }//end for i
                //now add the attraction values, but only for the query or hit sequences in my part of the selectnames array (assigned in recluster3d)
                for (int i = attnum; --i >= 0;) {
                    if ((tmphash.containsKey(String.valueOf(attvals[i].query))) && (((Integer) tmphash.get(String.valueOf(attvals[i].query))).intValue() == myi)) {
                        //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                        ClusterMethods.getattract3d(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt, attvalpow, attfactor);
                        movearr[attvals[i].query][0] += currmoveatt[0];
                        movearr[attvals[i].query][1] += currmoveatt[1];
                        movearr[attvals[i].query][2] += currmoveatt[2];
                        if ((tmphash.containsKey(String.valueOf(attvals[i].hit))) && (((Integer) tmphash.get(String.valueOf(attvals[i].hit))).intValue() == myi)) {
                            movearr[attvals[i].hit][0] -= currmoveatt[0];
                            movearr[attvals[i].hit][1] -= currmoveatt[1];
                            movearr[attvals[i].hit][2] -= currmoveatt[2];
                        }
                    } else if ((tmphash.containsKey(String.valueOf(attvals[i].hit))) && (((Integer) tmphash.get(String.valueOf(attvals[i].hit))).intValue() == myi)) {
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
        } else {//if no sequences were selected or all should be used
            //now get from where to where I should do my calculations
            int start, end;//,attstart,attend;
            if (myi == (cpu - 1)) {
                //special case, do everything from here to end to avoid rounding errors
                start = myi * (int) (elements / cpu);
                end = elements;
                //attstart=myi*(int)(attnum/cpu);
                //attend=attnum;
            } else {
                start = myi * (int) (elements / cpu);
                end = (myi + 1) * (int) (elements / cpu);
                //attstart=myi*(int)(attnum/cpu);
                //attend=(myi+1)*(int)(attnum/cpu);
            }
            if (parent.cluster2d == true) {
                //cluster all in 2D
                //System.out.println("cluster2D");
                for (int i = start; i < end; i++) {
                    //currmoveatt=getminattract(posarr[i],currmoveatt,minattract);
                    ClusterMethods.getminattract(posarr[i], currmoveatt, minattract);
                    movearr[i][0] += currmoveatt[0];
                    movearr[i][1] += currmoveatt[1];
                    for (j = elements; --j >= 0;) {
                        //currmoverep=getrepulse2d(posarr[i],posarr[j],currmoverep,repvalpow,repfactor,rand);
                        ClusterMethods.getrepulse2d(posarr[i], posarr[j], currmoverep, repvalpow, repfactor, ClusterMethods.rand);
                        movearr[i][0] += currmoverep[0];
                        movearr[i][1] += currmoverep[1];
                    }//end for j
                }//end for i
                for (int i = attnum; --i >= 0;) {
                    if (attvals[i].query >= start && attvals[i].query < end) {
                        //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                        ClusterMethods.getattract2d(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt, attvalpow, attfactor);
                        movearr[attvals[i].query][0] += currmoveatt[0];
                        movearr[attvals[i].query][1] += currmoveatt[1];
                    }
                    if (attvals[i].hit >= start && attvals[i].hit < end) {
                        //currmoveatt=getattract2d(posarr[attvals[i].hit],posarr[attvals[i].query],attvals[i].att,currmoveatt,attvalpow, attfactor);
                        ClusterMethods.getattract2d(posarr[attvals[i].hit], posarr[attvals[i].query], attvals[i].att, currmoveatt, attvalpow, attfactor);
                        movearr[attvals[i].hit][0] += currmoveatt[0];
                        movearr[attvals[i].hit][1] += currmoveatt[1];
                    }
                }//end for i
                //double totaldistsq;
                for (int i = start; i < end; i++) {
                    movearr[i][0] /= elements;
                    movearr[i][1] /= elements;
                    movearr[i][2] = 0;
                    totaldist = java.lang.Math.sqrt((movearr[i][0] * movearr[i][0]) + (movearr[i][1] * movearr[i][1]));
                    if (totaldist > maxmove) {
                        movearr[i][0] *= maxmove / totaldist;
                        movearr[i][1] *= maxmove / totaldist;
                    }
                }//end for i
            } else {
                //cluster all in 3D
                //System.out.println("cluster3D");
                for (int i = start; i < end; i++) {
                    //currmoveatt=getminattract(posarr[i],currmoveatt,minattract);
                    ClusterMethods.getminattract(posarr[i], currmoveatt, minattract);
                    movearr[i][0] += currmoveatt[0];
                    movearr[i][1] += currmoveatt[1];
                    movearr[i][2] += currmoveatt[2];
                    for (j = 0; j < elements; j++) {
                        //currmoverep=getrepulse3d(posarr[i],posarr[j],currmoverep,repvalpow,repfactor,rand);
                        ClusterMethods.getrepulse3d(posarr[i], posarr[j], currmoverep, repvalpow, repfactor, ClusterMethods.rand);
                        movearr[i][0] += currmoverep[0];
                        movearr[i][1] += currmoverep[1];
                        movearr[i][2] += currmoverep[2];
                    }//end for j
                }//end for i
                for (int i = attnum; --i >= 0;) {
                    if (attvals[i].query >= start && attvals[i].query < end) {
                        //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow, attfactor);
                        ClusterMethods.getattract3d(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt, attvalpow, attfactor);
                        movearr[attvals[i].query][0] += currmoveatt[0];
                        movearr[attvals[i].query][1] += currmoveatt[1];
                        movearr[attvals[i].query][2] += currmoveatt[2];
                    }
                    if (attvals[i].hit >= start && attvals[i].hit < end) {
                        //currmoveatt=getattract3d(posarr[attvals[i].hit],posarr[attvals[i].query],attvals[i].att,currmoveatt,attvalpow, attfactor);
                        ClusterMethods.getattract3d(posarr[attvals[i].hit], posarr[attvals[i].query], attvals[i].att, currmoveatt, attvalpow, attfactor);
                        movearr[attvals[i].hit][0] += currmoveatt[0];
                        movearr[attvals[i].hit][1] += currmoveatt[1];
                        movearr[attvals[i].hit][2] += currmoveatt[2];
                    }
                }//end for i
                //double totaldistsq;
                for (int i = start; i < end; i++) {
                    movearr[i][0] /= elements;
                    movearr[i][1] /= elements;
                    movearr[i][2] /= elements;
                    totaldist = java.lang.Math.sqrt((movearr[i][0] * movearr[i][0]) + (movearr[i][1] * movearr[i][1]) + (movearr[i][2] * movearr[i][2]));
                    if (totaldist > maxmove) {
                        movearr[i][0] *= maxmove / totaldist;
                        movearr[i][1] *= maxmove / totaldist;
                        movearr[i][2] *= maxmove / totaldist;
                    }
                }//end for i
            }//end clustering all in 3D
        }//end if not doselected
        this.done = true;
        synchronized (syncon) {
            syncon.notify();
        }// end synchronized parent
    }// end run
}//end class getmovethread

