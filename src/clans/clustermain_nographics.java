/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clans;

/**
 *
 * @author tancred
 */
public class clustermain_nographics {
    public clustermain_nographics(clusterdata data){
        this.data=data;
    }

    public clusterdata data=null;

    //the methods I want to use

    void initgraph(){
//        if(data.loadsaved!=null){
//            System.out.println("loading data from "+data.loadsaved);
//            clustermethods.loaddata(data);
//        }
        data.minpval=data.pval;
        data.mineval=data.eval;
        if(data.scval>=0){//in that case use a score cutoff
            data.minpval=data.scval;
            data.usescval=true;
        }else{
            data.usescval=false;
        }
        clustermethods.initgraph(data);
        
    }//end initgraph

    
    computethread mythread=new computethread(this);


    //--------------------------------------------------------------------------

    void loaddata(String inname){
        if(mythread!=null && mythread.stop!=true){
            System.err.println("Warning, you should stop the clustering thread before loading another file; stopping thread now");
            mythread.stop=true;
        }//else everything is OK
        data.loadsaved=inname;
        clustermethods.loaddata(data);
    }//end loaddata


    void startstopthread(){
        //does the iteration and the stop for a thread
        if(mythread.didrun==false || mythread.stop==true){//if this thread was never started or stopped
            //String tmpstr="";
            if(data.mymovearr==null){
                System.err.println("WARNING: No data currently available; please load some data!");
                return;
            }
            mythread.stop=true;
            mythread=new computethread(this);
            mythread.start();
        }else{//if this thread is running
            mythread.stop=true;
            //is unavailable until the thread has stopped running
            //the thread then sets the text to "resume" and re-enables the button.
        }
    }//end startstopthread



    //--------------------------------------------------------------------------
    //------------------------THREADS-------------------------------------------
    //--------------------------------------------------------------------------

    class computethread extends java.lang.Thread{

        public computethread(clustermain_nographics parent){
            this.parent=parent;
            this.stop=false;
            this.didrun=false;
        }

        boolean stop=true;
        boolean didrun=false;
        String tmpstr="";
        float tmpcool=1;
        clustermain_nographics parent;
        int roundsdone=0;

@Override
        public void run(){
            this.didrun=true;
            data.roundsdone=0;
            while (stop==false){
                data.rounds++;
                if(data.roundslimit!=-1){
                    data.roundsdone++;
                    //stopbutton.setText("STOP ("+roundsdone+"/"+roundslimit+")");
                    if(data.roundsdone>=data.roundslimit){
                        stop=true;
                        synchronized(parent){
                            parent.notify();
                        }
                    }
                }
                //parent.data.myposarr=clustermethods.recluster3d(parent.data);
                clustermethods.recluster3d(data);
                data.posarr=data.myposarr;
                tmpcool=(((float)((int)(data.currcool*100000)))/100000);
                if(tmpcool<=1e-5){
                    stop=true;
                }
            }// end while
        }// end run

    }// end class computethread

}
