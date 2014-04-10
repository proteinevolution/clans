package clans;

public class ClusteringWithoutGui {

	public ClusterData data=null;
    computethread mythread=new computethread(this);

	public ClusteringWithoutGui(ClusterData data){
        this.data=data;
        this.data.nographics = true;
        
        this.setup_cutoffs();
	}
	
	public void setup_cutoffs() {
		data.pvalue_threshold = data.pval;
        data.mineval = data.eval;
        
        if (data.scval >= 0) {//in that case use a score cutoff
            data.pvalue_threshold = data.scval;
            data.usescval = true;
        } else {
            data.usescval = false;
        }
        
        data.compute_attraction_values();
	}

    class computethread extends java.lang.Thread{
        //--------------------------------------------------------------------------
        //------------------------THREADS-------------------------------------------
        //--------------------------------------------------------------------------

        public computethread(ClusteringWithoutGui parent){
            this.parent=parent;
            this.stop=false;
            this.didrun=false;
        }

        boolean stop=true;
        boolean didrun=false;
        float tmpcool=1;
        ClusteringWithoutGui parent;

        @Override
        public void run(){
            this.didrun=true;
            data.roundsdone=0;
            while (stop==false){
                data.rounds++;
                if(data.roundslimit!=-1){
                    data.roundsdone++;
                  
                    if(data.roundsdone>=data.roundslimit){
                        stop=true;
                        synchronized(parent){
                            parent.notify();
                        }
                    }
                }

                ClusterMethods.iterateOneRound(data, false);
                
                data.posarr=data.positions;
                tmpcool=(((float)((int)(data.currcool*100000)))/100000);
                if(tmpcool<=1e-5){
                    stop=true;
                }
            }
        }

    }

    void initialize(){
    	data.initialize();
    	setup_cutoffs();
    }

    void startstopthread(){
        //does the iteration and the stop for a thread
        if(mythread.didrun==false || mythread.stop==true){//if this thread was never started or stopped
            //String tmpstr="";
            if(data.movements==null){
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
}