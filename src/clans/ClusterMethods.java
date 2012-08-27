package clans;

/**
 *
 * @author tancred
 */
public class ClusterMethods {
    
	static final java.util.Random rand=new java.util.Random(System.currentTimeMillis());
    
	static void setup_attraction_values_and_initialize(ClusterData data){
	    // compute "attraction" values for all sequence pairs from the hsp objects
		// and initialize the positions randomly
		
        data.rounds=0;
        data.currcool=1;

		data.mymovearr = null;
        data.lastmovearr = null;       

        data.elements = data.namearr.length;
	    
        if(data.elements == 0){ // no elements means nothing to do
	    	data.myposarr = new float[0][0];
	    	data.posarr=data.myposarr;
	    	return;
	    }
	    
	    data.myposarr=new float[data.elements][ClusterData.dimensions];
	    data.posarrtmp=new float[data.elements][ClusterData.dimensions];
	    data.drawarrtmp=new int[data.elements][ClusterData.dimensions];
	    data.mymovearr=new float[data.elements][ClusterData.dimensions];
	    data.lastmovearr=new float[data.elements][ClusterData.dimensions];
	    for(int i = 0; i < data.elements; i++){
	        data.mymovearr[i][0]=0;
	        data.mymovearr[i][1]=0;
	        data.mymovearr[i][2]=0;
	    	
	        data.lastmovearr[i][0]=0;
	        data.lastmovearr[i][1]=0;
	        data.lastmovearr[i][2]=0;
	    }
	    
	    // compute the "attraction values"
	    if(data.blasthits != null){
	        if(data.myattvals == null){
	            //synchronized(myattvals){//myattvals is null here; cannot sync on it
	            data.myattvals=compute_attraction_values(data.blasthits, data.minpval, data);
	            //}
	        }
	    }

	    data.myposarr=initialize_positions_randomly(data.myposarr);
	    data.posarr=data.myposarr;
	}


	static minattvals[] compute_attraction_values(minhsp[] indata,double minpval,ClusterData data){
	    if(indata==null){//possible (if alternate data source was loaded)
	        System.out.println("indata is null");
	        return data.myattvals;
	    }
	    System.out.println("indata is size:"+java.lang.reflect.Array.getLength(indata));
	    java.util.ArrayList<minattvals> tmpvec=new java.util.ArrayList<minattvals>();
	    int datanum=java.lang.reflect.Array.getLength(indata);
	    java.util.HashMap myhash=new java.util.HashMap(datanum);
	    float newatt;
	    String key;
	    float maxattval=0;
	    minattvals curratt=null;
	    data.maxvalfound=0;//init to zero, is assigned value in getattvalsimple or mult
	    //NOTE: this is not necessarily a symmetrical array. compute all values
	    //and then symmetrize computing the average values
	    int elements=data.elements;
	    if(data.rescalepvalues==false){
	        //make the attraction values
	        if(data.attvalsimple){
	            for(int i=datanum;--i>=0;){
	                if(indata[i].query<indata[i].hit){
	                    key=indata[i].query+"_"+indata[i].hit;
	                }else{
	                    key=indata[i].hit+"_"+indata[i].query;
	                }
	                if(myhash.containsKey(key)){
	                    curratt=(minattvals)myhash.get(key);
	                    if(curratt.att==-1){
	                        //in this case keep the -1
	                    }else{
	                        newatt=getattvalsimple(indata[i].val,elements,minpval,data);
	                        if(newatt==-1){
	                            curratt.att=-1;
	                        }else{
	                            newatt/=2;
	                            curratt.att+=newatt;
	                        }
	                    }
	                }else{
	                    //if I've never encountered this query-hit pair before
	                    curratt=new minattvals();
	                    if(indata[i].query<indata[i].hit){
	                        curratt.query=indata[i].query;
	                        curratt.hit=indata[i].hit;
	                    }else{
	                        curratt.hit=indata[i].query;
	                        curratt.query=indata[i].hit;
	                    }
	                    curratt.att=getattvalsimple(indata[i].val,elements,minpval,data);
	                    if(curratt.att!=-1){
	                        curratt.att/=2;
	                    }
	                    if(curratt.att!=0){
	                        myhash.put(key,curratt);
	                        tmpvec.add(curratt);
	                    }
	                }
	                if(curratt.att>maxattval){
	                    maxattval=curratt.att;
	                }
	            }//end for i
	        }else{
	            for(int i=0;i<datanum;i++){
	                if(indata[i].query<indata[i].hit){
	                    key=indata[i].query+"_"+indata[i].hit;
	                }else{
	                    key=indata[i].hit+"_"+indata[i].query;
	                }
	                if(myhash.containsKey(key)){
	                    curratt=(minattvals)myhash.get(key);
	                    if(curratt.att==-1){
	                        //in this case keep the -1
	                    }else{
	                        newatt=getattvalmult(indata[i].val,elements,minpval,data);
	                        if(newatt==-1){
	                            curratt.att=-1;
	                        }else{
	                            newatt/=2;
	                            curratt.att+=newatt;
	                        }
	                    }
	                }else{
	                    //if I've never encountered this query-hit pair before
	                    curratt=new minattvals();
	                    if(indata[i].query<indata[i].hit){
	                        curratt.query=indata[i].query;
	                        curratt.hit=indata[i].hit;
	                    }else{
	                        curratt.hit=indata[i].query;
	                        curratt.query=indata[i].hit;
	                    }
	                    curratt.att=getattvalmult(indata[i].val,elements,minpval,data);
	                    if(curratt.att !=-1){
	                        curratt.att/=2;
	                    }
	                    if(curratt.att!=0){
	                        tmpvec.add(curratt);
	                        myhash.put(key,curratt);
	                    }
	                }
	                if(curratt.att>maxattval){
	                    maxattval=curratt.att;
	                }
	            }//end for i
	        }
	        //divide all vals by maxattval (-->range: 0-1)
	        //standard, just divide all values by the maximum value
	        //note, this does NOT symmetrize the attractions
	        if(data.usescval==false){
	            for(int i=tmpvec.size()-1;i>=0;i--){
	                if(((minattvals)tmpvec.get(i)).att==-1){
	                    ((minattvals)tmpvec.get(i)).att=1;
	                    //System.out.println(((minattvals)tmpvec.elementAt(i)).query+" "+((minattvals)tmpvec.elementAt(i)).hit+" :"+((minattvals)tmpvec.elementAt(i)).att);
	                }else{
	                    ((minattvals)tmpvec.get(i)).att/=maxattval;
	                    //System.out.println(((minattvals)tmpvec.elementAt(i)).query+" "+((minattvals)tmpvec.elementAt(i)).hit+" :"+((minattvals)tmpvec.elementAt(i)).att);
	                }
	            }// end for i
	            //System.out.println("maxattval"+maxattval+" offset="+0);
	            data.p2attfactor=maxattval;
	            data.p2attoffset=0;
	        }else{//if using scval
	            data.p2attfactor=1;
	            data.p2attoffset=0;
	        }
	    }else{//if rescalepvaluecheckbox==true
	        float minattval=java.lang.Float.MAX_VALUE;
	        //rescale the attraction values to range from 0 to 1 (with the smallest positive non-zero value as zero.
	        if(data.attvalsimple){
	            for(int i=0;i<datanum;i++){
	                if(indata[i].query<indata[i].hit){
	                    key=indata[i].query+"_"+indata[i].hit;
	                }else{
	                    key=indata[i].hit+"_"+indata[i].query;
	                }
	                if(myhash.containsKey(key)){
	                    curratt=(minattvals)myhash.get(key);
	                    if(curratt.att==-1){
	                        //in this case keep the -1
	                    }else{
	                        newatt=getattvalsimple(indata[i].val,elements,minpval,data);
	                        if(newatt==-1){
	                            curratt.att=-1;
	                        }else{
	                            newatt/=2;
	                            curratt.att+=newatt;
	                        }
	                    }
	                }else{
	                    //if I've never encountered this query-hit pair before
	                    curratt=new minattvals();
	                    if(indata[i].query<indata[i].hit){
	                        curratt.query=indata[i].query;
	                        curratt.hit=indata[i].hit;
	                    }else{
	                        curratt.hit=indata[i].query;
	                        curratt.query=indata[i].hit;
	                    }
	                    curratt.att=getattvalsimple(indata[i].val,elements,minpval,data);
	                    if(curratt.att!=-1){
	                        curratt.att/=2;
	                    }
	                    if(curratt.att!=0){
	                        myhash.put(key,curratt);
	                        tmpvec.add(curratt);
	                    }
	                }
	                if(curratt.att>maxattval){
	                    maxattval=curratt.att;
	                }
	                if((curratt.att>0)&&(curratt.att<minattval)){
	                    minattval=curratt.att;
	                }
	            }//end for i
	        }else{
	            for(int i=0;i<datanum;i++){
	                if(indata[i].query<indata[i].hit){
	                    key=indata[i].query+"_"+indata[i].hit;
	                }else{
	                    key=indata[i].hit+"_"+indata[i].query;
	                }
	                if(myhash.containsKey(key)){
	                    curratt=(minattvals)myhash.get(key);
	                    if(curratt.att==-1){
	                        //in this case keep the -1
	                    }else{
	                        newatt=getattvalmult(indata[i].val,elements,minpval,data);
	                        if(newatt==-1){
	                            curratt.att=-1;
	                        }else{
	                            newatt/=2;
	                            curratt.att+=newatt;
	                        }
	                    }
	
	                }else{
	                    //if I've never encountered this query-hit pair before
	                    curratt=new minattvals();
	                    if(indata[i].query<indata[i].hit){
	                        curratt.query=indata[i].query;
	                        curratt.hit=indata[i].hit;
	                    }else{
	                        curratt.hit=indata[i].query;
	                        curratt.query=indata[i].hit;
	                    }
	                    curratt.att=getattvalmult(indata[i].val,elements,minpval,data);
	                    if(curratt.att!=-1){
	                        curratt.att/=2;
	                    }
	                    if(curratt.att!=0){
	                        myhash.put(key,curratt);
	                        tmpvec.add(curratt);
	                    }
	
	                }
	                if(curratt.att>maxattval){
	                    maxattval=curratt.att;
	                }
	                if((curratt.att>0)&&(curratt.att<minattval)){
	                    minattval=curratt.att;
	                }
	            }//end for i
	        }
	        //and divide all vals by maxattval and offset by minattval(-->range: 0-1)
	        float divval=maxattval-minattval;
	        for(int i=tmpvec.size()-1;i>=0;i--){
	            if(((minattvals)tmpvec.get(i)).att==-1){
	                ((minattvals)tmpvec.get(i)).att=1;
	            }else{
	                ((minattvals)tmpvec.get(i)).att=(((minattvals)tmpvec.get(i)).att-minattval)/divval;
	            }
	        }// end for i
	        //System.out.println("maxattval"+maxattval+" offset="+minattval);
	        data.p2attfactor=divval;
	        data.p2attoffset=minattval;
	    }
	    minattvals[] retarr=(minattvals[])tmpvec.toArray(new minattvals[0]);
	    System.out.println("attvals size="+java.lang.reflect.Array.getLength(retarr));
	    return retarr;
	}// end getattvals




	//--------------------------------------------------------------------------
	
	static float[][] initialize_positions_randomly(float[][] positions){
	    //seed the positions array with random numbers([-1 to 1[)
	
	    for(int i = 0; i < positions.length; i++){
	    	for(int j = 0; j < positions[j].length; j++){
	    		positions[i][j] = rand.nextFloat() * 2 - 1;
	        }
	    }
	    
	    return positions;
	}

    //--------------------------------------------------------------------------

     static void recluster3d(ClusterData data){
        //take the hsp objects from indata and compute "attraction" values for all sequence pairs
        //once you have those try to cluster the data in 2d by "energy minimization" approach.
        //iterative approch, might want to specify maximum number of iterations
        //use the positions array and the attracion/repulsion values to
        //compute movement vectors for each object
        //long time=System.currentTimeMillis();
        System.arraycopy(data.mymovearr, 0, data.lastmovearr, 0, data.elements);//move all values from mymovearr to lastmovearr
        float[][] mymovearr=data.mymovearr;
        for(int i=data.elements;--i>=0;){
            //lastmovearr[i][0]=mymovearr[i][0];
            //lastmovearr[i][1]=mymovearr[i][1];
            //lastmovearr[i][2]=mymovearr[i][2];
            mymovearr[i][0]=0;
            mymovearr[i][1]=0;
            mymovearr[i][2]=0;
        }
        int selnamenum=java.lang.reflect.Array.getLength(data.selectednames);
        int[] selectnames=new int[selnamenum];//(int[])selectednames.clone();//a copy of the array that won't change during this round of calculations
        System.arraycopy(data.selectednames,0, data.selectnames,0, data.selnamenum);
        String syncme="syncme";
        synchronized(data.myattvals){
            if(data.cpu==1){
                //mymovearr=getmovement(myposarr,myattvals,mymovearr);
                getmovement(data.myposarr,data.myattvals,data.mymovearr,data.selectnames,data);
                domove(data.myposarr,data.mymovearr,data.selectnames,data);
            }else{
                //compute the attraction values using multiple threads
                int selectnamesnum=java.lang.reflect.Array.getLength(data.selectednames);
                if((data.moveselectedonly)&&(selectnamesnum>0)&&(selectnamesnum<data.elements)){
                    selectnamesnum=java.lang.reflect.Array.getLength(selectnames);
                    java.util.HashMap <String,Integer>tmphash=new java.util.HashMap<String,Integer>();
                    int threadnum=selectnamesnum/data.cpu;
                    int tmpval;
                    for(int i=selectnamesnum-1;i>=0;i--){
                        tmpval=(int)i/threadnum;
                        if(tmpval>data.cpu-1){
                            tmpval=data.cpu-1;
                        }
                        tmphash.put(String.valueOf(selectnames[i]), Integer.valueOf(tmpval));
                    }//end for i
                    for(int i=0;i<data.cpu;i++){
                        //movethreads[i]=new getmovethread(myposarr,myattvals,mymovearr,i,cpu,tmphash,selectnames,this);//start cpu threads to write the attraction values to mymovearr
                        data.movethreads[i]=new getmovethread(data.myposarr,data.myattvals,data.mymovearr,i,data.cpu,tmphash,selectnames,syncme,data);//start cpu threads to write the attraction values to mymovearr
                        data.movethreads[i].start();
                    }
                }else{
                    for(int i=0;i<data.cpu;i++){
                        //movethreads[i]=new getmovethread(myposarr,myattvals,mymovearr,i,cpu,this);//start cpu threads to write the attraction values to mymovearr
                        data.movethreads[i]=new getmovethread(data.myposarr,data.myattvals,data.mymovearr,i,data.cpu,syncme,data);//start cpu threads to write the attraction values to mymovearr
                        data.movethreads[i].start();
                    }
                }
                //now wait for all threads to finish
                boolean alldone=false;
                try{
                    //synchronized(this){
                    synchronized(syncme){
                        while (alldone==false){
                            alldone=true;
                            for(int i=0;i<data.cpu;i++){
                                if(data.movethreads[i].done==false){
                                    alldone=false;
                                    break;
                                }
                            }
                            if(alldone==false){//if I still have to wait for a thread
                                //this.wait();//release lock on this and wait
                                syncme.wait();//release lock on this and wait
                            }
                        }
                    }//end while alldone
                }catch (InterruptedException e){
                    System.err.println("Interrupted wait in searchblast");
                    e.printStackTrace();
                }
                domove(data.myposarr,data.mymovearr,selectnames,data);
            }//end else cpu>1
        }//end sync myattvals
        //then move each according to movevec
        //System.out.println("calctime="+(System.currentTimeMillis()-time));
        //return data.myposarr;
    }// end recluster3d

    //--------------------------------------------------------------------------

    static void getmovement(float[][] posarr, minattvals[] attvals, float[][] movement,int[] selectnames, ClusterData data){
        //use the positions of all elements and their attraction/repulsion values to
        //calculate a movement vector for each (take into account the last movement).
        //repulsion doesn't have a specific value as all evalues below a certain point
        //are simply regarded as insignificant. therefore use a one formula for all
        //to compute the repulsive forces (see getrepulse)
        //inline all the getrepulse and getattract methods to reduce the number of method calls
        int i,j;
        int attnum=java.lang.reflect.Array.getLength(attvals);
        double[] currmoverep=new double[3];
        double[] currmoveatt=new double[3];
        float tmp;
        double totaldist=0;
        double totalmove=1;
        double distx,disty,distz;
        int selectnamesnum=java.lang.reflect.Array.getLength(selectnames);
        float repfac=data.repfactor;
        //float attfac=attfactor;
        //int attpow=attvalpow;
        int reppow=data.repvalpow;
        int hnum,qnum;
        float weight1=1,weight2=1;
        double minattract=data.minattract;
        int repvalpow=data.repvalpow;
        float repfactor=data.repfactor;
        int attvalpow=data.attvalpow;
        float attfactor=data.attfactor;
        int elements=data.elements;
        float[] weights=data.weights;
        if((data.moveselectedonly)&&(selectnamesnum>0)&&(selectnamesnum!=elements)){
            java.util.HashMap tmphash=new java.util.HashMap<Integer,Integer>((int)(selectnamesnum/0.8)+1,0.8f);
            Integer[] hashkeys=new Integer[elements];
            //hashkey[] hashkeys=new hashkey[elements];
            for(i=elements;--i>=0;){
                //hashkeys[i]=new hashkey(i);
                hashkeys[i] = Integer.valueOf(i);
            }
            if(data.cluster2d){
                //no point in inlining this bit; I have very few selected sequences and therefore the gain should be marginal
                //cluster only the selected sequences in 2D
                for(i=selectnamesnum;--i>=0;){
                    //currmoveatt=getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                    getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                    weight1=1;
                    if(weights!=null){
                        weight1=weights[selectnames[i]];
                    }
                    tmphash.put(hashkeys[selectnames[i]],null);
                    movement[selectnames[i]][0]+=currmoveatt[0]*weight1;
                    movement[selectnames[i]][1]+=currmoveatt[1]*weight1;
                    for(j=elements;--j>=0;){
                        if(j==selectnames[i]){
                            continue;
                        }
                        //currmoverep=getrepulse2d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow, repfactor, rand);
                        getrepulse2d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow,repfactor, rand);
                        weight2=1;
                        if(weights!=null){
                            weight2=weights[selectnames[j]];
                        }
                        movement[selectnames[i]][0]+=currmoverep[0]*weight2;
                        movement[selectnames[i]][1]+=currmoverep[1]*weight2;
                    }//end for j
                }//end for i
                for(i=attnum;--i>=0;){
                    if(tmphash.containsKey(hashkeys[attvals[i].query])){
                        //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        weight1=1;
                        weight2=1;
                        if(weights!=null){
                            weight1=weights[attvals[i].query];
                            weight2=weights[attvals[i].hit];
                        }
                        movement[attvals[i].query][0]+=currmoveatt[0]*weight2;
                        movement[attvals[i].query][1]+=currmoveatt[1]*weight2;
                        if(tmphash.containsKey(hashkeys[attvals[i].hit])){
                            movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                            movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                        }
                    }else if(tmphash.containsKey(hashkeys[attvals[i].hit])){
                        //currmoveatt=getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        weight1=1;
                        if(weights!=null){
                            weight1=weights[attvals[i].query];
                        }
                        movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                        movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                    }
                }//end for i
                for(i=selectnamesnum;--i>=0;){
                    movement[selectnames[i]][0]/=elements;
                    movement[selectnames[i]][1]/=elements;
                    movement[selectnames[i]][2]=0;
                    totaldist=java.lang.Math.sqrt((movement[selectnames[i]][0]*movement[selectnames[i]][0])+(movement[selectnames[i]][1]*movement[selectnames[i]][1]));
                    if(totaldist>data.maxmove){
                        tmp=(float)(data.maxmove/totaldist);
                        movement[selectnames[i]][0]*=tmp;
                        movement[selectnames[i]][1]*=tmp;
                        //movement[selectednames[i]][2]*=maxmove/totaldist;
                    }
                }//end for i
            }else{
                //cluster only the selected sequences in 3D
                for(i=selectnamesnum;--i>=0;){
                    //currmoveatt=getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                    getminattract(posarr[selectnames[i]],currmoveatt,minattract);
                    weight1=1;
                    if(weights!=null){
                        weight1=weights[selectnames[i]];
                    }
                    movement[selectnames[i]][0]+=currmoveatt[0]*weight1;
                    movement[selectnames[i]][1]+=currmoveatt[1]*weight1;
                    movement[selectnames[i]][2]+=currmoveatt[2]*weight1;
                    tmphash.put(hashkeys[selectnames[i]],null);
                    for(j=elements;--j>=0;){
                        if(j==selectnames[i]){
                            continue;
                        }
                        //currmoverep=getrepulse3d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow, repfactor,rand);
                        getrepulse3d(posarr[selectnames[i]],posarr[j],currmoverep,repvalpow, repfactor,rand);
                        weight2=1;
                        if(weights!=null){
                            weight2=weights[selectnames[j]];
                        }
                        movement[selectnames[i]][0]+=currmoverep[0]*weight2;
                        movement[selectnames[i]][1]+=currmoverep[1]*weight2;
                        movement[selectnames[i]][2]+=currmoverep[2]*weight2;
                    }//end for j
                }//end for i
                for(i=attnum;--i>=0;){
                    if(tmphash.containsKey(hashkeys[attvals[i].query])){
                        //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        //no point in inlining this bit; I reuse the results a second time
                        getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        weight1=1;
                        weight2=1;
                        if(weights!=null){
                            weight1=weights[attvals[i].query];
                            weight2=weights[attvals[i].hit];
                        }
                        movement[attvals[i].query][0]+=currmoveatt[0]*weight2;
                        movement[attvals[i].query][1]+=currmoveatt[1]*weight2;
                        movement[attvals[i].query][2]+=currmoveatt[2]*weight2;
                        if(tmphash.containsKey(hashkeys[attvals[i].hit])){
                            movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                            movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                            movement[attvals[i].hit][2]-=currmoveatt[2]*weight1;
                        }
                    }else if(tmphash.containsKey(hashkeys[attvals[i].hit])){
                        //currmoveatt=getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        getattract3d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        weight1=1;
                        if(weights!=null){
                            weight1=weights[attvals[i].query];
                        }
                        movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                        movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                        movement[attvals[i].hit][2]-=currmoveatt[2]*weight1;
                    }
                }//end for i
                for(i=selectnamesnum;--i>=0;){
                    movement[selectnames[i]][0]/=elements;
                    movement[selectnames[i]][1]/=elements;
                    movement[selectnames[i]][2]/=elements;
                    totaldist=java.lang.Math.sqrt((movement[selectnames[i]][0]*movement[selectnames[i]][0])+(movement[selectnames[i]][1]*movement[selectnames[i]][1])+(movement[selectnames[i]][2]*movement[selectnames[i]][2]));
                    if(totaldist>data.maxmove){
                        tmp=(float)(data.maxmove/totaldist);
                        movement[selectnames[i]][0]*=tmp;
                        movement[selectnames[i]][1]*=tmp;
                        movement[selectnames[i]][2]*=tmp;
                    }
                }//end for i
            }
        }else{//if no sequences were selected or all should be used
            if(data.cluster2d){
                //cluster all in 2D
                if(weights!=null){
                    for(i=elements;--i>=0;){
                        weight1=weights[i];
                        movement[i][0]-=posarr[i][0]*minattract*weight1;
                        movement[i][1]-=posarr[i][1]*minattract*weight1;
                        for(j=i+1;j<elements;j++){
                            weight2=weights[j];
                            distx=posarr[j][0]-posarr[i][0];
                            disty=posarr[j][1]-posarr[i][1];
                            if(distx==0 && disty==0){
                                //if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfactor*(rand.nextDouble()-0.5)*0.001;
                                //return movement;
                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty));
                                //here I scale the force**repvalpow
                                totalmove=1;
                                for(int r=repvalpow;--r>=0;){
                                    totalmove*=totaldist;
                                }
                                totalmove=(repfactor/totalmove);
                                movement[i][0]+=(-distx/totaldist)*totalmove*weight2;
                                movement[i][1]+=(-disty/totaldist)*totalmove*weight2;
                                movement[j][0]-=(-distx/totaldist)*totalmove*weight1;
                                movement[j][1]-=(-disty/totaldist)*totalmove*weight1;
                            }
                        }//end for j
                    }//end for i
                    for(i=attnum;--i>=0;){
                        getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        weight1=weights[attvals[i].query];
                        weight2=weights[attvals[i].hit];
                        movement[attvals[i].query][0]+=currmoveatt[0]*weight2;
                        movement[attvals[i].query][1]+=currmoveatt[1]*weight2;
                        movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                        movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                    }//end for i
                }else{//if weights==null, the use a default weighting of 1
                    for(i=elements;--i>=0;){
                        movement[i][0]-=posarr[i][0]*minattract;
                        movement[i][1]-=posarr[i][1]*minattract;
                        for(j=i+1;j<elements;j++){
                            distx=posarr[j][0]-posarr[i][0];
                            disty=posarr[j][1]-posarr[i][1];
                            if(distx==0 && disty==0){
                                //if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfactor*(rand.nextDouble()-0.5)*0.001;
                                //return movement;
                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty));
                                //here I scale the force**repvalpow
                                totalmove=1;
                                for(int r=repvalpow;--r>=0;){
                                    totalmove*=totaldist;
                                }
                                totalmove=(repfactor/totalmove);
                                movement[i][0]+=(-distx/totaldist)*totalmove;
                                movement[i][1]+=(-disty/totaldist)*totalmove;
                                movement[j][0]-=(-distx/totaldist)*totalmove;
                                movement[j][1]-=(-disty/totaldist)*totalmove;
                            }
                        }//end for j
                    }//end for i
                    for(i=attnum;--i>=0;){
                        getattract2d(posarr[attvals[i].query],posarr[attvals[i].hit],attvals[i].att,currmoveatt,attvalpow,attfactor);
                        movement[attvals[i].query][0]+=currmoveatt[0];
                        movement[attvals[i].query][1]+=currmoveatt[1];
                        movement[attvals[i].hit][0]-=currmoveatt[0];
                        movement[attvals[i].hit][1]-=currmoveatt[1];
                    }//end for i
                }//end else weights==1
                for(i=elements;--i>=0;){
                    movement[i][0]/=elements;
                    movement[i][1]/=elements;
                    movement[i][2]=0;
                    totaldist=java.lang.Math.sqrt((movement[i][0]*movement[i][0])+(movement[i][1]*movement[i][1]));
                    if(totaldist>data.maxmove){
                        tmp=(float)(data.maxmove/totaldist);
                        movement[i][0]*=tmp;
                        movement[i][1]*=tmp;
                    }
                }//end for i
            }else{
                //cluster all in 3D
                if(weights!=null){
                    for(i=elements;--i>=0;){
                        weight1=weights[i];
                        movement[i][0]-=posarr[i][0]*minattract*weight1;
                        movement[i][1]-=posarr[i][1]*minattract*weight1;
                        movement[i][2]-=posarr[i][2]*minattract*weight1;
                        for(j=i;++j<elements;){
                            weight2=weights[j];
                            distx=posarr[j][0]-posarr[i][0];
                            disty=posarr[j][1]-posarr[i][1];
                            distz=posarr[j][2]-posarr[i][2];
                            if(distx==0 && disty==0 && distz==0){
                                //if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][2]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][2]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                //return movement;
                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
                                //here I scale the force**repvalpow
                                totalmove=1;
                                for(int r=reppow;--r>=0;){
                                    totalmove*=totaldist;
                                }
                                totalmove=(repfac/totalmove);
                                movement[i][0]+=(-distx/totaldist)*totalmove*weight2;
                                movement[i][1]+=(-disty/totaldist)*totalmove*weight2;
                                movement[i][2]+=(-distz/totaldist)*totalmove*weight2;
                                movement[j][0]-=(-distx/totaldist)*totalmove*weight1;
                                movement[j][1]-=(-disty/totaldist)*totalmove*weight1;
                                movement[j][2]-=(-distz/totaldist)*totalmove*weight1;
                            }
                        }//end for j
                    }//end for i
                    for(i=attnum;--i>=0;){
                        hnum=attvals[i].hit;
                        qnum=attvals[i].query;
                        weight1=weights[hnum];
                        weight2=weights[hnum];
                        distx=posarr[hnum][0]-posarr[qnum][0];
                        disty=posarr[hnum][1]-posarr[qnum][1];
                        distz=posarr[hnum][2]-posarr[qnum][2];
                        totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
                        //scale totalmove with distance**attvalpow
                        totalmove=1;
                        for(int r=attvalpow;--r>=0;){
                            totalmove*=totaldist;
                        }
                        totalmove=totalmove*attvals[i].att*attfactor;
                        if(attvals[i].att<0){
                            //in cae of repulsion I want to react inversely to the attractive forces
                            totalmove=-1/totalmove;
                        }
                        if(totaldist!=0){
                            movement[qnum][0]+=(distx/totaldist)*totalmove*weight2;
                            movement[qnum][1]+=(disty/totaldist)*totalmove*weight2;
                            movement[qnum][2]+=(distz/totaldist)*totalmove*weight2;
                            movement[hnum][0]-=(distx/totaldist)*totalmove*weight1;
                            movement[hnum][1]-=(disty/totaldist)*totalmove*weight1;
                            movement[hnum][2]-=(distz/totaldist)*totalmove*weight1;
                        }//else{
                    }//end for i
                }else{//use weights==1
                    for(i=elements;--i>=0;){
                        movement[i][0]-=posarr[i][0]*minattract;
                        movement[i][1]-=posarr[i][1]*minattract;
                        movement[i][2]-=posarr[i][2]*minattract;
                        for(j=i;++j<elements;){
                            distx=posarr[j][0]-posarr[i][0];
                            disty=posarr[j][1]-posarr[i][1];
                            distz=posarr[j][2]-posarr[i][2];
                            if(distx==0 && disty==0 && distz==0){
                                //if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][2]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][2]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                //return movement;
                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
                                //here I scale the force**repvalpow
                                totalmove=1;
                                for(int r=reppow;--r>=0;){
                                    totalmove*=totaldist;
                                }
                                totalmove=(repfac/totalmove);
                                movement[i][0]+=(-distx/totaldist)*totalmove;
                                movement[i][1]+=(-disty/totaldist)*totalmove;
                                movement[i][2]+=(-distz/totaldist)*totalmove;
                                movement[j][0]-=(-distx/totaldist)*totalmove;
                                movement[j][1]-=(-disty/totaldist)*totalmove;
                                movement[j][2]-=(-distz/totaldist)*totalmove;
                            }
                        }//end for j
                    }//end for i
                    for(i=attnum;--i>=0;){
                        hnum=attvals[i].hit;
                        qnum=attvals[i].query;
                        distx=posarr[hnum][0]-posarr[qnum][0];
                        disty=posarr[hnum][1]-posarr[qnum][1];
                        distz=posarr[hnum][2]-posarr[qnum][2];
                        totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
                        //scale totalmove with distance**attvalpow
                        totalmove=1;
                        for(int r=attvalpow;--r>=0;){
                            totalmove*=totaldist;
                        }
                        totalmove=totalmove*attvals[i].att*attfactor;
                        if(attvals[i].att<0){
                            //in cae of repulsion I want to react inversely to the attractive forces
                            totalmove=-1/totalmove;
                        }
                        if(totaldist!=0){
                            movement[qnum][0]+=(distx/totaldist)*totalmove;
                            movement[qnum][1]+=(disty/totaldist)*totalmove;
                            movement[qnum][2]+=(distz/totaldist)*totalmove;
                            movement[hnum][0]-=(distx/totaldist)*totalmove;
                            movement[hnum][1]-=(disty/totaldist)*totalmove;
                            movement[hnum][2]-=(distz/totaldist)*totalmove;
                        }
                    }//end for i
                }//end else weights!=null
                for(i=elements;--i>=0;){
                    movement[i][0]/=elements;
                    movement[i][1]/=elements;
                    movement[i][2]/=elements;
                    totaldist=java.lang.Math.sqrt((movement[i][0]*movement[i][0])+(movement[i][1]*movement[i][1])+(movement[i][2]*movement[i][2]));
                    if(totaldist>data.maxmove){
                        tmp=(float)(data.maxmove/totaldist);
                        movement[i][0]*=tmp;
                        movement[i][1]*=tmp;
                        movement[i][2]*=tmp;
                    }
                }//end for i
            }//end clustering all in 3D
        }//end else cluster only the selected sequences
    }// end getmovement

    //--------------------------------------------------------------------------

    static void domove(float[][] posarr, float[][] movearr,int[] selectnames,ClusterData data){
        //move all objects according to their movement vector.
        data.currcool=data.currcool*data.cooling;
        //float totalpos;
        float multdamp=1-data.dampening;
        int selectnamesnum=java.lang.reflect.Array.getLength(selectnames);
        double currcool=data.currcool;
        float[][] lastmovearr=data.lastmovearr;
        if((data.moveselectedonly)&&(selectnamesnum>0)&&(selectnamesnum!=data.elements)){
            //System.out.println("domoveselected");
            if(data.cluster2d){
                for(int i=selectnamesnum;--i>=0;){
                    posarr[selectnames[i]][0]+=(currcool*((lastmovearr[selectnames[i]][0]*multdamp)+(movearr[selectnames[i]][0])));
                    posarr[selectnames[i]][1]+=(currcool*((lastmovearr[selectnames[i]][1]*multdamp)+(movearr[selectnames[i]][1])));
                    //posarr[selectednames[i]][2]+=(currcool*((lastmovearr[selectednames[i]][2]*(1-dampening))+(movearr[selectednames[i]][2])));
                }//end for i
            }else{//cluster in 3d
                for(int i=selectnamesnum;--i>=0;){
                    posarr[selectnames[i]][0]+=(currcool*((lastmovearr[selectnames[i]][0]*multdamp)+(movearr[selectnames[i]][0])));
                    posarr[selectnames[i]][1]+=(currcool*((lastmovearr[selectnames[i]][1]*multdamp)+(movearr[selectnames[i]][1])));
                    posarr[selectnames[i]][2]+=(currcool*((lastmovearr[selectnames[i]][2]*multdamp)+(movearr[selectnames[i]][2])));
                }//end for i
            }
        }else{
            if(data.cluster2d){
                for(int i=data.elements;--i>=0;){
                    posarr[i][0]+=(currcool*((lastmovearr[i][0]*multdamp)+(movearr[i][0])));
                    posarr[i][1]+=(currcool*((lastmovearr[i][1]*multdamp)+(movearr[i][1])));
                    //posarr[i][2]+=(currcool*((lastmovearr[i][2]*(1-dampening))+(movearr[i][2])));
                }// end for i
            }else{//cluster in 3d
                for(int i=data.elements;--i>=0;){
                    posarr[i][0]+=(currcool*((lastmovearr[i][0]*multdamp)+(movearr[i][0])));
                    posarr[i][1]+=(currcool*((lastmovearr[i][1]*multdamp)+(movearr[i][1])));
                    posarr[i][2]+=(currcool*((lastmovearr[i][2]*multdamp)+(movearr[i][2])));
                }// end for i
            }
        }
        //return posarr;
    }// end domove

    //--------------------------------------------------------------------------

    static void getminattract2d(float[] pos1,double[] movement,double minattract){
        //get minimum attraction for 2d only
        movement[0]=-pos1[0]*minattract;
        movement[1]=-pos1[1]*minattract;
        movement[2]=0;
        //return movement;
    }//end getminattract2d

    //--------------------------------------------------------------------------

    //static double[] getminattract(float[] pos1, double[] movement,double minattract){
    static void getminattract(float[] pos1, double[] movement,double minattract){
        //which way is pos1 going to move if it is attracted by the origin
        movement[0]=-pos1[0]*minattract;
        movement[1]=-pos1[1]*minattract;
        movement[2]=-pos1[2]*minattract;
        //return movement;
    }// end getminattract

    //--------------------------------------------------------------------------

    static void getrepulse2d(float[] pos1,float[] pos2,double[] movement, int repvalpow, float repfactor, java.util.Random rand){
        //get the repulsion in 2d only
        double distx=pos2[0]-pos1[0];
        double disty=pos2[1]-pos1[1];
        double totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty));
        if(totaldist==0){
            //if two points are at exactly the same position I need to add some random effect
            movement[0]=repfactor*(rand.nextDouble()-0.5)*0.001;
            movement[1]=repfactor*(rand.nextDouble()-0.5)*0.001;
            movement[2]=0;//repfactor*(rand.nextDouble()-0.5)*0.001;
            return;
        }
        //here I scale the force**repvalpow
        double totalmove=1;
        for(int i=repvalpow;--i>=0;){
            totalmove*=totaldist;
        }
        totalmove=repfactor/totalmove;
        movement[0]=(-distx/totaldist)*totalmove;
        movement[1]=(-disty/totaldist)*totalmove;
        movement[2]=0;
        //return movement;
    }//end getrepulse2d

    //--------------------------------------------------------------------------

    static void getattract2d(float[] pos1, float[] pos2, float attval, double[] movement, int attvalpow, float attfactor){
        //get the attractive forces for 2d only (forget Z-axis)
        //tmpattvals are between 0 and 1 (or ==2 for evalue==0) (o=no attraction, 1=max attraction)
        double distx=pos2[0]-pos1[0];
        double disty=pos2[1]-pos1[1];
        double totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty));
        //scale totalmove with distance**attvalpow
        double totalmove=1;
        for(int i=attvalpow;--i>=0;){
            totalmove*=totaldist;
        }
        totalmove=totalmove*attval*attfactor;
        if(attval<0){
            //in cae of repulsion I want to react inversely to the attractive forces
            totalmove=-1/totalmove;
        }
        if(totaldist!=0){
            movement[0]=(distx/totaldist)*totalmove;
            movement[1]=(disty/totaldist)*totalmove;
            movement[2]=0;
        }else{
            //repulsion values will differentially move them
            movement[0]=0;
            movement[1]=0;
            movement[2]=0;
        }
        //return movement;
    }//end getattract2d

    //--------------------------------------------------------------------------

    //static double[] getattract3d(float[] pos1, float[] pos2, float attval,double[] movement,int attvalpow, float attfactor){
    static void getattract3d(float[] pos1, float[] pos2, float attval,double[] movement,int attvalpow, float attfactor){
        //similar to getrepulse but this is an attractive force that scales with distance**2 (the further away, the greater)
        //which way is pos1 going to move, given pos2
        //tmpattvals are between 0 and 1 (or ==2 for evalue==0) (o=no attraction, 1=max attraction)
        double distx=pos2[0]-pos1[0];
        double disty=pos2[1]-pos1[1];
        double distz=pos2[2]-pos1[2];
        double totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
        //scale totalmove with distance**attvalpow
        double totalmove=1;
        for(int i=attvalpow;--i>=0;){
            totalmove*=totaldist;
        }
        totalmove=totalmove*attval*attfactor;
        if(attval<0){
            //in cae of repulsion I want to react inversely to the attractive forces
            totalmove=-1/totalmove;
        }
        if(totaldist!=0){
            movement[0]=(distx/totaldist)*totalmove;
            movement[1]=(disty/totaldist)*totalmove;
            movement[2]=(distz/totaldist)*totalmove;
        }else{
            //the repulsion values will differentially move them
            movement[0]=0;
            movement[1]=0;
            movement[2]=0;
        }
        //return movement;
    }// end getattract3d

    //--------------------------------------------------------------------------

    //static double[] getrepulse3d(float[] pos1,float[] pos2,double[] movement, int repvalpow, float repfactor, Random rand){
    static void getrepulse3d(float[] pos1,float[] pos2,double[] movement, int repvalpow, float repfactor, java.util.Random rand){
        //given these two objects, which way will object 1 move? force scales with 1/distance**2
        double distx=pos2[0]-pos1[0];
        double disty=pos2[1]-pos1[1];
        double distz=pos2[2]-pos1[2];
        double totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));
        if(totaldist==0){
            //if two points are at exactly the same position I need to add some random effect
            movement[0]=repfactor*(rand.nextDouble()-0.5)*0.001;
            movement[1]=repfactor*(rand.nextDouble()-0.5)*0.001;
            movement[2]=repfactor*(rand.nextDouble()-0.5)*0.001;
            return;
        }
        //here I scale the force**repvalpow
        double totalmove=1;
        for(int i=repvalpow;--i>=0;){
            totalmove*=totaldist;
        }
        totalmove=repfactor/totalmove;
        movement[0]=(-distx/totaldist)*totalmove;
        movement[1]=(-disty/totaldist)*totalmove;
        movement[2]=(-distz/totaldist)*totalmove;
        //return movement;
    }// end getrepulse3d

    static minattvals[] filter_attraction_values(minattvals[] attractions, double minpval){
        //use for filtering attraction values [-1<0<1]; -1 and 1 are max repulse/attract
        java.util.ArrayList <minattvals>retvec=new java.util.ArrayList<minattvals>();
        for(int i = 0; i < attractions.length; i++){
            if((attractions[i].att >= minpval) || (attractions[i].att <= -minpval)){
                retvec.add(attractions[i]);
            }
        }//end for i
        minattvals[] retarr=(minattvals[])retvec.toArray(new minattvals[0]);
        return retarr;
    }//end getattvals

    //--------------------------------------------------------------------------

    

    //--------------------------------------------------------------------------

    static float getattvalsimple(double[] invec,int dbsize,double minpval,ClusterData data){
        //System.out.println("simple");
        //this actually takes all data from vector and makes ONE number out of it
        //just use the BEST value!
        if(invec==null){
            return 0;
        }else if(java.lang.reflect.Array.getLength(invec)<1){//if I have no hits
            return 0;
        }
        double bestval=invec[0];//is a p-value (should be from 0 to 1)
        double currval;
        if(data.usescval){
            bestval=0;
            for(int i=java.lang.reflect.Array.getLength(invec)-1;i>=0;i--){
                currval=invec[i];
                if(currval>bestval){
                    bestval=currval;
                }
            }// end for i
            if(bestval<data.maxvalfound){//maxvalfound=worst accepted value (comes from P-values where larger=worse)
                data.maxvalfound=bestval;
            }
            currval=bestval;
            if(currval<minpval){//minpval also functions as minscoreval here
                //System.out.println(" currval="+currval+" is less than minpval="+minpval+" returning zero");
                return 0;
            }
        }else{
            for(int i=java.lang.reflect.Array.getLength(invec)-1;i>=1;i--){
                currval=invec[i];
                if(currval<bestval){
                    bestval=currval;
                }
            }// end for i
            if(bestval>data.maxvalfound){
                data.maxvalfound=bestval;
            }
            if(bestval==0){
                return -1;//this is identity
            }else if(bestval>1){//should never happen to p-values!
                return 0;
            }else if(bestval>minpval){//if this value is worse than permitted
                return 0;
            }
            //now all pvalues between 0 and 1
            currval=(-1*java.lang.Math.log(bestval));//ln10;//don't need it here as I convert all to relative attractions
        }
        return (float)(currval);
    }// end gettaval

    //--------------------------------------------------------------------------

    static float getattvalmult(double[] invec,int dbsize,double minpval,ClusterData data){
        //System.out.println("complex");
        //this actually takes all data from vector and makes ONE number out of it
        //new: multiply the pvalues of different hsp's
        if(invec==null){
            return 0;
        }else if(java.lang.reflect.Array.getLength(invec)<1){
            return 0;
        }
        double currval=invec[0];
        if(data.usescval){
            //then I am using score values (no logarithming at the end)
            currval=0;
            for(int i=java.lang.reflect.Array.getLength(invec);--i>=0;){//sum the values
                currval+=invec[i];
            }// end for i
            if(currval<data.maxvalfound){//maxvalfound=worst accepted value (comes from P-values where larger=worse)
                data.maxvalfound=currval;
            }
            if(currval<minpval){//minpval also functions as minscoreval here
                return 0;
            }
        }else{//then I am using P-values
            for(int i=java.lang.reflect.Array.getLength(invec);--i>=1;){
                currval*=invec[i];
            }// end for i
            if(currval>data.maxvalfound){
                data.maxvalfound=currval;
            }
            if(currval==0){
                return -1;//this is identity
            }else if(currval>1){//should never happen to p-values!
                return 0;
            }else if(currval>minpval){//if this value is worse than permitted
                return 0;
            }
            //now all pvalues between 0 and 1
            currval=(-1*java.lang.Math.log(currval));///ln10;
        }
        return (float)(currval);
    }// end gettaval

    //--------------------------------------------------------------------------

    static AminoAcidSequence[] remove_gaps_from_sequences(AminoAcidSequence[] inseqs){
        for(int i=0; i < inseqs.length; i++){
            inseqs[i].seq=inseqs[i].seq.replaceAll("-","");
        }
        return inseqs;
    }
}
