package clans.algorithms;

import java.util.Arrays;

import clans.model.ClusterData;
import clans.model.proteins.AminoAcidSequence;
import clans.model.proteins.MinimalAttractionValue;

public class ClusterMethods {

	public static final java.util.Random rand = new java.util.Random(System.currentTimeMillis());

	/**
	 * take the hsp objects from indata and compute "attraction" values for all sequence pairs once you have those try
	 * to cluster the data in 2d by "energy minimization" approach. iterative approch, might want to specify maximum
	 * number of iterations use the positions array and the attracion/repulsion values to compute movement vectors for
	 * each object long time=System.currentTimeMillis();
	 * 
	 * @param data
	 *            The data for which an interation should be computed.
	 * @param move_only_selected_sequences
	 *            If true, only the selected sequences are moved, if false all are moved. If no sequences or all are
	 *            selected this parameter is ignored.
	 */
	public static void iterateOneRound(ClusterData data, boolean move_only_selected_sequences) {
		
		shiftMovementsToLastMovements(data);
		
        synchronized(data.attractionValues){
			if (data.cpu == 1) {
				computeOnSingleCpu(data, move_only_selected_sequences);

			} else {
				computeOnManyCpus(data, move_only_selected_sequences);
            }
        }
    }
	
	/**
	 * Backup the movements calculated for the previous round as they are used, together with the ones from this round,
	 * in calculating the complete movements. Then resets the current movements to allow starting computations.
	 * 
	 * @param data
	 *            The data.
	 */
	private static void shiftMovementsToLastMovements(ClusterData data) {
		// backup the old movements to later use them along the new ones to arrive at this iteration's movements
        System.arraycopy(data.movements, 0, data.movementsLastIteration, 0, data.elements);

        float[][] mymovearr = data.movements;
        for (int i=0; i<mymovearr.length; i++) {
        	Arrays.fill(mymovearr[i], 0f);
        }
	}

	/**
	 * Starts the computations in a single thread (this one).
	 * 
	 * @param data
	 *            The data.
	 * @param move_only_selected_sequences
	 *            If true, only the selected sequences are moved, if false all are moved. If no sequences or all are
	 *            selected this parameter is ignored.
	 */
	protected static void computeOnSingleCpu(ClusterData data, boolean move_only_selected_sequences) {
		computeMovments(data, move_only_selected_sequences);
		applyMovement(data, move_only_selected_sequences);
	}
	
	/**
	 * Starts the computations in {@code data.cpu} many threads.
	 * 
	 * @param data
	 *            The data.
	 * @param move_only_selected_sequences
	 *            If true, only the selected sequences are moved, if false all are moved. If no sequences or all are
	 *            selected this parameter is ignored.
	 */
	protected static void computeOnManyCpus(ClusterData data, boolean move_only_selected_sequences) {
		final Object computationThreadLock = new Object();
	
        if (move_only_selected_sequences) {

			// assign each sequence to a thread
			java.util.HashMap<String, Integer> sequence_thread_mapping = new java.util.HashMap<String, Integer>();
			int sequences_per_threads = data.selectedSequencesIndicesStableCopy.length / data.cpu;
			int thread_for_this_sequence;

			for (int i = 0; i < data.selectedSequencesIndicesStableCopy.length; i++) {
				thread_for_this_sequence = i / sequences_per_threads;
				if (thread_for_this_sequence > data.cpu - 1) {
					thread_for_this_sequence = data.cpu - 1;
				}
				sequence_thread_mapping.put(String.valueOf(data.selectedSequencesIndicesStableCopy[i]), Integer.valueOf(thread_for_this_sequence));
			}

			for (int i = 0; i < data.cpu; i++) {
				data.movethreads[i] = new MovementComputerThread(data, i, computationThreadLock,
						sequence_thread_mapping);
				data.movethreads[i].start();
			}
			
        } else {
            for (int i = 0; i < data.cpu; i++) {
                data.movethreads[i] = new MovementComputerThread(data, i, computationThreadLock);
                data.movethreads[i].start();
            }
        }

        /**
		 * This method contains a wait() that will notice interrupt() on this thread and treat it as signal to
		 * stop after completing this round of computation.
		 */
		boolean stop_after_this_round = waitForThreadsToFinish(data.movethreads, computationThreadLock);

        applyMovement(data, move_only_selected_sequences);
        
		// now that the round is completed, escalate the interrupt to indicate "stop ASAP" to the caller
        if (stop_after_this_round) {
        	Thread.currentThread().interrupt();
        }
	}
	
	/**
	 * Waits for all child threads to finish their computations, then tells the caller whether computations should be
	 * stopped. Necessity to stop is triggered by an interrupt on this thread.
	 * 
	 * @param child_threads
	 *            The child threads.
	 * @param computationThreadLock
	 *            Synchronized on while waiting. This instance is used by the children to notify of completed
	 *            calculations.
	 * @return true if the caller is supposed to stop the computations ASAP, false if computations can continue.
	 */
	static boolean waitForThreadsToFinish(MovementComputerThread[] child_threads, final Object computationThreadLock) {
		boolean alldone = false;
		boolean stop_after_this_round = false;
		
		// this sync object is notified when a child thread finishes its computation which resumes at sync.wait() below
		synchronized (computationThreadLock) {
			
			while (!alldone) {
				alldone = true;

				for (MovementComputerThread child_thread: child_threads) {
					if (!child_thread.done) {
						alldone = false;
						break;
					}
				}

				if (alldone) {
					break;
				}
				
				try {
					/**
					 * wait on the sync object pauses this thread until a computation is finished in a child thread,
					 * which reawakens this thread by sync.notify()
					 */
					computationThreadLock.wait();

				} catch (InterruptedException e) {
					/**
					 * We want the InterruptedException to inform the thread to stop iterating after this round.
					 * However, when .wait() throws this Exception, the interrupted state of this thread is reset to
					 * false. As we want to finish this round, we later need to tell the thread to stop.
					 */
					stop_after_this_round = true;
				}
			}
		}
		
		return stop_after_this_round;
	}

    /**
     * use the positions of all elements and their attraction/repulsion values to calculate a movement vector for each
     * (take into account the last movement). repulsion doesn't have a specific value as all evalues below a certain
     * point are simply regarded as insignificant. therefore use a one formula for all to compute the repulsive forces
     * inline all the getrepulse and getattract methods to reduce the number of method calls
     * 
     * @param posarr
     * @param attvals
     * @param movement
     * @param selectedSequencesIndicesStableCopy
     * @param data
     */
	static void computeMovments(ClusterData data, boolean move_only_selected_sequences) {
    	
    	float[][] posarr = data.positions;
    	MinimalAttractionValue[] attvals = data.attractionValues;
    	float[][] movement = data.movements;
    	int[] selected_sequences_indices = data.selectedSequencesIndicesStableCopy; 

        int i,j;
        int attnum=attvals.length;
        double[] currmoverep=new double[3];
        double[] currmoveatt=new double[3];
        float tmp;
        double totaldist=0;
        double totalmove=1;
        double distx,disty,distz;
        float repfac=data.repfactor;
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
        
		if (move_only_selected_sequences) {
			java.util.HashMap<Integer, Integer> tmphash = new java.util.HashMap<Integer, Integer>(
					(int) (selected_sequences_indices.length / 0.8) + 1, 0.8f);

            Integer[] hashkeys=new Integer[elements];
            for(i=elements;--i>=0;){
                hashkeys[i] = Integer.valueOf(i);
            }
            
            if(data.cluster2d){
				// no point in inlining this bit; I have very few selected sequences and the gain should be marginal
				// cluster only the selected sequences in 2D
				for (i = selected_sequences_indices.length; --i >= 0;) {
                    addAttractionTowardsOrigin(posarr[selected_sequences_indices[i]],currmoveatt,minattract);

                    weight1=1;
                    if(weights!=null){
                        weight1=weights[selected_sequences_indices[i]];
                    }
                    
                    tmphash.put(hashkeys[selected_sequences_indices[i]],null);
                    movement[selected_sequences_indices[i]][0]+=currmoveatt[0]*weight1;
                    movement[selected_sequences_indices[i]][1]+=currmoveatt[1]*weight1;
                    for(j=elements;--j>=0;){
                        if(j==selected_sequences_indices[i]){
                            continue;
                        }

                        addRepulsion(posarr[selected_sequences_indices[i]], posarr[j], currmoverep, repvalpow, repfactor, rand,
                                data.cluster2d);
                        
                        weight2=1;
                        if(weights!=null){
                            weight2=weights[selected_sequences_indices[j]];
                        }
                        
                        movement[selected_sequences_indices[i]][0]+=currmoverep[0]*weight2;
                        movement[selected_sequences_indices[i]][1]+=currmoverep[1]*weight2;
                    }
                }
                
                for(i=attnum;--i>=0;){
                    if(tmphash.containsKey(hashkeys[attvals[i].query])){
                        addAttraction(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt,
                                attvalpow, attfactor, data.cluster2d);

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
                        addAttraction(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt,
                                attvalpow, attfactor, data.cluster2d);

                        weight1=1;
                        if(weights!=null){
                            weight1=weights[attvals[i].query];
                        }
                        movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                        movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                    }
                }
				for (i = selected_sequences_indices.length; --i >= 0;) {
                    movement[selected_sequences_indices[i]][0]/=elements;
                    movement[selected_sequences_indices[i]][1]/=elements;
                    movement[selected_sequences_indices[i]][2]=0;
                    totaldist=java.lang.Math.sqrt((movement[selected_sequences_indices[i]][0]*movement[selected_sequences_indices[i]][0])+(movement[selected_sequences_indices[i]][1]*movement[selected_sequences_indices[i]][1]));
                    if(totaldist>data.maxmove){
                        tmp=(float)(data.maxmove/totaldist);
                        movement[selected_sequences_indices[i]][0]*=tmp;
                        movement[selected_sequences_indices[i]][1]*=tmp;
                        //movement[selectednames[i]][2]*=maxmove/totaldist;
                    }
                }
            }else{
                // cluster only the selected sequences in 3D
				for (i = selected_sequences_indices.length; --i >= 0;) {
                    addAttractionTowardsOrigin(posarr[selected_sequences_indices[i]],currmoveatt,minattract);
                    weight1=1;
                    if(weights!=null){
                        weight1=weights[selected_sequences_indices[i]];
                    }
                    movement[selected_sequences_indices[i]][0]+=currmoveatt[0]*weight1;
                    movement[selected_sequences_indices[i]][1]+=currmoveatt[1]*weight1;
                    movement[selected_sequences_indices[i]][2]+=currmoveatt[2]*weight1;
                    tmphash.put(hashkeys[selected_sequences_indices[i]],null);
                    for(j=elements;--j>=0;){
                        if(j==selected_sequences_indices[i]){
                            continue;
                        }

                        addRepulsion(posarr[selected_sequences_indices[i]],posarr[j],currmoverep,repvalpow, repfactor,rand, data.cluster2d);
                        
                        weight2=1;
                        if(weights!=null){
                            weight2=weights[selected_sequences_indices[j]];
                        }
                        movement[selected_sequences_indices[i]][0]+=currmoverep[0]*weight2;
                        movement[selected_sequences_indices[i]][1]+=currmoverep[1]*weight2;
                        movement[selected_sequences_indices[i]][2]+=currmoverep[2]*weight2;
                    }
                }
                
                for(i=attnum;--i>=0;){
                    if(tmphash.containsKey(hashkeys[attvals[i].query])){
                        // no point in inlining this bit; I reuse the results a second time
                        addAttraction(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt,
                                attvalpow, attfactor, data.cluster2d);
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
                        addAttraction(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt,
                                attvalpow, attfactor, data.cluster2d);
                        weight1=1;
                        if(weights!=null){
                            weight1=weights[attvals[i].query];
                        }
                        movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                        movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                        movement[attvals[i].hit][2]-=currmoveatt[2]*weight1;
                    }
                }
				for (i = selected_sequences_indices.length; --i >= 0;) {
                    movement[selected_sequences_indices[i]][0]/=elements;
                    movement[selected_sequences_indices[i]][1]/=elements;
                    movement[selected_sequences_indices[i]][2]/=elements;
                    totaldist=java.lang.Math.sqrt((movement[selected_sequences_indices[i]][0]*movement[selected_sequences_indices[i]][0])+(movement[selected_sequences_indices[i]][1]*movement[selected_sequences_indices[i]][1])+(movement[selected_sequences_indices[i]][2]*movement[selected_sequences_indices[i]][2]));
                    if(totaldist>data.maxmove){
                        tmp=(float)(data.maxmove/totaldist);
                        movement[selected_sequences_indices[i]][0]*=tmp;
                        movement[selected_sequences_indices[i]][1]*=tmp;
                        movement[selected_sequences_indices[i]][2]*=tmp;
                    }
                }
            }
        
		} else { // move all sequences
            if(data.cluster2d){
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
                                // if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfactor*(rand.nextDouble()-0.5)*0.001;

                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty));

                                // here I scale the force**repvalpow
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
                        }
                    }
                    
                    for(i=attnum;--i>=0;){
                        addAttraction(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt,
                                attvalpow, attfactor, data.cluster2d);
                        
                        weight1=weights[attvals[i].query];
                        weight2=weights[attvals[i].hit];
                        movement[attvals[i].query][0]+=currmoveatt[0]*weight2;
                        movement[attvals[i].query][1]+=currmoveatt[1]*weight2;
                        movement[attvals[i].hit][0]-=currmoveatt[0]*weight1;
                        movement[attvals[i].hit][1]-=currmoveatt[1]*weight1;
                    }
                    
                } else { // weights==null, use default weight of 1
                    for(i=elements;--i>=0;){
                        movement[i][0]-=posarr[i][0]*minattract;
                        movement[i][1]-=posarr[i][1]*minattract;
                        for(j=i+1;j<elements;j++){
                            distx=posarr[j][0]-posarr[i][0];
                            disty=posarr[j][1]-posarr[i][1];

                            if(distx==0 && disty==0){
                                // if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfactor*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfactor*(rand.nextDouble()-0.5)*0.001;

                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty));
                                // here I scale the force**repvalpow
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
                        }
                    }
                    for(i=attnum;--i>=0;){
                        addAttraction(posarr[attvals[i].query], posarr[attvals[i].hit], attvals[i].att, currmoveatt,
                                attvalpow, attfactor, data.cluster2d);
                        
                        movement[attvals[i].query][0]+=currmoveatt[0];
                        movement[attvals[i].query][1]+=currmoveatt[1];
                        movement[attvals[i].hit][0]-=currmoveatt[0];
                        movement[attvals[i].hit][1]-=currmoveatt[1];
                    }
                }
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
                }
			} else { // cluster in 3D
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
                                // if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][2]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][2]-=repfac*(rand.nextDouble()-0.5)*0.001;
                            
                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));

                                // here I scale the force**repvalpow
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
                        }
                    }
                    for(i=attnum;--i>=0;){
                        hnum=attvals[i].hit;
                        qnum=attvals[i].query;
                        weight1=weights[hnum];
                        weight2=weights[hnum];
                        distx=posarr[hnum][0]-posarr[qnum][0];
                        disty=posarr[hnum][1]-posarr[qnum][1];
                        distz=posarr[hnum][2]-posarr[qnum][2];
                        totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));

                        // scale totalmove with distance**attvalpow
                        totalmove=1;
                        for(int r=attvalpow;--r>=0;){
                            totalmove*=totaldist;
                        }
                        
                        totalmove=totalmove*attvals[i].att*attfactor;
                        if(attvals[i].att<0){
                            // in case of repulsion I want to react inversely to the attractive forces
                            totalmove=-1/totalmove;
                        }
                        if(totaldist!=0){
                            movement[qnum][0]+=(distx/totaldist)*totalmove*weight2;
                            movement[qnum][1]+=(disty/totaldist)*totalmove*weight2;
                            movement[qnum][2]+=(distz/totaldist)*totalmove*weight2;
                            movement[hnum][0]-=(distx/totaldist)*totalmove*weight1;
                            movement[hnum][1]-=(disty/totaldist)*totalmove*weight1;
                            movement[hnum][2]-=(distz/totaldist)*totalmove*weight1;
                        }
                    }
				} else { // use weights==1
                    for(i=elements;--i>=0;){
                        movement[i][0]-=posarr[i][0]*minattract;
                        movement[i][1]-=posarr[i][1]*minattract;
                        movement[i][2]-=posarr[i][2]*minattract;
                        for(j=i;++j<elements;){
                            distx=posarr[j][0]-posarr[i][0];
                            disty=posarr[j][1]-posarr[i][1];
                            distz=posarr[j][2]-posarr[i][2];
                            
                            if(distx==0 && disty==0 && distz==0){
                                // if two points are at exactly the same position I need to add some random effect
                                movement[i][0]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][1]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[i][2]+=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][0]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][1]-=repfac*(rand.nextDouble()-0.5)*0.001;
                                movement[j][2]-=repfac*(rand.nextDouble()-0.5)*0.001;

                            }else{
                                totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));

                                // here I scale the force**repvalpow
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
                        }
                    }
                    for(i=attnum;--i>=0;){
                        hnum=attvals[i].hit;
                        qnum=attvals[i].query;
                        distx=posarr[hnum][0]-posarr[qnum][0];
                        disty=posarr[hnum][1]-posarr[qnum][1];
                        distz=posarr[hnum][2]-posarr[qnum][2];
                        totaldist=java.lang.Math.sqrt((distx*distx)+(disty*disty)+(distz*distz));

                        // scale totalmove with distance**attvalpow
                        totalmove=1;
                        for(int r=attvalpow;--r>=0;){
                            totalmove*=totaldist;
                        }
                        
                        totalmove=totalmove*attvals[i].att*attfactor;
                        if(attvals[i].att<0){
                            // in case of repulsion I want to react inversely to the attractive forces
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
                    }
                }
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
                }
            }
        }
    }

	/**
	 * Moves all points according to the computed movement vectors.
	 * 
	 * @param data
	 *            The data.
	 * @param move_only_selected_sequences
	 *            Whether only selected sequences should be moved.
	 */
	static void applyMovement(ClusterData data, boolean move_only_selected_sequences) {
    		
		float[][] posarr = data.positions;
		float[][] movearr = data.movements;
		int[] selected_sequences_indices = data.selectedSequencesIndicesStableCopy;

		// move all objects according to their movement vector.
		data.currcool = data.currcool * data.cooling;

		float multdamp = 1 - data.dampening;
		double currcool = data.currcool;
		float[][] lastmovearr = data.movementsLastIteration;

        if (move_only_selected_sequences) {
			if (data.cluster2d) {
				for (int i = selected_sequences_indices.length; --i >= 0;) {
					posarr[selected_sequences_indices[i]][0] += (currcool * ((lastmovearr[selected_sequences_indices[i]][0] * multdamp) + (movearr[selected_sequences_indices[i]][0])));
					posarr[selected_sequences_indices[i]][1] += (currcool * ((lastmovearr[selected_sequences_indices[i]][1] * multdamp) + (movearr[selected_sequences_indices[i]][1])));
					// posarr[selectednames[i]][2]+=(currcool*((lastmovearr[selectednames[i]][2]*(1-dampening))+(movearr[selectednames[i]][2])));
				}

			} else { // cluster in 3D
				for (int i = selected_sequences_indices.length; --i >= 0;) {
					posarr[selected_sequences_indices[i]][0] += (currcool * ((lastmovearr[selected_sequences_indices[i]][0] * multdamp) + (movearr[selected_sequences_indices[i]][0])));
					posarr[selected_sequences_indices[i]][1] += (currcool * ((lastmovearr[selected_sequences_indices[i]][1] * multdamp) + (movearr[selected_sequences_indices[i]][1])));
					posarr[selected_sequences_indices[i]][2] += (currcool * ((lastmovearr[selected_sequences_indices[i]][2] * multdamp) + (movearr[selected_sequences_indices[i]][2])));
				}
			}
  
		} else { // move all sequences
			if (data.cluster2d) {
                for(int i=data.elements;--i>=0;){
                    posarr[i][0]+=(currcool*((lastmovearr[i][0]*multdamp)+(movearr[i][0])));
                    posarr[i][1]+=(currcool*((lastmovearr[i][1]*multdamp)+(movearr[i][1])));
                    //posarr[i][2]+=(currcool*((lastmovearr[i][2]*(1-dampening))+(movearr[i][2])));
                }

			} else { // cluster in 3D
                for(int i=data.elements;--i>=0;){
                    posarr[i][0]+=(currcool*((lastmovearr[i][0]*multdamp)+(movearr[i][0])));
                    posarr[i][1]+=(currcool*((lastmovearr[i][1]*multdamp)+(movearr[i][1])));
                    posarr[i][2]+=(currcool*((lastmovearr[i][2]*multdamp)+(movearr[i][2])));
                }
            }
        }
    }
    
    /**
     * Adds the attractions of the origin of the coordinate system to the points' movement. 
     * 
     * @param position the position of the point
     * @param movement the direction of movement so far computed for this point
     * @param origin_attraction the degree of attraction the origin has on a point
     */
	static void addAttractionTowardsOrigin(float[] position, double[] movement, double origin_attraction) {
        movement[0] = -position[0] * origin_attraction;
        movement[1] = -position[1] * origin_attraction;
        movement[2] = -position[2] * origin_attraction;
    }

    /**
	 * Adds the attractive force that scales with distance**2. Determines which way pos1 is going to move, given pos2.
	 * 
	 * @param position_1
	 * @param position_2
	 * @param attraction_value
	 * @param movement
	 * @param attraction_exponent
	 * @param attraction_factor
	 * @param is_2d
	 */
	static void addAttraction(float[] position_1, float[] position_2, float attraction_value, double[] movement,
			int attraction_exponent, float attraction_factor, boolean is_2d) {

        double distx = position_2[0] - position_1[0];
        double disty = position_2[1] - position_1[1];
        double distz = position_2[2] - position_1[2];

        double totaldist = Math.sqrt((distx * distx) + (disty * disty) + (distz * distz));

        double totalmove = attraction_value * attraction_factor * Math.pow(totaldist, attraction_exponent);

        if (attraction_value < 0) { // in case of repulsion act inversely to the attractive forces
            totalmove = -1 / totalmove;
        }

        if (totaldist != 0) {
            movement[0] = (distx / totaldist) * totalmove;
            movement[1] = (disty / totaldist) * totalmove;
            movement[2] = (distz / totaldist) * totalmove;

        } else {
            // attraction does not work on these points if their at the same location
            movement[0] = 0;
            movement[1] = 0;
            movement[2] = 0;
        }
    }

    /**
	 * Adds the repulsive force between objects at position_1 and position_2. Force scales with 1/distance**2.
	 * 
	 * @param position_1
	 *            position of the first point
	 * @param position_2
	 *            position of the second point
	 * @param movement
	 * @param repulsion_exponent
	 * @param repulsion_factor
	 * @param rand
	 */
	static void addRepulsion(float[] position_1, float[] position_2, double[] movement, int repulsion_exponent,
			float repulsion_factor, java.util.Random rand, boolean is_2d) {
        double distx = position_2[0] - position_1[0];
        double disty = position_2[1] - position_1[1];
        double distz = position_2[2] - position_1[2];
        double totaldist = java.lang.Math.sqrt((distx * distx) + (disty * disty) + (distz * distz));
        
        if (totaldist == 0) {
            // if two points are at exactly the same position I need to add some random effect
            movement[0]=repulsion_factor*(rand.nextDouble()-0.5)*0.001;
            movement[1]=repulsion_factor*(rand.nextDouble()-0.5)*0.001;
            
            if (!is_2d) {
                movement[2]=repulsion_factor*(rand.nextDouble()-0.5)*0.001;
            }
            return;
        }
        
        // scale force**repvalpow
        double totalmove = repulsion_factor / Math.pow(totaldist, repulsion_exponent);

        movement[0] = (-distx / totaldist) * totalmove;
        movement[1] = (-disty / totaldist) * totalmove;
        movement[2] = (-distz / totaldist) * totalmove;
    }

	/**
	 * Use for filtering attraction values [-1<0<1]; -1 and 1 are max repulse/attract.
	 * 
	 * @param attractions
	 * @param minpval
	 * @return
	 */
	public static MinimalAttractionValue[] filterAttractionValues(MinimalAttractionValue[] attractions, double minpval){
        // use for filtering attraction values [-1<0<1]; -1 and 1 are max repulse/attract
        java.util.ArrayList <MinimalAttractionValue>retvec=new java.util.ArrayList<MinimalAttractionValue>();
        for(int i = 0; i < attractions.length; i++){
            if((attractions[i].att >= minpval) || (attractions[i].att <= -minpval)){
                retvec.add(attractions[i]);
            }
        }
        MinimalAttractionValue[] retarr=(MinimalAttractionValue[])retvec.toArray(new MinimalAttractionValue[0]);
        return retarr;
    }

	/**
	 * Computes simple attraction values. This actually takes all data from vector and makes ONE number out of it just
	 * use the BEST value!
	 * 
	 * @param invec
	 * @param dbsize
	 * @param minpval
	 * @param data
	 * @return
	 */
    public static float computeSimpleAttractionValue(double[] invec,int dbsize,double minpval,ClusterData data){

		if (invec == null) {
			return 0;
		} else if (invec.length < 1) { // if I have no hits
			return 0;
		}
       
		double bestval = invec[0]; // is a p-value (should be from 0 to 1)
        double currval;
        if(data.usescval){
            bestval=0;
            for(int i=invec.length-1;i>=0;i--){
                currval=invec[i];
                if(currval>bestval){
                    bestval=currval;
                }
            }
            if(bestval<data.maxvalfound){//maxvalfound=worst accepted value (comes from P-values where larger=worse)
                data.maxvalfound=bestval;
            }
            currval=bestval;
            if(currval<minpval){//minpval also functions as minscoreval here
                //System.out.println(" currval="+currval+" is less than minpval="+minpval+" returning zero");
                return 0;
            }
        }else{
            for(int i=invec.length-1;i>=1;i--){
                currval=invec[i];
                if(currval<bestval){
                    bestval=currval;
                }
            }
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
    }

	/**
	 * Returns a complex attraction value. This actually takes all data from vector and makes ONE number out of it.
	 * Multiply the pvalues of different HSPs.
	 * 
	 * @param invec
	 * @param dbsize
	 * @param minpval
	 * @param data
	 * @return
	 */
	public static float computeComplexAttractionValue(double[] invec, int dbsize, double minpval, ClusterData data) {

		if (invec == null) {
			return 0;
		} else if (invec.length < 1) {
			return 0;
		}

        double currval=invec[0];
        if(data.usescval){
			// then I am using score values (no logarithming at the end)
            currval=0;
            
			for (int i = invec.length; --i >= 0;) { // sum the values
                currval+=invec[i];
            }
			
			// maxvalfound=worst accepted value (comes from P-values where larger=worse)
			if (currval < data.maxvalfound) {
				data.maxvalfound = currval;
			}
			
			if (currval < minpval) { // minpval also functions as minscoreval here
                return 0;
            }
		} else { // then I am using P-values
            for(int i=invec.length;--i>=1;){
                currval*=invec[i];
            }
            if(currval>data.maxvalfound){
                data.maxvalfound=currval;
            }
            if(currval==0){
				return -1; // this is identity
			} else if (currval > 1) { // should never happen to p-values!
                return 0;
			} else if (currval > minpval) { // if this value is worse than permitted
                return 0;
            }
			// now all pvalues between 0 and 1
			currval = (-1 * java.lang.Math.log(currval)); // ln10;
        }
        return (float)(currval);
    }

	/**
	 * In-place removes all gaps from the input sequences
	 * 
	 * @param sequences
	 *            The input sequences
	 * @return The sequences (even though they're changed in-place)
	 */
	public static AminoAcidSequence[] removeGapsFromSequences(AminoAcidSequence[] sequences) {
		for (int i = 0; i < sequences.length; i++) {
			sequences[i].seq = sequences[i].seq.replaceAll("-", "");
		}
		return sequences;
	}
}