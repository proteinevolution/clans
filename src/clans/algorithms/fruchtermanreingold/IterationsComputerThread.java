package clans.algorithms.fruchtermanreingold;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import clans.gui.ProgramWindow;
import clans.io.FileHandling;
import clans.model.ClusterData;

/**
 * Class that handles the actual iteration computation.
 */
public class IterationsComputerThread extends java.lang.Thread {

	/**
	 * Creates a new thread instance for computing iterations of the algorithm with thead name "computation thread".
	 * <p>
	 * Note: Named threads are easier to debug! On a command line run command "jps -l" to get the programs process
	 * id (first column). Then run "jstack <process id>" to see all threads of your process.
	 * 
	 * @param parent
	 *            The GUI to which this thread belongs.
	 */
	public IterationsComputerThread(ProgramWindow parent, ClusterData data) {
		this.parent = parent;
		this.data = data;
		
		setName("IterationsComputerThread");
	}

	String tmpstr = "";
	float tmpcool = 1;
	final Object syncon = new Object(); // dummy object to sync on
	
	ProgramWindow parent;
	ClusterData data;

	/**
	 * Starts the computation by either computing itself (if {@code data.cpu}==1) or spawning {@code data.cpu} many
	 * child threads that do the work.
	 */
	@Override
	public void run() {

		data.roundsCompleted = 0;
		
		while (!Thread.currentThread().isInterrupted()) {
			data.rounds++;
			if (data.hasRoundsLimit()) {
				data.roundsCompleted++;

				if (data.roundsCompleted < data.getRoundsLimit()) {
					parent.updateStartStopResumeButtonLabel();

				} else {
					Thread.currentThread().interrupt();
					
					String completion_message = "completed the planned " + data.getRoundsLimit() + " rounds";
					if (parent.messageOverlayActive) {
						parent.messageOverlay.setCustomMessage(completion_message, null,
								parent.messageOverlay.getColorSuccess(), parent.messageOverlay.getDurationInfo(), true, false);
					} else {
						System.out.println(completion_message);
					}

					break;
				}
			}

			// first see whether the main window has the focus
			if (parent.draw_area.isFocusOwner()) {
				// mainwindow has the focus, then see whether I am done drawing
				try {
					while (parent.recalc == false) {
						Thread.sleep(100);
					}
				} catch (InterruptedException ie) {
					System.err.println("Interrupted sleep in computethread");
				}
			}

			if (data.changedvals) {
				parent.updateOptionValuesFromOptionsWindow();
				data.changedvals = false;
			}

			boolean optimize_only_selected_sequences;
			synchronized (parent.selectedSequencesLock) {
				/**
				 * Keep this from changing by user GUI interaction while we copy the selected names.
				 */
				optimize_only_selected_sequences = parent.optimizeOnlySelectedSequences();

				if (optimize_only_selected_sequences) {

					if ((data.selectedSequencesIndices.length == 0) || (data.selectedSequencesIndices.length == data.elements)) {
						// nothing or everything selected -> no special treatment necessary
						optimize_only_selected_sequences = false;

					} else {
						/**
						 * Create a copy of the selected entries that will remain stable during this round of
						 * calculations.
						 */
						data.selectedSequencesIndicesStableCopy = new int[data.selectedSequencesIndices.length];
						System.arraycopy(data.selectedSequencesIndices, 0, data.selectedSequencesIndicesStableCopy, 0, data.selectedSequencesIndices.length);
					}
				}
			}

			ClusterMethods.iterateOneRound(data, optimize_only_selected_sequences);

			if (parent.level == 0 && data.save_intermediate_results) {
				try {
					PrintWriter outwriter = new PrintWriter(new BufferedWriter(new FileWriter(
							data.getIntermediateResultfileName())));

					outwriter.println("sequences: " + data.positions.length);
					outwriter.println("values: ");
					outwriter.println("minattract " + data.minattract);
					outwriter.println("maxmove " + data.maxmove);
					outwriter.println("cooling " + data.cooling);
					outwriter.println("currcool " + data.currcool);
					outwriter.println("mineval " + data.mineval);
					outwriter.println("minpval " + data.pvalue_threshold);
					outwriter.println("repfactor " + data.repfactor);
					outwriter.println("attfactor " + data.attfactor);
					outwriter.println("dampening " + data.dampening);
					outwriter.println("rounds " + data.rounds);

					FileHandling.save_semicolon_delimited_positions(outwriter, data.positions);

					outwriter.flush();
					outwriter.close();

				} catch (IOException e) {
					System.err.println("unable to save positions to " + data.getIntermediateResultfileName());
					e.printStackTrace();
				}
			}

			data.posarr = data.positions;
			tmpcool = (((float) ((int) (data.currcool * 100000))) / 100000);
			if (parent.options_window != null) {
				parent.options_window.currcoolfield.setText(String.valueOf(tmpcool));
			}
			if (tmpcool <= 1e-5) {
				Thread.currentThread().interrupt();
				parent.modifyButtonStartStopResume("DONE (absolute zero)", true);
			}
			if (data.rounds % parent.skiprounds == 0 && data.nographics == false) {
				parent.recalc = false;
				parent.repaint();
			}
		}

		parent.repaint();
		
		synchronized (parent.computationLock) {
			parent.computationLock.notify();
		}
	}
}