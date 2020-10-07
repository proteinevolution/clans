package clans.gui;

import java.text.ParseException;
import java.util.*;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;

import javax.swing.*;
import javax.swing.Timer;

import clans.UnusedWeirdClassPreviouslyCalledSelectclass;
import clans.algorithms.ClusterDetection;
import clans.algorithms.ClusterDetectionBootstrapping;
import clans.algorithms.fruchtermanreingold.ClusterMethods;
import clans.algorithms.fruchtermanreingold.IterationsComputerThread;
import clans.io.BlastHitHandling;
import clans.io.ClusterDataLoadHelper;
import clans.io.FileHandling;
import clans.io.FileHandling2;
import clans.misc.MyMath;
import clans.misc.ReporterThread;
import clans.model.ClusterData;
import clans.model.SelectedSubsetHandling;
import clans.model.SequenceCluster;
import clans.model.SequenceGroup;
import clans.model.microarray.Replicates;
import clans.model.proteins.AminoAcidSequence;
import clans.model.proteins.HighScoringSegmentPair;
import clans.model.proteins.MinimalAttractionValue;
import clans.model.proteins.MinimalHsp;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.*;

/**
 * Creates and manages the CLANS graphical user interface (GUI).
 * 
 * TODO: Currently this class also contains computational code that should eventually be refactored to its own classes.
 */
public class ProgramWindow extends javax.swing.JFrame {

	private static final long serialVersionUID = -1310615823184297259L;

	/**
	 * The default labels for the start/stop/resume button so that they are consistent.
	 */
	enum ButtonStartStopMessage {
		START("Start run"), RESUME("Resume"), STOP("Stop");

		private String value;

		private ButtonStartStopMessage(String value) {
			this.value = value;
		}

		protected String get() {
			return value;
		}
	}

	/**
	 * Invokes an autosave every {@code data.getAutosaveIntervalMinutes())} minutes, if autosaving is enabled.
	 */
	private Timer autosaveTimer; // if autosaving is enabled this timer executes the autosave after the set interval
	/**
	 * If true, postpones the autosave setup to after the first successful save operation. This is necessary as
	 * sometimes, no filename is associated with {@data}, e.g. when coming directly from parsing BLAST result.
	 */
	private boolean autosaveAutoReenable;

	/**
	 * Runs the load/save operations in their own background thread while keeping the GUI responsive.
	 */
	private SwingWorker<Void, Integer> saveLoadWorker;

	/**
	 * If true, disables almost all mouse events as their listeners check this variable. This is used in @{code
	 * disableUserControls()} and {@code enableUserControls()} to block and reenable user input during, e.g., loading
	 * and saving.
	 */
	private boolean mouseEventsDisabled;
	/**
	 * Listens for escape key presses during long running operations and cancels them if the key is pressed. Currently
	 * used during loading and saving.
	 */
	private KeyListener cancelWorkInProgress;

	/**
	 * Creates the GUI.
	 * 
	 * @param data
	 *            The input data. Can be "empty" (i.e. without real data) to start the GUI without immediately loading a
	 *            file.
	 */
	public ProgramWindow(ClusterData data) {

		this.data = data;

		// closing events are handled in a WindowListener with custom windowClosing method
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);

		draw_area = new DrawArea(this, data);

		initializeComputationThread();

		initializeGuiComponents();

		consumeDataDuringConstruction();

		textfield_threshold_value.setText(String.valueOf(data.pvalue_threshold));
		textfield_info_min_blast_evalue.setText(String.valueOf(data.pvalue_threshold));

		// now initialize the stuff
		mousestart[0] = 0;
		mousestart[1] = 0;
		mousemove[0] = 0;
		mousemove[1] = 0;
		draw_area.init();

		// and now initialize a run
		resetGraph();

		updateSelectionButtonLabel();

		if (data.errbuff.length() > 0) {
			// If I have had errors up to this point
			new WindowErrorMessages(this, true, data.errbuff.toString()).setVisible(true);
		}

		if (data.getAbsoluteInputfileName() != null) {
			loadDataClansFormatInBackground(data.getAbsoluteInputfileName());
		}

		if (new File(data.getIntermediateResultfileName()).canRead()) {
			System.out.println("reading former data");
			FileHandling2.parse_intermediate_results(data);
			data.posarr = data.positions;
		}
	}

	/**
	 * Creates a new instance of the computation thread and names it "computation thread".
	 * <p>
	 * Note: Named threads are easier to debug! On a command line run command "jps -l" to get the programs process id
	 * (first column). Then run "jstack <process id>" to see all threads of your process.
	 */
	private void initializeComputationThread() {
		reporterThread = new ReporterThread(this, new IterationsComputerThread(this, data),
				reporterWaitForIteratorThreadStartLock, computationLock, "computationHasStopped", "ReporterThread");
	}

	/**
	 * Initializes the GUI elements. This is called from the constructor.
	 */
	private void initializeGuiComponents() {
		// graphpanelxx = new javax.swing.JPanel();
		drawbuttonpanel = new javax.swing.JPanel();

		initializeButton = new javax.swing.JButton();
		startStopResumeButton = new javax.swing.JButton();
		roundsCompletedTextfield = new javax.swing.JTextField();
		
		showSelectedButton = new javax.swing.JButton();
		selectionOrMovementModeButton = new javax.swing.JToggleButton();
		setThresholdButton = new javax.swing.JButton();
		textfield_threshold_value = new javax.swing.JTextField();
		textfield_info_min_blast_evalue = new javax.swing.JTextField();
		selectionButton = new javax.swing.JButton();
		showSequenceNamesCheckbox = new javax.swing.JCheckBox();
		showSequenceNumbersCheckbox = new javax.swing.JCheckBox();
		showConnectionsCheckbox = new javax.swing.JCheckBox();
		zoomOnSelectedButton = new javax.swing.JButton();
		resetZoomButton = new javax.swing.JButton();
		
		jMenuBar1 = new javax.swing.JMenuBar();
		menu_file = new javax.swing.JMenu();
		loadmenuitem = new javax.swing.JMenuItem();
		savemenuitem = new javax.swing.JMenuItem();
		saveasmenuitem = new javax.swing.JMenuItem();
		saveattvalsmenuitem = new javax.swing.JMenuItem();
		addseqsmenuitem = new javax.swing.JMenuItem();
		savemtxmenuitem = new javax.swing.JMenuItem();
		save2dmenuitem = new javax.swing.JMenuItem();
		printmenuitem = new javax.swing.JMenuItem();
		loadalternatemenuitem = new javax.swing.JMenuItem();
		loadtabsmenuitem = new javax.swing.JMenuItem();
		loadgroupsmenuitem = new javax.swing.JMenuItem();
		menu_misc = new javax.swing.JMenu();
		getseqsmenuitem = new javax.swing.JMenuItem();
		hidesingletonsmenuitem = new javax.swing.JMenuItem();
		remove_selected_sequences_menu_item = new javax.swing.JMenuItem();
		getchildmenuitem = new javax.swing.JMenuItem();
		getparentmenuitem = new javax.swing.JMenuItem();
		setrotmenuitem = new javax.swing.JMenuItem();
		attvalcompcheckbox = new javax.swing.JCheckBoxMenuItem();
		moveOnlySelectedSequencesCheckbox = new javax.swing.JCheckBoxMenuItem();
		cluster2dbutton = new javax.swing.JCheckBoxMenuItem();
		rescalepvaluescheckbox = new javax.swing.JCheckBoxMenuItem();
		messageOverlayActiveCheckbox = new javax.swing.JCheckBoxMenuItem();
		autosaveSetup = new javax.swing.JMenuItem();
		skipdrawingrounds = new javax.swing.JMenuItem();
		menu_draw = new javax.swing.JMenu();
		changefontmenuitem = new javax.swing.JMenuItem();
		getdotsizemenuitem = new javax.swing.JMenuItem();
		getovalsizemenuitem = new javax.swing.JMenuItem();
		changecolormenuitem = new javax.swing.JMenuItem();
		changefgcolormenuitem = new javax.swing.JMenuItem();
		changebgcolormenuitem = new javax.swing.JMenuItem();
		changeselectcolormenuitem = new javax.swing.JMenuItem();
		changenumbercolor = new javax.swing.JMenuItem();
		changeblastcolor = new javax.swing.JMenuItem();
		colorDotsBySequenceLengthCheckbox = new javax.swing.JCheckBoxMenuItem();
		showConnectionsFrustrationColoredCheckbox = new javax.swing.JCheckBoxMenuItem();
		showorigcheckbox = new javax.swing.JCheckBoxMenuItem();
		showInfoOnDrawAreaCheckbox = new javax.swing.JCheckBoxMenuItem();
		shownamesselectcheckbox = new javax.swing.JCheckBoxMenuItem();
		showBlastHitNumberCheckbox = new javax.swing.JCheckBoxMenuItem();
		zoommenuitem = new javax.swing.JMenuItem();
		centermenuitem = new javax.swing.JMenuItem();
		antialiasingcheckboxmenuitem = new javax.swing.JCheckBoxMenuItem();
		stereocheckboxmenuitem = new javax.swing.JCheckBoxMenuItem();
		stereoanglemenuitem = new javax.swing.JMenuItem();
		menu_windows = new javax.swing.JMenu();
		showoptionsmenuitem = new javax.swing.JMenuItem();
		sequencesitem = new javax.swing.JMenuItem();
		evalueitem = new javax.swing.JMenuItem();
		getblasthitsmenuitem = new javax.swing.JMenuItem();
		clustermenuitem = new javax.swing.JMenuItem();
		getseqsforselectedhits = new javax.swing.JMenuItem();
		seqscoloring = new javax.swing.JMenuItem();
		showseqsmenuitem = new javax.swing.JMenuItem();
		rotationmenuitem = new javax.swing.JMenuItem();
		affymenuitem = new javax.swing.JMenuItem();
		mapmanmenuitem = new javax.swing.JMenuItem();
		taxonomymenuitem = new javax.swing.JMenuItem();
		menu_help = new javax.swing.JMenu();
		aboutmenuitem = new javax.swing.JMenuItem();
		helpmenuitem = new javax.swing.JMenuItem();

		setTitle("3D-View");
		addWindowListener(new java.awt.event.WindowAdapter() {
			public void windowClosing(java.awt.event.WindowEvent evt) {
				openReallyExitSafetyDialog();
			}
		});

		final Dimension initial_size = new Dimension(640, 480);
		draw_area.setPreferredSize(initial_size);
		draw_area.setSize(initial_size);

		draw_area.addMouseMotionListener(new java.awt.event.MouseMotionAdapter() {
			public void mouseDragged(java.awt.event.MouseEvent evt) {
				graphpanelMouseDragged(evt);
			}
		});
		
		draw_area.addMouseListener(new java.awt.event.MouseAdapter() {
			public void mousePressed(java.awt.event.MouseEvent evt) {
				graphpanelMousePressed(evt);
			}

			public void mouseReleased(java.awt.event.MouseEvent evt) {
				graphpanelMouseReleased(evt);
			}
		});

		draw_area.addMouseWheelListener(new java.awt.event.MouseWheelListener() {
			public void mouseWheelMoved(java.awt.event.MouseWheelEvent evt) {
				zoomUsingMouseWheel(evt);
			}
		});
		
		draw_area.setLayout(new javax.swing.BoxLayout(draw_area, javax.swing.BoxLayout.LINE_AXIS));
		getContentPane().add(draw_area, java.awt.BorderLayout.CENTER);

		initializeButton.setText("Initialize");
		initializeButton.setToolTipText("Initialize the graph positions");
		initializeButton.setMnemonic(KeyEvent.VK_I);
		initializeButton.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				initializeGraphPositions();
			}
		});

		startStopResumeButton.setToolTipText("start/resume/stop the current run");
		startStopResumeButton.setMnemonic(KeyEvent.VK_S);
		startStopResumeButton.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				toggleComputationRunning();
			}
		});
		updateStartStopResumeButtonLabel();

		showSelectedButton.setText("Show selected");
		showSelectedButton.setMnemonic(KeyEvent.VK_O);
		showSelectedButton.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openShowSelectedSequencesWindow();
			}
		});

		selectionOrMovementModeButton.setToolTipText("Toggles between moving the coordinate system and selecting sequences");
		selectionOrMovementModeButton.setMnemonic(KeyEvent.VK_V);
		selectionOrMovementModeButton.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				updateSelectionOrMovementModeButtonLabel();
			}
		});
		updateSelectionOrMovementModeButtonLabel();

		roundsCompletedTextfield.setEditable(false);
		roundsCompletedTextfield.setText("0");
		roundsCompletedTextfield.setFont(getContentPane().getFont().deriveFont(Font.BOLD));
		roundsCompletedTextfield.setToolTipText("rounds completed since last initialization");
		roundsCompletedTextfield.setHorizontalAlignment(JTextField.CENTER);

		setThresholdButton.setText("Use p-values better than:");
		setThresholdButton.setMnemonic(KeyEvent.VK_B);
		setThresholdButton.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				setThresholdButtonPressed();
			}
		});

		textfield_threshold_value.setText("1");
		textfield_threshold_value.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				confirmedThresholdTextfieldValue();
			}
		});

		textfield_info_min_blast_evalue.setEditable(false);
		textfield_info_min_blast_evalue.setText("1");
		textfield_info_min_blast_evalue.setToolTipText("evalue limit used for blast");

		updateSelectionButtonLabel();
		selectionButton.setMnemonic(KeyEvent.VK_A);
		selectionButton.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				selectionButtonActivated();
			}
		});

		showSequenceNamesCheckbox.setText("show names");
		showSequenceNamesCheckbox.setToolTipText("show sequence names in graph");
		showSequenceNamesCheckbox.setMnemonic(KeyEvent.VK_N);
		showSequenceNamesCheckbox.addItemListener(new java.awt.event.ItemListener() {
			public void itemStateChanged(java.awt.event.ItemEvent evt) {
				requestRepaint();
			}
		});

		showSequenceNumbersCheckbox.setText("show numbers");
		showSequenceNumbersCheckbox.setToolTipText("show sequence numbers in graph");
		showSequenceNumbersCheckbox.setMnemonic(KeyEvent.VK_U);
		showSequenceNumbersCheckbox.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				requestRepaint();
			}
		});
		
		showConnectionsCheckbox.setText("show connections");
		showConnectionsCheckbox.setToolTipText("draw lines for all connections better than the selected cutoff");
		showConnectionsCheckbox.setMnemonic(KeyEvent.VK_T);
		showConnectionsCheckbox.addItemListener(new java.awt.event.ItemListener() {
			public void itemStateChanged(java.awt.event.ItemEvent evt) {
				requestRepaint();
			}
		});
		
		zoomOnSelectedButton.setMnemonic(KeyEvent.VK_Z);
		zoomOnSelectedButton.setToolTipText("Zoom to the selected sequences [>= 2 selected sequences].");
		zoomOnSelectedButton.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				zoomOnSelectedTriggered(true);
			}
		});
		
		resetZoomButton.setText("Reset Zoom");
		resetZoomButton.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				resetZoom();
			}
		});

		informAboutZoomChange();
		
		// we manage the buttons here so it becomes easier to change their layout
		drawbuttonpanel.setLayout(new java.awt.GridLayout(3, 5));

		// row 1
		drawbuttonpanel.add(initializeButton);
		drawbuttonpanel.add(startStopResumeButton);
		drawbuttonpanel.add(roundsCompletedTextfield);
		drawbuttonpanel.add(showSelectedButton);
		drawbuttonpanel.add(selectionOrMovementModeButton);
				
		// row 2
		drawbuttonpanel.add(setThresholdButton);
		drawbuttonpanel.add(textfield_threshold_value);
		drawbuttonpanel.add(textfield_info_min_blast_evalue);
		drawbuttonpanel.add(selectionButton);
		drawbuttonpanel.add(new JLabel("")); // empty cell to shift the following components
		
		// row 3
		drawbuttonpanel.add(showSequenceNamesCheckbox);
		drawbuttonpanel.add(showSequenceNumbersCheckbox);
		drawbuttonpanel.add(showConnectionsCheckbox);
		drawbuttonpanel.add(zoomOnSelectedButton);
		drawbuttonpanel.add(resetZoomButton);
		
		getContentPane().add(drawbuttonpanel, java.awt.BorderLayout.SOUTH);

		menu_file.setText("File");
		menu_file.setMnemonic(KeyEvent.VK_F);

		loadmenuitem.setText("Load Run");
		loadmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				loadDataClansFormat();
			}
		});
		menu_file.add(loadmenuitem);

		savemenuitem.setText("Save Run");
		savemenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				saveRun();
			}
		});
		menu_file.add(savemenuitem);

		saveasmenuitem.setText("Save Run As");
		saveasmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				saveRunAs();
			}
		});
		menu_file.add(saveasmenuitem);

		saveattvalsmenuitem.setText("Save attraction values to file");
		saveattvalsmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				saveAttractionValuesAs();
			}
		});
		menu_file.add(saveattvalsmenuitem);

		addseqsmenuitem.setText("Add Sequences");
		addseqsmenuitem.setToolTipText("This feature is currently limited to running CLANS in the command line");
		addseqsmenuitem.setEnabled(false);
		menu_file.add(addseqsmenuitem);

		savemtxmenuitem.setText("Save blast matrix p-values");
		savemtxmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				saveBlastMatrixPValuesAs();
			}
		});
		menu_file.add(savemtxmenuitem);

		save2dmenuitem.setText("Save 2d graph data");
		save2dmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				save2dGraphDataAs();
			}
		});
		menu_file.add(save2dmenuitem);

		printmenuitem.setText("Print view");
		printmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				printDrawArea();
			}
		});
		menu_file.add(printmenuitem);

		loadalternatemenuitem.setText("Load data in matrix format");
		loadalternatemenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				loadDataMatrixFormat();
			}
		});
		menu_file.add(loadalternatemenuitem);

		loadtabsmenuitem.setText("Load tabular data");
		loadtabsmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				loadDataTabularFormat();
			}
		});
		menu_file.add(loadtabsmenuitem);

		loadgroupsmenuitem.setText("Append sequence groups from file");
		loadgroupsmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				loadSequenceGroupsToAppend();
			}
		});
		menu_file.add(loadgroupsmenuitem);

		jMenuBar1.add(menu_file);

		menu_misc.setText("Misc");
		menu_misc.setMnemonic(KeyEvent.VK_M);

		getseqsmenuitem.setText("Extract selected sequences");
		getseqsmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openSaveSequenceAsFastaDialog();
			}
		});
		menu_misc.add(getseqsmenuitem);

		remove_selected_sequences_menu_item.setText("Remove selected sequences");
		remove_selected_sequences_menu_item.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				removeAllSelectedSequences();
			}
		});
		menu_misc.add(remove_selected_sequences_menu_item);

		hidesingletonsmenuitem.setText("Hide singletons");
		hidesingletonsmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				hideSingletonsMenuItemActionPerformed();
			}
		});
		menu_misc.add(hidesingletonsmenuitem);

		getchildmenuitem.setText("Use selected subset");
		getchildmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				useGraphSubset();
			}
		});
		menu_misc.add(getchildmenuitem);

		getparentmenuitem.setText("Use parent group (0)");
		getparentmenuitem.setEnabled(false);
		getparentmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				useGraphSuperset();
			}
		});
		menu_misc.add(getparentmenuitem);

		setrotmenuitem.setText("Set rotation values");
		setrotmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openRotationMatrixInputDialog();
			}
		});
		menu_misc.add(setrotmenuitem);

		attvalcompcheckbox.setSelected(true);
		attvalcompcheckbox.setText("Complex attraction");
		attvalcompcheckbox.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				updateUseComplexAttractionValuesState();
			}
		});
		menu_misc.add(attvalcompcheckbox);

		moveOnlySelectedSequencesCheckbox.setText("Optimize only selected sequences");
		moveOnlySelectedSequencesCheckbox.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				updateOptimizeOnlySelectedSequencesState();
			}
		});
		menu_misc.add(moveOnlySelectedSequencesCheckbox);

		cluster2dbutton.setText("Cluster in 2D");
		cluster2dbutton.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				toggle2d3dGraph();
			}
		});
		menu_misc.add(cluster2dbutton);

		rescalepvaluescheckbox.setText("Rescale attraction values");
		rescalepvaluescheckbox.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				updateRescaleAttractionValuesState();
			}
		});
		menu_misc.add(rescalepvaluescheckbox);

		autosaveSetup.setText("setup autosaving");
		autosaveSetup.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openAutosaveSetupDialog();
			}
		});
		menu_misc.add(autosaveSetup);

		messageOverlayActiveCheckbox.setText("activate message overlay");
		messageOverlayActive = true;
		messageOverlayActiveCheckbox.setSelected(true);
		messageOverlayActiveCheckbox.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				toggleMessageOverlayActive();
			}
		});
		menu_misc.add(messageOverlayActiveCheckbox);

		skipdrawingrounds.setText("Only draw every Nth round (speedup)");
		skipdrawingrounds.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				skipdrawingroundsActionPerformed();
			}
		});
		menu_misc.add(skipdrawingrounds);

		jMenuBar1.add(menu_misc);

		menu_draw.setText("Draw");
		menu_draw.setMnemonic(KeyEvent.VK_D);

		changefontmenuitem.setText("Change Font");
		changefontmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openChangeFontDialog();
			}
		});
		menu_draw.add(changefontmenuitem);

		getdotsizemenuitem.setText("Set dot size");
		getdotsizemenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openChaneDotSizeDialog();
			}
		});
		menu_draw.add(getdotsizemenuitem);

		getovalsizemenuitem.setText("Set selected circle size");
		getovalsizemenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openChangeSelectedCircleSizeDialog();
			}
		});
		menu_draw.add(getovalsizemenuitem);

		changecolormenuitem.setText("Change color (dot connections)");
		changecolormenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openChangeDotConnectionColorDialog();
			}
		});
		menu_draw.add(changecolormenuitem);

		changefgcolormenuitem.setText("Change color (Foreground)");
		changefgcolormenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openForegroundColorChangeDialog();
			}
		});
		menu_draw.add(changefgcolormenuitem);

		changebgcolormenuitem.setText("Change color (Background)");
		changebgcolormenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openBackgroundColorChangeDialog();
			}
		});
		menu_draw.add(changebgcolormenuitem);

		changeselectcolormenuitem.setText("Change color (Selecteds)");
		changeselectcolormenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openSelectedColorChangeDialog();
			}
		});
		menu_draw.add(changeselectcolormenuitem);

		changenumbercolor.setText("Change color (BLAST hit numbers)");
		changenumbercolor.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openNumberColorChangeDialog();
			}
		});
		menu_draw.add(changenumbercolor);

		changeblastcolor.setText("Change color (BLAST hit circles)");
		changeblastcolor.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openBlastColorChangeDialog();
			}
		});
		menu_draw.add(changeblastcolor);

		colorDotsBySequenceLengthCheckbox.setText("Color dots by sequence length (yellow=short, blue=long)");
		menu_draw.add(colorDotsBySequenceLengthCheckbox);

		showConnectionsFrustrationColoredCheckbox.setText("Color by edge \"frustration\" (red=too long, blue=too short)");
		showConnectionsFrustrationColoredCheckbox.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				requestRepaint();
			}
		});
		menu_draw.add(showConnectionsFrustrationColoredCheckbox);

		showorigcheckbox.setText("Show origin");
		menu_draw.add(showorigcheckbox);

		if (data == null) {
			showInfoOnDrawAreaCheckbox.setSelected(false);
		} else {
			showInfoOnDrawAreaCheckbox.setSelected(data.showinfo);
		}
		showInfoOnDrawAreaCheckbox.setText("Show info");
		showInfoOnDrawAreaCheckbox.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				requestRepaint();
			}
		});
		menu_draw.add(showInfoOnDrawAreaCheckbox);

		shownamesselectcheckbox.setText("Show names while selecting");
		menu_draw.add(shownamesselectcheckbox);

		showBlastHitNumberCheckbox.setText("Show hsp sequence numbers");
		menu_draw.add(showBlastHitNumberCheckbox);

		zoommenuitem.setText("Zoom");
		zoommenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openChangeZoomFactorDialog();
			}
		});
		menu_draw.add(zoommenuitem);

		centermenuitem.setText("Center graph");
		centermenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				centerGraph();
			}
		});
		menu_draw.add(centermenuitem);

		antialiasingcheckboxmenuitem.setText("Antialiasing (slow !)");
		antialiasingcheckboxmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				requestRepaint();
			}
		});
		menu_draw.add(antialiasingcheckboxmenuitem);

		stereocheckboxmenuitem.setText("Stereo");
		stereocheckboxmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				requestRepaint();
			}
		});
		menu_draw.add(stereocheckboxmenuitem);

		stereoanglemenuitem.setText("Change stereo angle (0-360)");
		stereoanglemenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openChangeStereoVisionAngleDialog();
			}
		});
		menu_draw.add(stereoanglemenuitem);

		jMenuBar1.add(menu_draw);

		menu_windows.setText("Windows");
		menu_windows.setMnemonic(KeyEvent.VK_W);

		showoptionsmenuitem.setText("Show options window");
		showoptionsmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openOptionsWindow();
			}
		});
		menu_windows.add(showoptionsmenuitem);

		sequencesitem.setText("Selecteds");
		sequencesitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openShowSelectedSequencesWindow();
			}
		});
		menu_windows.add(sequencesitem);

		evalueitem.setText("P-value plot");
		evalueitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openScoreDistributionPlot();
			}
		});
		menu_windows.add(evalueitem);

		getblasthitsmenuitem.setText("Show blast hits for sequence:");
		getblasthitsmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openShowBlastHitsForSequenceDialog();
			}
		});
		menu_windows.add(getblasthitsmenuitem);

		clustermenuitem.setText("Find clusters");
		clustermenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openFindClustersWindow();
			}
		});
		menu_windows.add(clustermenuitem);

		getseqsforselectedhits.setText("Get sequence with hits from/to selected");
		getseqsforselectedhits.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openGetSequencesConnectedToSelectedsWindow();
			}
		});
		menu_windows.add(getseqsforselectedhits);

		seqscoloring.setText("Edit Groups");
		seqscoloring.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openEditGroupsWindow();
			}
		});
		menu_windows.add(seqscoloring);

		showseqsmenuitem.setText("Show selected sequences as text (copy/pastable)");
		showseqsmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openShowCopyPasteableSequencesWindow();
			}
		});
		menu_windows.add(showseqsmenuitem);

		rotationmenuitem.setText("Rotation");
		rotationmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openRotationWindow();
			}
		});
		menu_windows.add(rotationmenuitem);

		affymenuitem.setText("Microarray_data");
		affymenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openMicroarrayWindow();
			}
		});
		menu_windows.add(affymenuitem);

		mapmanmenuitem.setText("Functional mapping");
		mapmanmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openFunctionalMappingWindow();
			}
		});
		menu_windows.add(mapmanmenuitem);

		taxonomymenuitem.setText("Taxonomy");
		taxonomymenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openTaxonomyWindow();
			}
		});
		menu_windows.add(taxonomymenuitem);

		jMenuBar1.add(menu_windows);

		menu_help.setText("Help");
		menu_help.setMnemonic(KeyEvent.VK_H);

		aboutmenuitem.setText("About");
		aboutmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openAboutWindow();
			}
		});
		menu_help.add(aboutmenuitem);

		helpmenuitem.setText("Help");
		helpmenuitem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				openHelpWindow();
			}
		});
		menu_help.add(helpmenuitem);

		jMenuBar1.add(menu_help);

		setJMenuBar(jMenuBar1);

		originalGlassPane = getGlassPane();
		activateMessageOverlay();
		addGuiMessageOverlayResizeListener();

		createLoadSaveCancelKeyListener();

		pack();
	}

	/**
	 * Adds listener to inform the message overlay of main window size changes.
	 */
	private void addGuiMessageOverlayResizeListener() {
		draw_area.addComponentListener(new ComponentAdapter() {
			@Override
			public void componentResized(ComponentEvent arg0) {
				if (messageOverlay != null) {
					messageOverlay.updateGlassSize(draw_area.getHeight());
				}
			}
		});
	}

	/**
	 * Creates the KeyListener that is added/removed when the escape key is activated/deactivated for canceling
	 * long-running operations. Currently used during file loading and saving.
	 */
	private void createLoadSaveCancelKeyListener() {
		cancelWorkInProgress = new KeyListener() {
			@Override
			public void keyTyped(KeyEvent e) {
			}

			@Override
			public void keyPressed(KeyEvent e) {
			}

			@Override
			public void keyReleased(KeyEvent e) {
				if (e.getKeyCode() == KeyEvent.VK_ESCAPE) {
					if (saveLoadWorker != null) {
						saveLoadWorker.cancel(true);
					}
				}
			}
		};
	}

	/**
	 * Consumes the data and sets some parts of the GUI up accordingly.
	 */
	private void consumeDataDuringConstruction() {
		data.nographics = false;

		data.mineval = data.eval;
		data.pvalue_threshold = data.pval;

		if (data.scval >= 0) { // in that case use a score cutoff
			data.usescval = true;

			data.pvalue_threshold = data.scval;
			setThresholdButton.setText("Use SC-vals better than");
			textfield_threshold_value.setText("0");
			evalueitem.setText("SC-value plot");
			attvalcompcheckbox.setSelected(false);

		} else {
			data.usescval = false;
		}

		data.seqlengths = new float[data.seqnum];
		float maxlength = 0;
		for (int i = 0; i < data.seqnum; i++) {
			data.seqlengths[i] = data.sequences[i].length();
			if (data.seqlengths[i] > maxlength) {
				maxlength = data.seqlengths[i];
			}
		}

		for (int i = 0; i < data.seqnum; i++) {
			data.seqlengths[i] /= maxlength;
		}
	}

	/**
	 * Opens file choice dialog and loads the selected tabular format file.
	 */
	private void loadDataTabularFormat() {
		boolean restart_computation = stopComputation(true);

		groupseqs = null;
		int returnVal = fc.showOpenDialog(this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			// update loadtabdata(fc.getSelectedFile().getAbsolutePath());
			if (myseqgroupwindow != null) {
				myseqgroupwindow.setVisible(false);
				myseqgroupwindow.dispose();
			}
			if (mymapfunctiondialog != null) {
				mymapfunctiondialog.setVisible(false);
				mymapfunctiondialog.dispose();
			}
			draw_area.repaint();
		}

		if (restart_computation) {
			startComputation();
		}
	}

	/**
	 * Opens a dialog that lets users enter a new stereo vision angle.
	 */
	private void openChangeStereoVisionAngleDialog() {
		String tmpstr = "";
		try {
			tmpstr = JOptionPane.showInputDialog(this, "Enter the new angle (int):",
					String.valueOf(draw_area.stereoangle));
			if (tmpstr != null) {
				draw_area.stereoangle = (int) (Float.parseFloat(tmpstr));
			}
		} catch (NumberFormatException ne) {
			javax.swing.JOptionPane.showMessageDialog(this, "ERROR, unable to parse integer from '" + tmpstr + "'");
		}
		draw_area.repaint();
	}

	/**
	 * Opens a windows that shows functional mapping data.
	 */
	private void openFunctionalMappingWindow() {
		if (mymapfunctiondialog != null) {
			mymapfunctiondialog.setVisible(false);
			mymapfunctiondialog.dispose();
		}
		mymapfunctiondialog = new WindowMicroarrayFunctionMapping(this);
		mymapfunctiondialog.setVisible(true);
	}

	/**
	 * Enables zoom with mousewheel+CTRL (coarse grained) or mousewheel+CTRL+SHIFT (fine grained)
	 * 
	 * @param evt
	 *            The mouse wheel event.
	 */
	private void zoomUsingMouseWheel(java.awt.event.MouseWheelEvent evt) {

		if (mouseEventsDisabled) {
			return;
		}

		if (!evt.isControlDown()) {
			return;
		}

		float zoom_factor = draw_area.getZoomFactor();

		if (evt.isShiftDown()) {
			zoom_factor += ((float) -evt.getWheelRotation()) / 100;
		} else {
			zoom_factor += ((float) -evt.getWheelRotation()) / 10;
		}
		
		if (zoom_factor < 0.1) {
			zoom_factor = 0.1f;
		}

		setZoom(zoom_factor);
	}

	/**
	 * Opens the Microarray Data window.
	 */
	private void openMicroarrayWindow() {
		if (myaffydialog != null) {
			myaffydialog.setVisible(false);
			myaffydialog.dispose();
		}
		myaffydialog = new WindowAffyMicroarrayData(this);
		myaffydialog.setVisible(true);
	}

	/**
	 * Opens a window that lets users specify rotations of the graph.
	 */
	private void openRotationWindow() {
		if (myrotationdialog != null) {
			myrotationdialog.setVisible(false);
			myrotationdialog.dispose();
		}
		
		myrotationdialog = new WindowRotatationAnimantionSetup(this);
		myrotationdialog.setVisible(true);
	}

	/**
	 * Loads additional sequence groups from a file and appends them to the current ones.
	 */
	private void loadSequenceGroupsToAppend() {
		
		int returnVal = fc.showOpenDialog(this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			data.append_groups_or_clusters_from_file(fc.getSelectedFile());
			if (myseqgroupwindow != null) {
				myseqgroupwindow.setVisible(false);
				myseqgroupwindow.dispose();
			}
			draw_area.repaint();
		}
	}

	/**
	 * Opens the options window.
	 */
	private void openOptionsWindow() {
		if (options_window != null) {
			options_window.setVisible(false);
			options_window.dispose();
		}
		options_window = new WindowOptions(this);
		options_window.setVisible(true);
	}

	/**
	 * Opens the "setup autosave" dialog and handles the user's interaction with it. The dialog lets users decide
	 * whether to save automatically and how often to do it.
	 */
	private void openAutosaveSetupDialog() {

		String user_input = (String) JOptionPane.showInputDialog(this,
				"autosave this file every N minutes (0 = disable autosave)", "setup autosaving for this file",
				JOptionPane.PLAIN_MESSAGE, null, null, data.getAutosaveIntervalMinutes());

		if (user_input == null) { // dialog was canceled -> no change
			return;
		}

		try {
			int new_interval = Integer.parseInt(user_input);

			if (new_interval < 0) { // time intervals cannot be negative
				throw new NumberFormatException();
			}

			data.setAutosaveIntervalMinutes(new_interval);
		} catch (NumberFormatException e) {
			javax.swing.JOptionPane.showMessageDialog(this, "unable to parse positive integer number from \""
					+ user_input + "\"");
			return;
		}

		if (data.getAutosaveIntervalMinutes() > 0) {
			enableAutosave(true);

		} else {
			disableAutosave(true);
			return;
		}
	}
	
	/**
	 * Enables autosaving if it is set up to occur. For loaded data without autosave parameters, a reasonably long
	 * autosave interval is automatically set up and the user is informed about it in a pop-up dialog.
	 * 
	 * @param show_dialog
	 *            If true, a dialogs informs about newly added autostart times for previously autostart unaware files.
	 *            This dialog is used only when the file is first loaded.
	 */
	private void initializeAutosave(boolean show_dialog) {
		if (!data.knowsAutosave()) {
			data.makeAutosaveAware();
			
			if (show_dialog) {
				javax.swing.JOptionPane.showMessageDialog(this, "Autosaving every " + data.getAutosaveIntervalMinutes()
						+ " minutes has been enabled for this file." + "\nTo change this go to Menu->Misc->"
						+ autosaveSetupLabel + ".", "autosaving enabled", JOptionPane.INFORMATION_MESSAGE);
			}
		}

		if (data.getAutosaveIntervalMinutes() > 0) {
			enableAutosave(false);
		}
		
	}

	private boolean isAutosaving() {
		return autosaveTimer != null && autosaveTimer.isRunning();
	}
	
	/**
	 * Enables autosaving if an autosave interval is set. Does nothing if the autosave interval is 0.
	 * 
	 * @param show_overlay_message
	 *            If true and if the message overlay is active, the user is informed about the change in an overlay
	 *            message.
	 */
	private void enableAutosave(boolean show_overlay_message) {

		if (data.getAutosaveIntervalMinutes() == 0) {
			return;
		}

		if (isAutosaving()) {
			autosaveTimer.stop();
		}

		// convert minutes to milliseconds for time
		autosaveTimer = new Timer(data.getAutosaveIntervalMinutes() * 60 * 1000, new ActionListener() {
			public void actionPerformed(ActionEvent ae) {
				autosave();
			}
		});
		autosaveTimer.start();

		if (messageOverlayActive && show_overlay_message) {
			messageOverlay.setCustomMessage("autosaving ENABLED", "every " + data.getAutosaveIntervalMinutes()
					+ " minutes", messageOverlay.getColorSuccess(), GuiMessageOverlay.Duration.INFO, true, false);
		}
	}

	/**
	 * Disables autosaving.
	 * 
	 * @param show_overlay_message
	 *            If true and if the message overlay is active, the user is informed about the change in an overlay
	 *            message.
	 */
	private void disableAutosave(boolean show_overlay_message) {

		data.setAutosaveIntervalMinutes(0);

		if (isAutosaving()) {
			autosaveTimer.stop();
		}
		autosaveTimer = null;
		

		if (messageOverlayActive && show_overlay_message) {
			messageOverlay.setCustomMessage("autosaving DISABLED", null, messageOverlay.getColorWorking(),
					GuiMessageOverlay.Duration.WARNING, true, false);
		}
	}
	
	/**
	 * Stops the autosave timer.
	 * 
	 * @return true if a timer exists and it was running, false else.
	 */
	private boolean pauseAutosave() {
		if (!isAutosaving()){
			return false;
		}
		
		autosaveTimer.stop();
		return true;
	}
	
	/**
	 * Restarts the autosave timer. This resets the time until the next autosave to the autosave interval.
	 */
	private void restartAutosave() {
		if (autosaveTimer != null) {
			autosaveTimer.restart();
		}
	}
	
	/**
	 * Toggles whether the message overlay is shown or disabled.
	 */
	private void toggleMessageOverlayActive() {
		
		if (messageOverlayActive) {
			deactivateMessageOverlay();
		} else {
			activateMessageOverlay();
		}
		
		messageOverlayActive = !messageOverlayActive;
	}
	
	/**
	 * Activates the message overlay.
	 * <p>
	 * TODO: remember: production versions MUST use GuiMessageOverlay instead of the debugging class
	 * GuiMessageOverlayLogged!
	 */
	private void activateMessageOverlay() {
		messageOverlay = new GuiMessageOverlay();
		messageOverlay.updateGlassSize(draw_area.getHeight());

		setGlassPane(messageOverlay);
	}
	
	/**
	 * Deactivates the message overlay.
	 */
	private void deactivateMessageOverlay() {
		setGlassPane(originalGlassPane);
		messageOverlay = null;
	}
	
	/**
	 * Changes the "enabled" state of all menus in the menubar.
	 * 
	 * @param enabled
	 *            If true all menus will be enabled, if false disabled.
	 */
	private void modifyMenusEnabled(boolean enabled) {
		for (int i=0; i < jMenuBar1.getMenuCount(); i++) {
			jMenuBar1.getMenu(i).setEnabled(enabled);
		}
	}

	/**
	 * Enables all menus in the menubar.
	 */
	private void enableMenus() {
		modifyMenusEnabled(true);
	}

	/**
	 * Disables all menus in the menubar.
	 */
	private void disableMenus() {
		modifyMenusEnabled(false);
	}

	/**
	 * Changes the "enabled" state of all buttons in the control bar (bottom of the CLANS window).
	 * 
	 * @param enabled
	 *            If true all buttons will be enabled, if false disabled.
	 */
	private void modifyControlButtonsEnabled(boolean enabled) {
		for (int i = 0; i < drawbuttonpanel.getComponentCount(); i++) {
			drawbuttonpanel.getComponent(i).setEnabled(enabled);
		}

		informAboutZoomChange();
	}

	/**
	 * Enables all buttons in the control bar (bottom of the CLANS window).
	 */
	private void enableControlButtons() {
		modifyControlButtonsEnabled(true);
	}

	/**
	 * Disables all buttons in the control bar (bottom of the CLANS window).
	 */
	private void disableControlButtons() {
		modifyControlButtonsEnabled(false);
	}

	/**
	 * Enables all mouse event handling in event listeners.
	 */
	private void enableMouseEventHandling() {
		mouseEventsDisabled = false;
	}
	
	/**
	 * Disables all mouse event handling in event listeners.
	 */
	private void disableMouseEventHandling() {
		mouseEventsDisabled = true;
	}
	
	/**
	 * Enables detecting escape key presses to cancel operations.
	 */
	private void enableCancelWorkInProgressKey() {
		addKeyListener(cancelWorkInProgress);
	}
	
	/**
	 * Disables detecting escape key presses to cancel operations.
	 */
	private void disableCancelWorkInProgressKey() {
		removeKeyListener(cancelWorkInProgress);
	}
	
	/**
	 * Enables user contols.
	 */
	private void enableUserControls() {
		enableMenus();
		
		enableControlButtons();
		updateStartStopResumeButtonLabel();
		updateZoomOnSelectedButtonLabel();
		
		enableMouseEventHandling();
		disableCancelWorkInProgressKey();
	}

	/**
	 * Disables user contols. This is useful to block the user from interacting with the GUI in long running operations
	 * like saving and loading, where he could break s.th. with his input.
	 */
	private void disableUserControls() {
		disableMenus();
		disableControlButtons();
		disableMouseEventHandling();
		enableCancelWorkInProgressKey();
	}

	/**
	 * Setup the number of rounds as chosen in the options window.
	 * 
	 * @param new_limit
	 *            The new rounds number limit.
	 * @return true if the options window needs to refresh its input field values from the model.
	 */
	private boolean setRoundsLimit(int new_limit) {
		if (new_limit > 1) {
			data.setRoundsLimit(new_limit);
			updateStartStopResumeButtonLabel();
			return false;
		}
		
		data.disableRoundsLimit();
		updateStartStopResumeButtonLabel();

		// disabling might change user input e.g. from -5 to -1, hence we need to update the option window values
		return true;
	}
	
	/**
	 * Opens an input dialog where users can choose to only repaint the GUI after every N rounds of calculation.
	 */
	private void skipdrawingroundsActionPerformed() {
		String tmpstr;
		tmpstr = javax.swing.JOptionPane.showInputDialog(this, "Draw each Nth round. N=", String.valueOf(skiprounds));
		try {
			if (tmpstr != null) {
				skiprounds = Integer.parseInt(tmpstr);
			}
		} catch (NumberFormatException ne) {
			javax.swing.JOptionPane.showMessageDialog(this, "ERROR, unable to parse integer from '" + tmpstr + "'");
		}
	}

	/**
	 * Opens the "Help" window.
	 */
	private void openHelpWindow() {
		new WindowHelp(this, true).setVisible(true);
	}

	/**
	 * Opens the "About" window.
	 */
	private void openAboutWindow() {
		new WindowAbout(this, true).setVisible(true);
	}

	/**
	 * Opens a window that offers the currently selected sequences in FASTA format as copy & pastable text.
	 */
	private void openShowCopyPasteableSequencesWindow() {
		if (data.selectedSequencesIndices.length < 1) {
			javax.swing.JOptionPane.showMessageDialog(this, "Please select some sequences");
			return;
		}
		
		StringBuffer outbuff = new StringBuffer();
		for (int i = 0; i < data.selectedSequencesIndices.length; i++) {
			outbuff.append(">" + data.sequence_names[data.selectedSequencesIndices[i]] + " " + data.selectedSequencesIndices[i] + "\n");
			outbuff.append(data.sequences[data.selectedSequencesIndices[i]].seq + "\n");
		}

		new WindowShowCopyPasteableSequences(new javax.swing.JFrame(), outbuff).setVisible(true);
	}

	/**
	 * Opens a dialog that lets users choose a new zoom factor.
	 */
	private void openChangeZoomFactorDialog() {
		String user_input = "";

		user_input = javax.swing.JOptionPane.showInputDialog(this, "New zoom factor (in percent)",
				String.valueOf((int) (draw_area.getZoomFactor() * 100)));

		if (user_input == null) {
			return;
		}

		float new_zoom_factor;
		try {
			new_zoom_factor = Integer.parseInt(user_input) / 100f;

			if (new_zoom_factor <= 0) {
				throw new NumberFormatException();
			}

		} catch (NumberFormatException ne) {
			javax.swing.JOptionPane.showMessageDialog(this, user_input
					+ " is not a positive integer number; zoom not changed");
			return;
		}

		setZoom(new_zoom_factor);
	}
		
	/**
	 * Returns whether the user wants to zoom in on the selected sequences.
	 */
	boolean isZoomingOnSelectedSequences() {
		return zoomOnSelected;
	}

	/**
	 * Updates the zoom of the view to the currently set zoom value while keeping the same area centered.
	 * 
	 * @param old_zoom
	 *            The previous zoom value.
	 */
	private void setZoom(float new_zoom_factor) throws IllegalStateException {
		draw_area.setZoom(new_zoom_factor);
		
		informAboutZoomChange();
	}
	
	/**
	 * Resets the zoom to none and centers the graph.
	 */
	private void resetZoom() {
		zoomOnSelected = false;

		draw_area.resetZoom();
		
		informAboutZoomChange();
	}

	/**
	 * Centers the graph while maintaining the current zoom level.
	 */
	private void centerGraph() {
		draw_area.centerGraph();
	}
	
	/**
	 * Opens a dialog to let the user change the GUI font.
	 */
	private void openChangeFontDialog() {
		draw_area.myfont = WindowFontChoice.getfont("Select Font", draw_area.myfont);
		repaint();
		draw_area.repaint();
	}
	
	/**
	 * Opens a window in which the user can create and manage custom groups of sequences. These group definitions are
	 * used to draw differently colored shapes in the draw area instead of simple dots with identical colors.
	 */
	private void openEditGroupsWindow() {

		if (myseqgroupwindow != null) {
			myseqgroupwindow.setVisible(false);
			myseqgroupwindow.dispose();
		}
		
		myseqgroupwindow = new WindowEditGroups(this);
		myseqgroupwindow.setVisible(true);
	}

	/**
	 * Opens a window that lets users see all sequences with connections to the currently selected sequences. Further
	 * options like the selection of those connected sequences are available to the user.
	 */
	private void openGetSequencesConnectedToSelectedsWindow() {
		int[] blasthitsarr = WindowSelectedsAndTheirMatches.getblasthits(data.attractionValues,
				data.selectedSequencesIndices, data.sequence_names);
		new WindowSelectedsAndTheirMatches(this, blasthitsarr, data.selectedSequencesIndices).setVisible(true);
	}

	/**
	 * Opens the "find clusters" window to let the user perform automatic cluster detection on the data.
	 */
	private void openFindClustersWindow() {
		// don't define optionsvec as it contains both strings and numbers in a defined order
		Vector<String> optionsvec = new Vector<String>();
		new WindowAutomaticClusterDetection(this, optionsvec).setVisible(true);
		if (optionsvec.size() == 0) {// if I canceled
			return;
		}
		boolean didbootstrap = false;
		String tmpstr = (String) optionsvec.remove(0);
		if (tmpstr.equals("convex")) {
			tmpstr = (String) optionsvec.remove(0);
			try {
				float sigmafac = Float.parseFloat(tmpstr);
				tmpstr = (String) optionsvec.remove(0);
				int minseqnum = Integer.parseInt(tmpstr);
				
				System.out.println("searching for convex clusters");
				Vector<SequenceCluster> clustervec = ClusterDetection.getConvex(data.attractionValues, sigmafac, minseqnum,
						data.elements);
				
				System.out.println("done searching for clusters; opening window");
				if (((String) optionsvec.remove(0)).equalsIgnoreCase("true")) {// if do bootstrap
					didbootstrap = true;
					tmpstr = (String) optionsvec.remove(0);
					int replicates = Integer.parseInt(tmpstr);
					tmpstr = (String) optionsvec.remove(0);
					float remove = Float.parseFloat(tmpstr);
					if (remove > 1) {
						remove /= 100;
					}
					
					if (ClusterDetectionBootstrapping.bootstrapconvex(data.attractionValues, clustervec, "convex", replicates,
							remove, sigmafac, minseqnum, data.elements) == false) {
						javax.swing.JOptionPane.showMessageDialog(this, "ERROR while bootstrapping");
						return;
					}
				}
				new WindowClusterDetectionResults(this, clustervec, "convex: " + minseqnum + ";" + sigmafac, didbootstrap)
						.setVisible(true);
			} catch (NumberFormatException ne) {
				javax.swing.JOptionPane.showMessageDialog(this, "Unable to parse float from " + tmpstr);
			}
		} else if (tmpstr.equals("linkage")) {
			tmpstr = (String) optionsvec.remove(0);
			try {
				int minlinkage = Integer.parseInt(tmpstr);
				tmpstr = (String) optionsvec.remove(0);
				int minseqnum = Integer.parseInt(tmpstr);

				System.out.println("searching for linkage clusters");
				Vector<SequenceCluster> clustervec = ClusterDetection.getlinkage(data.attractionValues, minlinkage, minseqnum,
						data.elements);
				
				System.out.println("done searching for clusters; opening window");
				if (((String) optionsvec.remove(0)).equalsIgnoreCase("true")) {// if do bootstrap
					didbootstrap = true;
					tmpstr = (String) optionsvec.remove(0);
					int replicates = Integer.parseInt(tmpstr);
					tmpstr = (String) optionsvec.remove(0);
					float remove = Float.parseFloat(tmpstr);
					if (remove > 1) {
						remove /= 100;
					}
					if (ClusterDetectionBootstrapping.bootstraplinkage(data.attractionValues, clustervec, "linkage", replicates,
							remove, minlinkage, minseqnum, data.elements) == false) {
						javax.swing.JOptionPane.showMessageDialog(this, "Error while bootstrapping");
						return;
					}
				}
				new WindowClusterDetectionResults(this, clustervec, "linkage: " + minseqnum + ";" + minlinkage, didbootstrap)
						.setVisible(true);
			} catch (NumberFormatException ne) {
				javax.swing.JOptionPane.showMessageDialog(this, "Unable to parse int from " + tmpstr);
			}
		} else if (tmpstr.equals("network")) {
			tmpstr = (String) optionsvec.remove(0);
			try {
				int minseqnum = Integer.parseInt(tmpstr);
				boolean dooffset = false;
				boolean globalaverage = false;
				tmpstr = (String) optionsvec.remove(0);
				if (tmpstr.equalsIgnoreCase("true")) {
					dooffset = true;
				}
				tmpstr = (String) optionsvec.remove(0);
				if (tmpstr.equalsIgnoreCase("true")) {
					globalaverage = true;
				}
				int maxrounds = Integer.parseInt(optionsvec.remove(optionsvec.size() - 1));

				System.out.println("searching for network clusters, maxrounds=" + maxrounds);
				Vector<SequenceCluster> clustervec = ClusterDetection.getnetwork(data.attractionValues, minseqnum, dooffset,
						globalaverage, data.elements, maxrounds);
				
				System.out.println("done searching for clusters; opening window");
				if (((String) optionsvec.remove(0)).equalsIgnoreCase("true")) {// if do bootstrap
					didbootstrap = true;
					tmpstr = (String) optionsvec.remove(0);
					int replicates = Integer.parseInt(tmpstr);
					tmpstr = (String) optionsvec.remove(0);
					float remove = Float.parseFloat(tmpstr);
					if (remove > 1) {
						remove /= 100;
					}
					
					if (ClusterDetectionBootstrapping.bootstrapnetwork(data.attractionValues, clustervec, "network", replicates,
							remove, minseqnum, dooffset, globalaverage, data.elements, maxrounds) == false) {
						javax.swing.JOptionPane.showMessageDialog(this, "Error while bootstrapping");
						return;
					}
				}
				new WindowClusterDetectionResults(this, clustervec, "network:" + minseqnum + ";" + dooffset + ";" + globalaverage,
						didbootstrap).setVisible(true);
			} catch (NumberFormatException ne) {
				javax.swing.JOptionPane.showMessageDialog(this, "Unable to parse int from " + tmpstr);
			}
		} else {
			javax.swing.JOptionPane.showMessageDialog(this, "Error in selecting clustering method: " + tmpstr);
		}
	}

	/**
	 * Opens file choice dialog and loads the selected matrix format file (only one value per pair).
	 */
	private void loadDataMatrixFormat() {
		boolean restart_computation = stopComputation(true);

		groupseqs = null;
		int returnVal = fc.showOpenDialog(this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			ClusterDataLoadHelper saveddata = ClusterData.load_matrix_file(fc.getSelectedFile());
			if (saveddata.file != null) {// if the data was read all right
				// repaint = "Error Loading Data";
				data.sequences = ClusterMethods.removeGapsFromSequences(saveddata.inaln);
				int seqs = data.sequences.length;

				// data.nameshash holds info about which name is which array number
				data.nameshash = new HashMap<String, Integer>((int) (seqs / 0.75) + 1, (float) 0.75);

				data.sequence_names = new String[seqs];
				data.positions = new float[seqs][3];
				Random rand = ClusterMethods.rand;
				for (int i = 0; i < seqs; i++) {
					data.sequence_names[i] = data.sequences[i].name.trim();
					data.sequences[i].name = new String("sequence" + i);
					data.nameshash.put(data.sequences[i].name, new Integer(i));
					data.positions[i][0] = rand.nextFloat();
					data.positions[i][1] = rand.nextFloat();
					data.positions[i][2] = rand.nextFloat();
				}
				data.blasthits = null;
				data.orgattvals = null;
				data.attractionValues = saveddata.attvals;
				data.attractionValues = saveddata.attvals;
				setThresholdButton.setText("Use Attraction values better than");
				textfield_threshold_value.setText("0");
				data.elements = data.sequence_names.length;
				// now symmetrize and normalize the attvals to range from -1 to +1
				float minval = 0;
				float maxval = 0;
				for (int i = data.attractionValues.length - 1; i >= 0; i--) {
					if (data.attractionValues[i].att > maxval) {
						maxval = data.attractionValues[i].att;
					} else if (data.attractionValues[i].att < minval) {
						minval = data.attractionValues[i].att;
					}
				}
				if (-minval > maxval) {// decide wether to divide by maxval or -minval
					maxval = -minval;
				}
				for (int i = data.attractionValues.length - 1; i >= 0; i--) {
					// now normalize the values
					data.attractionValues[i].att /= maxval;
				}
				data.p2attfactor = maxval;
				data.selectedSequencesIndices = new int[0];
				data.seqgroupsvec = saveddata.seqgroupsvec;
				data.posarr = data.positions;
				data.movementsLastIteration = new float[data.elements][ClusterData.dimensions];
				data.movements = new float[data.elements][ClusterData.dimensions];
				data.posarrtmp = new float[data.elements][ClusterData.dimensions];
				data.drawarrtmp = new int[data.elements][ClusterData.dimensions];
				data.resetDrawOrder();
				data.attvalsimple = true;
				// repaint = null;
				data.pvalue_threshold = 1;
				textfield_threshold_value.setText("1");
				textfield_info_min_blast_evalue.setText("1");
			} else {// if the data had errors
				JOptionPane.showMessageDialog(this, "Error reading data", "Error reading", JOptionPane.ERROR_MESSAGE);
				return;
			}
		}
		int seqnum = data.sequence_names.length;
		System.out.println("seqnum=" + seqnum);
		data.seqlengths = new float[seqnum];
		float maxlength = 0;
		for (int i = 0; i < seqnum; i++) {
			data.seqlengths[i] = data.sequences[i].seq.replaceAll("-", "").length();
			if (data.seqlengths[i] > maxlength) {
				maxlength = data.seqlengths[i];
			}
		}
		if (maxlength > 0) {
			for (int i = 0; i < seqnum; i++) {
				data.seqlengths[i] /= maxlength;
			}
		}

		if (restart_computation) {
			startComputation();
		}

		draw_area.repaint();
	}

	/**
	 * Toggles between 2 and 3 dimensions for the graph.
	 */
	private void toggle2d3dGraph() {

		// reset rotatation
		MyMath.setTo3x3IdentityMatrix(data.rotmtx);
		MyMath.setTo3x3IdentityMatrix(data.myrotmtx);
		MyMath.setTo3x3IdentityMatrix(draw_area.tmprotmtx);

		// get rid of all z-axis information
		for (int i = data.positions.length; --i >= 0;) {
			data.positions[i][2] = 0;
		}
		
		if (cluster2dbutton.isSelected()) {
			data.cluster2d = true;
		
		} else {
			data.cluster2d = false;
		}
		
		draw_area.repaint();
	}

	/**
	 * Opens a standard print dialog to print the contents of the draw area.
	 */
	private void printDrawArea() {
		boolean restart_computation = stopComputation(true);
		boolean restart_autosave = pauseAutosave();

		java.awt.print.PrinterJob printJob = java.awt.print.PrinterJob.getPrinterJob();
		printJob.setPrintable(draw_area);
		if (printJob.printDialog()) {
			try {
				printJob.print();
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		}

		if (restart_autosave) {
			restartAutosave();
		}
		
		if (restart_computation) {
			startComputation();
		}
	}

	/**
	 * Opens a color choice dialog that is safe in keeping the original color if canceled.
	 * 
	 * @param title
	 *            The title of the dialog.
	 * @param current_color
	 *            The current color, which is used as starting color.
	 * @return The picked color or {@code current_color} if the dialog is closed/canceled.
	 */
	protected java.awt.Color safe_change_color_dialog(String title, java.awt.Color current_color) {

		java.awt.Color new_color = null;

		try {
			new_color = JColorChooser.showDialog(this, title, current_color);
		} catch (java.awt.HeadlessException e) {
			System.err.println("HeadlessException!");
		}

		if (new_color == null) {
			return current_color;
		}
		return new_color;
	}

	/**
	 * Opens a color chooser dialog to change the color of selected sequence circles.
	 */
	private void openSelectedColorChangeDialog() {
		draw_area.selectedSequenceCircleColor = safe_change_color_dialog("Choose a new selected color", draw_area.selectedSequenceCircleColor);
		draw_area.repaint();
	}

	/**
	 * Opens a color chooser dialog to change the color of blast hit circles.
	 */
	private void openBlastColorChangeDialog() {
		draw_area.blastHitCircleColor = safe_change_color_dialog("Choose a new blast circle color", draw_area.blastHitCircleColor);
		draw_area.repaint();
	}

	/**
	 * Opens the "selected sequences" windows that displays the names of all selected sequences or of all sequences if
	 * none are selected.
	 */
	private void openShowSelectedSequencesWindow() {
		
		if (!containsData(true)) {
			return;
		}
		
		if (shownames != null) {
			shownames.setVisible(false);
			shownames.dispose();
		}
		
		shownames = new WindowShowSelectedSequences(data.sequence_names, this);
		shownames.setVisible(true);
	}

	/**
	 * Opens a color chooser dialog to change the foreground color.
	 */
	private void openForegroundColorChangeDialog() {
		draw_area.mainColor = safe_change_color_dialog("Choose a new foreground color", draw_area.mainColor);
		draw_area.repaint();
	}

	/**
	 * Opens a color chooser dialog to change the background color.
	 */
	private void openBackgroundColorChangeDialog() {
		draw_area.innerAreaBackgroundColor = safe_change_color_dialog("Choose a new background color", draw_area.innerAreaBackgroundColor);
		draw_area.repaint();
	}

	private void save2dGraphDataAs() {
		boolean restart_computation = stopComputation(true);

		int returnVal = fc.showSaveDialog(this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			File savefile = fc.getSelectedFile();
			try {
				PrintWriter outwrite = new PrintWriter(new BufferedWriter(new FileWriter(savefile)));

				int namenum = data.sequence_names.length;
				outwrite.println("ID\tNAME\tX\tY");

				for (int i = 0; i < namenum; i++) {
					outwrite.println(i + "\t" + data.sequence_names[i] + "\t" + data.posarrtmp[i][0]
							/ draw_area.drawWidth + "\t" + data.posarrtmp[i][1] / draw_area.drawHeight);
				}

				outwrite.close();
			} catch (IOException ioe) {
				javax.swing.JOptionPane.showMessageDialog(this, "IOERROR writing to '" + savefile.getName() + "'");
			}
		}

		if (restart_computation) {
			startComputation();
		}
	}

	/**
	 * Opens a color chooser dialog to change the background color.
	 */
	private void openNumberColorChangeDialog() {
		draw_area.blastHitNumberColor = safe_change_color_dialog("Select New Color", draw_area.blastHitNumberColor);
		draw_area.repaint();
	}

	/**
	 * Sets the threshold according to the threshold textfield content when the return key is pressed in the textfield.
	 */
	private void confirmedThresholdTextfieldValue() {
		double new_threshold;
		try {
			new_threshold = Double.parseDouble(textfield_threshold_value.getText());
		} catch (NumberFormatException e) {
			javax.swing.JOptionPane.showMessageDialog(this, "ERROR; unable to parse double from '"
					+ textfield_threshold_value.getText() + "'");
			return;
		}
		set_threshold(new_threshold);
	}

	/**
	 * Sets the threshold according to the threshold textfield content when the "use values better than" button is
	 * pressed.
	 */
	private void setThresholdButtonPressed() {
		confirmedThresholdTextfieldValue();
	}

	/**
	 * Sets the threshold used for clustering. In case the clustering is running, it is temporarily stopped then
	 * resumed.
	 * 
	 * @param threshold
	 *            The new threshold.
	 */
	private void set_threshold(double threshold) {

		if (!containsData(true)) {
			return;
		}

		boolean restart_computation = stopComputation(true);

		data.pvalue_threshold = threshold;

		if (data.blasthits != null) {
			synchronized (data.attractionValues) {
				data.compute_attraction_values();
			}
		} else if (data.attractionValues != null) {// remove all attvals below the specified value
			if (data.orgattvals == null) {
				data.orgattvals = data.attractionValues;
			}
			data.attractionValues = ClusterMethods.filterAttractionValues(data.orgattvals, data.pvalue_threshold);
		}

		data.resetDrawOrder();
		draw_area.repaint();

		if (restart_computation) {
			startComputation();
		}
	}

	private void saveBlastMatrixPValuesAs() {

		boolean restart_computation = stopComputation(true);

		int returnVal = fc.showSaveDialog(this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			File savefile = fc.getSelectedFile();
			int elementsl = data.sequence_names.length;
			try {
				PrintWriter outwrite = new PrintWriter(new BufferedWriter(new FileWriter(savefile)));

				outwrite.println("sequences=" + elementsl);
				outwrite.println("<seqs>");

				for (int i = 0; i < elementsl; i++) {
					outwrite.println(">" + data.sequence_names[i]);
				}

				outwrite.println();
				outwrite.println("</seqs>");

				if (data.seqgroupsvec.size() > 0) {

					outwrite.println("#user defined sequence groups");
					outwrite.println("<seqgroups>");

					SequenceGroup mygroup;
					for (int i = data.seqgroupsvec.size() - 1; i >= 0; i--) {
						mygroup = (SequenceGroup) data.seqgroupsvec.elementAt(i);
						outwrite.println("name=" + mygroup.name);
						outwrite.println("color=" + mygroup.color.getRed() + ";" + mygroup.color.getGreen() + ";"
								+ mygroup.color.getBlue());
						outwrite.print("numbers=");

						for (int j = mygroup.sequences.length - 1; j >= 0; j--) {
							outwrite.print(mygroup.sequences[j] + ";");
						}

						outwrite.println();
					}

					outwrite.println("</seqgroups>");
				}

				outwrite.println("<att>");
				MinimalAttractionValue[] myattvals = data.attractionValues;
				int datnum = myattvals.length;

				for (int i = 0; i < datnum; i++) {
					outwrite.println(myattvals[i].query + " " + myattvals[i].hit + " :" + myattvals[i].att);
				}

				outwrite.println("</att>");
				outwrite.close();
				System.out.println("done");
			} catch (IOException e) {
				javax.swing.JOptionPane.showMessageDialog(this, "IOERROR writing to '" + savefile.getAbsolutePath()
						+ "'");
			}
		}

		if (restart_computation) {
			startComputation();
		}
	}

	/**
	 * Opens a dialog that lets users change the size of all dots.
	 */
	private void openChaneDotSizeDialog() {
		String user_input = "";
		try {
			user_input = JOptionPane.showInputDialog(this, "Enter the new size (int):", String.valueOf(data.dotsize));
			if (user_input != null) {
				data.dotsize = (int) (Float.parseFloat(user_input));
			}
		} catch (NumberFormatException ne) {
			javax.swing.JOptionPane.showMessageDialog(this, "ERROR, unable to parse integer from '" + user_input + "'");
		}
		draw_area.repaint();
	}

	/**
	 * Opens a dialog that lets users enter specific rotation matrix values.
	 */
	private void openRotationMatrixInputDialog() {
		String tmpstr = "";
		try {
			tmpstr = JOptionPane.showInputDialog(this, "Enter the new rotation values: x , y , z (9 values total)");
			if (tmpstr == null) {
				return;
			}
			String[] tmparr = tmpstr.split(",");
			if (tmparr.length != 9) {
				JOptionPane.showMessageDialog(this, "You have to enter nine values separated by commas ','", "Error",
						JOptionPane.ERROR_MESSAGE);
				return;
			}
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					data.myrotmtx[i][j] = Double.parseDouble(tmparr[i * 3 + j]);
				}
			}
			data.rotmtx[0][0] = data.myrotmtx[0][0];
			data.rotmtx[0][1] = data.myrotmtx[0][1];
			data.rotmtx[0][2] = data.myrotmtx[0][2];
			data.rotmtx[1][0] = data.myrotmtx[1][0];
			data.rotmtx[1][1] = data.myrotmtx[1][1];
			data.rotmtx[1][2] = data.myrotmtx[1][2];
			data.rotmtx[2][0] = data.myrotmtx[2][0];
			data.rotmtx[2][1] = data.myrotmtx[2][1];
			data.rotmtx[2][2] = data.myrotmtx[2][2];
		} catch (NumberFormatException ne) {
			javax.swing.JOptionPane.showMessageDialog(this, "ERROR, unable to parse double from '" + tmpstr + "'");
		}
		draw_area.repaint();
	}

	/**
	 * Opens a dialog that lets users set the size for circles representing selected sequences.
	 */
	private void openChangeSelectedCircleSizeDialog() {
		// set the size for the circles for selected sequences
		String tmpstr = "";
		try {
			tmpstr = JOptionPane.showInputDialog(this, "Enter the new size(int):", String.valueOf(data.ovalsize));
			if (tmpstr != null) {
				data.ovalsize = Integer.parseInt(tmpstr);
			}
		} catch (NumberFormatException ne) {
			javax.swing.JOptionPane.showMessageDialog(this, "ERROR, unable to parse int from '" + tmpstr + "'");
		}
		draw_area.repaint();
	}

	/**
	 * Opens a dialog that lets users pick a sequence for which BLAST is run against the present database. The region
	 * and e-values of the HSPs are then shown mapped to this sequence.
	 * <p>
	 * i.e. region 1-200 hits cluster A, region 210-300 cluster b, ergo 2 domains
	 */
	private void openShowBlastHitsForSequenceDialog() {
		int referenceseqnum = WindowShowBlastHits.getrefseq(data.sequence_names);
		if (referenceseqnum == -1) {
			javax.swing.JOptionPane.showMessageDialog(this, "Please select a sequence");
			return;
		}

		// get the blast hits to this sequence
		HighScoringSegmentPair[] thishsp = BlastHitHandling.gethsps(referenceseqnum, data.sequences, data.cmd,
				data.formatdbpath, data.blastpath, data.addblastvbparam, data.referencedb, data.mineval,
				data.pvalue_threshold);
		
		// plot these blast hits on to the sequence
		WindowShowBlastHitsOnSequence myview = new WindowShowBlastHitsOnSequence(this, thishsp, referenceseqnum, data.sequence_names,
				data.sequences[referenceseqnum], data.nameshash);
		viewblasthitsvec.addElement(myview);
		
		myview.setVisible(true);
	}

	/**
	 * Opens a file chooser dialog that asks the user whether a selected file should be overwritten if it exists
	 * 
	 * @return The selected filename or "" if the dialog was canceled or closed.
	 */
	private String safe_output_file_chooser() {
		// this line makes using TAB possible to go from Yes to No in the FileChooser. Otherwise pressing Enter ALWAYS
		// selected YES, irrespective of the currently selected button
		UIManager.put("Button.defaultButtonFollowsFocus", Boolean.TRUE);

		JFileChooser jFileChooser = new JFileChooser(new File(System.getProperty("user.dir"))) {
			/**
			 * For existing files, ask whether to overwrite or not or cancel the save.
			 * 
			 * This code snippet stems from: http://stackoverflow.com/a/3729157/454402
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public void approveSelection() {
				File f = getSelectedFile();
				if (f.exists() && getDialogType() == SAVE_DIALOG) {
					int result = JOptionPane.showConfirmDialog(this, "The file exists, overwrite?", "Existing file",
							JOptionPane.YES_NO_CANCEL_OPTION);
					switch (result) {
					case JOptionPane.YES_OPTION:
						super.approveSelection();
						return;
					case JOptionPane.NO_OPTION:
						return;
					case JOptionPane.CLOSED_OPTION:
						return;
					case JOptionPane.CANCEL_OPTION:
						cancelSelection();
						return;
					}
				}
				super.approveSelection();
			}
		};

		// set input filename as default output filename
		if (data.getAbsoluteInputfileName() != null) {
			jFileChooser.setSelectedFile(new File(data.getAbsoluteInputfileName()));
		}

		int returnVal = jFileChooser.showSaveDialog(this);

		if (returnVal == JFileChooser.APPROVE_OPTION) {
			return jFileChooser.getSelectedFile().getPath();
		}

		return null; // signal "no file was selected" to caller
	}

	/**
	 * Saves the run to a CLANS format file and informs the user about progress and success in overlay and STDERR
	 * messages.
	 * 
	 * @param output_filename
	 *            The file to which the output is written.
	 */
	private void threaded_save_run(final String output_filename) {

		if (messageOverlayActive) {
			messageOverlay.setSaving();	
		} else {
			System.out.println("Note: saving can be canceled by pressing escape");
		}

		if (output_filename == null) {
			if (messageOverlayActive) {
				messageOverlay.setCanceled();
			}
			return;
		}
		
		final boolean restart_computation = stopComputation(true);
		final boolean restart_autosave = pauseAutosave();
		
		disableUserControls();

		if (saveLoadWorker != null) {
			saveLoadWorker.cancel(true); // kill potentially running old worker
		}
		
		saveLoadWorker = new SwingWorker<Void, Integer>() {

			@Override
			protected Void doInBackground() throws Exception {
				Thread.currentThread().setName("save in background");
				data.safer_save_to_file(output_filename, this);
				return null;
			}

			protected void done() {
				String error_message = "";
				boolean saving_successful = true;
				boolean saving_canceled = false;

				try {
					// we query the output of doInBackground only to get to the catch blocks in case of failures
					get();

				} catch (CancellationException e) {
					saving_successful = false;
					saving_canceled = true;
					error_message = "saving canceled";

				} catch (InterruptedException e) {
					saving_successful = false;
					error_message = e.getMessage();

				} catch (ExecutionException e) {
					saving_successful = false;
					
					Throwable exception = e.getCause();
					if (exception instanceof IOException) {
						error_message = exception.getMessage();

					} else if (exception instanceof IllegalStateException) {
						error_message = exception.getMessage();

					} else {
						error_message = "unexpected exception: " + e.getMessage();
					}
				}

				if (saving_successful) {

					if (autosaveAutoReenable) {
						initializeAutosave(false);
						autosaveAutoReenable = false;
					}

					setTitle("Clustering of " + data.getBaseInputfileName());

					if (messageOverlayActive) {
						messageOverlay.setCompleted();
					}

					if (restart_computation) {
						startComputation();
					}

				} else {
					System.err.println(error_message);

					if (messageOverlayActive) {
						if (saving_canceled) {
							messageOverlay.setCanceled();
						} else {
							messageOverlay.setFailed(error_message);
						}
					}
					
					if (restart_autosave) {
						restartAutosave();
					}
				}

				enableUserControls();
			}
		};

		saveLoadWorker.execute();
	}

	/**
	 * Gets the filename associated with the currently open data.
	 * 
	 * @return The filename or null if the data is not associated with a filename yet.
	 */
	private String getSaveFilename() {
		return data.getAbsoluteInputfileName();
	}

	/**
	 * Opens a file chooser dialog to let the user pick an output file. If the file exists, the user will be asked
	 * whether to use it or select a different one or cancel.
	 * 
	 * @return The filename or null if the user closed/canceled the input dialog.
	 */
	private String getSaveAsFilename() {
		return safe_output_file_chooser();
	}

	/**
	 * Saves the data in CLANS format to the filename associated with it. Used by menu entry "File->Save Run".
	 * <p>
	 * If the data is not associated with a filename (e.g. CLANS parsed BLAST results), this falls back to
	 * {@code saveRunAs} behavior of showing the user a file choice dialog.
	 */
	private void saveRun() {
		String output_filename = getSaveFilename();

		if (output_filename == null) {
			output_filename = getSaveAsFilename();

			if (output_filename == null) {
				return;
			}
		}

		threaded_save_run(output_filename);
	}

	/**
	 * Saves the data in CLANS format to a user chosen file. Used by menu entry "File->Save Run As". If the file
	 * selection dialog is canceled/closed, no save file is created.
	 */
	private void saveRunAs() {
		String output_filename = getSaveAsFilename();

		if (output_filename == null) {
			return;
		}

		threaded_save_run(output_filename);
	}

	/**
	 * Saves the data in CLANS format to its associated file. Used by the autosave feature to save data on a regular
	 * basis. If the data is not associated with a filename (e.g. CLANS parsed BLAST results), activating the autosave
	 * feature is postponed until a manual save is successfully performed.
	 */
	private void autosave() {
		String output_filename = getSaveFilename();

		if (output_filename == null) {
			autosaveAutoReenable = true;
			disableAutosave(false);
			return;
		}

		threaded_save_run(output_filename);
	}

	private void loadDataClansFormat() {

		groupseqs = null;

		int returnVal = fc.showOpenDialog(this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {

			String filename = fc.getSelectedFile().getAbsolutePath();
			loadDataClansFormatInBackground(filename);

			if (myseqgroupwindow != null) {
				myseqgroupwindow.setVisible(false);
				myseqgroupwindow.dispose();
			}
			if (mymapfunctiondialog != null) {
				mymapfunctiondialog.setVisible(false);
				mymapfunctiondialog.dispose();
			}
			draw_area.repaint();

		} else if (returnVal == JFileChooser.CANCEL_OPTION) {
			if (messageOverlayActive) {
				messageOverlay.setCanceled();
			}
		}
	}

	/**
	 * Permanently deletes sequences based on a selection mask.
	 * 
	 * @param keep_mask
	 *            An array with true at indices of sequences that should be kept and false at indices of sequences that
	 *            should be deleted.
	 */
	private void deleteSequences(boolean[] keep_mask) {
		int counter = 0;
		// shift for sequence ids in groups
		int[] sequence_index_shift = new int[keep_mask.length];
		for (int i = 0; i < keep_mask.length; i++) {
			sequence_index_shift[i] = i - counter;
			if (keep_mask[i]) {
				counter++;
			}
		}

		if (data.selectedSequencesIndices.length < 2) {
			return;
		}

		if (shownames != null) {
			shownames.setVisible(false);
			shownames.dispose();
			shownames = null;
		}

		if (myseqgroupwindow != null) {
			myseqgroupwindow.setVisible(false);
			myseqgroupwindow.dispose();
			myseqgroupwindow = null;
		}

		if (viewblasthitsvec.size() > 0) {
			for (int i = 0; i < viewblasthitsvec.size(); i++) {
				blastselectseqs = new int[0];
				viewblasthitsvec.elementAt(i).setVisible(false);
				viewblasthitsvec.elementAt(i).dispose();
			}
			viewblasthitsvec.setSize(0);
		}

		level++;
		getparentmenuitem.setEnabled(true);
		getparentmenuitem.setText("use parent group (" + level + ")");
		parentnameshash.addElement(data.nameshash);

		// holds info about which name is which array number
		data.nameshash = new HashMap<String, Integer>((int) (data.selectedSequencesIndices.length / 0.75) + 1, (float) 0.75);
		for (int i = 0; i < data.selectedSequencesIndices.length; i++) {
			data.nameshash.put(data.sequences[data.selectedSequencesIndices[i]].name, new Integer(i));
		}

		if (data.blasthits != null) {
			parentblasthits.addElement(data.blasthits);
			data.blasthits = SelectedSubsetHandling.get_blasthits(data.blasthits, data.selectedSequencesIndices);
		} else {
			if (data.orgattvals == null) {
				data.orgattvals = data.attractionValues;
			}
			parent_orgattvals.addElement(data.orgattvals);
			data.orgattvals = SelectedSubsetHandling.get_myattvals(data.orgattvals, data.selectedSequencesIndices);
			data.compute_attraction_values();
		}

		parentmovearr.addElement(data.movements);
		data.movements = SelectedSubsetHandling.get_mymovearr(data.movements, data.selectedSequencesIndices);
		if (data.movements.length > 0) {
			data.movementsLastIteration = new float[data.movements.length][data.movements[0].length];
		}

		parentposarr.addElement(data.positions);
		data.positions = SelectedSubsetHandling.get_myposarr(data.positions, data.selectedSequencesIndices);
		parentaln.addElement(data.sequences);
		data.sequences = SelectedSubsetHandling.get_sequences(data.sequences, data.selectedSequencesIndices);
		parentnamearr.addElement(data.sequence_names);
		data.sequence_names = SelectedSubsetHandling.get_namearr(data.sequence_names, data.selectedSequencesIndices);
		data.selectedSequencesIndices = new int[0];
		if (data.blasthits != null) {
			synchronized (data.attractionValues) {
				data.compute_attraction_values();
			}
		}

		// fix the sequence ids stored in groups
		SequenceGroup current_seqgroup;
		Vector<Integer> keep_group_sequences;
		for (int i = 0; i < data.seqgroupsvec.size(); i++) {
			current_seqgroup = data.seqgroupsvec.elementAt(i);

			// collect and fix the sequence ids to keep for this group
			keep_group_sequences = new Vector<Integer>();
			for (int j = 0; j < current_seqgroup.sequences.length; j++) {
				if (keep_mask[current_seqgroup.sequences[j]]) {
					keep_group_sequences.add(current_seqgroup.sequences[j]
							- sequence_index_shift[current_seqgroup.sequences[j]]);
				}
			}

			// copy the ids over to the seqgroup
			current_seqgroup.sequences = new int[keep_group_sequences.size()];
			for (int j = 0; j < keep_group_sequences.size(); j++) {
				current_seqgroup.sequences[j] = keep_group_sequences.elementAt(j);
			}

		}

		data.elements = data.sequence_names.length;
		data.posarr = data.positions;
		data.posarrtmp = new float[data.elements][ClusterData.dimensions];
		data.drawarrtmp = new int[data.elements][ClusterData.dimensions];

		data.resetDrawOrder();
		draw_area.repaint();
	}

	/**
	 * Permanently remove all currently selected sequences from the data.
	 */
	private void removeAllSelectedSequences() {

		boolean restart_computation = stopComputation(true);

		// selected_group_sequences == groupseqs
		groupseqs = null;

		int sequences = data.sequence_names.length;
		int[] sequences_to_remove = data.selectedSequencesIndices.clone();
		data.selectedSequencesIndices = new int[sequences - data.selectedSequencesIndices.length];
		boolean[] sequences_to_keep = new boolean[sequences];
		for (int i = 0; i < sequences; i++) {
			sequences_to_keep[i] = true;
		}

		for (int i = 0; i < sequences_to_remove.length; i++) {
			sequences_to_keep[sequences_to_remove[i]] = false;
		}

		int counter = 0;
		for (int i = 0; i < sequences; i++) {
			if (sequences_to_keep[i]) {
				data.selectedSequencesIndices[counter] = i;
				counter++;
			}
		}

		deleteSequences(sequences_to_keep);
		
		informAboutSequenceSelectionChange();

		if (restart_computation) {
			startComputation();
		}
	}

	/**
	 * Permanently removes all sequences without connections at the currently selected p-value cutoff.
	 */
	private void hideSingletonsMenuItemActionPerformed() {

		boolean restart_computation = stopComputation(true);

		groupseqs = null;
		int sequence_number = data.sequence_names.length;
		boolean[] sequences_to_keep = new boolean[sequence_number];

		for (int i = 0; i < data.attractionValues.length; i++) {
			if (data.attractionValues[i].att != 0) {
				sequences_to_keep[data.attractionValues[i].query] = true;
				sequences_to_keep[data.attractionValues[i].hit] = true;
			}
		}

		int counter = 0;
		// shift for sequence ids in groups
		int[] sequence_index_shift = new int[sequence_number];
		for (int i = 0; i < sequence_number; i++) {
			sequence_index_shift[i] = i - counter;
			if (sequences_to_keep[i]) {
				counter++;
			}
		}

		data.selectedSequencesIndices = new int[counter];
		counter = 0;
		for (int i = 0; i < sequence_number; i++) {
			if (sequences_to_keep[i]) {
				data.selectedSequencesIndices[counter] = i;
				counter++;
			}
		}

		deleteSequences(sequences_to_keep);

		if (restart_computation) {
			startComputation();
		}
	}

	/**
	 * Go to super-set of currently shown subset (use all sequences from the level before).
	 */
	private void useGraphSuperset() {

		boolean restart_computation = stopComputation(true);

		groupseqs = null;
		if (shownames != null) {
			shownames.setVisible(false);
			shownames.dispose();
			shownames = null;
		}
		if (myseqgroupwindow != null) {
			myseqgroupwindow.setVisible(false);
			myseqgroupwindow.dispose();
			myseqgroupwindow = null;
		}
		if (viewblasthitsvec.size() > 0) {
			for (int i = 0; i < viewblasthitsvec.size(); i++) {
				blastselectseqs = new int[0];
				viewblasthitsvec.elementAt(i).setVisible(false);
				viewblasthitsvec.elementAt(i).dispose();
			}
			viewblasthitsvec.setSize(0);
		}
		level--;
		if (level == 0) {
			getparentmenuitem.setEnabled(false);
		}
		getparentmenuitem.setText("use parent group (" + level + ")");
		if (data.blasthits != null) {
			data.blasthits = (MinimalHsp[]) parentblasthits.elementAt(level);
			parentblasthits.removeElementAt(level);
		} else {
			data.orgattvals = (MinimalAttractionValue[]) parent_orgattvals.elementAt(level);
			parentblasthits.removeElementAt(level);
			data.attractionValues = ClusterMethods.filterAttractionValues(data.orgattvals, data.pvalue_threshold);
		}
		data.movements = (float[][]) parentmovearr.elementAt(level);
		if (data.movements.length > 0) {
			data.movementsLastIteration = new float[data.movements.length][data.movements[0].length];
		}
		parentmovearr.removeElementAt(level);
		data.positions = (float[][]) parentposarr.elementAt(level);
		parentposarr.removeElementAt(level);
		data.sequences = (AminoAcidSequence[]) parentaln.elementAt(level);
		parentaln.removeElementAt(level);
		data.sequence_names = (String[]) parentnamearr.elementAt(level);
		parentnamearr.removeElementAt(level);
		data.nameshash = parentnameshash.elementAt(level);
		parentnameshash.removeElementAt(level);
		data.weights = (float[]) parentweights.elementAt(level);
		parentweights.removeElementAt(level);
		data.selectedSequencesIndices = new int[0];
		data.elements = data.sequence_names.length;
		if (data.blasthits != null) {
			synchronized (data.attractionValues) {
				data.compute_attraction_values();
			}
		}
		if (myseqgroupwindow != null) {
			myseqgroupwindow.setVisible(false);
			myseqgroupwindow.dispose();
		}
		if (parentseqgroups.size() > level) {
			data.seqgroupsvec = parentseqgroups.remove(level);
		}
		if (mymapfunctiondialog != null) {
			mymapfunctiondialog.makenameshash();
		}
		data.posarr = data.positions;
		data.posarrtmp = new float[data.elements][ClusterData.dimensions];
		data.drawarrtmp = new int[data.elements][ClusterData.dimensions];
		data.resetDrawOrder();
		draw_area.repaint();

		if (restart_computation) {
			startComputation();
		}
	}

	/**
	 * Go to sub-set of currently shown sequences (use all sequences from the level before).
	 */
	private void useGraphSubset() {

		boolean restart_computation = stopComputation(true);

		groupseqs = null;
		int selectednum = data.selectedSequencesIndices.length;
		if (selectednum < 2) {
			enableButtonStartStopResume();
			return;
		}
		if (shownames != null) {
			shownames.setVisible(false);
			shownames.dispose();
			shownames = null;
		}
		if (myseqgroupwindow != null) {
			myseqgroupwindow.setVisible(false);
			myseqgroupwindow.dispose();
			myseqgroupwindow = null;
		}
		if (myseqgroupwindow != null) {
			myseqgroupwindow.setVisible(false);
			myseqgroupwindow.dispose();
			data.showSequenceGroups = false;
		}
		selectvec.clear();
		if (viewblasthitsvec.size() > 0) {
			for (int i = 0; i < viewblasthitsvec.size(); i++) {
				blastselectseqs = new int[0];
				((WindowShowBlastHitsOnSequence) viewblasthitsvec.elementAt(i)).setVisible(false);
				((WindowShowBlastHitsOnSequence) viewblasthitsvec.elementAt(i)).dispose();
			}
			viewblasthitsvec.setSize(0);
		}
		level++;
		getparentmenuitem.setEnabled(true);
		getparentmenuitem.setText("use parent group (" + level + ")");
		parentnameshash.addElement(data.nameshash);
		data.nameshash = new HashMap<String, Integer>((int) (selectednum / 0.75) + 1, (float) 0.75);// holds info about
																									// which name is
																									// which array
																									// number
		for (int i = 0; i < selectednum; i++) {
			data.nameshash.put(data.sequences[data.selectedSequencesIndices[i]].name, new Integer(i));
		}
		if (data.blasthits != null) {
			parentblasthits.addElement(data.blasthits);
			data.blasthits = SelectedSubsetHandling.get_blasthits(data.blasthits, data.selectedSequencesIndices);
		}
		parentmovearr.addElement(data.movements);
		data.movements = SelectedSubsetHandling.get_mymovearr(data.movements, data.selectedSequencesIndices);
		parentposarr.addElement(data.positions);
		data.positions = SelectedSubsetHandling.get_myposarr(data.positions, data.selectedSequencesIndices);
		parentaln.addElement(data.sequences);
		data.sequences = SelectedSubsetHandling.get_sequences(data.sequences, data.selectedSequencesIndices);
		parentnamearr.addElement(data.sequence_names);
		data.sequence_names = SelectedSubsetHandling.get_namearr(data.sequence_names, data.selectedSequencesIndices);
		parentweights.addElement(data.weights);
		data.weights = SelectedSubsetHandling.get_weights(data.weights, data.selectedSequencesIndices);
		data.elements = data.sequence_names.length;
		if (data.blasthits == null) {
			if (data.orgattvals == null) {
				data.orgattvals = data.attractionValues;
			}
			parent_orgattvals.addElement(data.orgattvals);
			data.orgattvals = SelectedSubsetHandling.get_myattvals(data.orgattvals, data.selectedSequencesIndices);
			data.attractionValues = ClusterMethods.filterAttractionValues(data.orgattvals, data.pvalue_threshold);
		}
		parentseqgroups.addElement(data.seqgroupsvec);
		data.seqgroupsvec = new Vector<SequenceGroup>();
		data.selectedSequencesIndices = new int[0];

		synchronized (data.attractionValues) {
			data.compute_attraction_values();
		}

		data.posarr = data.positions;
		data.posarrtmp = new float[data.elements][ClusterData.dimensions];
		data.drawarrtmp = new int[data.elements][ClusterData.dimensions];
		data.resetDrawOrder();
		if (mymapfunctiondialog != null) {
			mymapfunctiondialog.makenameshash();
		}

		draw_area.repaint();

		if (restart_computation) {
			startComputation();
		}
	}

	/**
	 * Opens a windows with a score distribution plot. The plot differs for similarity data given as scores, p-values,
	 * or attraction values.
	 */
	private void openScoreDistributionPlot() {
		if (data.blasthits != null) {
			if (data.usescval) { // score mode
				WindowPValuePlot eplot = new WindowPValuePlot(data.blasthits, data.pvalue_threshold, true);
				eplot.setVisible(true);

			} else { // p-value mode
				WindowPValuePlot eplot = new WindowPValuePlot(data.blasthits, data.pvalue_threshold, false);
				eplot.setVisible(true);
			}
		} else { // attraction value mode
			WindowAttractionValuePlot attplot = new WindowAttractionValuePlot(data.attractionValues, data.pvalue_threshold);
			attplot.setVisible(true);
		}
	}

	/**
	 * Saves the currently selected sequences in FASTA format to a user-chosen file.
	 */
	private void openSaveSequenceAsFastaDialog() {

		boolean restart_computation = stopComputation(true);

		int returnVal = fc.showSaveDialog(this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			AminoAcidSequence[] selectedseqs = getselectedseqs();
			if (selectedseqs.length == 0) {
				javax.swing.JOptionPane.showMessageDialog(this, "Please select some sequences");
			} else {
				FileHandling.printfasta(selectedseqs, fc.getSelectedFile());
			}
		}

		if (restart_computation) {
			startComputation();
		}
	}

	/**
	 * Opens a dialog for changing connection colors.
	 */
	private void openChangeDotConnectionColorDialog() {
		new DialogChangeConnectionColors(this);
		draw_area.repaint();
	}

	/**
	 * Deselects everything if sequences are currently selected or selects all sequences if none are currently selected.
	 */
	private void selectionButtonActivated() {

		if (!containsData(true)) {
			return;
		}

		if (data.selectedSequencesIndices.length > 0) {
			deselectAllSequences();

		} else {// if nothing is selected, select all
			selectAllSequences();
		}
		
		if (shownames != null) {
			shownames.seqnamelist.setSelectedIndices(data.selectedSequencesIndices);
		}

		draw_area.repaint();
	}

	/**
	 * Zooms on selected sequences or resets zoom if no sequences are selected.
	 * 
	 * @param show_error_message
	 *            If true, an error message dialog is shown if no sequences are selected.
	 */
	private void zoomOnSelectedTriggered(boolean show_error_message) {

		if (!containsData(true)) {
			zoomOnSelected = false;
			return;
		}
		
		if (isZoomingOnSelectedSequences()) {
			zoomOnSelected = false;

		} else if (hasSelectedSequences()) {
			zoomOnSelected = true;

		} else {
			if (show_error_message) {
				JOptionPane.showMessageDialog(null, "Please select some sequences", "Message",
						JOptionPane.INFORMATION_MESSAGE);
			}

			zoomOnSelected = false;

		}
		
		draw_area.resetZoom();
		
		informAboutZoomChange();
	}

	
	/**
	 * Based on the GUI state, updates the select/move button label.
	 */
	private void updateSelectionOrMovementModeButtonLabel() {
		if (selectionOrMovementModeButton.isSelected()) {
			selectionOrMovementModeButton.setText("SELECT/move");
		} else {
			selectionOrMovementModeButton.setText("select/MOVE");
		}
	}
	
	/**
	 * @return true if in selection mode, false if in movement mode
	 */
	protected boolean inSelectionMode() {
		return selectionOrMovementModeButton.isSelected();
	}
	
	/**
	 * @return true if any sequences are selected, false else
	 */
	private boolean hasSelectedSequences() {
		return data != null && data.selectedSequencesIndices != null && data.selectedSequencesIndices.length > 0;
	}
	
	/**
	 * @return the number of selected sequences
	 */
	int getNumberOfSelectedSequences() {
		if (!hasSelectedSequences()) {
			return 0;
		}
		return data.selectedSequencesIndices.length;
	}

	/**
	 * Based on the GUI state, updates the "select all"/"clear selection" button label.
	 */
	void updateSelectionButtonLabel() {
		if (hasSelectedSequences()) {
			selectionButton.setText("Clear Selection");
		} else {
			selectionButton.setText("Select All");
		}
	}
	
	/**
	 * Based on the GUI state, updates the "Zoom On Selected"/"Show all" button label and enabled state.
	 */
	private void updateZoomOnSelectedButtonLabel() {
		if (zoomOnSelected) {
			zoomOnSelectedButton.setText("Show all");
			zoomOnSelectedButton.setEnabled(true);

		} else {
			zoomOnSelectedButton.setText("Zoom on selected");
			zoomOnSelectedButton.setEnabled(getNumberOfSelectedSequences() > 1);
			zoomOnSelectedButton.setSelected(getNumberOfSelectedSequences() > 1);
		}
	}
	
	/**
	 * Based on the GUI state, updates the "Reset Zoom" button label and enabled state.
	 */
	private void updateResetZoomButtonLabel() {
		if (isZoomingOnSelectedSequences()) {
			// this is independent of the zoom factor, hence containsData below is not enough to cover this case
			resetZoomButton.setEnabled(true);
		
		} else if (containsData(false)) {
			resetZoomButton.setEnabled(Math.abs(draw_area.getZoomFactor() - 1) > 1e-3);

		} else {
			resetZoomButton.setEnabled(false);
		}
	}

	/**
	 * Handles starting mouse dragging events in the draw area, i.e. the mouse button has just been pressed and the
	 * mouse started to move so it's not a click any more but a drag.
	 * 
	 * @param evt
	 *            The mouse event.
	 */
	private void graphpanelMousePressed(java.awt.event.MouseEvent evt) {
	
		if (mouseEventsDisabled) {
			return;
		}
	
		mouse_is_pressed = true;
		
		if (inSelectionMode()) {
			mousemove[0] = 0;
			mousemove[1] = 0;
			selectstart[0] = evt.getX() - draw_area.xTranslate;
			selectstart[1] = evt.getY() - draw_area.yTranslate;
			currmousepos[0] = selectstart[0];
			currmousepos[1] = selectstart[1];
		
		} else {
			if (evt.isAltDown() || evt.isControlDown() || evt.isMetaDown()) {// if I want to drag a sequence in 2d
				// move all selected sequences a certain amount
				moveseqs = true;
				mousestart[0] = evt.getX();
				mousestart[1] = evt.getY();
			} else if (evt.isShiftDown()) {
				mousestart[0] = evt.getX();
				mousestart[1] = evt.getY();
				translate[0] = draw_area.xTranslate;
				translate[1] = draw_area.yTranslate;
			} else {
				mousestart[0] = evt.getX();
				mousestart[1] = evt.getY();
				// if mouse is inside a certain area rotate x,y; outside rotate z
				if (stereocheckboxmenuitem.isSelected() == false) {
					if ((mousestart[0] < 0) || (mousestart[0] > draw_area.drawWidth)
							|| (mousestart[1] < 0) || (mousestart[1] > draw_area.drawHeight)) {
						draw_area.rotatez = true;
					} else {
						draw_area.rotatez = false;
					}
				} else {// if in stereo mode
					if ((mousestart[0] < 0) || (mousestart[0] > draw_area.drawWidth * 2)
							|| (mousestart[1] < 0) || (mousestart[1] > draw_area.drawHeight)) {
						draw_area.rotatez = true;
					} else {
						draw_area.rotatez = false;
					}
				}
				mousemove[0] = 0;
				mousemove[1] = 0;
			}
		}

		draw_area.repaint();
	}

	/**
	 * Handles mouse dragging events in the draw area during the dragging.
	 * 
	 * @param evt
	 *            The mouse event.
	 */
	private void graphpanelMouseDragged(java.awt.event.MouseEvent evt) {
	
		if (mouseEventsDisabled) {
			return;
		}
	
		if (inSelectionMode()) {
			if (shownamesselectcheckbox.isSelected()) {
				if (shownames == null) {
					shownames = new WindowShowSelectedSequences(data.sequence_names, this);
					shownames.setVisible(true);
				}
				int[] tmpreg = new int[4];
				tmpreg[0] = selectstart[0];
				tmpreg[1] = selectstart[1];
				tmpreg[2] = evt.getX() - draw_area.xTranslate;
				tmpreg[3] = evt.getY() - draw_area.yTranslate;
				updatetmpselected(tmpreg);
			}
			mousemove[0] = 0;
			mousemove[1] = 0;
			currmousepos[0] = evt.getX() - draw_area.xTranslate;
			currmousepos[1] = evt.getY() - draw_area.yTranslate;
		
		} else {
			if (evt.isShiftDown()) {
				draw_area.xTranslate = evt.getX() - mousestart[0] + translate[0];
				draw_area.yTranslate = evt.getY() - mousestart[1] + translate[1];
			} else if (moveseqs) {

			} else {
				mousemove[0] = evt.getX() - mousestart[0];
				mousemove[1] = evt.getY() - mousestart[1];
			}
		}

		draw_area.repaint();
	}

	/**
	 * Repaints the draw area.
	 */
	private void requestRepaint() {
		mousemove[0] = 0;
		mousemove[1] = 0;
		draw_area.repaint();
	}

	/**
	 * Starts stopped or stops running computations when the start/stop/resume button is activated.
	 */
	private void toggleComputationRunning() {
		if (!containsData(true)) {
			return;
		}

		if (isComputing()) {
			stopComputation(false);
		} else {
			startComputation();
		}

		draw_area.repaint();
	}

	private void initializeGraphPositions() {
		if (!containsData(true)) {
			return;
		}
		
		if (isComputing()) {
			stopComputation(true);
		}

		resetGraph();

		draw_area.repaint();
	}

	/**
	 * Handles complete mouse dragging events in the draw area, i.e. the mouse button has been released after dragging.
	 * 
	 * @param evt
	 */
	private void graphpanelMouseReleased(java.awt.event.MouseEvent evt) {
		
		if (mouseEventsDisabled) {
			return;
		}

		mousemove[0] = 0;
		mousemove[1] = 0;
		mouse_is_pressed = false;

		if (moveseqs) {
			moveseqs = false;
			int movex = evt.getX() - mousestart[0];
			int movey = evt.getY() - mousestart[1];
			moveselected(movex, movey);

		} else {
			if (inSelectionMode()) {
				int[] tmpreg = new int[4];
				tmpreg[0] = selectstart[0];
				tmpreg[1] = selectstart[1];
				tmpreg[2] = evt.getX() - draw_area.xTranslate;
				tmpreg[3] = evt.getY() - draw_area.yTranslate;
				if (evt.isAltDown() || evt.isControlDown() || evt.isShiftDown()) {
					// modifier keys indicate deselection instead of selection
					updateSelectedSequences(tmpreg, true);
				} else {
					updateSelectedSequences(tmpreg, false);
				}
				clusterconf = null;
			}

			// deepcopy the myrotmtx to rotmtx (do NOT assign reference)
			data.rotmtx[0][0] = data.myrotmtx[0][0];
			data.rotmtx[0][1] = data.myrotmtx[0][1];
			data.rotmtx[0][2] = data.myrotmtx[0][2];
			data.rotmtx[1][0] = data.myrotmtx[1][0];
			data.rotmtx[1][1] = data.myrotmtx[1][1];
			data.rotmtx[1][2] = data.myrotmtx[1][2];
			data.rotmtx[2][0] = data.myrotmtx[2][0];
			data.rotmtx[2][1] = data.myrotmtx[2][1];
			data.rotmtx[2][2] = data.myrotmtx[2][2];
			draw_area.tmprotmtx = new double[3][3];
		}
		draw_area.repaint();
	}

	/**
	 * Shows a prompt whether the user really wants to exit and do so if approved.
	 */
	private void openReallyExitSafetyDialog() {

		Object[] options = { "Yes", "No"};
		int result = JOptionPane.showOptionDialog(this, "Really exit CLANS?\nUnsaved changes will be lost!",
				"Really exit?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null,
				options, options[1]);

		switch (result) {
		case 0:
			System.exit(0);
		default:
			return;
		}
	}

	/**
	 * Opens a window that lets users add taxonomy information to their data.
	 * <p>
	 * This feature depends on files from the NCBI taxonomy database that users must supply themselves.
	 */
	private void openTaxonomyWindow() {
		if (taxonomydialog != null) {
			taxonomydialog.setVisible(false);
			taxonomydialog.dispose();
		}
		taxonomydialog = new WindowTaxonomy(this);
		taxonomydialog.setVisible(true);
	}

	/**
	 * Based on the GUI state, updates the internal state correctly to use complex or simple attraction values.
	 */
	private void updateUseComplexAttractionValuesState() {
		if (attvalcompcheckbox.isSelected()) {
			data.attvalsimple = false;
		} else {
			data.attvalsimple = true;
		}
	}

	/**
	 * Based on the GUI state, updates the internal state correctly to rescale or not attraction values.
	 */
	private void updateRescaleAttractionValuesState() {
		if (rescalepvaluescheckbox.isSelected()) {
			data.rescalepvalues = true;
		} else {
			data.rescalepvalues = false;
		}
	}

	/**
	 * If running, switches to "optimize only selected" mode with the next iteration, otherwise sets that mode for
	 * future calculations.
	 */
	private void updateOptimizeOnlySelectedSequencesState() {
		restartComputation();
	}

	/**
	 * Saves the attraction values to a user-chosen file.
	 */
	private void saveAttractionValuesAs() {

		boolean restart_computation = stopComputation(true);

		int returnVal = fc.showSaveDialog(this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			setTitle("Clustering of " + fc.getSelectedFile().getName());
			data.save_attraction_values_to_file(fc.getSelectedFile());
		}

		if (restart_computation) {
			startComputation();
		}
	}

	public int level = 0;
	Vector<float[][]> parentmovearr = new Vector<float[][]>();
	Vector<float[][]> parentposarr = new Vector<float[][]>();
	Vector<AminoAcidSequence[]> parentaln = new Vector<AminoAcidSequence[]>();
	Vector<String[]> parentnamearr = new Vector<String[]>();
	Vector<HashMap<String, Integer>> parentnameshash = new Vector<HashMap<String, Integer>>();

	Vector<MinimalHsp[]> parentblasthits = new Vector<MinimalHsp[]>();
	Vector<MinimalAttractionValue[]> parent_orgattvals = new Vector<MinimalAttractionValue[]>();
	Vector<Vector<SequenceGroup>> parentseqgroups = new Vector<Vector<SequenceGroup>>();
	Vector<float[]> parentweights = new Vector<float[]>();

	Vector<WindowShowBlastHitsOnSequence> viewblasthitsvec = new Vector<WindowShowBlastHitsOnSequence>();
	Vector<UnusedWeirdClassPreviouslyCalledSelectclass> selectvec = new Vector<UnusedWeirdClassPreviouslyCalledSelectclass>();

	boolean mouse_is_pressed = false;
	boolean moveseqs = false;// flag to set if I want to draw sequences in 3d-space
	public boolean recalc = true;// synchronize drawing and calculating
	public DrawArea draw_area;

	/**
	 * dataLock prevents the GUI from drawing data that is currently chaning, e.g. while adjusting to new data after
	 * loading.
	 */
	final public Object dataLock = new Object();
	
	/**
	 * Prevents the selection from changing during operations that need a stable selection.
	 */
	final public Object selectedSequencesLock = new Object();

	WindowShowSelectedSequences shownames;
	WindowEditGroups myseqgroupwindow;
	
	ReporterThread reporterThread;
	final Object reporterWaitForIteratorThreadStartLock = new Object();
	
	IterationsComputerThread iterationComputerThread;
	public final Object computationLock = new Object();

	public int skiprounds = 1;
	int[] mousemove = new int[2];
	int[] translate = new int[2];
	int[] mousestart = new int[2];
	int[] selectstart = new int[2];
	int[] currmousepos = new int[2];

	int[] groupseqs = null;// used by seqgroupsdialog to show which group is currently selected
	java.awt.Color groupseqscolor = new java.awt.Color(0, 0, 0);
	int[] blastselectseqs = new int[0];
	float[] clusterconf = null;// confidence values for the selected sequences

	static final JFileChooser fc = new JFileChooser(new File("."));
	//String repaint = null;

	public WindowOptions options_window = null;
	WindowRotatationAnimantionSetup myrotationdialog = null;// clans rotation values
	WindowAffyMicroarrayData myaffydialog = null;// loads/shows affymetrix data
	WindowMicroarrayFunctionMapping mymapfunctiondialog = null;// loads/shows metabolic/functional mapping
	WindowTaxonomy taxonomydialog = null;// show the NCBI taxonomic keywords
	String namesdmp_file = null;// "names.dmp";
	String nodesdmp_file = null;// "nodes.dmp";
	ArrayList<java.io.File> mapfiles = new ArrayList<java.io.File>();
	ArrayList<java.io.File> lookupfiles = new ArrayList<java.io.File>();
	Vector<Replicates> affyfiles = null;
	public boolean usefoldchange = false;
	public boolean avgfoldchange = false;
	boolean dotsfirst = false;

	ClusterData data = null;

	private javax.swing.JMenu menu_misc;
	private javax.swing.JMenuItem aboutmenuitem;
	private javax.swing.JMenuItem addseqsmenuitem;
	private javax.swing.JMenuItem affymenuitem;
	javax.swing.JCheckBoxMenuItem antialiasingcheckboxmenuitem;
	private javax.swing.JCheckBoxMenuItem attvalcompcheckbox;
	private javax.swing.JMenuItem centermenuitem;
	private javax.swing.JMenuItem changebgcolormenuitem;
	private javax.swing.JMenuItem changeblastcolor;
	private javax.swing.JMenuItem changecolormenuitem;
	private javax.swing.JMenuItem changefgcolormenuitem;
	private javax.swing.JMenuItem changefontmenuitem;
	private javax.swing.JMenuItem changenumbercolor;
	private javax.swing.JMenuItem changeselectcolormenuitem;
	public javax.swing.JButton selectionButton;
	javax.swing.JCheckBoxMenuItem cluster2dbutton;
	private javax.swing.JMenuItem clustermenuitem;
	javax.swing.JCheckBoxMenuItem showConnectionsFrustrationColoredCheckbox;
	private javax.swing.JPanel drawbuttonpanel;
	private javax.swing.JMenu menu_draw;
	private javax.swing.JMenuItem evalueitem;
	private javax.swing.JMenu menu_file;
	private javax.swing.JMenuItem getblasthitsmenuitem;
	private javax.swing.JMenuItem getchildmenuitem;
	private javax.swing.JMenuItem getdotsizemenuitem;
	private javax.swing.JMenuItem getovalsizemenuitem;
	private javax.swing.JMenuItem getparentmenuitem;
	private javax.swing.JMenuItem getseqsforselectedhits;
	private javax.swing.JMenuItem getseqsmenuitem;
	public GuiMessageOverlay messageOverlay;
	private Component originalGlassPane;
	public boolean messageOverlayActive;
	//private javax.swing.JPanel graphpanelxx;
	private javax.swing.JMenu menu_help;
	private javax.swing.JMenuItem helpmenuitem;
	private javax.swing.JMenuItem remove_selected_sequences_menu_item;
	private javax.swing.JMenuItem hidesingletonsmenuitem;
	private javax.swing.JMenuBar jMenuBar1;
	javax.swing.JCheckBoxMenuItem colorDotsBySequenceLengthCheckbox;
	private javax.swing.JMenuItem loadalternatemenuitem;
	private javax.swing.JMenuItem loadgroupsmenuitem;
	private javax.swing.JMenuItem loadmenuitem;
	private javax.swing.JMenuItem loadtabsmenuitem;
	private javax.swing.JMenuItem mapmanmenuitem;
	private javax.swing.JButton setThresholdButton;
	private javax.swing.JTextField textfield_threshold_value;
	private javax.swing.JTextField textfield_info_min_blast_evalue;
	private javax.swing.JTextField roundsCompletedTextfield;
	

	/**
	 * The selection status of this item is forwarded to the computation threads to let them optimize all sequences or
	 * just a subset.
	 */
	private javax.swing.JCheckBoxMenuItem moveOnlySelectedSequencesCheckbox;
	
	private javax.swing.JMenuItem printmenuitem;
	private javax.swing.JCheckBoxMenuItem rescalepvaluescheckbox;
	private javax.swing.JCheckBoxMenuItem messageOverlayActiveCheckbox;
	private javax.swing.JMenuItem rotationmenuitem;
	private javax.swing.JMenuItem save2dmenuitem;
	private javax.swing.JMenuItem saveattvalsmenuitem;
	private javax.swing.JMenuItem savemenuitem;
	private javax.swing.JMenuItem saveasmenuitem;
	private javax.swing.JMenuItem savemtxmenuitem;
	/**
	 * If selected, we are in SELECT mode, if not selected in MOVE mode.
	 */
	javax.swing.JToggleButton selectionOrMovementModeButton;
	private javax.swing.JMenuItem seqscoloring;
	private javax.swing.JMenuItem sequencesitem;
	private javax.swing.JMenuItem setrotmenuitem;
	javax.swing.JCheckBoxMenuItem showBlastHitNumberCheckbox;
	javax.swing.JCheckBox showConnectionsCheckbox;
	javax.swing.JCheckBoxMenuItem showInfoOnDrawAreaCheckbox;
	javax.swing.JCheckBox showSequenceNamesCheckbox;
	private javax.swing.JCheckBoxMenuItem shownamesselectcheckbox;
	javax.swing.JCheckBox showSequenceNumbersCheckbox;
	private javax.swing.JMenuItem showoptionsmenuitem;
	javax.swing.JCheckBoxMenuItem showorigcheckbox;
	private javax.swing.JButton showSelectedButton;
	private javax.swing.JMenuItem showseqsmenuitem;
	private final String autosaveSetupLabel = "setup autosaving";
	private javax.swing.JMenuItem autosaveSetup;
	private javax.swing.JMenuItem skipdrawingrounds;
	private javax.swing.JButton initializeButton;
	private javax.swing.JMenuItem stereoanglemenuitem;
	javax.swing.JCheckBoxMenuItem stereocheckboxmenuitem;
	
	private javax.swing.JButton startStopResumeButton;
	private boolean reenableButtonAfterStoppingComputation;
	
	private javax.swing.JMenuItem taxonomymenuitem;
	private javax.swing.JMenu menu_windows;
	
	private javax.swing.JButton zoomOnSelectedButton;
	private boolean zoomOnSelected = false;
	
	private javax.swing.JButton resetZoomButton;

	
	private javax.swing.JMenuItem zoommenuitem;
	
	/**
	 * Based on the GUI state, updates the start/stop/resume button label.
	 */
	public void updateStartStopResumeButtonLabel() {
		
		if (data == null) {
			modifyButtonStartStopResume("no data", false);
			return;
		}

		Object message = null;

		if (data.hasRoundsLimit()) {

			if (isComputing()) {
				message = String
						.format("STOP (%" + String.valueOf(data.getRoundsLimit()).length() + "s/"
								+ data.getRoundsLimit() + ")", data.roundsCompleted);

			} else { // reached previous round limit or just set one up
				message = ButtonStartStopMessage.START + " for " + data.getRoundsLimit() + " rounds";
			}

		} else { // without round limit
			if (isComputing()) {
				message = ButtonStartStopMessage.STOP;

			} else if (data.rounds == 0) {
				message = ButtonStartStopMessage.START;

			} else {
				message = ButtonStartStopMessage.RESUME;
			}
		}

		modifyButtonStartStopResume(message, true);
	}
	
	/**
	 * Changes label and availability state of the start/resume/stop button.
	 * 
	 * @param label
	 *            The new label or null if the old label should be kept.
	 * @param enabled
	 *            True if the button should be enabled, false if disabled.
	 */
	public void modifyButtonStartStopResume(Object label, boolean enabled) {
		if (label != null) {
			if (label instanceof ButtonStartStopMessage) {
				startStopResumeButton.setText(((ButtonStartStopMessage) label).get());
			} else {
				startStopResumeButton.setText(label.toString());
			}
		}

		startStopResumeButton.setEnabled(enabled);
	}
	
	/**
	 * Enables the start/stop/resume button.
	 */
	private void enableButtonStartStopResume() {
		startStopResumeButton.setEnabled(true);
	}

	/**
	 * Disables the start/stop/resume button.
	 */
	private void disableButtonStartStopResume() {
		startStopResumeButton.setEnabled(false);
	}

	/**
	 * Checks if we are currently computing new iterations of the clustering.
	 * 
	 * @return true iff the clustering is running
	 */
	protected boolean isComputing() {
		return reporterThread != null && reporterThread.hasRunningIterationsComputerThread();
	}
	
	/**
	 * Starts the computation.
	 * 
	 * @return true iff the computation was already running before the call.
	 */
	private boolean startComputation() {
		if (isComputing()) {
			return true;
		}

		if (!containsData(true)) {
			return false;
		}
		
		reenableButtonAfterStoppingComputation = true;

		updateOptionValuesFromOptionsWindow();

		initializeComputationThread(); // thread cannot be restarted
		
		/**
		 * we need to wait for the iterator thread to start, otherwise the start/stop/resume button will not get a
		 * "STOP" label, as the threads status is checked there.
		 */
		synchronized (reporterWaitForIteratorThreadStartLock) {
			reporterThread.start();
			
			while(!reporterThread.hasRunningIterationsComputerThread()) {
				
				try {
					reporterWaitForIteratorThreadStartLock.wait();
				
				} catch (InterruptedException e) {
					System.err.println("fatal error: GUI was interrupted while waiting for iterator thread!");
					e.printStackTrace();
					System.exit(11);
				}
			}
		}

		updateStartStopResumeButtonLabel();

		return false;
	}

	/**
	 * Stops the computation.
	 * 
	 * @param keep_button_disabled
	 *            if true, the start/stop/resume button stays disabled after returning, leaving it to the following code
	 *            to re-enable it. Useful if the button needs to stay disabled during an operation.
	 * @return true if the computation was running.
	 */
	private boolean stopComputation(boolean keep_button_disabled) {
		if (!isComputing()) {
			return false;
		}
		
		reenableButtonAfterStoppingComputation = !keep_button_disabled;

		disableButtonStartStopResume();

		reporterThread.stopChild();
		reporterThread.joinChild();

		return true;
	}
	
	private void restartComputation() {
		if (stopComputation(false)) {
			startComputation();
		}
	}
	
	/**
	 * Inform the GUI that the iteration computations have stopped and housekeeping can start. This method is called by
	 * the {@code reporterThread} once a started {@code iterationComputerThread} has died.
	 * <p>
	 * Note: Further housekeeping after stopping can be done here easily.
	 */
	public void computationHasStopped() {
		// WARNING: if this method is ever renamed, its name also has to be changed in initializeComputationThread 
		if (reenableButtonAfterStoppingComputation) {
			updateStartStopResumeButtonLabel();
		}
	}

	/**
	 * Returns whether all or just the selected subset of sequences should be optimized during calculations.
	 * 
	 * @return true if only selected sequences should be optimized, false if all should be optimized.
	 */
	public boolean optimizeOnlySelectedSequences() {
		return moveOnlySelectedSequencesCheckbox.isSelected();
	}
	
	/**
	 * Resets the graph to random sequence positions. This is used when and equal to just having loaded data off BLAST
	 * results.
	 */
	void resetGraph() {
		data.changedvals = false;

		updateOptionValuesFromOptionsWindow();

		if (!isComputing()) {
			draw_area.repaint();
		}

		data.initialize();
		updateRoundsCompletedTextfield();

		if (options_window != null) {
			options_window.currcoolfield.setText(String.valueOf(data.currcool));
		}

		updateStartStopResumeButtonLabel();

		mousemove[0] = 0;
		mousemove[1] = 0;

		resetZoom();
		centerGraph();
	}

	/**
	 * Starts a new thread that loads a CLANS file in the background to keep the GUI responsive. Once the loading is
	 * done, the GUI is briefly blocked from redrawing (read: it is frozen) to inject the new values into it.
	 * 
	 * @param filename
	 *            The CLANS file to be loaded.
	 */
	private void loadDataClansFormatInBackground(final String filename) {

		if (messageOverlayActive) {
			messageOverlay.setLoading();
		} else {
			System.out.println("Note: loading can be canceled by pressing escape");
		}
		
		final boolean restart_computation = stopComputation(true);
		final boolean restart_autosave = pauseAutosave();
		
		disableUserControls();
	
		if (saveLoadWorker != null) {
			saveLoadWorker.cancel(true); // kill potentially running old worker
		}
		
		final ClusterDataLoadHelper[] new_data = new ClusterDataLoadHelper[1];
		final long start_time = System.currentTimeMillis();
		
		saveLoadWorker = new SwingWorker<Void, Integer>() {

			@Override
			protected Void doInBackground() throws Exception {
				Thread.currentThread().setName("load in background");
				new_data[0] = ClusterData.load_run_from_file(new java.io.File(filename), saveLoadWorker);
				return null;
			}

			@Override
			protected void done() {
				String error_message = "";
				boolean loading_successful = true;
				boolean loading_canceled = false;
				
				try {
					// we query the output of doInBackground only to get to the catch blocks in case of failures
					get();

				} catch (CancellationException e) {
					loading_successful = false;
					loading_canceled = true;
					error_message = "loading canceled";

				} catch (InterruptedException e) {
					loading_successful = false;
					error_message = e.getMessage();

				} catch (ExecutionException e) {
					loading_successful = false;

					Throwable exception = e.getCause();
					if (exception instanceof FileNotFoundException) {
						error_message = exception.getMessage();

					} else if (exception instanceof ParseException) {
						error_message = "line " + ((ParseException) exception).getErrorOffset() + ": "
								+ exception.getMessage();

					} else if (exception instanceof IOException) {
						error_message = exception.getMessage();

					} else {
						error_message = "unexpected exception: " + e.getMessage();
					}
				}

				// completion_time cannot be recorded here as loading_successful triggers further actions
				long completion_time;
				
				if (loading_successful) {

					synchronized (dataLock) {
						data.injectLoadedDataIntoExisting(filename, new_data[0]);
					}
					injectLoadedSetupIntoGui();
					
					// completion_time is slightly off as we must capture it before initializeAutosave to avoid wait
					// time during user interaction with a dialog opened in that method
					completion_time = System.currentTimeMillis();
					
					initializeAutosave(true);

					if (messageOverlayActive) {
						messageOverlay.setCompleted("autsaving every " + data.getAutosaveIntervalMinutes()
								+ " minutes");
					}

					if (restart_computation) {
						startComputation();
					}

				} else {
					completion_time = System.currentTimeMillis();
					
					System.err.println(error_message);

					if (messageOverlayActive) {
						if (loading_canceled) {
							messageOverlay.setCanceled();
						} else {
							messageOverlay.setFailed(error_message);
						}
					}
					
					if (restart_autosave) {
						restartAutosave();
					}
				}

				System.out.println("loading took " + String.format("%s", (completion_time - start_time) / 1000f)
						+ " seconds.\n");

				enableUserControls();
			}
		};

		saveLoadWorker.execute();
	}

	/**
	 * Loads clans run from file.
	 * 
	 * @param filename
	 *            The input file.
	 * @return true iff loading succeeded.
	 * @throws IOException
	 * @throws ParseException
	 * @throws FileNotFoundException
	 */
	void injectLoadedSetupIntoGui() {

		textfield_info_min_blast_evalue.setText(String.valueOf(data.maxvalfound));
		if (data.blasthits == null) {
			setThresholdButton.setText("Use Attraction values better than");
			textfield_threshold_value.setText("0");
			savemtxmenuitem.setText("Save Attraction values as matrix");
			evalueitem.setText("Attraction value plot");
		} else {
			if (data.usescval) {
				setThresholdButton.setText("Use SC-vals better than");
				textfield_threshold_value.setText("0");
				evalueitem.setText("SC-value plot");
			} else {
				setThresholdButton.setText("Use P-values better than");
				textfield_threshold_value.setText("1");
				evalueitem.setText("P-value plot");
			}
		}

		if (data.complexatt == true) {
			attvalcompcheckbox.setSelected(true);
		} else {
			attvalcompcheckbox.setSelected(false);
		}
		if (data.cluster2d) {
			cluster2dbutton.setSelected(true);
		} else {
			cluster2dbutton.setSelected(false);
		}
		if (data.showinfo) {
			showInfoOnDrawAreaCheckbox.setSelected(true);
		} else {
			showInfoOnDrawAreaCheckbox.setSelected(false);
		}
		textfield_threshold_value.setText(String.valueOf(data.pvalue_threshold));

		if (options_window != null) {
			options_window.refreshValuesFromModel();
		}

		textfield_info_min_blast_evalue.setText(String.valueOf(data.maxvalfound));
		updateRoundsCompletedTextfield();
		
		setTitle("Clustering of " + data.getBaseInputfileName());

		resetZoom();
	}

	boolean contains_data() {
		return data.movements != null;
	}

	boolean containsData(boolean showinfo) {
		if (contains_data()) {
			return true;
		}

		if (showinfo) {
			javax.swing.JOptionPane.showMessageDialog(this, "you have to load data first", "no data loaded",
					javax.swing.JOptionPane.ERROR_MESSAGE);
		}
		return false;
	}

	public void initaddedseqs(MinimalHsp[] blastvec, AminoAcidSequence[] allaln, String[] allnamearr,
			HashMap<String, Integer> allnameshash, int[] newnumarr, float[][] allposarr, float maxmove, double pval,
			boolean useselectedonly) {
		// initialize the necessary variable for the case where new sequences are added to an existing run.
		data.sequences = allaln;
		data.positions = allposarr;
		data.blasthits = blastvec;
		data.maxmove = maxmove;
		data.pvalue_threshold = pval;
		textfield_threshold_value.setText(String.valueOf(data.pvalue_threshold));
		textfield_info_min_blast_evalue.setText(String.valueOf(data.pvalue_threshold));
		data.selectedSequencesIndices = newnumarr;
		data.nameshash = allnameshash;
		data.sequence_names = allnamearr;
		data.elements = data.sequence_names.length;
		data.posarr = data.positions;
		data.movementsLastIteration = new float[data.elements][ClusterData.dimensions];
		data.movements = new float[data.elements][ClusterData.dimensions];
		data.posarrtmp = new float[data.elements][ClusterData.dimensions];
		data.drawarrtmp = new int[data.elements][ClusterData.dimensions];
		data.resetDrawOrder();

		data.compute_attraction_values();

		moveOnlySelectedSequencesCheckbox.setSelected(useselectedonly);
		int seqnum = data.sequences.length;
		System.out.println("seqnum=" + seqnum);
		data.seqlengths = new float[seqnum];
		float maxlength = 0;
		for (int i = 0; i < seqnum; i++) {
			data.seqlengths[i] = data.sequences[i].seq.replaceAll("-", "").length();
			if (data.seqlengths[i] > maxlength) {
				maxlength = data.seqlengths[i];
			}
		}
		for (int i = 0; i < seqnum; i++) {
			data.seqlengths[i] /= maxlength;
		}

	}

	// --------------------------------------------------------------------------

	AminoAcidSequence[] getselectedseqs() {
		// get the currently selected sequences
		ArrayList<AminoAcidSequence> tmpvec = new ArrayList<AminoAcidSequence>();
		int sequences = data.selectedSequencesIndices.length;
		AminoAcidSequence curraaseq;
		for (int i = 0; i < sequences; i++) {
			curraaseq = new AminoAcidSequence();
			curraaseq.name = data.sequence_names[data.selectedSequencesIndices[i]];
			curraaseq.seq = data.sequences[data.selectedSequencesIndices[i]].seq;
			tmpvec.add(curraaseq);
		}
		AminoAcidSequence[] retarr = (AminoAcidSequence[]) tmpvec.toArray(new AminoAcidSequence[0]);
		return retarr;
	}

	/**
	 * Updates the data with the values of all available fields of the options window.
	 */
	public void updateOptionValuesFromOptionsWindow() {
		if (options_window == null) {
			return;
		}

		String tmpstr = "";
		boolean refresh_values_from_model = false;
		try {
			tmpstr = options_window.attfield.getText();
			data.attfactor = Float.parseFloat(tmpstr);

			tmpstr = options_window.repfield.getText();
			data.repfactor = Float.parseFloat(tmpstr);

			tmpstr = options_window.dampfield.getText();
			data.dampening = Float.parseFloat(tmpstr);

			tmpstr = options_window.coolfield.getText();
			data.cooling = Double.parseDouble(tmpstr);

			tmpstr = options_window.currcoolfield.getText();
			data.currcool = Double.parseDouble(tmpstr);

			tmpstr = options_window.minattfield.getText();
			data.minattract = Double.parseDouble(tmpstr);

			tmpstr = options_window.maxmovefield.getText();
			data.maxmove = Float.parseFloat(tmpstr);

			tmpstr = options_window.attvalpowtextfield.getText();
			data.attvalpow = Integer.parseInt(tmpstr);

			tmpstr = options_window.repvalpowtextfield.getText();
			data.repvalpow = Integer.parseInt(tmpstr);

			tmpstr = options_window.roundstextfield.getText();
			int new_rounds_limit = Integer.parseInt(tmpstr);
			refresh_values_from_model = setRoundsLimit(new_rounds_limit);
					
		} catch (NumberFormatException e) {
			javax.swing.JOptionPane.showMessageDialog(this, "ERROR, unable to parse number from '" + tmpstr + "'");
		}

		if (refresh_values_from_model) {
			options_window.refreshValuesFromModel();
		}
	}

	/**
	 * Updates the selection given a selected area of the window. *
	 * 
	 * @param selected_region_corners
	 *            the selected region as {min X, max X, min Y, max Y}
	 * @param deselect
	 *            if true, selected points in the region will be deselected, otherwise the selection of all points will
	 *            be toggled
	 */
	void updateSelectedSequences(int[] selected_region_corners, boolean deselect) {

		ArrayList<Integer> tmpselect = new ArrayList<Integer>();

		int minx = selected_region_corners[0];
		int maxx = selected_region_corners[2];
		
		int miny = selected_region_corners[1];
		int maxy = selected_region_corners[3];
		
		int tmpint;
		if (minx > maxx) {
			tmpint = maxx;
			maxx = minx;
			minx = tmpint;
		}
		if (miny > maxy) {
			tmpint = maxy;
			maxy = miny;
			miny = tmpint;
		}
		
		// now get all sequences that are within those 4 coords
		float[][] currpos = data.posarrtmp;
		for (int i = 0; i < data.elements; i++) {
			// see if this object is within the selected region
			if ((currpos[i][0] > minx) && (currpos[i][0] < maxx) && (currpos[i][1] > miny) && (currpos[i][1] < maxy)) {
				tmpselect.add(new Integer(i));
			}
		}
		
		// now select those that are not present in selectednames and deselect those that are
		int tmparrsize = tmpselect.size();
		Integer[] tmparr = (Integer[]) tmpselect.toArray(new Integer[0]);
		tmpselect.clear();
		boolean isnew;
		
		synchronized (selectedSequencesLock) {

			for (int i = 0; i < tmparrsize; i++) {
				tmpint = tmparr[i].intValue();
				isnew = true;

				for (int j = 0; j < data.selectedSequencesIndices.length; j++) {
					if (data.selectedSequencesIndices[j] == tmpint) { // if this element was already selected
						isnew = false;// do not add it to the new Vector
						data.selectedSequencesIndices[j] = -1;// mark this as unselected
					}
				}

				if (isnew && (deselect == false)) {// if this element should be selected
					tmpselect.add(tmparr[i]);
				}
			}

			for (int i = 0; i < data.selectedSequencesIndices.length; i++) {
				// now add those that were selected before
				if (data.selectedSequencesIndices[i] > -1) {
					tmpselect.add(new Integer(data.selectedSequencesIndices[i]));
				}
			}

			// now assign the new selecteds to selectelements and update shownamedialog
			setSelectedSequences(tmpselect);
		}

		if (shownames != null) {
			shownames.setselected(data.selectedSequencesIndices, !shownames.showall);
		}
	}

	/**
	 * Informs the relevant GUI elements about a change in the selected sequences.
	 */
	private void informAboutSequenceSelectionChange() {
		updateSelectionButtonLabel();
		updateZoomOnSelectedButtonLabel();

		if (myseqgroupwindow != null) { // update the highlighting of groups with selected sequences
			myseqgroupwindow.repaint();
		}
	}

	/**
	 * Informs the relevant GUI elements about a completed clustering iteration.
	 */
	public void informAboutCompletedIteration() {
		updateRoundsCompletedTextfield();
	}

	/**
	 * Informs the relevant GUI elements about a zoom level/type change.
	 */
	private void informAboutZoomChange() {
		updateZoomOnSelectedButtonLabel();
		updateResetZoomButtonLabel();
	}
	
	/**
	 * Updates the display of completed rounds from the model.
	 */
	private void updateRoundsCompletedTextfield() {
		if (containsData(false)) {
			roundsCompletedTextfield.setText("Round: " + Integer.toString(data.rounds));

		} else {
			roundsCompletedTextfield.setText("Round: 0");
		}
	}
	
	/**
	 * Sets the selected entries and does some housekeeping for that.
	 * 
	 * @param new_selecteds
	 *            The indices of the selected entries.
	 */
	protected void setSelectedSequences(ArrayList<Integer> new_selecteds) {
		
		synchronized (selectedSequencesLock) {
			data.selectedSequencesIndices = new int[new_selecteds.size()];
			for (int i = 0; i < new_selecteds.size(); i++) {
				data.selectedSequencesIndices[i] = new_selecteds.get(i).intValue();
			}
		}
		
		informAboutSequenceSelectionChange();

		// to get to an overview from zoom on selected mode when 
		if (zoomOnSelected && data.selectedSequencesIndices.length < 2) {
			resetZoom();
		}
	}

	protected void selectAllSequences() {
		ArrayList<Integer> all_indices = new ArrayList<Integer>();
		for (int i = 0; i < data.sequences.length; i++) {
			all_indices.add(i);
		}

		setSelectedSequences(all_indices);
	}
	
	/**
	 * Deselects everything.
	 */
	protected void deselectAllSequences() {
		setSelectedSequences(new ArrayList<Integer>());
	}

	void updatetmpselected(int[] tmpreg) {
		// as input this has 4 coordinates
		// what this does is to take all sequences that are within those coordinates
		// and select/deselect them. It also needs to update the selection in shownamedialog
		ArrayList<Integer> tmpselect = new ArrayList<Integer>();
		int minx = tmpreg[0];
		int maxx = tmpreg[2];
		int miny = tmpreg[1];
		int maxy = tmpreg[3];
		int tmpint;
		if (minx > maxx) {
			tmpint = maxx;
			maxx = minx;
			minx = tmpint;
		}
		if (miny > maxy) {
			tmpint = maxy;
			maxy = miny;
			miny = tmpint;
		}
		// now get all sequences that are within those 4 coords
		float[][] currpos = data.posarrtmp;
		for (int i = 0; i < data.elements; i++) {
			// see if this object is within the selected region
			if ((currpos[i][0] > minx) && (currpos[i][0] < maxx) && (currpos[i][1] > miny) && (currpos[i][1] < maxy)) {
				tmpselect.add(new Integer(i));
			}
		}
			// now select those that are not present in selectednames and deselect those that are
		int j;
		int selectedelements = data.selectedSequencesIndices.length;
		int[] tmpdeselected = new int[selectedelements];
		boolean isnew;
		for (int i = 0; i < tmpselect.size(); i++) {
			tmpint = ((Integer) tmpselect.get(i)).intValue();
			isnew = true;
			for (j = 0; j < selectedelements; j++) {
				if (data.selectedSequencesIndices[j] == tmpint) {
					isnew = false;
					tmpdeselected[j] = 1;
					break;
				}
			}
			if (isnew == false) {
				tmpselect.remove(i);
				i--;
			}
		}
			// now add the already selected names
		for (int i = selectedelements - 1; i >= 0; i--) {
			// now add those that were selected before
			if (tmpdeselected[i] != 1) {
				tmpselect.add(0, new Integer(data.selectedSequencesIndices[i]));
			}
		}
			// now assign the new selecteds to tmpselectelements and update shownamedialog
		int[] tmpselectednames = new int[tmpselect.size()];
		for (int i = tmpselect.size() - 1; i >= 0; i--) {
			tmpselectednames[i] = ((Integer) tmpselect.get(i)).intValue();
		}
		if (shownames != null) {
			shownames.setselected(tmpselectednames, true);
		}
	}

	void moveselected(int xint, int yint) {
		// move selected sequences by movex and movy in 3d space
		// use the rotmtx to compute the x,y and z coord to move them by
		float x = (float) ((float) xint / draw_area.xscale);
		float y = (float) ((float) yint / draw_area.yscale);
		double[][] invrot = getinvrot(data.rotmtx);
		float[] move = new float[3];
		move[0] = (float) (invrot[0][0] * x + invrot[0][1] * y);
		move[1] = (float) (invrot[1][0] * x + invrot[1][1] * y);
		move[2] = (float) (invrot[2][0] * x + invrot[2][1] * y);
		// now move each of the selected seqs by that amount
		int selectedelements = data.selectedSequencesIndices.length;
		// for each of the selected sequences shift it's position in posarr by x,y,z coords.
		// selectednames is a int[] with the selected sequences indices
		for (int i = 0; i < selectedelements; i++) {
			data.positions[data.selectedSequencesIndices[i]][0] += move[0];
			data.positions[data.selectedSequencesIndices[i]][1] += move[1];
			data.positions[data.selectedSequencesIndices[i]][2] += move[2];
		}
	}

	static double[][] getinvrot(double[][] m) {
		// get the inverse of a 3x3 matrix inmtx
		double[][] retmtx = {
				{ ((m[1][1] * m[2][2]) - (m[2][1] * m[1][2])), -((m[0][1] * m[2][2]) - (m[2][1] * m[0][2])),
						((m[0][1] * m[1][2]) - (m[1][1] * m[0][2])) },
				{ -((m[1][0] * m[2][2]) - (m[2][0] * m[1][2])), ((m[0][0] * m[2][2]) - (m[2][0] * m[0][2])),
						-((m[0][0] * m[1][2]) - (m[1][0] * m[0][2])) },
				{ ((m[1][0] * m[2][1]) - (m[2][0] * m[1][1])), -((m[0][0] * m[2][1]) - (m[2][0] * m[0][1])),
						((m[0][0] * m[1][1]) - (m[1][0] * m[0][1])) } };
		return retmtx;
	}
}
