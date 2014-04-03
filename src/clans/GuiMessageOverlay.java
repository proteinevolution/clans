package clans;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.Timer;
import javax.swing.UIManager;

public class GuiMessageOverlay {

	/**
	 * Enum of the different states of the message overlay, e.g. OFF, LOADING, SAVING_COMPLETED.
	 */
	class States {
		private static final int OFF = 0;

		private static final int LOADING = 11;
		private static final int LOADING_COMPLETED = 12;
		private static final int LOADING_CANCELED = 13;
		private static final int LOADING_FAILED = 14;

		private static final int SAVING = 21;
		private static final int SAVING_COMPLETED = 22;
		private static final int SAVING_CANCELED = 23;
		private static final int SAVING_FAILED = 24;
	}

	/**
	 * Enum of the messages for different {@code States}.
	 */
	class Messages {
		private static final String LOADING = "loading file";
		private static final String LOADING_COMPLETED = "loaded file successfully";
		private static final String LOADING_CANCELED = "canceled loading file";
		private static final String LOADING_FAILED = "LOADING FILE FAILED!";

		private static final String SAVING = "saving file";
		private static final String SAVING_COMPLETED = "saved file successfully";
		private static final String SAVING_CANCELED = "canceled saving file";
		private static final String SAVING_FAILED = "SAVING FILE FAILED!";
	}

	/**
	 * Enum of the durations for different message types, e.g. INFO or ERROR.
	 */
	class Durations {
		private static final int INFO = 3000; // in milliseconds
		private static final int ERROR = 3000; // in milliseconds
		private static final int INFINITE = Integer.MAX_VALUE;
	}

	private JPanel parent; // the parent is used to adjust the glass when the parent is resized.

	private JPanel glass; // the transparent container for the labels
	private Component top_spacer; // moves the labels down
	private JLabel mainLabel; // shows the main message
	private JLabel detailsLabel; // shows messageDetails

	private String messageDetails; // the details of the message, e.g. what went wrong during loading

	private Timer timer; // draw only every X milliseconds instead of using up one CPU for drawing continuously.

	private int mode; // the current mode
	private int duration; // the duration for the current mode
	private long completionTime = 0; // time of last completed operation to trigger message display

	private Color colorWorkInProgress; // color for work-in-progress operations
	private Color colorCompleted; // color for completed operations
	private Color colorFailed; // color for failed operations

	public GuiMessageOverlay(JPanel glass, JPanel parent) {

		this.glass = glass;
		this.parent = parent;

		glass.setLayout(new BoxLayout(glass, BoxLayout.PAGE_AXIS));

		mainLabel = new JLabel();
		mainLabel.setFont(UIManager.getFont("Label.font").deriveFont(Font.BOLD).deriveFont(20f));
		mainLabel.setAlignmentX(Component.CENTER_ALIGNMENT);

		detailsLabel = new JLabel();
		detailsLabel.setFont(UIManager.getFont("Label.font").deriveFont(Font.BOLD).deriveFont(16f));
		detailsLabel.setAlignmentX(Component.CENTER_ALIGNMENT);

		// the initial values for top_spacer ARE NOT USED as we overwrite them in updateGlassSize() below
		top_spacer = Box.createRigidArea(new Dimension(0, 300));
		glass.add(top_spacer);
		glass.add(mainLabel);
		glass.add(Box.createRigidArea(new Dimension(0, 3)));
		glass.add(detailsLabel);
		updateGlassSize();

		colorWorkInProgress = new Color(255, 128, 0, 0);
		colorCompleted = new Color(0, 220, 0, 0);
		colorFailed = new Color(220, 0, 0, 0);

		resetMessageDetails();

		int timer_delay = 50; // in milliseconds
		setupTimer(timer_delay);

		setOff();
	}

	/**
	 * Sets up Timer to only redraw every X milliseconds instead of continuously, which would use up 1 CPU entirely just
	 * for the overlay.
	 */
	private void setupTimer(int delay) {
		timer = new Timer(delay, new ActionListener() {
			public void actionPerformed(ActionEvent ae) {
				invoke();
			}
		});
		timer.setInitialDelay(0);
	}

	/**
	 * Resets detailed message to none.
	 */
	private void resetMessageDetails() {
		messageDetails = "";
	}

	/**
	 * Switches the overlay off and resets its state. Used once in constructor.
	 */
	private void setOff() {
		if (timer != null) {
			timer.stop();
		}

		mode = States.OFF;
		duration = Durations.INFINITE;

		resetMessageDetails();

		glass.setVisible(false);
	}

	/**
	 * Changes the state to "loading has completed".
	 */
	private void setLoadingCompleted() {
		mode = States.LOADING_COMPLETED;
		duration = Durations.INFO;
		resetMessageDetails();
	}

	/**
	 * Changes the state to "loading was canceled".
	 */
	private void setLoadingCanceled() {
		mode = States.LOADING_CANCELED;
		duration = Durations.INFO;
		resetMessageDetails();
	}

	/**
	 * Changes the state to "loading has failed".
	 */
	private void setLoadingFailed() {
		mode = States.LOADING_FAILED;
		duration = Durations.ERROR;
	}

	/**
	 * Changes the state to "saving has completed".
	 */
	private void setSavingCompleted() {
		mode = States.SAVING_COMPLETED;
		duration = Durations.INFO;
		resetMessageDetails();
	}

	/**
	 * Changes the state to "saving was canceled".
	 */
	private void setSavingCanceled() {
		mode = States.SAVING_CANCELED;
		duration = Durations.INFO;
		resetMessageDetails();
	}

	/**
	 * Changes the state to "saving has failed".
	 */
	private void setSavingFailed() {
		mode = States.SAVING_FAILED;
		duration = Durations.ERROR;
	}

	/**
	 * Displays the message (and the message details, if any) in the given color.
	 * 
	 * @param message
	 *            Main message.
	 * @param color
	 *            Text color.
	 */
	private void showMessage(String message, Color color) {
		mainLabel.setForeground(color);
		mainLabel.setText(message);

		detailsLabel.setForeground(color);
		detailsLabel.setText(messageDetails);

		glass.setVisible(true);
		glass.repaint();
	}

	/**
	 * Displays the message (and the message details, if any) and fades the color accoring to the fraction of time of
	 * the complete {@code duration} for the message.
	 * 
	 * @param message
	 *            Main message.
	 * @param color
	 *            Text color.
	 */
	private void showFadingMessage(String message, Color color) {
		float fade_fraction = getFadingFraction();

		if (fade_fraction < 0) { // fading is over -> deactivate
			setOff();
			return;
		}

		showMessage(message, fadeColor(color, fade_fraction));
	}

	/**
	 * Computes the fraction of time already elapsed w.r.t. the complete duration or signals "duration exeeded" (see
	 * return).
	 * 
	 * @return Fraction of duration already passed as value in [0, 1] or -1 if duration exceeded
	 */
	private float getFadingFraction() {
		long time_since_start = System.currentTimeMillis() - completionTime;
		if (time_since_start > duration) {
			return -1f;
		}
		return (float) time_since_start / duration;
	}

	/**
	 * Fades a {@code base_color} by adjusting its alpha value according to {@code fade_fraction}.
	 * 
	 * @param base_color
	 *            The color that will be faded.
	 * @param fade_fraction
	 *            The percentage of fading as value in [0, 1].
	 * @return
	 */
	private Color fadeColor(Color base_color, float fade_fraction) {
		return new Color(base_color.getRed(), base_color.getGreen(), base_color.getBlue(),
				(int) (255 - 255 * fade_fraction));
	}

	/**
	 * Adjusts the label placement when the glass's parent is resized. Used once in constructor.
	 */
	protected void updateGlassSize() {
		glass.remove(0); // top_spacer is the first entry

		top_spacer = Box.createRigidArea(new Dimension(0, this.parent.getHeight() / 3 * 2));

		glass.add(top_spacer, 0);

		glass.validate();
	}

	/**
	 * Returns whether the overlay is currently visible.
	 * @return True if and only if the overlay is visible.
	 */
	protected boolean isVisible() {
		return glass.isVisible();
	}

	/**
	 * Shows the overlay with its messages if there are any.
	 */
	protected void show() {
		if (mode != States.OFF && !timer.isRunning()) {
			timer.start();
		}
	}

	/**
	 * Changes the state to "currently loading".
	 */
	protected void setLoading() {
		mode = States.LOADING;
		duration = Durations.INFINITE;
		resetMessageDetails();
		show();
	}

	/**
	 * Changes the state to "currently saving".
	 */
	protected void setSaving() {
		mode = States.SAVING;
		duration = Durations.INFINITE;
		resetMessageDetails();
		show();
	}

	/**
	 * Changes the state to "loading/saving has completed" depending on previous state.
	 */
	protected void setCompleted() {
		completionTime = System.currentTimeMillis();
		resetMessageDetails();

		switch (mode) {
		case States.LOADING:
			setLoadingCompleted();
			break;

		case States.SAVING:
			setSavingCompleted();
			break;

		default:
			setOff();
			break;
		}

		show();
	}

	/**
	 * Changes the state to "loading/saving was canceled" depending on previous state.
	 */
	protected void setCanceled() {
		completionTime = System.currentTimeMillis();
		resetMessageDetails();

		switch (mode) {
		case States.LOADING:
			setLoadingCanceled();
			break;

		case States.SAVING:
			setSavingCanceled();
			break;

		default:
			setOff();
			break;
		}

		show();
	}

	/**
	 * Changes the state to "loading/saving has failed" depending on previous state.
	 */
	protected void setFailed() {
		completionTime = System.currentTimeMillis();
		// don't call resetMessageDetails(); here this method is called by setFailed(String)

		switch (mode) {
		case States.LOADING:
			setLoadingFailed();
			break;

		case States.SAVING:
			setSavingFailed();
			break;

		default:
			setOff();
			break;
		}

		show();
	}

	/**
	 * Changes the state to "loading/saving has failed" depending on previous state and adds a more detailed message.
	 * 
	 * @param details
	 *            A message stating details, e.g. file name or reason for failure.
	 */
	protected void setFailed(String details) {
		if (details.length() > 70) { // shorten long messages
			details = details.substring(0,  65) + "[...]";
		}
		messageDetails = details;
		setFailed();
	}

	/**
	 * Checks if and which type of message is waiting to be displayed and initiates its display. Does nothing it no
	 * messages is about to be displayed.
	 */
	protected void invoke() {
		switch (mode) {

		case States.OFF:
			return;

		case States.LOADING:
			showMessage(Messages.LOADING, colorWorkInProgress);
			return;

		case States.LOADING_COMPLETED:
			showFadingMessage(Messages.LOADING_COMPLETED, colorCompleted);
			return;

		case States.LOADING_CANCELED:
			showFadingMessage(Messages.LOADING_CANCELED, colorFailed);
			return;

		case States.LOADING_FAILED:
			showFadingMessage(Messages.LOADING_FAILED, colorFailed);
			return;

		case States.SAVING:
			showMessage(Messages.SAVING, colorWorkInProgress);
			return;

		case States.SAVING_COMPLETED:
			showFadingMessage(Messages.SAVING_COMPLETED, colorCompleted);
			return;

		case States.SAVING_CANCELED:
			showFadingMessage(Messages.SAVING_CANCELED, colorFailed);
			return;

		case States.SAVING_FAILED:
			showFadingMessage(Messages.SAVING_FAILED, colorFailed);
			return;

		default:
			setOff(); // safety off-switch
			return;
		}
	}
}
