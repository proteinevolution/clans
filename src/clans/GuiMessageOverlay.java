package clans;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.Timer;
import javax.swing.UIManager;

public class GuiMessageOverlay extends JComponent {
	private static final long serialVersionUID = 818220626180736376L;

	/**
	 * States of the message overlay, e.g. OFF, LOADING, SAVING_COMPLETED.
	 */
	enum State {
		OFF,
		
		LOADING,
		LOADING_COMPLETED,
		LOADING_CANCELED,
		LOADING_FAILED,
		
		SAVING,
		SAVING_COMPLETED,
		SAVING_CANCELED,
		SAVING_FAILED;
	}

	/**
	 * Messages for different {@code State}s.
	 */
	enum Message {
		OFF("I'm an idle overlay"),
		
		LOADING("loading file, please wait"),
		LOADING_COMPLETED("loaded file successfully"),
		LOADING_CANCELED("canceled loading file"),
		LOADING_FAILED("LOADING FILE FAILED!"),

		SAVING("saving file, please wait"),
		SAVING_COMPLETED("saved file successfully"),
		SAVING_CANCELED("canceled saving file"),
		SAVING_FAILED("SAVING FILE FAILED!");

		private String value;

		private Message(String value) {
			this.value = value;
		}
		
		protected String get() {
			return value;
		}
	}

	/**
	 * Durations in milliseconds for different message types, e.g. INFO or ERROR.
	 */
	enum Duration {
		
		FADING(2000),
		INFO(2500),
		ERROR(4000),
		INFINITE(Integer.MAX_VALUE);
		
		private int value;
		
		private Duration(int value) {
			this.value = value;
		}
		
		protected int get() {
			return value;
		}
	}

	private JPanel parent; // the parent is used to adjust the glass when the parent is resized.

	private Component top_spacer; // moves the labels down
	private JLabel mainLabel; // shows the main message
	private JLabel detailsLabel; // shows messageDetails
	private Color colorLabelBackground;

	private String messageDetails; // the details of the message, e.g. what went wrong during loading

	private Timer timer; // draw only every X milliseconds instead of using up one CPU for drawing continuously.

	protected State state; // the current state
	private Duration duration; // the duration for the current mode
	private long completionTime = 0; // time of last completed operation to trigger message display

	private String dots; // the progress dots on loading and saving
	private int dotTicks;
	final private int millisecondsPerDot = 350;
	
	private float alpha;
	
	private Color colorWorkInProgress; // color for work-in-progress operations
	private Color colorCompleted; // color for completed operations
	private Color colorFailed; // color for failed operations

	public GuiMessageOverlay(JPanel parent) {
		
		this.parent = parent;

		this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		
		colorLabelBackground = new Color(128, 128, 128);
		
		mainLabel = new JLabel();
		mainLabel.setFont(UIManager.getFont("Label.font").deriveFont(Font.BOLD).deriveFont(20f));
		mainLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		mainLabel.setOpaque(true);
		mainLabel.setBackground(colorLabelBackground);
		
		detailsLabel = new JLabel();
		detailsLabel.setFont(UIManager.getFont("Label.font").deriveFont(Font.BOLD).deriveFont(16f));
		detailsLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		detailsLabel.setOpaque(true);
		detailsLabel.setBackground(colorLabelBackground);

		// the initial values for top_spacer ARE NOT USED as we overwrite them in updateGlassSize() below
		top_spacer = Box.createRigidArea(new Dimension(0, 300));
		this.add(top_spacer);
		this.add(mainLabel);
		this.add(Box.createRigidArea(new Dimension(0, 3)));
		this.add(detailsLabel);
		updateGlassSize();

		colorWorkInProgress = new Color(255, 128, 0);
		colorCompleted = new Color(128, 255, 0);
		colorFailed = new Color(255, 102, 137);

		resetMessageDetails();

		int timer_delay = 25; // in milliseconds
		setupTimer(timer_delay);
	
		setOff();
	}

	protected void setState(State state) {
		this.state = state;
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
	 * Sets the correct number of dots for recurring calls to indicate work in progress.
	 */
	private void updateDots() {
		if (state == State.LOADING || state == State.SAVING) {
			dotTicks += 1;
			
			if (dotTicks * timer.getDelay() > millisecondsPerDot) {
				if (dots.length() < 3) {
					dots += ".";
				} else {
					dots = "";
				}
			
				dotTicks = 0;
			} 
		}
	}

	/**
	 * Switches the overlay off and resets its state. Used once in constructor.
	 */
	private void setOff() {
		if (timer != null) {
			timer.stop();
		}

		alpha = 1;
		
		setState(State.OFF);
		duration = Duration.INFINITE;
		dots = "";

		resetMessageDetails();

		setVisible(false);
	}

	/**
	 * Changes the state to "loading has completed".
	 */
	private void setLoadingCompleted() {
		setState(State.LOADING_COMPLETED);
		duration = Duration.INFO;
		resetMessageDetails();
	}

	/**
	 * Changes the state to "loading was canceled".
	 */
	private void setLoadingCanceled() {
		setState(State.LOADING_CANCELED);
		duration = Duration.INFO;
		resetMessageDetails();
	}

	/**
	 * Changes the state to "loading has failed".
	 */
	private void setLoadingFailed() {
		setState(State.LOADING_FAILED);
		duration = Duration.ERROR;
	}

	/**
	 * Changes the state to "saving has completed".
	 */
	private void setSavingCompleted() {
		setState(State.SAVING_COMPLETED);
		duration = Duration.INFO;
		resetMessageDetails();
	}

	/**
	 * Changes the state to "saving was canceled".
	 */
	private void setSavingCanceled() {
		setState(State.SAVING_CANCELED);
		duration = Duration.INFO;
		resetMessageDetails();
	}

	/**
	 * Changes the state to "saving has failed".
	 */
	private void setSavingFailed() {
		setState(State.SAVING_FAILED);
		duration = Duration.ERROR;
	}

	/**
	 * Displays the message in the given color. Message details are added if available.
	 * 
	 * @param message
	 *            The main message.
	 * @param color
	 *            Text color for the messages.
	 * @param backgroundColor
	 *            Background color for the messages.
	 */
	private void showMessage(Message message, Color color, Color backgroundColor) {
		showMessage(message, color, backgroundColor, "");
	}

	/**
	 * Displays the message in the given color. Message details and progress dots are added if available.
	 * 
	 * @param message
	 *            The main message.
	 * @param color
	 *            Text color for the messages.
	 * @param backgroundColor
	 *            Background color for the messages.
	 * @param workInProgressDots
	 *            Dots to indicate work-in-progress, e.g. during loading
	 */
	protected void showMessage(Message message, Color color, Color backgroundColor, String workInProgressDots) {

		mainLabel.setForeground(color);
		mainLabel.setBackground(colorLabelBackground);
		mainLabel.setText(message.get() + workInProgressDots);

		detailsLabel.setForeground(color);
		detailsLabel.setBackground(colorLabelBackground);
		detailsLabel.setText(messageDetails);

		repaint();
		mainLabel.repaint();
		detailsLabel.repaint();
		setVisible(true);
	}

	/**
	 * Displays the message (and the message details, if any) and fades the color according to the fraction of time of
	 * the complete {@code duration} for the message.
	 * 
	 * @param message
	 *            Main message.
	 * @param color
	 *            Text color.
	 */
	private void showFadingMessage(Message message, Color color) {
		float fade_fraction = getFadingFraction();

		if (fade_fraction < 0) { // fading is over -> deactivate
			setOff();
			return;
		}
		
		alpha = 1 - fade_fraction; // fade by reducing visibility for the whole glass pane; see paintComponent below.
		showMessage(message, color, colorLabelBackground);
	}

	/**
	 * Computes the fraction of time already elapsed w.r.t. the complete duration or signals "duration exeeded" (see
	 * return).
	 * 
	 * @return Fraction of duration already passed as value in [0, 1] or -1 if duration exceeded
	 */
	private float getFadingFraction() {
		long time_since_start = System.currentTimeMillis() - completionTime;
		if (time_since_start > duration.get()) {
			return -1f;
		}
		
		float timespan_without_fading = duration.get() - Duration.FADING.get();
		
		if (time_since_start < timespan_without_fading) { // make the initial delay non-transparent
			return 0;
		}

		// fade during duration - Durations.DELAY and make the subtraction negative number proof
		float fade_fraction = (float) Math.max(time_since_start - timespan_without_fading, 0) / Duration.FADING.get();
		return fade_fraction;
	}

	/**
	 * Adjusts the label placement when the glass's parent is resized. Used once in constructor.
	 */
	protected void updateGlassSize() {
		remove(0); // top_spacer is the first entry

		top_spacer = Box.createRigidArea(new Dimension(0, parent.getHeight() / 3 * 2));

		add(top_spacer, 0);

		validate();
	}


	/**
	 * Shows the overlay with its messages if there are any.
	 */
	protected void activate_overlay() {
		if (state != State.OFF && !timer.isRunning()) {
			timer.start();
		}
	}

	/**
	 * Changes the state to "currently loading".
	 */
	protected void setLoading() {
		setState(State.LOADING);
		duration = Duration.INFINITE;
		resetMessageDetails();
		activate_overlay();
	}

	/**
	 * Changes the state to "currently saving".
	 */
	protected void setSaving() {
		setState(State.SAVING);
		duration = Duration.INFINITE;
		resetMessageDetails();
		activate_overlay();
	}

	/**
	 * Changes the state to "loading/saving has completed" depending on previous state.
	 */
	protected void setCompleted() {
		completionTime = System.currentTimeMillis();
		resetMessageDetails();

		switch (state) {
		case LOADING:
			setLoadingCompleted();
			break;

		case SAVING:
			setSavingCompleted();
			break;

		default:
			setOff();
			break;
		}
	}

	/**
	 * Changes the state to "loading/saving was canceled" depending on previous state.
	 */
	protected void setCanceled() {
		completionTime = System.currentTimeMillis();
		resetMessageDetails();

		switch (state) {
		case LOADING:
			setLoadingCanceled();
			break;

		case SAVING:
			setSavingCanceled();
			break;

		default:
			setOff();
			break;
		}
	}

	/**
	 * Changes the state to "loading/saving has failed" depending on previous state.
	 */
	protected void setFailed() {
		completionTime = System.currentTimeMillis();
		// don't call resetMessageDetails(); here this method is called by setFailed(String)

		switch (state) {
		case LOADING:
			setLoadingFailed();
			break;

		case SAVING:
			setSavingFailed();
			break;

		default:
			setOff();
			break;
		}
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
		switch (state) {

		case OFF:
			return;

		case LOADING:
			updateDots();
			showMessage(Message.LOADING, colorWorkInProgress, colorLabelBackground, dots);
			return;

		case LOADING_COMPLETED:
			showFadingMessage(Message.LOADING_COMPLETED, colorCompleted);
			return;

		case LOADING_CANCELED:
			showFadingMessage(Message.LOADING_CANCELED, colorFailed);
			return;

		case LOADING_FAILED:
			showFadingMessage(Message.LOADING_FAILED, colorFailed);
			return;

		case SAVING:
			updateDots();
			showMessage(Message.SAVING, colorWorkInProgress, colorLabelBackground, dots);
			return;

		case SAVING_COMPLETED:
			showFadingMessage(Message.SAVING_COMPLETED, colorCompleted);
			return;

		case SAVING_CANCELED:
			showFadingMessage(Message.SAVING_CANCELED, colorFailed);
			return;

		case SAVING_FAILED:
			showFadingMessage(Message.SAVING_FAILED, colorFailed);
			return;

		default:
			setOff(); // safety off-switch
			return;
		}
	}
	
	protected void paintComponent(Graphics g) {
		Graphics2D g2d = (Graphics2D) g;
		// if the parent enables AA, we can do it here to... does not seem necessary, though
//		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, alpha));
	}
}

/**
 * This class wraps {@code GuiMessageOverlay} for debugging and shows various information in stderr.
 * NEVER USE THIS CLASS IN A PRODUCTION RELEASE!!!
 */
class GuiMessageOverlayLogged extends GuiMessageOverlay {
	private static final long serialVersionUID = -8396089731118502477L;

	private int message_count;
	
	public GuiMessageOverlayLogged(JPanel parent) {
		super(parent);
		
		message_count = 0;
	}

	/**
	 * Shows state changes.
	 */
	protected void setState(State state) {
		super.setState(state);
		System.err.println("state is now: " + Message.valueOf(state.toString()));
		message_count = 0;
	}
	
	/**
	 * Shows number of message refreshs.
	 */
	protected void showMessage(Message message, Color color, Color backgroundColor, String workInProgressDots) {
		super.showMessage(message, color, backgroundColor, workInProgressDots);
		
		message_count++;
		System.err.println("\tmessage refresh No. " + message_count);
	}
}
