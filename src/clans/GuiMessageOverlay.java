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

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.Timer;
import javax.swing.UIManager;

public class GuiMessageOverlay extends JComponent {
	private static final long serialVersionUID = 818220626180736376L;

	/**
	 * States of the message overlay, e.g. OFF, LOADING, SAVING_COMPLETED.
	 */
	enum State {
		OFF,
		CUSTOM,
		
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
		OFF("MESSAGE OVERLAY OFF"),
		
		CUSTOM("CUSTOM MESSAGE SET"),
		
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
		WARNING(3500),
		ERROR(4500),
		INFINITE(Integer.MAX_VALUE);
		
		private int value;
		
		private Duration(int value) {
			this.value = value;
		}
		
		protected int get() {
			return value;
		}
	}

	private Component top_spacer; // moves the labels down
	private JLabel mainLabel; // shows the main message
	private JLabel detailsLabel; // shows messageDetails
	private int labelInsetTop;
	private int labelInsetBottom;
	private int labelInsetLeftRight;
	private Color colorLabelBackground;

	private Timer timer; // draw only every X milliseconds instead of using up one CPU for drawing continuously.


	private long startTime; // time of last state change
	private State state; // the current state
	private String currentMessage; // the details of the message, e.g. what went wrong during loading
	private String messageDetails; // the details of the message, e.g. what went wrong during loading
	private Color currentColor; // the details of the message, e.g. what went wrong during loading
	private Duration currentDuration; // the duration for the current mode
	private boolean currentIsFading; // whether the current message fades or abruptly vanishes
	private boolean currentHasProgressDots; // whether the current message gets appended progress dots

	private String currentDots; // the progress dots on loading and saving
	private int dotTicks;
	final private int millisecondsPerDot = 350;
	
	private float alpha;
	
	private Color colorWorkInProgress; // color for work-in-progress operations
	private Color colorCompleted; // color for completed operations
	private Color colorFailed; // color for failed operations

	public GuiMessageOverlay() {

		this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		
		colorLabelBackground = new Color(128, 128, 128);
		
		labelInsetTop = 4;
		labelInsetBottom = 2;
		labelInsetLeftRight = 4;
				
		mainLabel = new JLabel();
		mainLabel.setFont(UIManager.getFont("Label.font").deriveFont(Font.BOLD).deriveFont(20f));
		mainLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		mainLabel.setOpaque(true);
		mainLabel.setBackground(colorLabelBackground);
		mainLabel.setBorder(BorderFactory.createEmptyBorder(labelInsetTop, labelInsetLeftRight, labelInsetBottom,
				labelInsetLeftRight));
		
		detailsLabel = new JLabel();
		detailsLabel.setFont(UIManager.getFont("Label.font").deriveFont(Font.BOLD).deriveFont(16f));
		detailsLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		detailsLabel.setOpaque(true);
		detailsLabel.setBackground(colorLabelBackground);

		// the initial values for top_spacer ARE NOT USED as we overwrite them in updateGlassSize() below
		top_spacer = Box.createRigidArea(new Dimension(0, 300));
		this.add(top_spacer);
		this.add(mainLabel);
//		this.add(detailsLabel);

		colorWorkInProgress = new Color(255, 128, 0);
		colorCompleted = new Color(128, 255, 0);
		colorFailed = new Color(255, 102, 137);

		int timer_delay = 100; // in milliseconds
		setupTimer(timer_delay);
	
		setOff();
	}

	protected void setState(State state) {
		this.state = state;
		startTime = System.currentTimeMillis();
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
	 * Checks if and which type of message is waiting to be displayed and initiates its display. Does nothing it no
	 * messages is about to be displayed. Disables the message overlay if message duration is over.
	 */
	private void invoke() {

		// time since start is also relevant for fading, hence we forward it for faded message (see switch below)
		long time_since_start = System.currentTimeMillis() - startTime;
		if (time_since_start > currentDuration.get()) {
			setOff();
			return;
		}

		if (state == State.OFF) {
			return;
		}

		resetTransparency();

		if (currentIsFading) {
			showFadingMessage(time_since_start);
		} else {
			showMessage();
		}
	}
	
	/**
	 * Sets the main message.
	 * @param message The main message text.
	 */
	private void setMessage(String message){
		this.currentMessage = message;
	}
	
	/**
	 * Resets message to none.
	 */
	private void resetMessage() {
		currentMessage = "";
	}
	
	/**
	 * Sets message details.
	 */
	private void setMessageDetails(String details) {
		messageDetails = details;
		detailsLabel.setBorder(BorderFactory.createEmptyBorder(labelInsetTop, labelInsetLeftRight, labelInsetBottom,
				labelInsetLeftRight));
		this.add(detailsLabel);
	}
	
	/**
	 * Resets detailed message to none.
	 */
	private void resetMessageDetails() {
		messageDetails = "";
		this.remove(detailsLabel);
	}

	private void setColor(Color color) {
		this.currentColor = color;
	}
	
	private void resetColor() {
		this.currentColor = colorCompleted;
	}
	
	private void setDuration(Duration duration) {
		this.currentDuration = duration;
	}
	
	private void resetDuration() {
		this.currentDuration = Duration.INFINITE;
	}
	
	private void enableFading() {
		this.currentIsFading = true;
	}
	
	private void disableFading() {
		this.currentIsFading = false;
	}
	
	private void enableProgressDots() {
		currentDots = "";
		this.currentHasProgressDots= true;
	}
	
	private void disableProgressDots() {
		currentDots = "";
		this.currentHasProgressDots = false;
	}

	/**
	 * Resets the transparency of the glass pane.
	 */
	private void resetTransparency() {
		alpha = 1;
	}
	
	/**
	 * Sets the correct number of dots for recurring calls to indicate work in progress.
	 */
	private void updateDots() {
		if (state == State.LOADING || state == State.SAVING) {
			dotTicks += 1;
			
			if (dotTicks * timer.getDelay() > millisecondsPerDot) {
				if (currentDots.length() < 3) {
					currentDots += ".";
				} else {
					currentDots = "";
				}
			
				dotTicks = 0;
			} 
		}
	}
	
	/**
	 * Switches the overlay off and resets its state. Used once in constructor.
	 */
	private void setOff() {
		setVisible(false);
		
		if (timer != null) {
			timer.stop();
		}

		resetMessage();
		resetMessageDetails();
		resetColor();
		resetDuration();
		
		disableFading();
		disableProgressDots();
		resetTransparency();
		
		setState(State.OFF);
	}

	/**
	 * Changes the state to "loading has completed".
	 */
	private void setLoadingCompleted() {
		setupStandardMessage(State.LOADING_COMPLETED);
	}

	/**
	 * Changes the state to "loading was canceled".
	 */
	private void setLoadingCanceled() {
		setupStandardMessage(State.LOADING_CANCELED);
	}

	/**
	 * Changes the state to "loading has failed".
	 */
	private void setLoadingFailed() {
		setupStandardMessage(State.LOADING_FAILED);
	}

	/**
	 * Changes the state to "saving has completed".
	 */
	private void setSavingCompleted() {
		setupStandardMessage(State.SAVING_COMPLETED);
	}

	/**
	 * Changes the state to "saving was canceled".
	 */
	private void setSavingCanceled() {
		setupStandardMessage(State.SAVING_CANCELED);
	}

	/**
	 * Changes the state to "saving has failed".
	 */
	private void setSavingFailed() {
		setupStandardMessage(State.SAVING_FAILED);
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
	protected void showMessage() {
		mainLabel.setForeground(currentColor);
		
		
		if (currentHasProgressDots) {
			updateDots();
			mainLabel.setText(currentMessage + currentDots);
		
		} else {
			mainLabel.setText(currentMessage);
		}
		
		detailsLabel.setForeground(currentColor);
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
	private void showFadingMessage(long time_since_start) {
		alpha = 1 - getFadingValue(time_since_start); // fade by reducing visibility for the whole glass pane; see paintComponent below.
		showMessage();
	}

	/**
	 * Computes the fading value. After an initial delay of (Message duration - {@code Duration.FADING}) seconds in
	 * which no fading (==0) is returned, the fraction of time already elapsed w.r.t. the complete message duration is
	 * returned. The speed of fading from 0 -> 1 values is linear.
	 * 
	 * @return Fraction of duration already passed as value in [0, 1]
	 */
	private float getFadingValue(long time_since_start) {
		
		float timespan_without_fading = currentDuration.get() - Duration.FADING.get();
		
		if (time_since_start < timespan_without_fading) { // make the initial delay period non-transparent
			return 0;
		}

		// fade during duration - Durations.DELAY and make the subtraction negative number proof
		float fade_fraction = (float) Math.max(time_since_start - timespan_without_fading, 0) / Duration.FADING.get();
		return fade_fraction;
	}

	/**
	 * Adjusts the label placement when the glass's parent is resized. Used once in constructor.
	 */
	protected void updateGlassSize(int new_height) {
		remove(0); // top_spacer is the first entry

		top_spacer = Box.createRigidArea(new Dimension(0, new_height / 3 * 2));

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

	private void setupMessage(State state, String main_message, String message_details, Color color,
			Duration duration, boolean is_fading, boolean has_progress_dots) {

		setState(state);
		
		setMessage(main_message);
		
		if (message_details == null) {
			resetMessageDetails();
		} else {
			setMessageDetails(message_details);
		}
		
		setColor(color);
		
		setDuration(duration);
		
		if (is_fading) {
			enableFading();
		}
		
		if (is_fading) {
			enableFading();
		} else {
			disableFading();
		}
		
		if (has_progress_dots) {
			enableProgressDots();
		} else {
			disableProgressDots();
		}
		
		activate_overlay();
	}
	
	private void setupProgressMessage(State state, String main_message, String message_details) {
		setupMessage(state, main_message, message_details, colorWorkInProgress, Duration.INFINITE, false, true);
	}
	
	private void setupFadingMessage(State state, String main_message, String message_details, Color color,
			Duration duration) {
		setupMessage(state, main_message, message_details, color, duration, true, false);
	}

	private void setupStandardMessage(State state) {
		String main_message = Message.valueOf(state.toString()).get();

		switch (state) {
		case OFF:
			return;
			
		case CUSTOM:
			// custom message are set up by their own method
			return;
			
		case LOADING:
		case SAVING:
			setupProgressMessage(state, main_message, null);
			return;
			
		case LOADING_CANCELED:
		case SAVING_CANCELED:
			setupFadingMessage(state, main_message, null, colorWorkInProgress, Duration.INFO);
			return;
			
		case LOADING_COMPLETED:
		case SAVING_COMPLETED:
			setupFadingMessage(state, main_message, null, colorCompleted, Duration.INFO);
			return;

		case LOADING_FAILED:
		case SAVING_FAILED:
			setupFadingMessage(state, main_message, null, colorFailed, Duration.ERROR);
			return;

		default:
			break;
		}
	}

	protected void paintComponent(Graphics g) {
		Graphics2D g2d = (Graphics2D) g;
		// if the parent enables AA, we can do it here to... does not seem necessary, though
//		g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, alpha));
	}
	
	protected Color getColorCompleted() {
		return colorCompleted;
	}
	
	protected Color getColorDefault() {
		return colorWorkInProgress;
	}
	
	protected Color getColorError() {
		return colorFailed;
	}
	
	
	protected void setCustomMessage(String main_message, String message_details, Color color, Duration duration,
			boolean is_fading, boolean with_progress_dots) {
		setupMessage(State.CUSTOM, main_message, message_details, color, duration, is_fading, with_progress_dots);
	}
	
	/**
	 * Changes the state to "currently loading".
	 */
	protected void setLoading() {
		setupStandardMessage(State.LOADING);
	}

	/**
	 * Changes the state to "currently saving".
	 */
	protected void setSaving() {
		setupStandardMessage(State.SAVING);
	}

	/**
	 * Changes the state to "loading/saving has completed" depending on previous state.
	 */
	protected void setCompleted() {
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
	 * 
	 * Changes the state to "loading/saving was canceled" depending on previous state.
	 */
	protected void setCanceled() {
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
		setMessageDetails(details);
		setFailed();
	}
}

/**
 * This class wraps {@code GuiMessageOverlay} for debugging and shows various information in stderr.
 * NEVER USE THIS CLASS IN A PRODUCTION RELEASE!!!
 */
class GuiMessageOverlayLogged extends GuiMessageOverlay {
	private static final long serialVersionUID = -8396089731118502477L;

	private int message_count;
	private long last_message_time;
	
	public GuiMessageOverlayLogged() {
		super();
		message_count = 0;
		last_message_time = 0;
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
	 * Shows statistic on message frequency.
	 */
	protected void showMessage() {
		super.showMessage();
		
		long current_time = System.currentTimeMillis();
		System.err.println("\ttime since last message: " + (current_time - last_message_time));
		last_message_time = current_time;
		
		message_count++;
		System.err.println("\tmessage refresh No. " + message_count);
	}
}
